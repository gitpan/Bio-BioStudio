#!/usr/bin/env perl

use Bio::BioStudio;
use Bio::BioStudio::ConfigData;
use Bio::GeneDesign::Basic qw(_compare_sequences);
use Bio::Tools::Run::StandAloneBlastPlus;
use Digest::MD5;
use Time::Format qw(%time);
use Getopt::Long;
use Pod::Usage;
use Carp;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_PCRTagger_$VERSION";

local $| = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'      => \$p{CHROMOSOME},
      'EDITOR=s'          => \$p{EDITOR},
      'MEMO=s'            => \$p{MEMO},
      'SCALE=s'           => \$p{SCALE},
      'SCOPE=s'           => \$p{SCOPE},
      'STARTPOS=i'        => \$p{STARTPOS},
      'STOPPOS=i'         => \$p{STOPPOS},
      'MINTAGMELT=i'      => \$p{MINTAGMELT},
      'MAXTAGMELT=i'      => \$p{MAXTAGMELT},
      'MINTAGLEN=i'       => \$p{MINTAGLEN},
      'MAXTAGLEN=i'       => \$p{MAXTAGLEN},
      'MINAMPLEN=i'       => \$p{MINAMPLEN},
      'MAXAMPLEN=i'       => \$p{MAXAMPLEN},
      'MAXAMPOLAP=i'      => \$p{MAXAMPOLAP},
      'MINPERDIFF=i'      => \$p{MINPERDIFF},
      'MINORFLEN=i'       => \$p{MINORFLEN},
      'FIVEPRIMESTART=i'  => \$p{FIVEPRIMESTART},
      'MINRSCUVAL=f'      => \$p{MINRSCUVAL},
      'OUTPUT=s'          => \$p{OUTPUT},
			'help'              => \$p{HELP}
);
$p{BSVERSION} = $bsversion;

################################################################################
################################ SANITY  CHECK #################################
################################################################################
pod2usage(-verbose=>99, -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE")
  if ($p{HELP});

my $bs = Bio::BioStudio::ConfigData->config('blast_support');
if ($bs ne 'Y')
{
  die "Cannot make PCR tags: blast support is not enabled\n";
}
  
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n"  unless ($p{CHROMOSOME});
my $chr    = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $chrseq = $chr->sequence;
my $chrlen = length $chrseq;

$p{SCALE}            = $p{SCALE}            || "chrom";
$p{SCOPE}            = $p{SCOPE}            || "chrom";
$p{OUTPUT}           = $p{OUTPUT}           || "txt";
$p{MINTAGMELT}       = $p{MINTAGMELT}       || 57.9;
$p{MAXTAGMELT}       = $p{MAXTAGMELT}       || 60.9;
$p{MINPERDIFF}       = $p{MINPERDIFF}       || 33;
$p{MINTAGLEN}        = $p{MINTAGLEN}        || 19;
$p{MAXTAGLEN}        = $p{MAXTAGLEN}        || 28;
$p{MINAMPLEN}        = $p{MINAMPLEN}        || 200;
$p{MAXAMPLEN}        = $p{MAXAMPLEN}        || 500;
$p{MAXAMPOLAP}       = $p{MAXAMPOLAP}       || 25;
$p{MINORFLEN}        = $p{MINORFLEN}        || 501;
$p{FIVEPRIMEBUFFER}  = $p{FIVEPRIMEBUFFER}  || 101;
$p{THREEPRIMEBUFFER} = $p{THREEPRIMEBUFFER} || 0;
$p{MINRSCUVAL}       = $p{MINRSCUVAL}       || .04;
$p{STARTPOS}         = $p{STARTPOS}         || 1;
$p{STOPPOS}          = $p{STOPPOS}          || $chrlen;
$p{ORF_TAG_INC}      = $p{ORF_TAG_INC}      || 1000;

$p{FIVEPRIMEBUFFER}++ while ($p{FIVEPRIMEBUFFER} % 3 != 0);
$p{THREEPRIMEBUFFER}++ while ($p{THREEPRIMEBUFFER} % 3 != 0);

unless ($p{EDITOR} && $p{MEMO})
{
  print "\n ERROR: Both an editor's id and a memo must be supplied.\n\n";
}

if ($p{SCOPE} eq "seg" && ($p{STOPPOS} <= $p{STARTPOS}))
{
  die "\n ERROR: The start and stop coordinates do not parse.\n\n";
}

if ( $p{MAXTAGMELT} < $p{MINTAGMELT}
  || $p{MAXTAGMELT} >= 90 || $p{MINTAGMELT} <= 20)
{
  die "\n ERROR: The tag melting parameters do not parse.\n\n";
}

if ($p{MAXTAGLEN} < $p{MINTAGLEN} || $p{MINTAGLEN} <= 0 ||
  ((($p{MAXTAGLEN} - 1) % 3 != 0 ) || (($p{MINTAGLEN} - 1) % 3 != 0)))
{
  die "\n ERROR: The tag length parameters do not parse. Tag lengths must be "
    . "multiples of 3 plus 1.\n\n";
}

if ($p{MAXAMPLEN} < $p{MINAMPLEN} || $p{MINAMPLEN} <= 0)
{
  die "\n ERROR: The amplicon length parameters do not parse.\n\n";
}

if ($p{MINPERDIFF} > 66 || $p{MINPERDIFF} <= 20 )
{
  die "\n ERROR: The minimum percent difference does not parse.\n\n";
}

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $newchr = $chr->iterate();
my $mask = Bio::BioStudio::Mask->new(-sequence => $newchr->seqobj);
my $GD     = $newchr->GD;

my $genome = $BS->gather_latest_genome($chr->species);
my $factory = _make_BLAST_factory($genome);

my %report = ();

my @CDSes = $newchr->db->features(
  -seq_id     => $newchr->seq_id,
  -types      => 'CDS',
);
$mask->add_to_mask(\@CDSes);

my @introns = $newchr->db->features(
  -seq_id     => $newchr->seq_id,
  -types      => 'intron',
);
$mask->add_to_mask(\@introns);
my %intronlist = map {$_->display_name => $_} @introns;

my @genes = $newchr->db->features(
  -seq_id     => $newchr->seq_id,
  -types      => 'gene',
  -start      => $p{STARTPOS},
  -end        => $p{STOPPOS},
  -range_type => 'contains'
);

#Single family codons can't be first or last codons
#codons that don't share their siblings' first two bases can't be first codons
my %fams;
my %dicodons;
my $codon_t = $GD->codontable;
foreach my $cod (keys %{$codon_t})
{
  my $aa = $codon_t->{$cod};
  $fams{$aa} = [] if (! exists $fams{$aa});
  $dicodons{$aa} = {} if (! exists $dicodons{$aa});
  push @{$fams{$aa}}, $cod if ($GD->rscutable->{$cod} > $p{MINRSCUVAL});
  $dicodons{$aa}->{substr($cod, 0, 2)}++;
}
my %badfirstaas = map {$_ => 1} grep { scalar(@{$fams{$_}}) == 1 } keys %fams;
foreach my $aa (grep {scalar (@{$fams{$_}}) > 1} keys %dicodons)
{
  my $flag = 0;
  foreach (keys %{$dicodons{$aa}})
  {
    $flag++ if ($dicodons{$aa}->{$_} == 1);
  }
  $badfirstaas{$aa}++ if ($flag != 0);
}
#print "bads: ", keys %badfirstaas, "\n\n";


#Find tags, pick tags, implement tags, check for errors
print "Parsing genes for tags:\n";
foreach my $gene (@genes)
{
  my $gid  = $gene->Tag_load_id;
  my $cDNA = $newchr->make_cDNA($gene);
  my $dir = $gene->strand;
  
  #filter genes that are too small
  if (length $cDNA < $p{MINORFLEN})
  {
    $report{$gid} = 'too short';
    next;
  }
  my ($gstart, $gstop) = ($gene->start, $gene->end);
  if ($dir != 1)
  {
    $gstart += $p{THREEPRIMEBUFFER};
    $gstop  -= $p{FIVEPRIMEBUFFER};
  }
  else
  {
    $gstart += $p{FIVEPRIMEBUFFER};
    $gstop  -= $p{THREEPRIMEBUFFER};
  }
  my $rstart = $gstart;
  
  my $geneseq = $gene->seq->seq;
  my $aaseq = $GD->translate(-sequence => $cDNA);
    
  #Extract all possible tags from the gene
  #
  my @pretags = ();
  my $count = $p{MAXTAGLEN};
  while ($rstart + $count < $gstop - $p{MINTAGLEN})
  {
    #Make the wide oligo and the oligo
    my $woligo = substr ($chrseq, $rstart - 1, $count + 2);
    my $wolen = length($woligo);

    my $tstart = $dir == 1  ? $rstart + 2 : $rstart;
    my $oligo = substr ($chrseq, $tstart - 1, $count);
    my $olen = length($oligo);
    
    #Exclude oligos that begin or end in unswappable codons
    my $peptide = $GD->translate(-sequence => $woligo, -frame => $dir);
    my $firstres = substr($peptide,  0, 1);
    my $lastres  = substr($peptide, -1, 1);
    if (exists $badfirstaas{$firstres} || exists $badfirstaas{$lastres})
    {
      ($count, $rstart) = counterset($count, $rstart);
#      print "\t\t\tFAIL aa $firstres-$lastres\n";
      next;
    }
    
    #Exclude oligos that don't meet melting standard
    my $currTm = $GD->melt(-sequence => $oligo, -nearest_neighbor => 1);
    if ($currTm  < $p{MINTAGMELT} || $currTm > $p{MAXTAGMELT})
    {
      if ($currTm < $p{MINTAGMELT})
      {
        $rstart = $rstart + 3;
        $count = $p{MAXTAGLEN};
      }
      elsif ($currTm > $p{MAXTAGMELT})
      {
        ($count, $rstart) = counterset($count, $rstart)
      }
#      print "\t\t\tFAIL Tm $currTm\n";
      next;
    }
        
    #Exclude oligos that aren't entirely in exons or lapped by other genes
    my @ofeats = $mask->features_in_range($rstart, $count);
    my @ifeats = grep { exists $intronlist{$_} } @ofeats;
    if (scalar @ifeats || scalar @ofeats > 1)
    {
      ($count, $rstart) = counterset($count, $rstart);
#      print "\t\t\tFAIL mask\n";
      next;
    }

    #Recode oligos
    $woligo = $GD->complement($woligo, 1) if ($dir == -1);
    my $wmdoligo = $GD->codon_juggle(
      -sequence => $woligo,
      -algorithm => 'most_different_sequence'
    );
    
    #Check that the first two bases are the same
    if (substr($wmdoligo, 0, 2) ne substr($woligo, 0, 2))
    {
      my $codon = substr($woligo, 0, 3);
      my $aa = $GD->codontable->{$codon};
      my $di = substr($codon, 0, 2);
      my $possibles = $GD->reversecodontable->{$aa};
      my @choices = grep {substr($_, 0, 2) eq $di} @{$possibles};
      @choices = grep {substr($_, -1) ne substr($codon, -1)} @choices;
      substr($wmdoligo, 0, 3, $choices[0]);
    }
    $woligo = $GD->complement($woligo, 1) if ($dir == -1);
    $wmdoligo = $GD->complement($wmdoligo, 1) if ($dir == -1);
    my $mdoligo = $dir == 1
      ? substr($wmdoligo, 2, $wolen - 2)
      : substr($wmdoligo, 0, $wolen - 2);
    
    
   
    #Exclude oligos whose recodes don't meet percent difference standard
    my $comps = _compare_sequences($oligo, $mdoligo);
    if ( $comps->{'P'} < $p{MINPERDIFF})
    {
      ($count, $rstart) = counterset($count, $rstart);
#      print "\t\t\tFAIL diff " . $comps->{'P'} . "\n";
      next;
    }
    
    #Exclude oligos whose recodes don't meet melting standard
    my $MDTm = $GD->melt(-sequence => $mdoligo, -nearest_neighbor => 1);
    if ($MDTm  < $p{MINTAGMELT} || $MDTm > $p{MAXTAGMELT})
    {
      ($count, $rstart) = counterset($count, $rstart);
#      print "\t\t\tFAIL mdTm $MDTm\n";
      next;
    }
    
    my $offset = $tstart - $gstart + 1;
    my $tagid = $gid . "_" . $offset;
    my $tag = Bio::SeqFeature::Generic->new(
      -start          => $offset,
      -end            => $offset + $olen - 1,
      -display_name   => $tagid,
      -primary_tag    => "tag",
      -tag            => {
        wtseq         => $oligo,
        newseq        => $mdoligo,
        wtpos         => $tstart,
        ingene        => $gid,
        difference    => $comps->{'P'},
        translation   => $peptide
      }
    );
    push @pretags, $tag;

    #print "\t$rstart\t$tstart\n";
    #print "\t\t$woligo\n\t\t$oligo\n";
    #print "\t\t$wmdoligo\n\t\t$mdoligo\n";
    #print "\t\t $peptide\n";
    
    ($count, $rstart) = counterset($count, $rstart);
  }
  
  if (scalar @pretags < 2)
  {
    $report{$gid} = 'no more than one good tag preBLAST';
    next;
  }
  
  #BLAST all tags against the latest version of the genome
  # only keep tags that hit wtseq once and newseq never
  my @tagobjs = ();
  foreach my $tag (@pretags)
  {
    my $id = $tag->display_name;
    my $wtseq  = join q{}, $tag->get_tag_values('wtseq');
    my $newseq = join q{}, $tag->get_tag_values('newseq');
    my $wtobj  = Bio::Seq->new(-seq => $wtseq,  id => 'wt_' . $id);
    my $newobj = Bio::Seq->new(-seq => $newseq, id => 'md_' . $id);
    push @tagobjs, $wtobj, $newobj;
  }
  
  $factory->run(
    -method           => "blastn",
    -query            => \@tagobjs,
    -method_args => [
      -word_size      => 17,
      -perc_identity  => 70
  ]);
  my %hits;
  while (my $result = $factory->next_result)
  {
    my $name = $result->query_name();
    while( my $hit = $result->next_hit())
    {
      while( my $hsp = $hit->next_hsp())
      {
        $hits{$name}++;
      }
    }
  }
  
  my @posttags;
  foreach my $tag (sort {$a->start <=> $b->start} @pretags)
  {
    my $id = $tag->display_name;
    my $wthits = $hits{'wt_' . $id} || 0;
    my $mdhits = $hits{'md_' . $id} || 0;
    next if ($wthits > 1 || $mdhits > 0);
    push @posttags, $tag;
  }
  if (scalar @posttags < 2)
  {
    $report{$gid} = 'no more than one good tag postBLAST';
    next;
  }
  if($posttags[0]->start + $p{MINAMPLEN} > $posttags[-1]->end)
  {
    $report{$gid} = 'bad tag range postBLAST';
    next;
  }

  #choose tag pairs
  # 
  my $tagtarget = (length($cDNA) - $p{MINORFLEN});
  $tagtarget = $tagtarget / $p{ORF_TAG_INC};
  $tagtarget = ($tagtarget % $p{ORF_TAG_INC}) + 1;
  my $tagcount = 0;
  my @chosen;
  my %usedtags = ();
  my $tmask = Bio::BioStudio::Mask->new(-sequence => $gene);
  my $amask = Bio::BioStudio::Mask->new(-sequence => $gene);

  foreach my $utag (sort sortdiff @posttags)
  {
    my $ustart = $utag->start;
    my $partner = undef;
    my $utagid = $utag->display_name;
    my $len = $utag->end - $ustart + 1;
    my $mcount = $tmask->count_features_in_range($ustart, $len);
    next if ($mcount > 0);
    
    my %possibles;
    my @pool = grep {
      ! exists $usedtags{$_} &&
      (  ($_->end - $ustart <= $p{MAXAMPLEN} && $_->end - $ustart >= $p{MINAMPLEN})
      || ($utag->end - $_->start <= $p{MAXAMPLEN} && $utag->end - $_->start >= $p{MINAMPLEN}))
      } @posttags;
    foreach my $dtag (@pool)
    {
      my $dstart = $dtag->start;
      my $dtagid = $dtag->display_name;
      my $tstart = $ustart < $dstart  ? $ustart  : $dstart;
      my $tend   = $ustart < $dstart  ? $dtag->end    : $utag->end;
      my $pampsize = $tend - $tstart + 1;
      
      #nominal size of the amplicon, with and without intron removal
      my @feats = $mask->features_in_range($tstart, $pampsize);
      my @ifeats = grep { exists $intronlist{$_} } @feats;
      my $intronbp = 0;
      foreach my $ifeat (@ifeats)
      {
        $intronbp += $intronlist{$ifeat}->length();
      }
      my $ampsize = $pampsize - $intronbp;
      next if ($ampsize > $p{MAXAMPLEN} || $ampsize < $p{MINAMPLEN});
     
      #would this tag overlap a previously chosen tag
      my $dtaglen = $dtag->end - $dstart + 1;
      my $dcount = $tmask->count_features_in_range($dstart, $dtaglen);
      next if ($dcount > 0);
     
      #how badly would would overlap an existing amplicon
      my %occs = $amask->occlusion($tstart, $pampsize);
      my $occflag = 0;
      my $maxolap = 0;
      foreach my $occ (keys %occs)
      {
        my $olap = $occs{$occ};
        $occflag++ if ($olap > ($p{MAXAMPOLAP} / 100));
        $maxolap = $olap if ($olap > $maxolap);
      }
      next if ($occflag > 0);
      
      my $diffd = join q{}, $dtag->get_tag_values('difference');
      $possibles{$dtagid} = [$maxolap, $diffd, $dtag];
    }
    if (! scalar(keys %possibles))
    {
      next;
    }
    my @downstreams = sort {  $possibles{$a}->[0] <=> $possibles{$b}->[0]
                           || $possibles{$a}->[1] <=> $possibles{$b}->[1] }
                      keys %possibles;
    my @downbuds = map {$possibles{$_}->[2]} @downstreams;
    
    #BLAST the pair to check amplification uniqueness.
    #
    my %phits = ();
    my %objs = ();
    # The utag wt and new sequences, forwards
    my $uid = $utag->display_name;
    my $uwtid = 'rcwt_' . $uid;
    my $uwtseq  = join q{}, $utag->get_tag_values('wtseq');
    $phits{$uwtid} = [];
    $objs{$uwtid} = Bio::Seq->new(-seq => $uwtseq,  id => $uwtid);
    my $unewid = 'rcmd_' . $uid;
    my $unewseq = join q{}, $utag->get_tag_values('newseq');
    $phits{$unewid} = [];
    $objs{$unewid} = Bio::Seq->new(-seq => $unewseq, id => $unewid);
    # The dtags wt and new sequences, reverses
    my @dtagids = ();
    foreach my $dtag (@downbuds)
    {
      my $did = $dtag->display_name;
      push @dtagids, $did;
      my $dwtid = 'rcwt_' . $did;
      my $dwtseq  = join q{}, $dtag->get_tag_values('wtseq');
      $phits{$dwtid} = [];
      $objs{$dwtid} = Bio::Seq->new(-seq => $dwtseq,  id => $dwtid);
      my $dnewid = 'rcmd_' . $did;
      my $dnewseq = join q{}, $dtag->get_tag_values('newseq');
      $phits{$dnewid} = [];
      $objs{$dnewid} = Bio::Seq->new(-seq => $dnewseq, id => $dnewid);
    }
    
    my @pairobjs = values %objs;
    $factory->run(
      -method               => 'blastn',
      -query                => \@pairobjs,
      -method_args => [
        -word_size          => 4,
        -gapextend          => 2,
        -gapopen            => 1,
        -penalty            => -1,
        -perc_identity      => 70,
        -best_hit_overhang  => .25
    ]);
    while (my $result = $factory->next_result)
    {
      my $name = $result->query_name();
      my $qlen = length $objs{$name}->seq;
      while( my $hit = $result->next_hit())
      {
        while( my $hsp = $hit->next_hsp())
        {
          next if $hsp->percent_identity < 85;
          next if ($hsp->length < .7 * $qlen);
          push @{$phits{$name}}, [$hit, $hsp];
        }
      }
    }
    my @Fwts = @{$phits{'rcwt_' . $uid}};
    my @Fmds = @{$phits{'rcmd_' . $uid}};
    foreach my $dtagid (@dtagids)
    {
      my ($wtnogo, $mdnogo) = (0, 0);
      my @Rwts = @{$phits{'rcwt_' . $dtagid}};
      foreach my $Rwpair (@Rwts)
      {
        my ($Rwhit, $Rwhsp) = @{$Rwpair};
        foreach my $Fwpair (@Fwts)
        {
          my ($Fwhit, $Fwhsp) = @{$Fwpair};
          #Pass if hits are on different chromosomes
          next if ($Rwhit->name ne $Fwhit->name);
          #Pass if hits are on the same strand
          next if ($Rwhsp->strand eq $Fwhsp->strand);
          $wtnogo++;
        }
      }
      my @Rmds = @{$phits{'rcmd_' . $dtagid}};
      foreach my $Rmpair (@Rmds)
      {
        my ($Rmhit, $Rmhsp) = @{$Rmpair};
        foreach my $Fmpair (@Fmds)
        {
          my ($Fmhit, $Fmhsp) = @{$Fmpair};
          #Pass if hits are on different chromosomes
          next if ($Rmhit->name ne $Fmhit->name);
          #Pass if hits are on the same strand
          next if ($Rmhsp->strand eq $Fmhsp->strand);
          $mdnogo++;
        }
      }
      if ($mdnogo + $wtnogo == 0)
      {
        $partner = $possibles{$dtagid}->[2];
        last;
      }
    }
    if ($partner)
    {
      $tagcount++;
      $usedtags{$utag}++;
      $usedtags{$partner}++;
      $tmask->add_to_mask([$utag, $partner]);
      my ($tstart, $tend) = (undef, undef);
      if ($ustart < $partner->start)
      {
        $tstart = $ustart;
        $tend = $partner->end;
        push @chosen, [$utag, $partner];
      }
      else
      {
        $tstart = $partner->start;
        $tend = $utag->end;
        push @chosen, [$partner, $utag];
      }
      my $amp = Bio::SeqFeature::Generic->new(
        -start          => $tstart,
        -end            => $tend,
        -display_name   => "anonamp",
        -primary_tag    => "amplicon",
      );
      $amask->add_to_mask([$amp]);
    }
    last if ($tagcount == $tagtarget);
  }

  if (scalar(@chosen) == 0)
  {
    $report{$gid} = 'no good tag pairs';
    next;
  }
  if ($tagcount < $tagtarget)
  {
    $report{$gid} .= "$tagcount of $tagtarget pairs chosen;";
  }
  else
  {
    $report{$gid} .= "$tagtarget pairs chosen;";
  }
  
  #Add tags and amplicons to database, make sequence changes
  my $tcount = 1;
  foreach my $tagpair (@chosen)
  {
    my ($gutag, $gdtag) = @{$tagpair};
    my $adjust = $gstart - 1;
    $gutag->start($gutag->start() + $adjust);
    $gdtag->start($gdtag->start() + $adjust);
    $gutag->end($gutag->end() + $adjust);
    $gdtag->end($gdtag->end() + $adjust);
    $newchr->add_feature(-feature => $gutag);
    $newchr->add_feature(-feature => $gdtag);
    
    my $utagid = $gutag->display_name;
    my $dtagid = $gdtag->display_name;
    my $aid = $gid . "_amp" . $tcount;
    my $comment  = "PCR_product $aid annotated ";
    $comment .= "(tags $utagid and $dtagid added)\n";
    
    my $atts = {ingene => $gid, uptag => $utagid, dntag => $dtagid};
    my $amp = Bio::SeqFeature::Generic->new(
      -start          => $gutag->start,
      -end            => $gdtag->end,
      -primary_tag    => "PCR_product",
      -display_name   => $aid,
    );
    $newchr->add_feature(
        -feature => $amp,
        -comments => [$comment],
        -attributes => $atts
    );
    $tcount++;
  }
  
  #Do error checking
  my $glen = $gene->end - $gstart + 1;
  my $newgeneseq = substr($newchr->sequence, $gstart - 1, $glen);
  my $oldgeneseq = substr($chrseq, $gstart - 1, $glen);
  if ($newgeneseq eq $oldgeneseq)
  {
    $report{$gid} .= " No change in sequence;";
  }
	my $newcDNA = $newchr->make_cDNA($gene);
  my $newpep = $GD->translate(-sequence => $newcDNA);
  my $oldpep = $GD->translate(-sequence => $cDNA);
  if ($newpep ne $oldpep)
  {
    $report{$gid} .= " Change in amino acid sequence;";
  }
  print "$gid: ", scalar(@chosen), " / $tagtarget\n";
}

print "\n\n";

$factory->cleanup;

foreach my $gid (sort keys %report)
{
  print "$gid : $report{$gid}\n";
}

#Tell chromosome to write itself
$newchr->write_chromosome();

exit;

sub counterset
{
  my ($count, $start) = @_;
  $start = $start + 3 if ($count < $p{MINTAGLEN} + 2);
  $count = $count >= $p{MINTAGLEN} + 2 ? $count - 3 : $p{MAXTAGLEN};
  return ($count, $start);
}

sub sortdiff
{
  my $diffa = join q{}, $a->get_tag_values('difference');
  my $diffb = join q{}, $b->get_tag_values('difference');
  return $diffb <=> $diffa || $a->start <=> $b->start;
}

sub _make_BLAST_factory
{
  my ($dbdata) = @_;

  my $bdir = Bio::BioStudio::ConfigData->config('tmp_path');
  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
    -db_dir  => $bdir,
    -db_data => $dbdata,
    -create  => 1
  );
  return $factory;
}
__END__


=head1 NAME

  BS_PCRTagger.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility creates unique tags for open reading frames to aid the analysis
    of synthetic content in a nascent synthetic genome. Each tag in a gene has
    a wildtype and a synthetic version that correspond to the same offset in the
    gene; each tag can be paired with another to form gene specific amplicons
    which are also specific to either wildtype or synthetic sequence, depending
    on which tags are used.
   
  To pick tags for a chromosome, each open reading frame over I<MINORFLEN> base
   pairs long will be slightly recoded to contain a set of PCR tags. The
   locations and sequences of these tags are carefully chosen to maximize the
   selectivity of the tags for either wild type or synthetic sequence. Each wild
   type or synthetic tag and its reverse complement are unique in the entire
   wild type genome; this is accomplished by creating a BLAST database for the
   entire wild type genome and BLASTing each potential tag against it (this
   requires that a complete wild type genome is available in the BioStudio
   repository). Pairs of tags are selected in such a way that they will not
   amplify any other genomic sequence under 1000 bases long. Each synthetic
   counterpart to a wild type tag is recoded with GeneDesign's “most different”
   algorithm to guarantee maximum nucleotide sequence difference while
   maintaining identical protein sequence and, hopefully, minimizing any effect
   on gene expression. The synthetic tags are all at least I<MINPERDIFF> percent
   recoded from the wild type tags. Each tag is positioned in such a way that
   the first and last nucleotides correspond to the wobble of a codon that can
   be edited to change its wobble without changing its amino acid.  This usually
   automatically excludes methionine or tryptophan, but it can exclude others
   when a I<MINRSCUVAL> filter is in place. The wobble restriction ensures that
   the synthetic and wild type counterparts have different 5’ and 3’
   nucleotides, minimizing the chances that they (and their complements) will
   cross-prime. This means that tags will be between I<MINTAGLEN> and
   I<MAXTAGLEN> base pairs long, where I<TAGLEN> is a multiple of 3 plus 1. All
   tags have melting temperature between I<MINTAGMELT> and I<MAXTAGMELT> so they
   can be used in a single set of PCR conditions.
  
  Tag pairs are chosen to form amplicons specific for each ORF, with at least
   one amplicon chosen per kilobase of ORF. Each amplicon is between
   I<MINAMPLEN> and I<MAXAMPLEN> base pairs long, ensuring that they will all
   fall within an easily identifiable range on an agarose gel. No amplicon will
   be chosen within the first I<FIVEPRIMESTART> base pairs of an ORF to avoid
   disrupting unknown regulatory features. Amplicons are forbidden from
   overlapping each other by more than I<MAXAMPOLAP> percent.

=head1 ARGUMENTS

Required arguments:

  -C, --CHROMOSOME : The chromosome to be modified
  -E, --EDITOR : The person responsible for the edits
  -ME, --MEMO : Justification for the edits
 
Optional arguments:

  -SCA, --SCALE : [genome, chrom (def)] Which version number to increment
  -SCO, --SCOPE : [seg, chrom (def)] How much sequence will the edit affect.
                  seg requires STARTPOS and STOPPOS.
  -STA, --STARTPOS : The first base for analysis;ignored unless SCOPE = seg
  -STO, --STOPPOS  : The last base for analysis;ignored unless SCOPE = seg
  --MINTAGMELT : (default 58) Minimum melting temperature for tags
  --MAXTAGMELT : (default 60) Maximum melting temperature for tags
  --MINPERDIFF : (default 33) Minimum base pair difference between synthetic and
                 wildtype versions of a tag
  --MINTAGLEN  : (default 19) Minimum length for tags. Must be a multiple of 3,
                 plus 1
  --MAXTAGLEN  : (default 28) Maximum length for tags. Must be a multiple of 3,
                 plus 1
  --MINAMPLEN  : (default 200) Minimum span for a pair of tags
  --MAXAMPLEN  : (default 500) Maximum span for a pair of tags
  --MAXAMPOLAP : (default 25) Maximum percentage of overlap allowed between
                 different tag pairs
  --MINORFLEN  : (default 501) Minimum size of gene for tagging eligibility
  --FIVEPRIMESTART : (default 101) The first base in a gene eligible for a tag
  --MINRSCUVAL : (default 0.06) The minimum RSCU value for any replacement codon
                 in a tag
  --OUTPUT    : [html, txt (def)] Format of reporting and output.
  -h, --help : Display this message
 
=cut