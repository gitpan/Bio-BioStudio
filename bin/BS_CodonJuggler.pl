#!/usr/bin/env perl

use Bio::BioStudio;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_CodonJuggler_$VERSION";

local $| = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'   => \$p{CHROMOSOME},
      'EDITOR=s'       => \$p{EDITOR},
      'MEMO=s'         => \$p{MEMO},
      'SCALE=s'        => \$p{SCALE},
      'SCOPE=s'        => \$p{SCOPE},
      'STARTPOS=i'     => \$p{STARTPOS},
      'STOPPOS=i'      => \$p{STOPPOS},
      'FROM=s'         => \$p{FROM},
      'TO=s'           => \$p{TO},
      'DUBWHACK'       => \$p{DUBWHACK},
      'VERWHACK'       => \$p{VERWHACK},
      'ALLWHACK'       => \$p{ALLWHACK},
      'OUTPUT'         => \$p{OUTPUT},
			'help'           => \$p{HELP}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named" unless ($p{CHROMOSOME});
my $chr = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $oldchrseq = $chr->sequence();

$p{OUTPUT} = $p{OUTPUT} || 'txt';

print "BSERROR: Both an editor's id and a memo must be supplied.\n\n"
  unless ($p{EDITOR} && $p{MEMO});

die "BSERROR: Two codons must be supplied." unless ($p{FROM} && $p{TO});
($p{FROM}, $p{TO}) = (uc $p{FROM}, uc $p{TO});
die "BSERROR: Both codons must be three bases long."
  if (length($p{FROM}) != 3 || length($p{TO}) != 3);

my $codon_t = $chr->GD->codontable;
unless (exists $codon_t->{$p{FROM}} && exists $codon_t->{$p{TO}})
{
  die "BSERROR: At least one of the specified codons doesn't exist in " .
    " the codon table. Make sure the codons contain only the characters " .
    "A, T, C, and G.\n";
}

if ($p{FROM} eq $p{TO})
{
  die "BSERROR: The two codons you selected are the same; no editing will " .
    "be done.\n";
}

$p{SCALE} = $p{SCALE} || "chrom";
$p{SCOPE} = $p{SCOPE} || "chrom";

if ($p{SCOPE} && $p{SCOPE} eq "seg" &&
  ((! $p{STARTPOS} || ! $p{STOPPOS}) || $p{STOPPOS} <= $p{STARTPOS}))
{
  die "BSERROR: The start and stop coordinates do not parse\n";
}

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $newchr = $chr->iterate();
my $GD = $newchr->GD;
my $chrseq = $newchr->sequence();
my $chrlen = length $chrseq;

my $changes  = $chr->allowable_codon_changes($p{FROM}, $p{TO});
$p{SWAPTYPE} = $GD->codon_change_type(-from => $p{FROM}, -to => $p{TO});
$p{MORF} = $GD->complement(-sequence => $p{FROM}, -reverse => 1);
$p{OT} = $GD->complement(-sequence => $p{TO}, -reverse => 1);

my @pgenes = $newchr->db->features(
  -seq_id     => $chr->seq_id,
  -range_type => 'contains',
  -types      => 'gene',
  -start      => $p{STARTPOS},
  -end        => $p{STOPPOS}
);

my @genes;
foreach my $pgene (@pgenes)
{
  if ($pgene->has_tag('orf_classification'))
  {
    next if $pgene->Tag_orf_classification eq 'RNA';
  }
  my @subfeats = $newchr->flatten_subfeats($pgene);
  my @CDSes = grep {$_->primary_tag eq 'CDS'} @subfeats;
  next if (! scalar @CDSes);
  push @genes, $pgene;
}

my %index_to_type;
my %state;
my $mask = Bio::BioStudio::Mask->new(-sequence => $newchr->seqobj);
my %gmasks = ();
foreach my $gene (@genes)
{
  my $gid = $gene->primary_id;
  $state{$gene->display_name} = [0, 0, q{}, $gid];
  my $offset = $gene->start() - 1;
  my $gmask = Bio::BioStudio::Mask->new(-sequence => $gene, -offset => $offset);
  my @subfeats = $newchr->flatten_subfeats($gene);
  $gmask->add_to_mask(\@subfeats);
  $mask->add_to_mask(\@subfeats);
  $gmasks{$gid} = $gmask;
  $index_to_type{$gid} = $gene->primary_tag;
  $index_to_type{$_->primary_id} = $_->primary_tag foreach (@subfeats);
}

foreach my $gene (@genes)
{
  my $gname = $gene->display_name;
  my $orient = $gene->strand();
  my $gstart = $gene->start();
  my $gend = $gene->end();
  my $glen = $gend - $gstart + 1;
  my $gid = $gene->primary_id;
  my $gmask = $gmasks{$gid};
  my %subs = map {$_->primary_id => 1} $newchr->flatten_subfeats($gene);
  $subs{$gid}++;
  
  my $x = $gstart;
  while ($x <= $gend - 2)
  {
    ##Adjust for introns
    my %posa = $gmask->what_overlaps($x);
    my @introns = grep {$index_to_type{$_} eq 'intron'} values %posa;
    while (scalar @introns)
    {
      $x++;
      %posa = $gmask->what_overlaps($x);
      @introns = grep {$index_to_type{$_} eq 'intron'} values %posa;
    }
    my $basea = $x;
    $x++;
    
    my %posb = $gmask->what_overlaps($x);
    @introns = grep {$index_to_type{$_} eq 'intron'} values %posb;
    while (scalar @introns)
    {
      $x++;
      %posb = $gmask->what_overlaps($x);
      @introns = grep {$index_to_type{$_} eq 'intron'} values %posb;
    }
    my $baseb = $x;
    $x++;
    
    my %posc = $gmask->what_overlaps($x);
    @introns = grep {$index_to_type{$_} eq 'intron'} values %posc;
    while (scalar @introns)
    {
      $x++;
      %posc = $gmask->what_overlaps($x);
      @introns = grep {$index_to_type{$_} eq 'intron'} values %posc;
    }
    my $basec = $x;
    $x++;
    
    ##Codon sequence
    my $cod = substr($chrseq, $basea - 1, 1)
            . substr($chrseq, $baseb - 1, 1)
            . substr($chrseq, $basec - 1, 1);
    next if ($orient < 0 && $cod ne $p{MORF});
    next if ($orient >= 0 && $cod ne $p{FROM});
    my $maskbit = $mask->count_features_in_range($basea, 3);
    
    my $newcod = $orient > 0  ? $p{TO}  : $p{OT};
    my $newaa = $codon_t->{$p{TO}};
    
    my $oldcod = $orient > 0 ? $p{FROM} : $p{MORF};
    my $oldaa = $codon_t->{$p{FROM}};
    
    my $oset = $basea - $gstart + 1;
    my $codname = $gname . "_" . $oldaa . $oset . $newaa;

    my $codon = Bio::SeqFeature::Generic->new(
      -start        => $basea,
      -end          => $basec,
      -primary_tag  => $p{SWAPTYPE},
      -display_name => $codname,
      -tag          => {
        wtseq       => $oldcod,
        newseq      => $newcod,
        parent_id   => $gname
      }
    );
    
    ## No overlap involved
    #
    if ($maskbit == 1)
    {
      $state{$gname}->[0]++;
      my $comment = "$p{SWAPTYPE} " . $codname . " added\n";
      $newchr->add_feature(-feature => $codon, -comments => [$comment]);
      next;
    }
    
    ## Overlap involved:
    #
    my @featlaps = $mask->features_in_range($basea, 3);
    my %laps = map {$_ => 1} grep {! exists $subs{$_}} @featlaps;
    my @theseolaps = keys %laps;
    my $lapcount = scalar @theseolaps;

    # If there are more than one overlapping features, no change can be made
    # If there are zero, we have run into an internal error
    if ($lapcount > 2)
    {
      $state{$gname}->[2] .= "cannot change codon at $basea (too occluded) ";
      next;
    }
    elsif ($lapcount < 1)
    {
      $state{$gname}->[2] .= "cannot change codon at $basea (INTERNAL ERROR) ";
      next;
    }
    
    # If the only overlapping feature is an intron, no change can be made
    if ($index_to_type{$featlaps[0]} eq 'intron')
    {
      $state{$gname}->[2] .= "cannot change codon at $basea (overlaps an intron)";
      next;
    }
    
    # Now the overlap must be a CDS; determine the frame
    my $lapfeat = $chr->db->fetch($theseolaps[0]);
    my $lapname = $lapfeat->display_name;
    my $lapstart = $lapfeat->start;
    my $laporient = $lapfeat->strand;
    my $laplen = $lapfeat->end - $lapstart + 1;
    my $CDSseq = substr($chrseq, $lapstart - 1, $laplen);
    my $frame = ($basea - $lapstart) % 3;
    my $length = $basec - ($basea - $frame) + 1;
    $length++ while ($length % 3 != 0);
    my $syncodstart = $basea - $lapstart - $frame;
    my $inframe = substr($CDSseq, $syncodstart, $length);
    my $lappep = $laporient == 1
      ?  $GD->translate($inframe)
      :  $GD->translate($GD->complement($inframe, 1));
    my $orientswit = $lapfeat->strand eq $orient  ?  0  :  1;
    my %allowed = %{$changes->{$orientswit}};
    my $allowflag = exists $allowed{$lappep}  ? 1 : 0;
    
    # If the change doesn't preserve overlapping sequence and no exceptions were
    # provided, skip
    my $lapstatus = $lapfeat->Tag_orf_classification;
    my $genstatus = $gene->Tag_orf_classification;
    my $dubwhack = $lapstatus eq 'Dubious' && $genstatus ne 'Dubious' ? 1 : 0;
    my $verwhack = $genstatus ne 'Dubious'  ? 1 : 0;
    
    if ($allowflag == 0 && (! $p{ALLWHACK}
      || ! ($dubwhack && $p{DUBWHACK}) || ! ($verwhack && $p{VERWHACK})))
    {
      $state{$gname}->[2] .= "cannot change codon at $basea (overlaps $lapname)";
      next;
    }

    # Add the codon for this gene
    $state{$gname}->[0]++;
    my $comment = "$p{SWAPTYPE} " . $codname . " added\n";
    $newchr->add_feature(-feature => $codon, -comments => [$comment]);
    $chrseq = $newchr->sequence();
    
    # Add the codon for the overlapping gene
    
    my $newCDSseq = substr($chrseq, $lapstart - 1, $laplen);
    my $newframe = substr($newCDSseq, $syncodstart, $length);
    
    for (my $offset = 0; $offset <= $length - 3; $offset += 3)
    {
      my $oldsyncod = substr($inframe, $offset, 3);
      my $newsyncod = substr($newframe, $offset, 3);
      my $pos = $syncodstart + $offset + $lapstart + 1;
      if ($oldsyncod ne $newsyncod)
      {
        my $newLcod = $newsyncod;
        if ($laporient == -1)
        {
          $newLcod = $GD->complement(-sequence => $newLcod, -reverse => 1);
        }
        my $newLaa = $codon_t->{$newLcod};
    
        my $oldLcod = $oldsyncod;
        if ($laporient == -1)
        {
          $oldLcod = $GD->complement(-sequence => $oldLcod, -reverse => 1);
        }
        my $oldLaa = $codon_t->{$oldLcod};

        my $swaptype = $GD->codon_change_type(-from => $oldLcod, -to =>$newLcod);
        $newLcod = $newsyncod if ($laporient == -1);
        $oldLcod = $oldsyncod if ($laporient == -1);
        
        my $Loffset = $pos - $lapstart + 1;
        my $parent = $lapfeat->has_tag('parent_id') ? $lapfeat->Tag_parent_id : $lapname;
        my $codLname = $parent . "_" . $oldLaa . $Loffset . $newLaa;
        my $modnote = " $codLname Modified to accommodate $p{FROM} to $p{TO} change in $gname";
        my $comment2 = "$swaptype $codLname added\n";

        my $Lcodon = Bio::SeqFeature::Generic->new(
          -start        => $pos-1,
          -end          => $pos+1,
          -primary_tag  => $swaptype,
          -display_name => $codLname,
          -tag          => {
            wtseq       => $oldLcod,
            newseq      => $newLcod,
            parent_id   => $parent,
            Note        => $modnote
          }
        );
        $newchr->add_feature(-feature => $Lcodon, -comments => [$comment2]);
        $chrseq = $newchr->sequence();
        $state{$parent}->[1]++;
        $state{$parent}->[2] .= $modnote;
      }
    }
  }
}

#Do error checking
$chrseq = $newchr->sequence();
foreach my $gene (@genes)
{
  my $gname = $gene->display_name;
  my $gstart = $gene->start();
  my $gend = $gene->end();
  my $glen = $gend - $gstart + 1;
  
  my $newseq = substr($chrseq, $gstart - 1, $glen);
  my $oldseq = substr($oldchrseq, $gstart - 1, $glen);
  if ($newseq eq $oldseq && ($state{$gname}->[0] || $state{$gname}->[1]))
  {
    $state{$gname}->[2] .= " No change in sequence;";
  }
  my $cDNA = $chr->make_cDNA($gene);
	my $newcDNA = $newchr->make_cDNA($gene);
  my $newpep = $GD->translate(-sequence => $newcDNA);
  my $oldpep = $GD->translate(-sequence => $cDNA);
  if ($newpep ne $oldpep)
  {
    $state{$gname}->[2] .= " Change in amino acid sequence;";
  }
}

#Do reporting
print "\nReport:\n";
foreach my $gname (sort keys %state)
{
  my $gid = $state{$gname}->[3];
  my @results = @{$state{$gname}};
  next unless($results[0] || $results[1] || $results[2]);
  print $gname, " : ";
  if ($results[0])
  {
    my $plural = $results[0] > 1  ? 's' : q{};
    print "$results[0] $p{FROM} codon$plural changed; ";
  }
  if ($results[1])
  {
    my $plural = $results[1] > 1  ? 's' : q{};
    print "$results[1] codon$plural changed; ";
  }
  print "$results[2]" if ($results[2]);
  print "\n";
}

$newchr->write_chromosome();

exit;

__END__

=head1 NAME

  BS_CodonJuggler.pl
 
=head1 VERSION

    Version 2.00

=head1 DESCRIPTION

  This utility switches any one codon to any other. By default, it will not make
    any change to a gene that will cause a translation change in an overlapping
    gene; this behavior can be overridden with the -D and -V flags, which will
    allow the utility to make nonsynonymous changes to ORFs marked dubious and
    verified, respectively.
   
	If a stop codon is changed to a different stop codon, the change will be
    marked "stop_retained_variant". Otherwise synonymous changes are marked
    "synonymous_codon". If a stop is changed to a non stop, it is a "stop_lost".
    If a non stop is changed to a stop it is "stop_gained"; any other change is
    a "non_synonymous_codon".

=head1 USAGE

Required arguments:

  -C,  --CHROMOSOME : The chromosome to be modified
  -E,  --EDITOR : The person responsible for the edits
  -M,  --MEMO   : Justification for the edits
  -F,  --FROM   : The codon to be replaced
  -T,  --TO     : The codon to be introduced

Optional arguments:

  -SCA, --SCALE : [genome, chrom (def)] Which version number to increment
  -SCO, --SCOPE : [seg, chrom (def)] How much sequence the edit will affect.
                  seg requires STARTPOS and STOPPOS.
  -STA, --STARTPOS : The first base for editing; ignored unless SCOPE = seg
  -STO, --STOPPOS  : The last base for editing; ignored unless SCOPE = seg
  -D,   --DUBWHACK : Allow nonsynonymous changes to dubious ORFs on behalf of
                     non dubious ORFs
  -V,   --VERWHACK : Allow nonsynonymous changes to verified ORFs on behalf of
                     non dubious ORFs
  -A,   --ALLWHACK : Allow even nonsynonymous changes to all ORFs
  -O,   --OUTPUT   : [html, txt (def)] Format of reporting and output.
  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013, BioStudio developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of Johns Hopkins, the Joint Genome Institute, the Lawrence Berkeley
National Laboratory, the Department of Energy, and the BioStudio developers may
not be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut