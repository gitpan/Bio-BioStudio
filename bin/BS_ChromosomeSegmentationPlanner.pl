#!/usr/bin/perl

use Bio::BioStudio;
use Bio::BioStudio::Mask;
use Bio::BioStudio::RestrictionEnzyme;
use Bio::BioStudio::RestrictionEnzyme::Store;
use Bio::BioStudio::Megachunk;
use Bio::BioStudio::Chunk;
use Getopt::Long;
use Pod::Usage;
use English qw(-no_match_vars);
use Carp;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_ChromosomeSegmentationPlanner_$VERSION";

local $OUTPUT_AUTOFLUSH = 1;

my $debug = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'    => \$p{CHROMOSOME},
      'WTCHR=s'         => \$p{WTCHR},
      'RESET=s'         => \$p{RESET},
      'MARKERS=s'       => \$p{MARKERS},
      'LASTMARKER=s'    => \$p{LASTMARKER},
      'SCOPE=s'         => \$p{SCOPE},
      'STARTPOS=i'      => \$p{STARTPOS},
      'STOPPOS=i'       => \$p{STOPPOS},
      'ISSMIN=i'        => \$p{ISSMIN},
      'ISSMAX=i'        => \$p{ISSMAX},
      'CHUNKLENMIN=i'   => \$p{CHUNKLENMIN},
      'CHUNKLENMAX=i'   => \$p{CHUNKLENMAX},
      'CHUNKNUM=i'      => \$p{CHUNKNUM},
      'CHUNKNUMMIN=i'   => \$p{CHUNKNUMMIN},
      'CHUNKNUMMAX=i'   => \$p{CHUNKNUMMAX},
      'CHUNKOLAP=i'     => \$p{CHUNKOLAP},
      'FPUTRPADDING=i'  => \$p{FPUTRPADDING},
      'TPUTRPADDING=i'  => \$p{TPUTRPADDING},
      'FIRSTLETTER=s'   => \$p{FIRSTLETTER},
      'REDOREDB'        => \$p{REDOREDB},
      'help'            => \$p{HELP}
);
if ($p{HELP})
{
  pod2usage(
    -verbose => 99,
    -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE");
}

################################################################################
############################### SANITY CHECKING ################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n" unless ($p{CHROMOSOME});
my $oldchr = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $GD = $oldchr->GD;
my $chrlen = $oldchr->len();

$p{RESET}        = $p{RESET}        || "nonpal_and_IIB";
my $RES = $GD->set_restriction_enzymes(-enzyme_set => $p{RESET});
my @enzes = values %{$RES};

die "\n ERROR: No target chromosome was named.\n\n" unless ($p{WTCHR});
my $wtchr = $BS->set_chromosome(-chromosome => $p{WTCHR});

$p{CHUNKLENMIN}  = $p{CHUNKLENMIN}  || 6000;
$p{CHUNKLENMAX}  = $p{CHUNKLENMAX}  || 9920;
die "\n ERROR: The chunk length parameters do not parse.\n"    
  if ($p{CHUNKLENMAX} < $p{CHUNKLENMIN});

$p{ISSMIN}       = $p{ISSMIN}       || 900;
$p{ISSMAX}       = $p{ISSMAX}       || 1500;
die "\n ERROR: The ISS length parameters do not parse.\n" 
  if ($p{ISSMAX} < $p{ISSMIN});

$p{SCOPE}        = $p{SCOPE}        || 'chrom';
$p{STARTPOS}     = $p{STARTPOS}     || 1;
$p{STOPPOS}      = $p{STOPPOS}      || $chrlen;
die "\n ERROR: The start and stop coordinates do not parse.\n"         
  if ($p{SCOPE} eq "seg" && ($p{STOPPOS} <= $p{STARTPOS}));

$p{CHUNKNUM}     = $p{CHUNKNUM}     || 4;
$p{CHUNKNUMMIN}  = $p{CHUNKNUMMIN}  || 3;
$p{CHUNKNUMMAX}  = $p{CHUNKNUMMAX}  || 5;
$p{CHUNKOLAP}    = $p{CHUNKOLAP}    || 40;
$p{FPUTRPADDING} = $p{FPUTRPADDING} || 500;
$p{TPUTRPADDING} = $p{TPUTRPADDING} || 100;
$p{FIRSTLETTER}  = $p{FIRSTLETTER}  || 'A';
$p{REDOREDB}     = $p{REDOREDB}     || 0;
$p{REDOREDB}     = $p{REDOREDB}     || 0;


die "\n BSERROR: No markers were named.\n\n" unless ($p{MARKERS});
my %markers = ();
my $BS_MARKERS = {};
$BS_MARKERS = $BS->fetch_custom_markers();
%markers = map {$_ => -1} split(",", $p{MARKERS});
foreach my $marker (keys %markers)
{
  die "\n BSERROR: $marker not recognized as a marker.\n"
    unless (exists $BS_MARKERS->{$marker});
}
$p{MARKERHSH} = \%markers;
die "\n ERROR: last marker $p{LASTMARKER} not recognized as a marker.\n"
  if ($p{LASTMARKER} && ! exists $BS_MARKERS->{$p{LASTMARKER}});
  
## Order the markers and determine which is the smallest
my ($x, $size, $smallest) = (0, undef, undef);
$p{MARKERORDER} = {};
if ($p{LASTMARKER} && exists ($p{MARKERHSH}->{$p{LASTMARKER}}))
{
  $p{MARKERHSH}->{$p{LASTMARKER}} = $x;
  $p{MARKERORDER}->{$x} = $p{LASTMARKER};
  $x++;
}
foreach my $marker (keys %{$p{MARKERHSH}})
{
  next if ($p{MARKERHSH}->{$marker} != -1);
  $p{MARKERORDER}->{$x} = $marker ;
  my $markersize = length($BS_MARKERS->{$marker}->sequence);
  $smallest = $marker if (! $size || $markersize < $size);
  $x++;
}
$p{SMALLESTMARKER} = $smallest;
$p{SMALLESTMARKERLEN} = length($BS_MARKERS->{$p{SMALLESTMARKER}}->sequence);
$p{MARKERCOUNT} = scalar(keys %{$p{MARKERHSH}});
 
################################################################################
################################# CONFIGURING ##################################
################################################################################
$p{BSVERSION} = $bsversion;
print "CONFIGURING...\n";
my $DBLIST = $BS->db_list();

#Create a global restriction enzyme database, if one does not exist
my $REDBNAME = $oldchr->name . '_RED';
if ($p{REDOREDB} || ! exists $DBLIST->{$REDBNAME})
{
  print "CREATING GLOBAL MARKUP DATABASE...\n";
  local $SIG{CHLD} = 'DEFAULT';
  system(
    'BS_GlobalREMarkup.pl',
    '--CHROMOSOME', $oldchr->name(),
    '--RESET ',  $p{RESET},
    '--CHUNKLENMIN ',  $p{CHUNKLENMIN},
    '--CHUNKLENMAX ',  $p{CHUNKLENMAX},
    '--CULLINC 2000',
    '--REDOREDB 1'
  ) == 0 || die "Cannot execute BS_GlobalREMarkup.pl: $OS_ERROR $CHILD_ERROR\n";
}

my $REDB = Bio::BioStudio::RestrictionEnzyme::Store->new(
  -name => $REDBNAME,
  -enzyme_definitions => $RES
);

my $planchr = $oldchr->iterate(-chrver => 1, -tag => '_PLAN');
my $pchrseq = $planchr->sequence();
my $pchrlen = length($pchrseq);

## Initialize parameters and constants
my $MAX_MEGA_CHUNK = $p{CHUNKNUM} * $p{CHUNKLENMAX};
my $MIN_MEGA_CHUNK = $p{CHUNKNUM} * $p{CHUNKLENMAX};

## Gather genes and check for essentiality.  Mask the UTRs of essential
## and fast growth genes by modifying their start and stop sites.
## Find the right telomere and make the junction just before it
my @UTCS = $planchr->db->features(-seq_id => $planchr->seq_id, -type => 'UTC');
@UTCS = sort {$b->start <=> $a->start} @UTCS;
my $rutc = $UTCS[0];
my $chrend = scalar(@UTCS)  ? $rutc->start -1 : $pchrlen;

my @genes = $planchr->db->get_features_by_type('gene');
my $mask = $planchr->type_mask('gene');

my @efgenes = ();
foreach my $gene (@genes)
{
  my $status = $gene->has_tag('essential_status')
    ? $gene->Tag_essential_status
    : 'Nonessential';
  if ($status ne 'Nonessential')
  {
    my $newstart = $gene->strand == 1
        ? $gene->start - $p{FPUTRPADDING}
        : $gene->start - $p{TPUTRPADDING};
    $newstart = 1 if ($newstart < 1);
    $gene->start($newstart);
    my $newstop = $gene->strand == 1
        ? $gene->stop + $p{TPUTRPADDING}
        : $gene->stop + $p{FPUTRPADDING};
    $newstop = $chrend if ($newstop > $chrend);
    $gene->stop($newstop);
    push @efgenes, $gene;
  }
}
my @CDSes = $planchr->db->get_features_by_type('CDS');
my %ESSENTIALS;
foreach my $CDS (@CDSes)
{
  my $status = $CDS->Tag_essential_status || 'Nonessential';
  $ESSENTIALS{$CDS->display_name} = $status;
}

# Mask the regions that are not viable ISS locations:  Essential / Fast-Growth
#   genes, regions smaller than the smallest possible ISS site,
#   and regions with no non-essential genes in them
my $altmask = Bio::BioStudio::Mask->new(-sequence => $planchr->seqobj);
$altmask->add_to_mask(\@efgenes);
my @intergenics = @{$altmask->find_deserts()};
foreach my $range (@intergenics)
{
  my ($start, $end) = ($range->start, $range->end);
  if ($end - $start < $p{ISSMIN})
  {
    my $newfeat = Bio::SeqFeature::Generic->new(
      -start        => $start,
      -end          => $end,
      -display_name => "no_ISS_$start",
      -primary_tag  => "forbid"
    );
    $altmask->add_to_mask([$newfeat]);
  }
}
my @ess = sort {abs($b->start - $b->end) <=> abs($a->start - $a->end)}
          @{$altmask->find_deserts()};

#fbound and tbound are the 5' and 3' bounds of segmentation, where we start and
# stop. default is the first and last base of the chromosome, respectively.
# lastpos is where we actually start counting from.
my $tbound = $p{STOPPOS}  ? $p{STOPPOS}  : $chrend;
$tbound = $chrend if ($tbound > $chrend);
my $lastpos = $tbound;
my $fbound = $p{STARTPOS} + ($p{CHUNKLENMAX} - (0.5*$p{ISSMIN}));
$lastpos = $p{LASTMARKER}
  ? $tbound + length($BS_MARKERS->{$p{LASTMARKER}}->sequence)
  : $tbound + length($BS_MARKERS->{$p{MARKERORDER}->{0}}->sequence);
$fbound = $fbound - $p{SMALLESTMARKERLEN};

my $markercount = 0;
my $firstmarker = 0;
my $foundcand = undef;
my $doISS = 0;
my $mchrollback = 0;
my $redomch = undef;
my $excisor = undef;
my $lastISSseq = undef;
my %usedsites;
my %deadends;

################################################################################
################################# SEGMENTING  ##################################
################################################################################
print "SEGMENTING...\n";
my @mchs;

while (! $lastpos || ($lastpos-$fbound) >= $p{CHUNKLENMIN} * 2)
{
  my $mch = $redomch  ? $redomch  : Bio::BioStudio::Megachunk->new();
  $lastpos = $mch->end  ? $mch->end : $lastpos;
 
  $mch->excisor($excisor) unless ($mch->excisor);
 
  #decide where the start and stop of this megachunk will be.
  # frange and trange are the 5' and 3' boundaries of the megachunk.
  my $frange = $lastpos - $MAX_MEGA_CHUNK >= $fbound 
    ? $lastpos - $MAX_MEGA_CHUNK
    : $fbound;
  $mch->frange($frange);
  $mch->trange($lastpos - $MIN_MEGA_CHUNK);
 
  #eligible regions have no masking - that is, no essential or fast growth genes
  # sort by distance from the last position and choose a start inside frange
  unless ($mch->start)
  {
    my @eligibles = grep { $_ >= $mch->frange && $_ <= $mch->trange}
                    sort {abs($b-$lastpos) <=> abs($a-$lastpos)} @ess;
    my $mchstart  = $eligibles[0] ? $eligibles[0] : $mch->frange;
    $mchstart = $fbound if ($lastpos < $MIN_MEGA_CHUNK);
    $mch->start($mchstart);
    $mch->end($lastpos);
  }
 
  #pick a marker for the megachunk and decide which marker will be next.
  # If lastmarker is specified but not in the marker group, take care
  $markercount = $mch->markercount  ? $mch->markercount : $markercount;
  unless ($mch->marker)
  {
    $mch->markercount($markercount);
    $mch->firstmarker($firstmarker);
    my ($markername, $omarkername);
    if ($p{LASTMARKER} && $firstmarker == 0
      && (! exists $p{MARKERHSH}->{$p{LASTMARKER}}))
    {
      $markername = $p{LASTMARKER};
      $omarkername = $p{MARKERORDER}->{$markercount};
      ($firstmarker, $markercount) = (1, 0);
    }
    else
    {
      $markername = $p{MARKERORDER}->{$markercount % $p{MARKERCOUNT}};
      $omarkername = $p{MARKERORDER}->{($markercount+1) % $p{MARKERCOUNT}};
    }
    $mch->marker($BS_MARKERS->{$markername});
    $mch->omarker($BS_MARKERS->{$omarkername});
  }
  
  #If there was a previous megachunk, its 5' enzyme becomes the 3' enzyme of
  # this megachunk. Mark it, its exclusions, and its creations, as used
  my @borders = ();
  if (! $mch->prevenz && $foundcand)
  {
    %usedsites = ();
    $mch->prevenz($foundcand);
    push @borders, $foundcand;
    %usedsites = ($foundcand->id => 1);
    $usedsites{$_}++ foreach (@{$RES->{$foundcand->id}->exclude});
    if ($foundcand->creates)
    {
      $usedsites{$_}++ foreach (@{$foundcand->creates});
    }
  }
 
my $mchlen = $mch->end - $mch->start + 1;
print "Making a megachunk " . $mch->start . ".." . $mch->end . " ~$mchlen bp ";
print "with $markercount " . $mch->marker->name . "\n";
 
  $mch->chunks([]) unless($mch->chunks);
  my %usedhangs = ();
  my $chunknum = $mch->chunknum ? $mch->chunknum : 1;
  my $firsterr = 0;
  my $redoch = $mchrollback ? pop @{$mch->chunks}  : undef;
 
  #Start to make chunks.
  while (abs($lastpos - $mch->start) > $p{CHUNKLENMIN} || $chunknum <= $p{CHUNKNUMMAX})
  {
    my $ch = $redoch  ? $redoch : Bio::BioStudio::Chunk->new();
    $ch->number($chunknum);
   
    $ch->prevcand($foundcand) unless ($ch->prevcand);
    $lastpos = $ch->prevcand ? $ch->prevcand->end  : $lastpos;
    $mch->lastlastpos($lastpos);
   
    #Keep track of which overhangs and which sites are off limits
    my %tempsites = %usedsites;
    my %temphangs = %usedhangs;
    $ch->used_enzymes(\%tempsites) unless ($ch->used_enzymes);
    $ch->used_overhangs(\%temphangs) unless ($ch->used_overhangs);
   
    #
    my $issflag = $doISS;
    $issflag++ if ($chunknum == $p{CHUNKNUM} && $firsterr == 0)
                || $chunknum == $p{CHUNKNUMMAX};
    $issflag-- if ($mch->start <= $fbound);

    #Get the list of viable candidates for enzyme border
    # If this is the first (3' most) chunk, allow for the ISS sequence
    my $isslen = $chunknum == 1
          ? $p{CHUNKLENMAX} - (0.5*$p{ISSMIN})
          : 0;
    $isslen -= length($mch->marker->sequence) if ($chunknum == 1);

    #Create the candidate list
    my @candlist;
    unless ($ch->enzlist)
    {
      my $candlistref = ProposeCandidates($lastpos, \%p, $isslen);
     
      #Filter candidates by position and sort by score, with lower scores first.
      # If this is the 3' most chunk of the megachunk, don't pick REs that can't
      # be removed from the marker. If an ISS sequence is 3' of the chunk, don't
      # pick REs that are present in the wildtype sequence of the ISS sequence.
      @candlist = @{$candlistref};
      @candlist = grep {abs($_->start - $lastpos) >= $p{CHUNKLENMIN}} @candlist;
      @candlist = grep {abs($_->end   - $lastpos) <= $p{CHUNKLENMAX}} @candlist;
      @candlist = grep {$_->start > $mch->start} @candlist;
      @candlist = sort {$a->score <=> $b->score} @candlist;
      if ($chunknum == 1)
      {
        @candlist = grep {! exists($mch->marker->static_enzymes->{$_->id})} @candlist;
        if ($lastISSseq)
        {
          my $ihsh = $GD->restriction_status(-sequence => $lastISSseq);
          @candlist = grep {$ihsh->{$_->id} == 0} @candlist;
        }
      }
      my @cutlist = @candlist;
      $ch->bkupenzlist(\@cutlist);
    }
    else
    {
      @candlist = @{$ch->enzlist};
    }
 
    #Take a survey of how often overhang sequences appear among candidates
    # We will want to use overhangs that appear the least first
    my %ohangsurvey;
    foreach my $cand (@candlist)
    {
      my $enztype = $cand->type;
      $ohangsurvey{"$_.$enztype"}++ foreach (keys %{$cand->overhangs});
    }

my @nonos = keys %{$ch->used_enzymes};
my @nohgs = keys %{$ch->used_overhangs};
print "\t ch $chunknum \t", scalar(@candlist), " possibilities\tflag $issflag";
    print "\t firsterr $firsterr \t nonos: @nonos / @nohgs";
print "\t(" . $ch->prevcand->id . " @ $lastpos)" if ($ch->prevcand);
print "\n";
       
    #Pick a candidate
    %usedhangs = %{$ch->used_overhangs};
    %usedsites = %{$ch->used_enzymes};
    $foundcand = NextCandidate(\@candlist, \%usedhangs, \%ohangsurvey,
      \%usedsites, $issflag, $ch->prevcand, $excisor, $mch->omarker, \%p);
   
    #If no candidate can be picked, we will either have to move forward or back
    if (! $foundcand)
    {
      #if we've at least made the minimum number of chunks and this is the first
      # error, we will try relaxing the standards and moving up.
      if ($chunknum >= $p{CHUNKNUMMIN} && $firsterr == 0)
      {
        $firsterr++;
        $ch->enzlist($ch->bkupenzlist);
        $redoch = $ch;
        next;
      }
      #Otherwise, we need to reduce the number of possible chunks and rechoose.
      # we will restore a previous candidate and its workspace
      elsif ($chunknum > 1)
      {
        $redoch = pop @{$mch->chunks};
        $chunknum--;
        next;
      }
      #If we've backed all the way up to chunk 1, we need to redo the previous
      # megachunk entirely.
      elsif ($chunknum == 1)
      {
         last;
      }
      else
      {
        die "sorry, boss - ran into a hole!\n\n";
      }
    }
    $ch->enzlist(\@candlist);
    $ch->enzyme($foundcand);
    if (exists $deadends{$foundcand->name}
            && $deadends{$foundcand->name} == $chunknum)
    {
      print "Been down this road... ";
      last;
    }
    if ($chunknum == 1 && ! $redoch)
    {
      $mch->firstchunk($foundcand) unless($mch->firstchunk);
    }
   
##DEBUG UPDATE
print "\t  picked ", $foundcand->name, q{ };
print $foundcand->phang, q{ }, $foundcand->score;
my @mustmovers = ();
@mustmovers = @{$foundcand->movers} if $foundcand->movers;
print "(must move @mustmovers)" if $foundcand->movers;
print "\n";
   
    #Finish picking.  Reset used sites, push the candidate onto the stack.
    push @{$mch->chunks}, $ch;
    $redoch = undef;
    %usedsites = ();
    $usedsites{$_}++ foreach (@{$foundcand->exclude});
    if ($foundcand->creates)
    {
      $usedsites{$_}++ foreach (@{$foundcand->creates});
    }
    $lastpos = $foundcand->end;
    last if $chunknum == $p{CHUNKNUM} && $firsterr == 0;
    last if $chunknum == $p{CHUNKNUMMAX};
    last if ($lastpos < $p{CHUNKLENMIN});
    $chunknum++;
  #End chunk picking
  }
 
  #If fewer chunks are picked than allowed, we have to redo the whole megachunk.
  if ($chunknum < $p{CHUNKNUMMIN} && $lastpos > $p{CHUNKLENMIN})
  {
    $mchrollback = 1;
    $deadends{$mch->firstchunk->name} = 1 if ($mch->firstchunk);
    $redomch = pop @mchs;
    $redomch->pexcisor(undef);
    next;
  }
  if (! $mch->pexcisor && $lastpos > $p{CHUNKLENMIN})
  {
    my $rtest = $altmask->count_features_in_range($lastpos, $p{ISSMIN});
    die ("wtf\n\n$rtest\n\n")  if ($rtest > 0);
    my $plist = ProposeExcisors($lastpos, $mch->omarker, \%p);
    my @partnerlist = sort {$a->score <=> $b->score} @{$plist};
    $excisor = $partnerlist[0];
    die ("can't find IIB site!!!\n") unless $excisor;
my $rISSbound = $lastpos + $p{ISSMIN};
print "\t  next excisor should be ", $excisor->id, q{ };
print $excisor->start, q{ }, $excisor->score, " @ $rISSbound\n";
    $lastISSseq = wtISS($lastpos, $excisor->start, $wtchr);
    $lastpos = $excisor->start;
    $mch->pexcisor($excisor);
  }
  push @mchs, $mch;
  $redomch = undef;
  $mchrollback = 0;
  $markercount++;
  print "\n";
#End megachunk picking
}

################################################################################
########################### COMMITTING TO DATABASE  ############################
################################################################################
print "SUMMARIZING...\n";
my $letter = $p{FIRSTLETTER};
@mchs = reverse @mchs;
foreach my $mch (@mchs)
{
  my $megachunkname = $p{NEWNAME} . ".$letter";
  my @chlist = sort {$a->start <=> $b->start} map {$_->enzyme} @{$mch->chunks};
  print "\t", $_->id, q{ }, $_->presence, q{ }, $_->start, "\n" foreach (@chlist);
  my $excisor = $mch->excisor ? $mch->excisor : undef;
  my $mcstart = $letter eq $p{FIRSTLETTER}  ? $fbound : $chlist[0]->start;
  my $mcend = $excisor ? $excisor->start : $tbound;
  my $d = $planchr->db->new_feature(
      -index        => 1,
      -seq_id       => $planchr->seq_id,
      -start        => $mcstart,
      -end          => $mcend,
      -primary_tag  => "megachunk",
      -display_name => $megachunkname,
      -source       => "BIO",
      -attributes   => {
        "load_id"   => $megachunkname,
        "intro"     => $p{INTRO},
        "marker"    => $mch->marker->name,
        "bsversion" => $p{BSVERSION}}
  );
  if ($excisor)
  {
    $d->add_tag_value("excisor", $excisor->id);
    $d->update();
  }
  print $d->Tag_load_id, " (", $d->start, "..", $d->end, ") ", $d->Tag_marker;
  print q{ }, $d->Tag_excisor, q{ }, $excisor->start if ($excisor);
  print "\n";
 
  my $chunknum = 1;
  my $lastenz = undef;
  my $chlim = scalar(@chlist);
  $chlim++ if ($letter eq $p{FIRSTLETTER});
  while ($chunknum <= $chlim)
  {
    my ($chnum, $nchnum) = ($chunknum, $chunknum+1);
    my $chunkname = $megachunkname . $chnum;
    my $nchunkname = $megachunkname . $nchnum;
    my $left = $chunknum == 1 && $letter ne $p{FIRSTLETTER}  ? shift @chlist : undef;
    my $right = shift @chlist;
    if ($left)
    {
      my $l = $planchr->db->new_feature(
        -index        => 1,
        -seq_id       => $planchr->seq_id,
        -start        => $left->start,
        -end          => $left->end,
        -score        => sprintf("%.3f", $left->score),
        -display_name => $left->name,
        -primary_tag  => "enzyme_recognition_site",
        -source       => "BIO",
        -strand       => $left->strand,
        -attributes   => {
          "presence"    => $left->presence,
          "infeat"      => $left->featureid,
          "enzyme"      => $left->id,
          "load_id"     => $left->name,
          "wanthang"    => $left->phang,
          "peptide"     => $left->peptide,
          "ohangoffset" => $left->offset,
          "megachunk"   => $megachunkname,
          "chunks"      => $chunkname}
      );
      if ($left->movers)
      {
        $l->add_tag_value("remove", $_) foreach (@{$left->movers});
      }
      $l->update();
#      print $l->Tag_enzyme, q{ }, $l->start, "\n";
    }
    if ($right)
    {
      my $r = $planchr->db->new_feature(
        -index        => 1,
        -seq_id       => $planchr->seq_id,
        -start        => $right->start,
        -end          => $right->end,
        -score        => sprintf("%.3f", $right->score),
        -display_name => $right->name,
        -primary_tag  => "enzyme_recognition_site",
        -source       => "BIO",
        -strand       => $right->strand,
        -attributes   => {
          "presence"    => $right->presence,
          "infeat"      => $right->featureid,
          "enzyme"      => $right->id,
          "load_id"     => $right->name,
          "wanthang"    => $right->phang,
          "peptide"     => $right->peptide,
          "ohangoffset" => $right->offset,
          "megachunk"   => $megachunkname,
          "chunks"      => [$chunkname, $nchunkname]}
      );
      if ($right->movers)
      {
        $r->add_tag_value("remove", $_) foreach (@{$right->movers});
      }
      $r->update();
#      print $r->Tag_enzyme, q{ }, $r->start, "\n";
    }
   
    #Make a chunk
    my $cstart = $left  ? $left->start  : $mcstart;
       $cstart = $lastenz ? $lastenz->start  : $mcstart;
    my $cend = $right ? $right->end : $mcend;
    my $c = $planchr->db->new_feature(
          -seq_id       => $planchr->seq_id,
          -start        => $cstart,
          -end          => $cend,
          -primary_tag  => "chunk",
          -display_name => $chunkname,
          -source       => "BIO",
          -index        => 1,
          -attributes   => {
            "load_id"     => $chunkname,
            "intro"       => $p{INTRO},
            "bsversion"   => $p{BSVERSION},
            "megachunk"   => $megachunkname}
    );
    $c->add_tag_value("enzymes", $left->id) if ($left);
    $c->add_tag_value("enzymes", $lastenz->id) if ($lastenz);
    $c->add_tag_value("enzymes", $right->id) if ($right);
    $c->update;
    print "\t$chunkname (", $c->start, "..", $c->end, ")\n";
    $lastenz = $right ? $right  : undef;
    $chunknum++;
  }
  print "\n\n";
  $letter++;
  $letter = substr($letter, -1) x length($letter) if (length($letter) > 1);
}

if ($p{LASTENZ} && $p{FINISHAT} && $p{LASTLETTER})
{
  my $oletter = $p{LASTLETTER};
  my $mchname = $p{NEWNAME} . ".$oletter";
  my $chname = $mchname . $p{FINISHAT};
  my $leftname = $p{LASTENZ};
  my @prepos = split("_", $leftname);
  my $pos = $prepos[1];
  my $res = $REDB->search( -left => $pos-16, -right => $pos+16);
  my @rearr = grep {$_->name eq $leftname} @{$res};
  die ("oops, no lastenz?!?\n") unless (scalar @rearr == 1);
  my $c = $planchr->db->new_feature(
    -index        => 1,
    -seq_id       => $planchr->seq_id,
    -start        => $rearr[0]->start,
    -end          => $excisor->start,
    -primary_tag  => "chunk",
    -display_name => $chname,
    -source       => "BIO",
    -attributes   => {
      "load_id"     => $chname,
      "intro"       => $p{INTRO},
      "bsversion"   => $p{BSVERSION},
      "megachunk"   => $mchname}
  );
  print "\t$chname (", $c->start, "..", $c->end, ")\n";
}

my @NEWGFF;
  my $seq_stream = $planchr->db->get_seq_stream()  or die "failed to get_seq_stream()";
  my @featarr;
  while (my $seq = $seq_stream->next_seq)
  {
    push @featarr, $seq;
  }
  @featarr = sort {$a->start <=> $b->start || (($b->end - $b->start) <=> ($a->end - $a->start))} @featarr;
  push @NEWGFF, gff3_string($_) . "\n" foreach (@featarr);

  push @NEWGFF, @{print_as_fasta($planchr->sequence, $planchr->seq_id, 1)};
  my $newfile = $planchr->path_to_GFF;
  open (my $OUT, '>', "$newfile") || die "can't open test output, $newfile $OS_ERROR\n";
  print $OUT @NEWGFF;
  close $OUT;
  
exit;

################################################################################
################################# SUBROUTINES ##################################
################################################################################
sub suitable_overhang
{
  my ($used_overhangs_ref, $overhang, $type) = @_;
  return 0 if (! $overhang);
  my $gnahrevo = $GD->complement(-sequence => $overhang, -reverse => 1);
  my %hsh = %{$used_overhangs_ref};
  return 0 if (  $overhang eq $gnahrevo  );
  return 0 if (  exists( $hsh{$overhang} )  && exists($hsh{$overhang}->{$type}));
  return 0 if (  exists( $hsh{$gnahrevo} )  && exists($hsh{$gnahrevo}->{$type}));
  return 1;
}

sub wtISS
{
  my ($left, $right, $wtchr) = @_;
  ($left, $right) = ($right, $left) if ($left > $right);
  my $size = abs($right - $left + 1);
 
  my @lefngenes  = sort {abs($a->end - $left)   <=> abs($b->end - $left)  } @genes;
  my @rigngenes  = sort {abs($a->start - $right) <=> abs($b->start - $right)} @genes;
  my ($lefgene, $riggene) = ($lefngenes[0], $rigngenes[0]);
 
  my @wtlefgenes = $wtchr->db->features(
      -seq_id => $wtchr->seq_id,
      -type => "gene",
      -name => $lefgene->display_name);
  my @wtriggenes = $wtchr->db->features(
      -seq_id => $wtchr->seq_id,
      -type => "gene",
      -name => $riggene->display_name);
  my ($wtlefgene, $wtriggene) = ($wtlefgenes[0], $wtriggenes[0]);
 
  my $lefendoffset   = $lefgene->end - $left;
  my $rigstartoffset = $right - $riggene->start;
 
  my $wtstart = $wtlefgene->end - $lefendoffset + 1;
  my $wtend   = $wtriggene->start + $rigstartoffset - 1;
 
  my $wtISS = lc substr($wtchr->sequence, $wtstart - 1, $wtend - $wtstart + 1);
  my $checkseq = substr($wtISS, $size/2);
 
  return $checkseq;
}

#Given a position along the chromosome, nominates IIB sites
sub ProposeExcisors
{
  my ($lastpos, $marker, $p) = @_;
  my @IIBs = grep {$RES->{$_}->class eq "IIB"} keys %{$RES};
  if ($p->{MARKERS})
  {
    @IIBs  = grep {! exists($marker->static_enzymes->{$_})} @IIBs;
  }
  my %lookfors = map {$_ => 1} @IIBs;
  my $right = $lastpos + $p->{ISSMIN};
  my $left = $right - $p->{CHUNKLENMAX};
  my $rawgrab = $REDB->search( -left => $left, -right => $right );
  my @pool = grep {exists $lookfors{$_->id}} @{$rawgrab};
  my %igenics = map {$_->id => 1} grep {$_->presence eq "intergenic"} @pool;
  @IIBs = grep {! exists $igenics{$_}} @IIBs;
  my $wtISS = wtISS($lastpos, $right, $wtchr);
  my @finalists;
  foreach my $enz (@IIBs)
  {
    my $checkhsh = $RES->{$enz}->positions($wtISS);
    next if (scalar(keys %{$checkhsh}) != 0);
    my ($ess_flag, $fgw_flag, $lap_flag) = (0, 0, 0);
    my $size = $RES->{$enz}->len();
    my @movers;
    my @exonics = grep {$_->presence eq "existing" && $_->id eq $enz} @pool;
    my @igenics = grep {$_->presence eq "intergenic" && $_->id eq $enz} @pool;
    next if (scalar(@igenics));
    foreach my $ex (@exonics)
    {
      my $egene = $ex->featureid;
      my $emaskbit = $mask->count_features_in_range($ex->start - 1, $size);
      $lap_flag++ if ($emaskbit > 1);
      $ess_flag++ if ($ESSENTIALS{$egene} eq "Essential");
      $fgw_flag++ if ($ESSENTIALS{$egene} eq "fast_growth");
      push @movers, $ex->name;
    }
    next if ($lap_flag != 0);
    my $score = 0;
    $score   += log($RES->{$enz}->score);
    $score   += 0.1 * scalar(@movers);
    $score   += 0.5 * $fgw_flag;
    $score   += 1.0 * $ess_flag;
    my $name = $p->{MARKERS}  ? $marker->name : "ISS";
    my $newfeat = Bio::BioStudio::RestrictionEnzyme->new(
      -enzyme => $RES->{$enz},
      -name => $enz . "_" . $name,
      -score => $score,
      -presence => "a",
      -start => $right
    );
    $newfeat->movers(\@movers) if (scalar(@movers));
    push @finalists, $newfeat;
  }
  return \@finalists;
}

#uses globals $RES, and $ESSENTIALS
#Given a position along the chromosome, nominates low penalty restriction sites
sub ProposeCandidates
{
  my ($lastpos, $p, $flag) = @_;
  my ($min, $max) = $flag ==  0
      ? ($p->{CHUNKLENMIN}, $p->{CHUNKLENMAX})
      : ($p->{CHUNKLENMIN}, $flag);
  my $rcpos = $lastpos - $min > 0 ?  $lastpos - $min  : 1;
  my $lcpos = $lastpos - $max > 0 ?  $lastpos - $max  : 1;
  my $lef = $lcpos - $p->{CHUNKLENMAX} > 0 ? $lcpos - $p->{CHUNKLENMAX} : 1;
  my $rig = $flag ? $flag + $lastpos : $lastpos;
  print "CHECKING  $lcpos, $rcpos, $lef, $rig, $lastpos, $flag\n" if ($debug);
  my $pool = $REDB->search(
        -left => $lef - $p->{CHUNKOLAP},
        -right => $rig + $p->{CHUNKOLAP});
  
  my @candidates = grep {$_->start >= $lcpos && $_->end <= $rcpos} @{$pool};
  my @finalists;
  foreach my $re (@candidates)
  {
    next if ($re->eligible);
    my $size = $re->end - $re->start + 1;
    my $maskbit = $mask->count_features_in_range($re->start - 1, $size);
    #Drop if exonic in exon overlap
    next if ($re->presence ne "intergenic" && $maskbit > 1);
    my @buds = grep {abs($_->start - $re->start) <= $p->{CHUNKLENMAX}
                && $_->presence  ne "potential" && $_->id eq $re->id
                && $_->name ne $re->name}
               @{$pool};
    my @igenics = grep {$_->presence eq "intergenic"} @buds;
    #Drop if too many intergenics around
    next if (scalar(@igenics));
    my $gene = $re->presence ne "intergenic"  ? $re->featureid  : undef;
    my @exonics = grep {$_->presence eq "existing"} @buds;
    my $ess_flag = $gene && $ESSENTIALS{$gene} eq "Essential"   ? 1 : 0;
    my $fgw_flag = $gene && $ESSENTIALS{$gene} eq "fast_growth" ? 1 : 0;
    my $lap_flag = 0;
    my @movers;
    foreach my $ex (@exonics)
    {
      next if ($ex->start > $rig + $p->{CHUNKOLAP});
      my $egene = $ex->featureid;
      my $emaskbit = $mask->count_features_in_range($ex->start - 1, $size);
      $lap_flag++ if ($emaskbit > 1);
      last if $lap_flag;
      $ess_flag++ if ($ESSENTIALS{$egene} eq "Essential");
      $fgw_flag++ if ($ESSENTIALS{$egene} eq "fast_growth");
      push @movers, $ex->name;
    }
    #Drop if it requires modification in exon overlap
    next if ($lap_flag != 0);
    my $distance = abs($re->start - $lcpos);
   
    # Score is a function of distance from largest possible mark
    my $score = $distance <= 5000 ?  0  :  0.0008  * $distance - 4;
    # Plus the log of the price per unit
    $score   += log($RES->{$re->id}->score);
    # Plus one tenth point for each orf modified
    $score   += 0.1 * scalar(@movers);
    # Plus one half point for each fast growth orf modified
    $score   += 0.5 * $fgw_flag;
    # Plus one point for each essential orf modified
    $score   += 1.0 * $ess_flag;
    # Ignore if score is too high
    next if ($score > 1);
    $re->score($score);
    $re->movers(\@movers) if (scalar(@movers));
    push @finalists, $re;
  }
  return \@finalists;
}

sub NextCandidate
{
  my ($candlist, $usedhangref, $hangsurvey, $usedsiteref,
      $flag, $prevenz, $excisor, $omarker, $p) = @_;
  my $foundcand = undef;
  while (! $foundcand && scalar(@{$candlist}))
  {
    my $candidate = shift @{$candlist};
    my $enz = $candidate->id;
    next if (exists $usedsiteref->{$candidate->id});
    if ($flag)
    {
      my $ISSseq = $altmask->count_features_in_range($candidate->end, $p->{ISSMIN});
      if ($ISSseq > 0)
      {
        print "\t\t\t\t discard: candidate $enz exists in an ISS inviable location...\n" if ($debug);
        next;
      }
      my $wtISSseq = wtISS($candidate->start, $candidate->start + $p->{ISSMIN}, $wtchr);
      my $ihsh = $candidate->positions($wtISSseq);
      if (scalar(keys %{$ihsh}))
      {
        print "\t\t\t\t discard: candidate $enz occurs in the targeting wildtype sequence...\n" if ($debug);
        next;
      }
      if ($p->{MARKER} && exists($omarker->static_enzymes->{$candidate->id}))
      {
        print "\t\t\t\t discard: candidate $enz occurs in the intergenic region of the target marker...\n" if ($debug);
        next;
      }
    }
    my $enztype = $candidate->type;
    my %hanghsh = %{$candidate->overhangs};
    my @changs = keys %hanghsh;
    @changs = sort {$hangsurvey->{"$a.$enztype"} <=> $hangsurvey->{"$b.$enztype"}} @changs;
    @changs = grep { suitable_overhang( $usedhangref, $_, $enztype ) && $_} @changs;
    
    unless (scalar @changs)
    {
      print "\t\t\t\t discard: candidate $enz has no suitable overhangs...\n" if ($debug);
      next;
    }
    my $hang = undef;
    if (! $prevenz || $candidate->presence eq "intergenic")
    {
      $hang = $changs[0];
    }
    else
    {
      foreach my $phang (@changs)
      {
        my ($oldseq, $newseq) = NewSequence($candidate, $phang, $p->{CODON_TABLE});
        unless($newseq)
        {
          next;
        }
        my $oSTATUS = $GD->restriction_status(-sequence => $oldseq);
        my $nSTATUS = $GD->restriction_status(-sequence => $newseq);
        next if ($nSTATUS->{$prevenz->id} != $oSTATUS->{$prevenz->id} && $prevenz->id ne $candidate->id);
        next if ($nSTATUS->{$candidate->id} != 1);
        next if ($excisor && $nSTATUS->{$excisor->id} != 0);
        $hang = $phang;
        my @createds = grep {$nSTATUS->{$_} > $oSTATUS->{$_} && $_ ne $candidate->id} keys %{$RES};
        $candidate->creates(\@createds) if (scalar @createds);
        last;
      }
    }
    unless ($hang)
    {
      print "\t\t\t\t discard: candidate $enz would introduce an untenable enzyme as a side effect...\n" if ($debug);
      next;
    }
    $candidate->phang($hang);
    my $gnah = $GD->complement($hang, 1);
    $usedhangref->{$hang} = {} unless (exists $usedhangref->{$hang});
    $usedhangref->{$gnah} = {} unless (exists $usedhangref->{$gnah});
    $usedhangref->{$hang}->{$enztype}++;
    $usedhangref->{$gnah}->{$enztype}++;
    $foundcand = $candidate;
  }
  return $foundcand ? $foundcand  : undef;
}

sub NewSequence
{
  my ($candidate, $ohang) = @_;
  print "\t\t\t determining what the chosen overhang for " . $candidate->id . " will look like in place...\n" if ($debug);
  my $start = $candidate->start;
  my $end = $candidate->end;
  my $strand = $candidate->strand;
  my $site = $candidate->recseq;
  my $cut = $candidate->cutseq;
  my ($lef, $rig) = (0, 0);
  ($lef, $rig) = ($1, $2) if ($cut =~ $candidate->class_regexes->{IIA});
  ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
  $rig = 0 if ($rig < 0);
  my @CDSes = $planchr->db->get_features_by_name($candidate->featureid);
  my $CDS = $CDSes[0];
  my $offset = $candidate->offset || 0;
  if ($CDS->strand == -1)
  {
    $end ++ while (($CDS->end - $end) % 3 != 0);
    $start -- while (($end - $start + 1) /3  < length($candidate->peptide));
  }
  else
  {
    $start -- while (($start - $CDS->start) % 3 != 0);
    $end ++ while (($end - $start + 1) /3  < length($candidate->peptide));
  }
  my $testseq = substr($pchrseq, $start - 1, $end - $start + 1);
  $testseq = $GD->complement($testseq, 1) if ($CDS->strand == -1);
  my $aa = $GD->translate(-sequence => $testseq);
  $site = $GD->complement($site, 1) if ($strand == -1);
  my $alnmat = $GD->pattern_aligner(-sequence => $testseq,
                                    -pattern => $site,
                                    -peptide => $aa);
  return (undef, undef) unless ($alnmat);
  unless ($candidate->type eq "b")
  {
    substr($alnmat, $offset, length($ohang)) = $ohang;
  }
  my $newpatt = $GD->pattern_adder(-sequence => $testseq, -pattern => $alnmat);
  $newpatt = $GD->complement($newpatt, 1) if ($CDS->strand == -1);
  $testseq = $GD->complement($testseq, 1) if ($CDS->strand == -1);
  my $pattlen = length($newpatt);
  my $contextlen = 60 + $pattlen;
  $start -=20;
  my $context = substr($pchrseq, $start -1, $contextlen);
  my $oldtext = $context;
  substr($context, 20, $pattlen) = $newpatt;
  return ($oldtext, $context);
}

__END__

=head1 NAME

  BS_ChromosomeSegmentationPlanner.pl
 
=head1 VERSION

    Version 2.00

=head1 DESCRIPTION

  This utility proposes a set of restriction enzyme recognition site changes to
    a chromosome. The goal is to make the chromosome assemblable from
    multikilobase pieces called chunks, which in groups of roughly CHUNKNUM form
    larger pieces called megachunks. Megachunks end in special regions called
    InterSiteSequences (ISS), which consist of synthetic sequence, followed by a
    marker, followed by wild-type sequence, followed by a type IIB restriction
    enzyme recognition site.  The wild-type sequence targets the megachunk to
    its target chromosome for homologous recombination (thus wild type here can
    be any homologous chromosome, as long as gene order is the same). Megachunks
    alternate markers to allow a simple selection for successful integration.
    Markers should be defined in config/markers.
 
  No edits will be made, this utility merely proposes a plan for you to vet.
    If the plan meets your requirements, run BS_ChromosomeSegmenter to implement
    it.
 
=head1 ARGUMENTS

Required arguments:

  -C,  --CHROMOSOME : The chromosome to be modified
  --RESET  : Which list of restriction enzymes to use (defined in GeneDesign)
  -W,  --WTCHR  : The chromosome that will receive chunks (usually wildtype)
  -M,  --MARKERS : Comma separated list which will be alternately inserted
      into megachunk ISS sequences (must be defined in config/markers)
 
Optional arguments:

  -SC,  --SCOPE : [seg, chrom (def)] How much sequence will the edit
                  affect. seg requires STARTPOS and STOPPOS.
  -STA, --STARTPOS : The first base for analysis;ignored unless SCOPE = seg
  -STO, --STOPPOS  : The last base for analysis;ignored unless SCOPE = seg
  --CHUNKLENMIN : The minimum size for chunks (default 6000 bp)
  --CHUNKLENMAX : The maximum size for chunks (default 9920 bp)
  --CHUNKNUM : The target number of chunks per megachunk (default 4)
  --CHUNKNUMMIN : The minimum number of chunks per megachunk (default 3)
  --CHUNKNUMMAX : The maximum number of chunks per megachunk (default 5)
  --CHUNKOLAP : The number of bases each chunk must overlap (default 40)
  --ISSMIN : Minimum size of the homologous intersite sequence (default 900)
  --ISSMAX : Maximum size of the homologous intersite sequence (default 1500)
  --FIRSTLETTER : The first letter to assign a megachunk (default 'A')
  -FP, --FPUTRPADDING : No edit zone upstream of the five prime end of
            essential/fast growth genes when no UTR is annotated (default 500)
  -TP, --TPUTRPADDING : No edit zone downstream of the three prime end of
            essential/fast growth genes when no UTR is annotated (default 100)
  -L, --LASTMARKER : Which marker should be the last marker inserted (must be
            defined in config/markers)
  --REDOREDB : Force recreation of the restriction enzyme database. Useful if
            the chromosome sequence or the chunk size parameters change between
            uses of the Segmentation Planner
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