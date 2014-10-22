#!/usr/bin/env perl

use Bio::BioStudio;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_ChromosomeSplicer_$VERSION";

local $| = 1;

my %ACTIONS = (segmentflank => 1, featflank => 1, featins => 1);

my %p;
GetOptions (
      'CHROMOSOME=s'      => \$p{CHROMOSOME},
      'EDITOR=s'          => \$p{EDITOR},
      'MEMO=s'            => \$p{MEMO},
      'SCALE=s'           => \$p{SCALE},
      'ACTION=s'          => \$p{ACTION},
      'STARTPOS=i'        => \$p{STARTPOS},
      'STOPPOS=i'         => \$p{STOPPOS},
      'FEATURE=s'         => \$p{FEATURE},
      'INSERT=s'          => \$p{INSERT},
      'DISTANCE=i'        => \$p{DISTANCE},
      'DIRECTION=i'       => \$p{DIRECTION},
      'NAME=s'            => \$p{NAME},
      'OUTPUT=s'          => \$p{OUTPUT},
      'DESTROY'           => \$p{DESTROY},
    	'help'              => \$p{HELP}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n"  unless ($p{CHROMOSOME});
my $chr    = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $chrseq = $chr->sequence;
my $chrlen = length $chrseq;

my $BS_FEATS = $BS->fetch_custom_features();

$p{SCALE}     = $p{SCALE}     || 'chrom';
$p{OUTPUT}    = $p{OUTPUT}    || 'txt';
$p{STARTPOS}  = $p{STARTPOS}  || 1;
$p{STOPPOS}   = $p{STOPPOS}   || $chrlen;
$p{DISTANCE}  = $p{DISTANCE}  || 0;
$p{DESTROY}   = $p{DESTROY}   || undef;

unless ($p{EDITOR} && $p{MEMO})
{
  print "\n ERROR: Both an editor's id and a memo must be supplied.\n\n";
}

unless ($p{ACTION})
{
  print "\n BSERROR: No action was specified.\n\n";
  pod2usage(-verbose=>99, -sections=>"ARGUMENTS");
}
unless (exists $ACTIONS{$p{ACTION}})
{
  die "\n BSERROR: Unrecognized action requested.\n";
}
if ($p{ACTION} eq "featflank" && ! $p{FEATURE})
{
  print "\n BSERROR: A feature based action was requested, ";
  print "but no feature was specified.\n";
  die();
}
if ($p{ACTION} eq "featflank" && (! $p{DISTANCE} || ! $p{DIRECTION}))
{
  print "\n BSERROR: Feature flank requested but distance or direction ";
  print "not specified.\n";
  die();
}

if ($p{DIRECTION} &&
   ! ($p{DIRECTION} == 3 || $p{DIRECTION} == 5 || $p{DIRECTION} == 35))
{
  die "\n BSERROR: Feature flank requested but direction not 3, 5, or 35.\n";
}
if ($p{ACTION} eq "featins" && ! $p{NAME})
{
  die "\n BSERROR: An insertion was requested but not named (use -N).\n";
}
if ($p{ACTION} =~ /ins/ && ! $p{INSERT})
{
  print "\n BSERROR: A insertion action was requested but no custom feature ";
  print "was specified for insertion.\n";
  die();
}
if ($p{INSERT} && ! exists $BS_FEATS->{$p{INSERT}})
{
  die "\n BSERROR: Unrecognized custom feature requested for insertion.\n";
}
if ($p{STOPPOS} <= $p{STARTPOS})
{
  die "\n BSERROR: The start and stop coordinates do not parse.\n";
}

################################################################################
################################## SPLICING  ###################################
################################################################################
my $newchr = $chr->iterate();
$p{REPORT} = {};
my @changes;

if ($p{ACTION} eq 'featins')
{
  my $bsfeat = $BS_FEATS->{$p{INSERT}}->clone();
  my $insname = $p{NAME} || $bsfeat->display_name . '_' . $newchr->signature;
  $insname =~ s{\s}{\_}msxg;
  $bsfeat->display_name($insname);

  my $newfeat = eval
  {
    $newchr->insert_feature(
      -feature => $bsfeat,
      -position => $p{STARTPOS},
      -destroy  => $p{DESTROY}
    );
  };
  my $e;
  if ($e = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
  {
    print "Can't insert $insname; would destroy " . $e->error . "\n";
  }
  elsif ($newfeat)
  {
    push @changes, $newfeat;
  }
}

elsif ($p{ACTION} eq 'segmentflank')
{
  my $lbsfeat = $BS_FEATS->{$p{INSERT}}->clone();
  my $linsname = $p{NAME} || $lbsfeat->display_name . '_' . $newchr->signature;
  $linsname .= '_5';
  $linsname =~ s{\s}{\_}msxg;
  $lbsfeat->display_name($linsname);
  my $movelen = length($lbsfeat->seq->seq);
  my $newLfeat = eval
  {
    $newchr->insert_feature(
      -feature => $lbsfeat,
      -position => $p{STARTPOS},
      -destroy  => $p{DESTROY}
    );
  };
  my $e1;
  if ($e1 = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
  {
    print "Can't 5' flank segment; would destroy " . $e1->error . "\n";
  }

  my $rbsfeat = $BS_FEATS->{$p{INSERT}}->clone();
  my $rinsname = $p{NAME} || $rbsfeat->display_name . '_' . $newchr->signature;
  $rinsname .= '_3';
  $rinsname =~ s{\s}{\_}msxg;
  $rbsfeat->display_name($rinsname);

  my $newRfeat = eval
  {
    $newchr->insert_feature(
      -feature => $rbsfeat,
      -position => $p{STOPPOS} + $movelen,
      -destroy  => $p{DESTROY}
    );
  };
  my $e2;
  if ($e2 = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
  {
    print "Can't 3' flank segment; would destroy " . $e2->error . "\n";
  }
  
  if ($newLfeat && $newRfeat)
  {
    push @changes, $newLfeat, $newRfeat;
  }
    
}

elsif ($p{ACTION} eq 'featflank')
{
  my @directs   = ();
  my %movehash  = ();
  push @directs, 5 if ($p{DIRECTION} =~ /5/);
  push @directs, 3 if ($p{DIRECTION} =~ /3/);
  my @targets  = $newchr->db->features(
    -seqid      => $newchr->seq_id,
    -start      => $p{STARTPOS},
    -end        => $p{STOPPOS},
    -range_type => 'contains',
    -type       => $p{FEATURE}
  );
  if (! scalar @targets)
  {
    print "\n ERROR: There are no features of type $p{FEATURE} in the ";
    print "interval $p{STARTPOS}..$p{STOPPOS}.\n";
    die();
  }
  foreach my $prefeat (@targets)
  {
    my $feature = $newchr->db->fetch($prefeat->primary_id);
    my $featid = $feature->display_name;
    my $movelen = 0;
    foreach my $direct (@directs)
    {
      my $insfeat = $BS_FEATS->{$p{INSERT}}->clone();
      my $insname = $p{NAME} || $insfeat->display_name . '_' . $newchr->signature;
      $insname .= "_$featid" . "_$direct";
      $insname =~ s{\s}{\_}msxg;
      $insfeat->display_name($insname);
      
      my ($start, $stop) = (undef, undef);
      if (  ($direct == 3 && $feature->strand  > -1)
         || ($direct == 5 && $feature->strand == -1))
      {
        $start = $feature->end() + $p{DISTANCE} + 1;
      }
      elsif ( ($direct == 3 && $feature->strand == -1)
           || ($direct == 5 && $feature->strand  > -1))
      {
        $start = $feature->start() - $p{DISTANCE};
      }
      my $newfeat = eval
      {
        $newchr->insert_feature(
          -feature => $insfeat,
          -position => $start,
          -destroy  => $p{DESTROY}
        );
      };
      my $e;
      if ($e = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
      {
        print "Can't $direct' flank $featid; would destroy " . $e->error . "\n";
      }
      elsif ($newfeat)
      {
        push @changes, $newfeat;
      }
    }
  }
}

################################################################################
################################### WRITING ####################################
################################################################################

if (scalar @changes)
{
  $newchr->write_chromosome();
}
else
{
  print "No changes made - no new version generated.\n";
}

exit;

__END__

=head1 NAME

  BS_ChromosomeSplicer.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility inserts features into a chromosome. Custom features should be
    defined in config/features.txt. You must define a sequence segment for the
    edits;it can be the entire chromosome.
   
  The utility can be run in three modes:
    segmentflank - In this mode, a custom feature will be inserted at both
      ends of the defined sequence segment. The insert will interrupt any
      feature that overlaps the segment edges.
    featflank - In this mode, all instances of a feature in the defined sequence
      segment will be flanked (3', 5', or both) at a defined distance with an
      instance of a custom feature. This insertion will NOT interrupt any
      features;if a target feature cannot be flanked without interrupting
      another chromosomal feature, you will be informed.
    featins - In this mode, a custom feature will be inserted at the five prime
      end of the defined sequence segment. You must give the feature a unique
      name.

=head1 ARGUMENTS

Required arguments:

  -C,   --CHROMOSOME : The chromosome to be modified
  -E,   --EDITOR : The person responsible for the edits
  -M,   --MEMO : Justification for the edits
  -A,   --ACTION : The action to take.
            segmentflank : flank the defined segment with custom features
                    requires -I
            featflank : flank a feature type in the defined segment with
                    custom features;requires -I, -F, -DIS, -DIR
            featins : insert a custom feature at the five prime end of the
                    defined segment;requires -I, -N

Optional arguments:

  -SC,  --SCALE : [genome, chrom (def)] Which version number to increment
  -STA, --STARTPOS : The first base eligible for editing, defaults to 1
  -STO, --STOPPOS : The last base eligible for editing, defaults to chr length
  -DES, --DESTROY : Whether or not other features should be disrupted by
            insertions; default 0 for no
  -F,   --FEATURE : The feature to be targeted
  -I,   --INSERT : The custom feature to be inserted;must be an entry in
            config/features.txt
  -DIS, --DISTANCE : The distance in base pairs from a feature for an
            insertion; required for action featflank
  -DIR, --DIRECTION : [3, 5, 35] Should custom features be inserted 3p, 5p, or
            both 3p and 5p of a feature;required for action featflank
  -N,   --NAME : What to name the inserted feature;required for action featins
  -OUT, --OUTPUT : [html, txt (def)] Format for reporting and output
  -h,   --help : Display this message
 
=cut
