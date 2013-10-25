#!/usr/bin/env perl

use Bio::BioStudio;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_ChromosomeCutter_$VERSION";

local $| = 1;

my @actionnames = qw(seqdel seqdelprp featdel featdelprp listdel);
my %ACTIONS = map {$_ => 1} @actionnames;

my %p;
GetOptions (
  'CHROMOSOME=s'   => \$p{CHROMOSOME},
  'EDITOR=s'       => \$p{EDITOR},
  'MEMO=s'         => \$p{MEMO},
  'SCALE=s'        => \$p{SCALE},
  'ACTION=s'       => \$p{ACTION},
  'STARTPOS=i'     => \$p{STARTPOS},
  'STOPPOS=i'      => \$p{STOPPOS},
  'TYPE=s'         => \$p{TYPE},
  'FEATURES=s'     => \$p{FEATURES},
  'INSERT=s'       => \$p{INSERT},
  'OUTPUT=s'       => \$p{OUTPUT},
  'DESTROY'        => \$p{DESTROY},
	'help'           => \$p{HELP}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
############################### SANITY CHECKING ################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n"  unless ($p{CHROMOSOME});
my $chr    = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $chrseq = $chr->sequence;
my $chrlen = length $chrseq;

my $BS_FEATS = $BS->fetch_custom_features();

$p{STARTPOS} = $p{STARTPOS} || 1;
$p{STOPPOS} = $p{STOPPOS} || $chrlen;
$p{SCALE}  = $p{SCALE}  || 'chrom';
$p{OUTPUT} = $p{OUTPUT} || 'html';

unless ($p{EDITOR} && $p{MEMO})
{
  print "\n ERROR: Both an editor's id and a memo must be supplied.\n\n";
}

unless ($p{ACTION})
{
  print "\n ERROR: No action was specified.\n\n";
  pod2usage(-verbose=>99, -sections=>"ARGUMENTS");
}
if (! exists $ACTIONS{$p{ACTION}})
{
  die "\n ERROR: Unrecognized action requested.\n";
}
if ($p{ACTION} =~ /feat/ && ! $p{TYPE})
{
  print "\n ERROR: A feature type based action was requested, ";
  print "but no feature type was specified.\n";
  die;
}
if ($p{ACTION} eq "listdel" && ! $p{FEATURES})
{
  print "\n ERROR: A list based action was requested, ";
  print "but no features were specified.\n";
  die;
}
if ($p{INSERT} && ! exists $BS_FEATS->{$p{INSERT}})
{
  die "\n ERROR: Unrecognized custom feature requested for insertion.\n";
}
if ($p{STOPPOS} <= $p{STARTPOS})
{
  die "\n ERROR: The start and stop coordinates do not parse.\n";
}

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $newchr = $chr->iterate();
$p{REPORT} = {};
my @changes;


################################################################################
################################### LISTDEL ####################################
################################################################################
if ($p{ACTION} eq "listdel")
{
  my $pretlist = $p{FEATURES};
  my @tlist = split(",", $pretlist);
  my @glist = ();
  foreach my $seek (@tlist)
  {
    my @res = $newchr->db->features(
      -seqid      => $newchr->seq_id,
      -start      => $p{STARTPOS},
      -end        => $p{STOPPOS},
      -name       => $seek,
      -range_type => "contains",
    );
    if (scalar(@res) == 0)
    {
      print "Can't find a feature called $seek - skipping.\n";
      next;
    }
    elsif (scalar(@res) > 1)
    {
      print "$seek is not specific enough as a feature name - skipping.\n";
      next;
    }
    push @glist, @res;
  }

  my %movehash = ();
  foreach my $prefeat (@glist)
  {
    my $featid   = $prefeat->primary_id;
    my $feature  = $newchr->db->fetch($featid);
    my $start = $feature->start;
    my $delfeat = eval
    {
      $newchr->delete_region(
        -start => $start,
        -stop  => $feature->end
      );
    };
    my $e1;
    if ($e1 = Bio::BioStudio::Exception::DeleteFeature->caught())
    {
      print "Can't delete region: " . $e1->error . "\n";
    }
    elsif ($delfeat)
    {
      push @changes, $delfeat;
    }
  
    if ($p{INSERT})
    {
      my $bsfeat = $BS_FEATS->{$p{INSERT}}->clone();
      my $insname = $bsfeat->display_name . '_' . $delfeat->display_name;
      $insname =~ s{\s}{\_}msxg;
      $bsfeat->display_name($insname);
      my $newfeat = eval
      {
        $newchr->insert_feature(
          -feature  => $bsfeat,
          -position => $start,
          -destroy  => $p{DESTROY}
        );
      };
      my $e2;
      if ($e2 = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
      {
        print "Can't insert $insname; would destroy " . $e2->error . "\n";
      }
      elsif ($newfeat)
      {
        push @changes, $newfeat;
      }
    }
  }
}

################################################################################
################################### FEATDEL ####################################
################################################################################
if ($p{ACTION} eq "featdel")
{
  my @targets  = $newchr->db->features(
    -seqid      => $newchr->seq_id,
    -start      => $p{STARTPOS},
    -end        => $p{STOPPOS},
    -type       => $p{TYPE},
    -range_type => "contains",
  );
  foreach my $prefeat (@targets)
  {
    my $featid   = $prefeat->primary_id;
    my $feature  = $newchr->db->fetch($featid);
    my $start = $feature->start;
    my $delfeat = eval
    {
      $newchr->delete_region(
        -start => $start,
        -stop  => $feature->end
      );
    };
    my $e1;
    if ($e1 = Bio::BioStudio::Exception::DeleteFeature->caught())
    {
      print "Can't delete region: " . $e1->error . "\n";
    }
    elsif ($delfeat)
    {
      push @changes, $delfeat;
    }
  
    if ($p{INSERT})
    {
      my $bsfeat = $BS_FEATS->{$p{INSERT}}->clone();
      my $insname = $bsfeat->display_name . '_' . $delfeat->display_name;
      $insname =~ s{\s}{\_}msxg;
      $bsfeat->display_name($insname);
      my $newfeat = eval
      {
        $newchr->insert_feature(
          -feature => $bsfeat,
          -position => $start,
          -destroy  => $p{DESTROY}
        );
      };
      my $e2;
      if ($e2 = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
      {
        print "Can't insert $insname; would destroy " . $e2->error . "\n";
      }
      elsif ($newfeat)
      {
        push @changes, $newfeat;
      }
    }
  }
}

################################################################################
################################# FEATDELPRP  ##################################
################################################################################
elsif ($p{ACTION} eq "featdelprp")
{
  my @targets  = $newchr->db->features(
    -seqid      => $newchr->seq_id,
    -start      => $p{STARTPOS},
    -end        => $p{STOPPOS},
    -range_type => "contains",
    -type       => $p{TYPE}
  );
  if (! scalar @targets)
  {
    die "There are no targets in range.\n";
  }
  foreach my $feature (@targets)
  {
    print "working on $feature\n";
    my $featname = $feature->display_name();
    my $pdel = Bio::SeqFeature::Generic->new(
      -start          => $feature->start,
      -end            => $feature->end,
      -primary_tag    => 'proposed_deletion',
      -display_name   => 'pdel_' . $featname
    );

    my $newfeat = eval
    {
      $newchr->add_feature(-feature => $pdel);
    };
    my $e;
    if ($e = Bio::BioStudio::Exception::AddFeature->caught())
    {
      print "Can't add $newfeat: " . $e->error . "\n";
    }
    elsif ($newfeat)
    {
      push @changes, $newfeat;
    }
  }
}

################################################################################
################################### SEQDEL  ####################################
################################################################################
elsif ($p{ACTION} eq "seqdel")
{
  my $delfeat = eval
  {
    $newchr->delete_region(
      -start => $p{STARTPOS},
      -stop  => $p{STOPPOS}
    );
  };
  my $e1;
  if ($e1 = Bio::BioStudio::Exception::DeleteFeature->caught())
  {
    print "Can't delete region: " . $e1->error . "\n";
  }
  elsif ($delfeat)
  {
    push @changes, $delfeat;
  }
  
  if ($p{INSERT})
  {
    my $bsfeat = $BS_FEATS->{$p{INSERT}}->clone();
    my $insname = $bsfeat->display_name . '_' . $delfeat->display_name;
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
    my $e2;
    if ($e2 = Bio::BioStudio::Exception::PreserveExsistingFeature->caught())
    {
      print "Can't insert $insname; would destroy " . $e2->error . "\n";
    }
    elsif ($newfeat)
    {
      push @changes, $newfeat;
    }
  }
}

################################################################################
################################## SEQDELPRP ###################################
################################################################################
elsif ($p{ACTION} eq "seqdelprp")
{
  my $segname = "$p{STARTPOS}..$p{STOPPOS}";
  my $pdel = Bio::SeqFeature::Generic->new(
    -start          => $p{STARTPOS},
    -end            => $p{STOPPOS},
    -primary_tag    => 'proposed_deletion',
    -display_name   => 'pdel_' . $segname
  );

  my $newfeat = eval
  {
    $newchr->add_feature(-feature => $pdel);
  };
  my $e;
  if ($e = Bio::BioStudio::Exception::AddFeature->caught())
  {
    print "Can't add $newfeat: " . $e->error . "\n";
  }
  elsif ($newfeat)
  {
    push @changes, $newfeat;
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

  BS_ChromosomeCutter.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility removes features from a chromosome and offers the chance to
    replace them with custom features. Custom features should be defined in
    config/features.txt. You must define a sequence segment for the edits;it
    can be the entire chromosome.
   
  The utility can be run in four modes:
    seqdel - The defined sequence segment will be deleted. You can have a
      custom feature inserted in its place.
    seqdelprp - A deletion proposal feature that spans the defined sequence
      segment will be created.  No sequence editing will take place.
    featdel - All features of the target type in the defined sequence segment
      will be deleted. You can have a custom feature inserted in their place.
    featdelprp - Deletion proposal features that span each feature of the target
      type in the defined sequence segment will be created. No sequence editing
      will take place.
    listdel - You provide a list of feature names; all of those features will be
      deleted.

=head1 ARGUMENTS

Required arguments:

  -C,   --CHROMOSOME : The chromosome to be modified
  -E,   --EDITOR : The person responsible for the edits
  -M,   --MEMO : Justification for the edits
  -A,   --ACTION : The action to take.
                seqdel : delete the defined segment
                seqdelprp : propose the defined segment for deletion
                featdel : delete features of type -T from the defined segment
                          requires -T
                featdelprp : propose features of type -T for deletion
                          requires -T
                listdel : delete the features listed in -F
                          requires -F

Optional arguments:

  -STA, --STARTPOS : The first base eligible for editing (default 1)
  -STO, --STOPPOS : The last base eligible for editing (default chr length)
  -SC,  --SCALE : [genome, chrom (def)] Which version number to increment
  -T,   --TYPE : The type of feature to be targeted
  -DES, --DESTROY : Whether or not other features should be disrupted by
            insertions; default 0 for no
  -F,   --FEATURES : A comma separated list of feature names to be deleted
  -I,   --INSERT : The feature to be inserted in each deletion;
                   must be an entry in config/features.txt
  -OU,  --OUTPUT : html or txt
  -h,   --help : Display this message

=cut
