#!/usr/bin/env perl

use Bio::BioStudio;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_SeqChecker_$VERSION";

local $| = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'  => \$p{CHROMOSOME},
      'OUTPUT=s'      => \$p{OUTPUT},
    	'help'          => \$p{HELP}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $BS = Bio::BioStudio->new();

$p{OUTPUT} = $p{OUTPUT} || 'txt';

die "BSERROR: No chromosome was named" unless ($p{CHROMOSOME});
my $chr = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});

################################################################################
################################# CONFIGURING ##################################
################################################################################
$p{BSVERSION} = $bsversion;
my $BS_FEATS = $BS->fetch_custom_features();
my %BSKINDS = map {$_->primary_tag => 1} values %{$BS_FEATS};


################################################################################
############################### ERROR  CHECKING ################################
################################################################################
print "Checking ", $chr->name(), "\n";
my $GD = $chr->GD();
my $chrseq = $chr->sequence();
my @features = $chr->db->features();
foreach my $feature (@features)
{
  my $fname = $feature->display_name;
  my $trueseq = $feature->seq->seq;
  if ($feature->has_tag('newseq'))
  {
    my $annseq = $feature->Tag_newseq;
    if ($annseq ne $trueseq)
    {
      print "WARNING: $fname has bad newseq tag!\n";
      print "\t$annseq tag vs $trueseq actual\n";
    }
  }
  if (exists $BSKINDS{$feature->primary_tag})
  {
    foreach my $bsfeature (values %{$BS_FEATS})
    {
      if ($feature =~ $bsfeature->display_name)
      {
        my $shouldseq = $bsfeature->sequence;
        my $isseq = $feature->seq->seq;
        $isseq = $GD->complement($isseq, 1) if ($feature->strand == -1);
        if ($shouldseq ne $isseq)
        {
          print "WARNING: $feature sequence looks weird!\n";
          print "\t$shouldseq != $isseq actual\n";
        }
      }
    }
  }
}
print "\n";

exit;

__END__

=head1 NAME

  BS_SeqChecker.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility checks the sequences of the features in a chromosome against
  their annotations and against the expected sequences found in BioStudio
  configuration.

=head1 ARGUMENTS

Required arguments:

  -C,   --CHROMOSOME : The chromosome to be checked

Optional arguments:

  -O,   --OUTPUT : [html, txt (def)] Format for reporting and output
  -h,   --help : Display this message
 
=cut
