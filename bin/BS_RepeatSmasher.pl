#!/usr/bin/env perl

use Bio::BioStudio;
use Bio::SeqFeature::Generic;
use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_RepeatSmasher_$VERSION";

local $| = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'   => \$p{CHROMOSOME},
      'EDITOR=s'       => \$p{EDITOR},
      'MEMO=s'         => \$p{MEMO},
      'SCALE=s'        => \$p{SCALE},
      'STARTPOS=i'     => \$p{STARTPOS},
      'STOPPOS=i'      => \$p{STOPPOS},
      'OUTPUT=s'       => \$p{OUTPUT},
			'help'           => \$p{HELP}
);
pod2usage(-verbose=>99, -sections=>"DESCRIPTION|ARGUMENTS") if ($p{HELP});

################################################################################
############################### SANITY CHECKING ################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n" unless ($p{CHROMOSOME});
my $chr    = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});
my $chrseq = $chr->sequence;
my $chrlen = length $chrseq;
my $GD = $chr->GD();
$p{STARTPOS} = $p{STARTPOS} || 1;
$p{STOPPOS} = $p{STOPPOS} || $chrlen;
$p{SCALE}  = $p{SCALE}  || 'chrom';
$p{OUTPUT} = $p{OUTPUT} || 'html';

unless ($p{EDITOR} && $p{MEMO})
{
  print "\n ERROR: Both an editor's id and a memo must be supplied.\n\n";
}

if ((! $p{STARTPOS} || ! $p{STOPPOS}) || $p{STOPPOS} <= $p{STARTPOS})
{
  die "\n ERROR: The start and stop coordinates do not parse.\n";
}

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $REPORT = {};
my $newchr = $chr->iterate();

################################################################################
############################### REPEAT SMASHING ################################
################################################################################
my @genes  = $newchr->db->features(
  -seqid      => $newchr->seq_id(),
  -start      => $p{STARTPOS},
  -end        => $p{STOPPOS},
  -range_type => 'contains',
  -type       => 'gene',
);

die "ERROR: There are no genes wholly contained in the range $p{STARTPOS}.."
  . "$p{STOPPOS}\n\n" unless (scalar(@genes));
 
foreach my $gene (@genes)
{
  my $genename = $gene->display_name;
  print "working on gene $genename\n";
  my $newseq = $GD->repeat_smash($gene);
  $newseq = $GD->complement($newseq, 1) if ($gene->strand == -1);
  $newchr->modify_feature(
    -feature      => $gene,
    -new_sequence => $newseq,
    -tags         => {Note => 'repeat smashed'},
    -preserve_overlapping_features => 1
  );
  $REPORT->{$genename} = "Repeatsmashed!"
}

#Summarize
print "\n\n";
print "Report:\n";
foreach my $featname (sort keys %{$REPORT})
{
  print "$featname : ", $REPORT->{$featname}, "\n";
}

################################################################################
################################### WRITING ####################################
################################################################################
if (scalar keys %{$REPORT})
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

  BS_RepeatSmasher.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility uses GeneDesign functions to minimize internal local alignments
    of gene features. The BioStudio twist is that it replaces internal features
    like PCRTags and codon swaps after the "repeat smashing" is done (GeneDesign
    has no respect for annotation). Given a chromosome and a range that contains
    at least one whole gene, this utility will generate a new version with all
    of the genes inside that range repeat smashed.

=head1 ARGUMENTS

Required arguments:

  -C,   --CHROMOSOME : The chromosome to be modified
  -E,   --EDITOR : The person responsible for the edits
  -M,   --MEMO : Justification for the edits
  -STA, --STARTPOS : The first base eligible for editing
  -STO, --STOPPOS : The last base eligible for editing

Optional arguments:

  -SC,  --SCALE : [genome, chrom (def)] Which version number to increment
  -OU,  --OUTPUT : html or txt
  -h,   --help : Display this message

=cut
