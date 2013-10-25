#!/usr/bin/perl

use Bio::BioStudio;
use Getopt::Long;
use Pod::Usage;
use CGI qw(:standard);

use strict;
use warnings;

my $VERSION = '2.00';
my $bsversion = "BS_PCRTagDumper_$VERSION";

local $| = 1;

my %p;
GetOptions (
      'CHROMOSOME=s'  => \$p{CHROMOSOME},
      'SCOPE=s'       => \$p{SCOPE},
      'STARTPOS=i'    => \$p{STARTPOS},
      'STOPPOS=i'     => \$p{STOPPOS},
      'OUTPUT=s'      => \$p{OUTPUT},
      'help'          => \$p{HELP}
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $BS = Bio::BioStudio->new();

$p{SCALE}  = $p{SCALE}  || 'chrom';
$p{OUTPUT} = $p{OUTPUT} || 'html';
$p{BSVERSION} = $bsversion;

die "BSERROR: No chromosome was named.\n"  unless ($p{CHROMOSOME});
my $chr    = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});

if ($p{SCOPE} && $p{SCOPE} eq "seg" &&
  ((! $p{STARTPOS} || ! $p{STOPPOS}) || $p{STOPPOS} <= $p{STARTPOS}))
{
  die "\n ERROR: The start and stop coordinates do not parse.\n";
}

################################################################################
################################# CONFIGURING ##################################
################################################################################
my $PCRPRODUCT = qr{\_amp(\d+)v(\d+)}msx;
my $GD = $chr->GD();
my $db = $chr->db();

my $dna = $chr->sequence();

my ($start, $end) = $p{SCOPE} eq "seg" 
                  ? ($p{STARTPOS}, $p{STOPPOS})
                  : (1, length($dna));
                  
my @amps = $db->features(
  -range_type => 'contains',
  -types      => 'PCR_product',
  -start      => $start,
  -end        => $end
);

print "PCR products in $p{CHROMOSOME} from $start to $end";
print <<"END",
All complete amplicons inside the range $start..$end are listed here in 5'-3'
orientation. Percent difference and melting temperatures are calculated as
averages.

CAUTION: if the annotated tag sequence does not match the current sequence, an
asterisk will appear by the amplicon number; see the bottom of the page for the
actual sequence.
END

print "   ORF      Amp# \t  5'-3' forward Wild Type   \t  ";
print "5'-3' reverse Wild Type   \t  5'-3' forward Synthetic   \t  ";
print "5'-3' reverse Synthetic   \tSize \t%Diff\tWT Tm\tSyn Tm";
print "\t  Links" if ($BS->{gbrowse} && $p{OUTPUT} eq "html");
print "\n\n";
my @takenotes;

foreach my $amplicon (sort {$a->start <=> $b->start} @amps)
{
  my $warning = 0;
  my $genename = $amplicon->Tag_ingene;
  my $number = 1;
  $number = $1 if ($amplicon =~ $PCRPRODUCT);
  my $intro = $amplicon->Tag_intro;
  my $wtsrc = $chr->species() . $chr->seq_id() . "_$intro";
  my @uptags = $db->features(-name => $amplicon->Tag_uptag);
  my @dntags = $db->features(-name => $amplicon->Tag_dntag);
  my ($uptag, $dntag)   = ($uptags[0], $dntags[0]);
  my ($fwtseq, $rwtseq) = ($uptag->Tag_wtseq, $dntag->Tag_wtseq);
  my ($fmdseq, $rmdseq) = ($uptag->Tag_newseq, $dntag->Tag_newseq);
  my ($fdiff, $rdiff)   = ($uptag->Tag_difference, $dntag->Tag_difference);
  my ($floc, $rloc)     = ($uptag->location(), $dntag->location());
 
  my $sp1 = q{ } x (28 - length($fwtseq));
  my $sp2 = q{ } x (28 - length($rwtseq));
  my $floclen = $floc->end - $floc->start + 1;
  my $rloclen = $rloc->end - $rloc->start + 1;
  my $f_check = substr($dna, $floc->start - $start, $floclen);
  my $r_check = substr($dna, $rloc->end - $rloclen - $start + 1, $rloclen);
 
  if ($rmdseq ne $r_check)
  {
    $warning++;
    my $warnmsg = "* The current sequence for the reverse synthetic primer ";
    $warnmsg = "of amplicon $number in $genename is $r_check\n";
    push @takenotes, $warnmsg;
  }
  if ($fmdseq ne $f_check)
  {
    $warning++;
    my $warnmsg = "* The current sequence for the forward synthetic primer ";
    $warnmsg = "of amplicon $number in $genename is $f_check\n";
    push @takenotes, $warnmsg;
  }
  my $wtstart = $uptag->Tag_wtpos;
  my $wtend = $wtstart + ($amplicon->end - $amplicon->start + 1) - 1;
  my $disclaimer = $warning > 0  ?  q{*}  :  q{};
  my $diff = int(($rdiff + $fdiff) / 2 + .5);
  my $size = $amplicon->stop - $amplicon->start + 1;
  my $wt_Tm = int( ( $GD->melt($fwtseq) + $GD->melt($rwtseq) ) / 2 + 0.5);
  my $md_Tm = int( ( $GD->melt($fmdseq) + $GD->melt($rmdseq) ) / 2 + 0.5);
  print $genename, q{      }, $number, q{ }, $disclaimer, "\t";
  print $fwtseq, "$sp1\t", $GD->complement($rwtseq, 1), "$sp2\t";
  print $fmdseq, "$sp1\t", $GD->complement($rmdseq, 1), "$sp2\t";
  print $size, "  \t  ", $diff, " \t ", $wt_Tm, " \t    ", $md_Tm;
  if ($BS->{gbrowse} && $p{OUTPUT} eq "html")
  {
    my $wthref  = "http://$BS->{this_server}/cgi-bin/gb2/gbrowse/$wtsrc/?start";
       $wthref .= "=$wtstart;stop=$wtend;ref=$p{SEQID};";
    my $wtlink = a({href => $wthref, -target => '_blank'}, 'wt');
    my $synhref = link_to_feature($chr, $amplicon);
    my $synlink = a({href => $synhref, -target => '_blank'}, 'syn');
    print " \t  ", $wtlink, q{ }, $synlink;
  }
  print "\n";
}

print "\n\n";
print "$_\n" foreach (@takenotes);
print "\n";

exit;

__END__

=head1 NAME

  BS_PCRTagDumper.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility creates a list of PCR Tags from a chromosome.  It will alert when
   the sequence for a synthetic tag is not what was expected; this usually means
   that a subsequent edit modified the sequence without considerately updating
   the tags annotation.

=head1 ARGUMENTS

Required arguments:

  -C, --CHROMOSOME : The chromosome to be parsed

Optional arguments:

  -SC,  --SCOPE : [seg, chrom (def)] How much sequence to parse for tags.
                  seg requires STARTPOS and STOPPOS.
  -STA, --STARTPOS : The first base for parsing; ignored unless SCOPE = seg
  -STO, --STOPPOS  : The last base for parsing; ignored unless SCOPE = seg
  -OU,  --OUTPUT   : [html, txt (def)] Format of the output
  -h,   --help : Display this message
 
=cut
