#!/usr/bin/env perl

use Bio::BioStudio;
use Bio::BioStudio::RestrictionEnzyme::Store;
use Bio::BioStudio::RestrictionEnzyme::Seek qw(:BS);
use English qw(-no_match_vars);
use Getopt::Long;
use Pod::Usage;
use Carp;

use strict;
use warnings;

my $VERSION = '2.00';

local $OUTPUT_AUTOFLUSH = 1;

my %p;
GetOptions (
			'CHROMOSOME=s'  => \$p{CHROMOSOME},
			'RESET=s'       => \$p{RESET},
			'CULLINC=i'     => \$p{CULLINC},
			'CHUNKLENMIN=i' => \$p{CHUNKLENMIN},
			'CHUNKLENMAX=i' => \$p{CHUNKLENMAX},
			'REDOREDB'      => \$p{REDOREDB},
			'help'          => \$p{HELP}
);
if ($p{HELP})
{
  pod2usage(
    -verbose => 99,
    -sections=>"NAME|VERSION|DESCRIPTION|ARGUMENTS|USAGE");
}

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $BS = Bio::BioStudio->new();

die "BSERROR: No chromosome was named.\n"  unless ($p{CHROMOSOME});
my $chr = $BS->set_chromosome(-chromosome => $p{CHROMOSOME});

$p{CULLINC} = $p{CULLINC} || 2000;
$p{REDOREDB} = $p{REDOREDB} || 0;
$p{RESET} = $p{RESET} || "nonpal_short";
$p{CHUNKLENMIN} = $p{CHUNKLENMIN} || 5000;
$p{CHUNKLENMAX} = $p{CHUNKLENMAX} || 10000;
  
################################################################################
################################# CONFIGURING ##################################
################################################################################
my $GD = $chr->GD;
my $db = $chr->database();
my $chrseq = $chr->sequence();
my $chrlen = length($chrseq);

print "Loading restriction enzymes... ";
my $RES = $GD->set_restriction_enzymes(-enzyme_set => $p{RESET});
my @enzes = values %{$RES};
my @ordlen = sort {$b->len <=> $a->len} @enzes;
my $maxenzlen = $ordlen[0]->len;
print scalar(@enzes), " restriction enzymes loaded (up to $maxenzlen bp).\n";

## find genes, find intergenic regions and add to database (unless they exist)
## Load the database for the original chromosome
print "Gene finding... ";
my @genes  = $db->get_features_by_type('gene');
my $mask = Bio::BioStudio::Mask->new(-sequence => $chr->seqobj);
$mask->add_to_mask(\@genes);

my @exons   = $db->get_features_by_type('CDS');
print scalar(@exons), " exons found\n";

################################################################################
################################# DATABASING  ##################################
################################################################################
my $DBLIST = $BS->db_list();
my $redbname = $chr->name() . "_RED";

if ($p{REDOREDB} || ! exists $DBLIST->{$redbname})
{
  ## Create and populate the suffix tree of restriction enzyme recognition sites
  my $tree = $GD->build_suffix_tree(-input => \@enzes, -peptide => 1);
  print "Suffix tree created.\n\n";

  my $filename = $BS->tmp_path() . $redbname . ".out";
  chmod 0777, $filename;
  open my $OUT, '>', $filename || croak ("Can't write to $filename: $OS_ERROR");
  my ($x, $y) = (0, 0);
  
  ## Process exons
  foreach my $exon (@exons)
  {
    print "\tWorking on $exon...\n";
    my $start = $exon->start;
    my $hits = find_enzymes_in_CDS($chr, $tree, $exon);
    foreach my $enz (@{$hits})
    {
      print $OUT $enz->line_report(q{.}, "\n");
      $y++;
    }
    $x++;
  }
  
  ## Process igens, catch lost enzymes at the interfaces
  my @igens = $chr->make_intergenic_features();
  push @igens, $db->get_features_by_type('intron');
  foreach my $igen (@igens)
  {
    print "\tWorking on $igen...\n";
    my $hits = find_enzymes_in_igen($chr, $igen, $maxenzlen);
    foreach my $enz (@{$hits})
    {
      print $OUT $enz->line_report(q{.}, "\n");
      $y++;
    }
    $x++;
  }

  close $OUT;
  print "Processed $x features, found $y restriction enzyme sites.\n\n";

  ## Fastload the MySQL database
  print "Loading database...\n";
  my $newREDB = Bio::BioStudio::RestrictionEnzyme::Store->new
  (
    -name => $redbname,
    -enzyme_definitions => $RES,
    -create => 1,
    -file => $filename
  );
  $newREDB->dumpfile($filename);
  $newREDB->load();
}

################################################################################
################################### CULLING ####################################
################################################################################
my $REDB = Bio::BioStudio::RestrictionEnzyme::Store->new(
  -name               => $redbname,
  -enzyme_definitions => $RES
);

## Gather genes and check for essentiality.

my %ESSENTIALS;
foreach my $CDS (@exons)
{
  my @parents = $db->get_feature_by_name($CDS->Tag_parent_id);
  my $parent = $parents[0];
  if ($parent && $parent->has_tag('essential_status'))
  {
    $ESSENTIALS{$CDS->Tag_load_id} = $parent->Tag_essential_status;
  }
  else
  {
    $ESSENTIALS{$CDS->Tag_load_id} = 'Nonessential';
  }
}

##Begin culling
# Given a position along the chromosome, find sites that will never be eligible
# for landmark status. If the site is "p", remove it from the db.  If the site
# is "e" or "i", mark it ineligible in the db so it can be skipped in future
# analyses.
print "Screening database...\n";
my $position = 1;
my $drcount = 0;
my $igcount = 0;
my $ticker = 0;
while ($position <= $chrlen)
{
  my (@culls, @ignores) = ((), ());
  print "@ $position bp\n" if ($ticker % 10 == 0);
  my $lcpos = $position - $p{CULLINC};
  $lcpos = 1 if ($lcpos < 0);
  my $rcpos = $position + $p{CULLINC};
  my $lef = $position - ($p{CHUNKLENMIN} + 2000);
  $lef = 1 if ($lef <= 0);
  my $rig = $position + ($p{CHUNKLENMIN} + 2000);
  #print "\tlooking for REs between $lcpos ($lef) and ($rig) $rcpos...\n";
  my $pool = $REDB->search(-left => $lef, -right => $rig);
  #print "\tpooled ", scalar @{$pool}, " enzymes...\n";
  my @candidates = grep {$_->start >= $lcpos && $_->end <= $rcpos} @{$pool};
  #print "\t\tvetting ", scalar @candidates, " enzymes...\n";
  foreach my $re (@candidates)
  {
    next if ($re->eligible && $re->eligible eq "no");
    my $size = $re->end - $re->start + 1;
    my $maskbit = $mask->count_features($re->start - 1);
    if ($maskbit > 1) #DROP: exonic in exon overlap
    {
      push @culls, $re->dbid if ($re->presence eq "potential");
      push @ignores, $re->dbid if ($re->presence ne "potential");
      next;
    }
    my @buds = grep {abs($_->start - $re->start) <= $p{CHUNKLENMIN}
                  && $_->presence ne "potential" && $_->id eq $re->id
                  && $_->name ne $re->name} @{$pool};
    my @igenics = grep {$_->presence eq "intergenic"} @buds;
    if (scalar(@igenics)) #DROP: too many intergenics around
    {
      push @culls, $re->dbid if ($re->presence eq "potential");
      push @ignores, $re->dbid if ($re->presence ne "potential");
      next;
    }
    my $gene = $re->presence ne "intergenic"  ? $re->featureid  : q{};
    my @exonics = grep {$_->presence eq "existing"} @buds;
    my $ess_flag = exists($ESSENTIALS{$gene}) && $ESSENTIALS{$gene} eq "Essential"   ? 1 : 0;
    my $fgw_flag = exists($ESSENTIALS{$gene}) && $ESSENTIALS{$gene} eq "fast_growth" ? 1 : 0;
    my $lap_flag = 0;
    my @movers;
    foreach my $ex (@exonics)
    {
      my $egene = $ex->featureid;
      my $emaskbit = $mask->count_features($ex->start - 1);
      $lap_flag++ if ($emaskbit != 1);
      $ess_flag++ if (exists ($ESSENTIALS{$gene}) && $ESSENTIALS{$egene} eq "Essential");
      $fgw_flag++ if (exists ($ESSENTIALS{$gene}) && $ESSENTIALS{$egene} eq "fast_growth");
      push @movers, $ex->name;
    }
    if ($lap_flag != 0) #DROP: modification in exon overlap
    {
      push @culls, $re->dbid if ($re->presence eq "potential");
      push @ignores, $re->dbid if ($re->presence ne "potential");
      next;
    }
    #my $distance = abs($position - $re->{START});
    #my $score = $distance <= $p{CULLINC} ?	0	:	0.002  * $distance - 4;
    # Plus the log of the price per unit
    my $score   = log($re->score);
    # Plus one tenth point for each orf modified
    $score   += 0.1 * scalar(@movers);
    # Plus one half point for each fast growth orf modified
    $score   += 0.5 * $fgw_flag;
    # Plus one point for each essential orf modified
    $score   += 1.0 * $ess_flag;

    if ($score > 1) #DROP: score too high (>1)
    {
      push @culls, $re->dbid if ($re->presence eq "potential");
      push @ignores, $re->dbid if ($re->presence ne "potential");
      next;
    }
  }
  
  if (scalar(@culls))
  {
    $REDB->cull(\@culls);
    $drcount += scalar(@culls);
  }
  
  if (scalar(@ignores))
  {
    $REDB->screen(\@ignores);
    $igcount += scalar(@ignores);
  }
  
  $position += $p{CULLINC};
  print "\t$drcount culls, $igcount ignores so far...\n" if ($ticker % 10 == 0);
  $ticker++;
}
print "dropped $drcount sites and marked $igcount sites ignorable.\n\n";

exit;

__END__

=head1 NAME

  BS_GlobalREMarkup.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  This utility creates an exhaustive database of restriction enzyme recognition
   sites along a chromsome, both existing and potential.
  Every intergenic sequence is parsed for existing recognition sites. Those that
   are found are marked (i)mmutable; A suffix tree is created from all possible
   6-frame translations of restriction enzyme recognition sites, such that each
   node in the tree is an amino acid string that may be reverse translated to be
   a recognition site. Every exonic sequence in the chromosome is then searched
   with the suffix tree both for (e)xisting recognition sites and for sites
   where a (p)otential recognition site could be introduced without changing the
   protein sequence of the gene. As long as they occur within protein coding
   genes, existing and potential recognition sites may be manipulated to yield
   any of several different overhangs, all of which are computed by the
   algorithm.
  A score is assigned to every restriction enzyme site. The score is a function
   of the log of the enzyme's price per unit, plus one tenth for each orf that
   must be modified to make the enzyme unique within the range of the chunk
   size, plus one half for each fast growth orf modified, plus one for each
   essential orf modified.
  Every extant site is indexed so that later in the design process, when
    potential sites are considered, the number of sites that must be modified is
    a known contribution to the cost. However, all potential sites that could
    not be made unique under any circumstances are culled from the database,
    with a window size that defaults to 2000 base pairs but may be defined by
    the user.

=head1 ARGUMENTS

Required arguments:

  -CHR, --CHROMOSOME : The chromosome to be modified

Optional arguments:

  --REDOREDB : Whether or not to reform the RE database from scratch
  -CU,  --CULLINC : The width of the culling window
  --RESET : Which list of restriction enzymes to use (default nonpal)
  --CHUNKLENMIN : The minimum size of chunks to be designed (default 5000)
  --CHUNKLENMAX : The maximum size of chunks to be designed (default 10000)
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