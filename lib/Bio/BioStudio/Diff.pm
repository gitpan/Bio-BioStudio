#
# BioStudio module for comparing normalized SeqFeatures
#

=head1 NAME

Bio::BioStudio::Diff

=head1 VERSION

Version 1.05

=head1 DESCRIPTION

BioStudio functions for comparing Bio::DB::SeqFeatures. These are used to
"diff" two chromosomes, to figure out changes in annotation and sequence.

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::BioStudio::Diff;

use Exporter;
use Bio::GeneDesign::Codons qw(translate);
use Bio::GeneDesign::Basic qw(compare_sequences);
use Bio::BioStudio::BLAST qw(bl2seq_blastn bl2seq_blastp);
use Text::Diff;

use strict;
use warnings;

our $VERSION = '1.05';
$| = 1;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
  compare_feature_sources
  compare_feature_types
  compare_feature_lengths
  compare_feature_orientations
  compare_feature_scores
  compare_feature_phases
  compare_feature_attributes
  compare_feature_annotations
  compare_feature_sequences
  compare_feature_translations
  compare_comments
);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

my %CODES = (
      1  => "deleted feature",
      2  => "added feature",
      3  => "lost subfeature",
      4  => "gained subfeature",
      5  => "lost sequence",
      6  => "gained sequence",
      7  => "change in translation",
      8  => "change in sequence",
      9  => "lost attributes",
      10 => "gained attributes",
      11 => "changed annotation",
      12 => "changed subfeature"
);

=head2 compare_feature_sources

=cut

sub compare_feature_sources
{
  my ($feat1, $feat2) = @_;
  return undef if ($feat1->source() eq $feat2->source());
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $feat1->source_tag(),
            -newatt => $feat2->source_tag(),
            -comment => "source",
            -code => 11
  );
}

=head2 compare_feature_types

=cut

sub compare_feature_types
{
  my ($feat1, $feat2) = @_;
  return undef if ($feat1->primary_tag() eq $feat2->primary_tag());
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $feat1->primary_tag(),
            -newatt => $feat2->primary_tag(),
            -comment => "primary_tag",
            -code => 11
  );
}

=head2 compare_feature_lengths

=cut

sub compare_feature_lengths
{
  my ($feat1, $feat2) = @_;
  my $len1 = $feat1->end - $feat1->start + 1;
  my $len2 = $feat2->end - $feat2->start + 1;
  my $difference = abs($len1 - $len2);
  return undef if ($difference == 0);
  my ($code, $baseloss, $basegain) = (0, 0, 0);
  if ($len1 > $len2)
  {
    $code = 5;
    $baseloss = $difference;
  }
  else
  {
    $code = 6;
    $basegain = $difference;
  }
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $len1,
            -newatt => $len2,
            -code => $code,
            -baseloss => $baseloss,
            -basegain => $basegain
  );
}

=head2 compare_feature_orientations

=cut

sub compare_feature_orientations
{
  my ($feat1, $feat2) = @_;
  return undef if ($feat1->strand() eq $feat2->strand());
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $feat1->strand(),
            -newatt => $feat2->strand(),
            -comment => "strand",
            -code => 11
  );
}

=head2 compare_feature_scores

=cut

sub compare_feature_scores
{
  my ($feat1, $feat2) = @_;
  my $f1score = $feat1->score ? $feat1->score : 0;
  my $f2score = $feat2->score ? $feat2->score : 0;
  return undef if ($f1score eq $f2score);
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $f1score,
            -newatt => $f2score,
            -comment => "score",
            -code => 11
  );
}

=head2 compare_feature_phases

=cut

sub compare_feature_phases
{
  my ($feat1, $feat2) = @_;
  my $f1phase = $feat1->phase ? $feat1->phase : 0;
  my $f2phase = $feat2->phase ? $feat2->phase : 0;
  return undef if ($f1phase eq $f2phase);
  return Bio::BioStudio::Diff::Difference->new(
            -oldfeat => $feat1,
            -newfeat => $feat2,
            -oldatt => $f1phase,
            -newatt => $f2phase,
            -comment => "phase",
            -code => 11
  );
}

=head2 compare_feature_attributes

=cut

sub compare_feature_attributes
{
  my ($feat1, $feat2) = @_;
  my @attChanges;

  ##Check to see if source has changed
  my $sourcechange = compare_feature_sources($feat1, $feat2);
  push @attChanges, $sourcechange if ($sourcechange);

  ##Check to see if type has changed
  my $typechange = compare_feature_types($feat1, $feat2);
  push @attChanges, $typechange if ($typechange);

  ##Check to see if length has changed
  my $lengthchange = compare_feature_lengths($feat1, $feat2);
  push @attChanges, $lengthchange if ($lengthchange);

  ##Check to see if score has changed
  my $scorechange = compare_feature_scores($feat1, $feat2);
  push @attChanges, $scorechange if ($scorechange);

  ##Check to see if orientation has changed
  my $orientchange = compare_feature_orientations($feat1, $feat2);
  push @attChanges, $orientchange if ($orientchange);

  ##Check to see if phase has changed
  my $phasechange = compare_feature_phases($feat1, $feat2);
  push @attChanges, $phasechange if ($phasechange);

  if (scalar (@attChanges))
  {
    return \@attChanges;
  }
  else
  {
    return undef;
  }
}

=head2 compare_feature_annotations

=cut

sub compare_feature_annotations
{
  my ($feat1, $feat2) = @_;
  my @annChanges;

  my @tags1 = $feat1->get_all_tags();
  my @tags2 = $feat2->get_all_tags();
  foreach my $tag (grep {$feat2->has_tag($_)} @tags1)
  {
    my %vals1 = map {$_ => 1} $feat1->get_tag_values($tag);
    my %vals2 = map {$_ => 1} $feat2->get_tag_values($tag);
    foreach my $val (grep {! exists $vals2{$_}} keys %vals1)
    {
      push @annChanges, Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -oldatt => "$tag = $val",
                -code => 9
      );
    }
    foreach my $val (grep {! exists $vals1{$_}} keys %vals2)
    {
      push @annChanges, Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -newatt => "$tag = $val",
                -code => 10
      );
    }
  }
  foreach my $tag (grep {! $feat2->has_tag($_)} @tags1)
  {
    push @annChanges, Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -oldatt => "$tag = " . join(", ", $feat1->get_tag_values($tag)),
                -code => 9
    );
  }
  foreach my $tag (grep {! $feat1->has_tag($_)} @tags2)
  {
    push @annChanges, Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -oldatt => "$tag = " . join(", ", $feat2->get_tag_values($tag)),
                -code => 10
    );
  }

  if (scalar (@annChanges))
  {
    return \@annChanges;
  }
  else
  {
    return undef;
  }
}

=head2 compare_feature_sequences

=cut

sub compare_feature_sequences
{
  my ($feat1, $feat2, $blast) = @_;

  my ($old, $var) = ($feat1->seq->seq, $feat2->seq->seq);
  return undef if ($old eq $var);
  my ($aligns, $basechange) = (undef, 0);
  if ($blast)
  {
    if ($feat2->has_tag("newseq"))
    {
      my $cthsh = compare_sequences($feat2->Tag_wtseq, $feat2->Tag_newseq);
      $basechange = $cthsh->{D};
    }
    elsif ($feat1->primary_tag() eq "chromosome")
    {
      ##WHOLE CHROMOSOME ALIGNMENT
    }
    else
    {
      my $alnz = bl2seq_blastn($feat1, $feat2, $blast);
      if ($alnz)
      {
        $aligns = $alnz;
        my $mismatchcount = 0;
        my $hitlen = 0;
        while (my $hit = $aligns->next_hit())
        {
          if (scalar($hit->hsps()))
          {
            $hitlen += $hit->length();
            $mismatchcount += ($hit->length() - $hit->matches('id'));
          }
        }
        $mismatchcount += (length($var) - $hitlen);
        $basechange = $mismatchcount;
      }
    }
  }
  return Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -oldatt => $old,
                -newatt => $var,
                -code => 8,
                -aligns => $aligns,
                -basechange => $basechange
  );
}

=head2 compare_feature_translations

=cut

sub compare_feature_translations
{
  my ($feat1, $feat2, $CODON_TABLE, $blast) = @_;
  my $orf1 = $feat1->seq->seq;
	my $orf2 = $feat2->seq->seq;
  my $f1phase = $feat1->phase ? $feat1->phase : 0;
  my $f2phase = $feat2->phase ? $feat2->phase : 0;
	my $old = translate($orf1, $f1phase + 1, $CODON_TABLE);
	my $var = translate($orf2, $f2phase + 1, $CODON_TABLE);
  return undef if ($old eq $var);

  if ((!$old && length($orf1) >=3) || (!$var && length($orf2) >=3))
  {
    print "original $feat1 " if (!$old);
    print "variant $feat2 " if (!$var);
    print "has no translation, strangely!!\n";
    return undef;
  }
  my $aligns = undef;
  if ($blast)
  {
    $aligns = bl2seq_blastp($feat1, $feat2, $CODON_TABLE, $blast);
  }
  return Bio::BioStudio::Diff::Difference->new(
                -oldfeat => $feat1,
                -newfeat => $feat2,
                -oldatt => $old,
                -newatt => $var,
                -code => 7,
                -aligns => $aligns
  );
}

=head2 compare_comments

=cut

sub compare_comments
{
	my ($ref1, $ref2) = @_;
	my $diff = diff $ref1, $ref2;
	return $diff;
}

1;

__END__

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, BioStudio developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.

* Neither the name of the Johns Hopkins nor the names of the developers may be 
used to endorse or promote products derived from this software without specific 
prior written permission.

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
