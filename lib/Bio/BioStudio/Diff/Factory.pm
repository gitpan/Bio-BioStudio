#
# BioStudio module for comparing normalized SeqFeatures
#

=head1 NAME

Bio::BioStudio::Diff::Factory

=head1 VERSION

Version 1.05

=head1 DESCRIPTION

An object that allows the comparison of  two Bio::DB::SeqFeature::Store objects.

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::BioStudio::Diff::Factory;

use Bio::BioStudio::Basic qw(flatten_subfeats);
use Bio::BioStudio::Diff qw(:all);
use Bio::BioStudio::Diff::Difference;
use Bio::GeneDesign::Basic qw(compare_sequences);
use List::Util qw(first);
use strict;

$| = 1;

use base qw(Bio::Root::Root);
=head1 CONSTRUCTOR METHOD

=head2 new

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($olddb, $newdb, $checktrx, $aligntrx, $checkseq, $alignseq, $codon_table,
    $blast_factory) =
     $self->_rearrange([qw(OLDDB
                           NEWDB
                           CHECKTRANSLATION
                           ALIGNTRANSLATION
                           CHECKSEQUENCE
                           ALIGNSEQUENCE
                           CODON_TABLE
                           BLAST_FACTORY)], @args);
                           
   $self->throw("\"Old\" database not supplied") unless ($olddb);
   $self->olddb($olddb);

   $self->throw("\"New\" database not supplied") unless ($newdb);
   $self->newdb($newdb);
   
   $self->throw("Align flags set but no blast_factory provided")
      if (($aligntrx || $alignseq) && ! $blast_factory);
      
   $self->throw("Translation flags set but no codon_table provided")
      if (($aligntrx || $checktrx) && ! $codon_table);
                  
   $checktrx && $self->checktrx($checktrx);
   $aligntrx && $self->aligntrx($aligntrx);
   $checkseq && $self->checkseq($checkseq);
   $alignseq && $self->alignseq($alignseq);
   $codon_table && $self->codon_table($codon_table);
   $blast_factory && $self->blast_factory($blast_factory);
   
   my (%feats1, %feats2) = ((), ());
   #PARENTS ONLY! Subfeatures will be dealt with downstream
   my $iterator1 = $olddb->get_seq_stream();
   while (my $feature = $iterator1->next_seq)
   {
     next if ($feature->has_tag("parent_id"));
     $feats1{$feature->Tag_load_id} = $feature;
   }

   my $iterator2 = $newdb->get_seq_stream();
   while (my $feature = $iterator2->next_seq)
   {
     next if ($feature->has_tag("parent_id"));
     $feats2{$feature->Tag_load_id} = $feature;
   }

   my @onlyolds = grep {! exists $feats2{$_} } keys %feats1;
   $self->only_old_features(\@onlyolds);
    
   my @onlynews = grep {! exists $feats1{$_} } keys %feats2;
   $self->only_new_features(\@onlynews);   
   
   my @commons = grep  {  exists $feats2{$_} } keys %feats1; 
   $self->common_features(\@commons);

   return $self;
}

=head1 COMPARISON METHODS

=head2 compare_dbs

=cut

sub compare_dbs
{
  my ($fact) = @_;
  
  my @D;
  
  #Common Features
  foreach my $featid (@{$fact->common_features()})
  {
    my $oldfeat = first {1} $fact->olddb->get_features_by_name($featid);
    my $newfeat = first {1} $fact->newdb->get_features_by_name($featid);
    $fact->throw("Can't find $featid in olddb!") unless ($oldfeat);
    $fact->throw("Can't find $featid in newdb!") unless ($newfeat);
    push @D, $fact->compare_features($oldfeat, $newfeat);
  }

  #Deleted Features
  foreach my $featid (@{$fact->only_old_features()})
  {
    my $feat = first {1} $fact->olddb->get_features_by_name($featid);
    $fact->throw("Can't find $featid in olddb!") unless ($feat);
    my $baseloss = $feat->end - $feat->start + 1;
    push @D, Bio::BioStudio::Diff::Difference->new(
        -oldfeat => $feat,
        -code => 1,
        -baseloss => $baseloss
    );
  }

  #Inserted Features
  foreach my $featid (@{$fact->only_new_features()})
  {
    my $feat = first {1} $fact->newdb->get_features_by_name($featid);
    $fact->throw("Can't find $featid in newdb!") unless ($feat);
    my $basechange = 0;
    my $basegain = $feat->end - $feat->start + 1;
    if ($feat->has_tag("newseq"))
    {
      my $cthsh = compare_sequences($feat->Tag_wtseq, $feat->Tag_newseq);
      $basechange = $cthsh->{D};
    }
    push @D, Bio::BioStudio::Diff::Difference->new(
        -newfeat => $feat,
        -code => 2,
        -basegain => $basegain,
        -basechange => $basechange
    );
  }
  
  return @D; 
}

=head2 compare_features

=cut

sub compare_features
{
  my ($fact, $feat1, $feat2) = @_;
  my @Changes = ();

  ##Check for subfeatures and compare
  if (scalar($feat1->get_SeqFeatures) || scalar($feat2->get_SeqFeatures))
  {
    my @allsubs1 = flatten_subfeats($feat1);
    my @allsubs2 = flatten_subfeats($feat2);
    my %types = map {$_->primary_tag() => 1} @allsubs1, @allsubs2;
    foreach my $type (keys %types)
    {
      my @subtype1s = grep {$_->primary_tag eq $type} @allsubs1;
      my %subs1 = map {$_->Tag_load_id => $_} @subtype1s;
      my @subtype2s = grep {$_->primary_tag eq $type} @allsubs2;
      my %subs2 = map {$_->Tag_load_id => $_} @subtype2s;

      #Compare all sub features of type $type
      foreach my $sub1id (grep {exists $subs2{$_}} keys %subs1)
      {
        my $subfeat1 = delete $subs1{$sub1id};
        my $subfeat2 = delete $subs2{$sub1id};
        foreach my $sobj ($fact->compare_features($subfeat1, $subfeat2))
        {
          if ($sobj->code() >= 4)
          {
            push @Changes, Bio::BioStudio::Diff::Difference->new(
                      -oldfeat => $feat1,
                      -newfeat => $feat2,
                      -oldsubfeat => $subfeat1,
                      -newsubfeat => $subfeat2,
                      -code => 12,
                      -subdiff => $sobj);
          }
        }
      }
      #Deleted subfeatures of type $type
      foreach (keys %subs1)
      {
        my $subfeat = $subs1{$_};
        my $baseloss = $subfeat->end - $subfeat->start + 1;

        push @Changes, Bio::BioStudio::Diff::Difference->new(
                  -oldfeat => $feat1,
                  -oldsubfeat => $subfeat,
                  -newfeat => $feat2,
                  -baseloss => $baseloss,
                  -code => 3);
      }
      #Inserted subfeatures of type $type
      foreach (keys %subs2)
      {
        my $subfeat = $subs2{$_};
        my $basechange = 0;
        if ($subfeat->has_tag("newseq"))
        {
          my $ct = compare_sequences($subfeat->Tag_wtseq, $subfeat->Tag_newseq);
          $basechange = $ct->{D};
        }
        push @Changes, Bio::BioStudio::Diff::Difference->new(
                  -oldfeat => $feat1,
                  -newfeat => $feat2,
                  -newsubfeat => $subfeat,
                  -code => 4,
                  -basechange => $basechange);
      }
    }
  }

  ##Check to see if attributes have changed
  my $attarray = compare_feature_attributes($feat1, $feat2);
  push @Changes, @$attarray if ($attarray);

  ##Check to see if annotations have changed
  my $annarray = compare_feature_annotations($feat1, $feat2);
  push @Changes, @$annarray if ($annarray);

  ##If the translation switch is on, see if translation has changed (CDS ONLY)
  if ($feat1->primary_tag() eq "CDS" && ($fact->checktrx || $fact->aligntrx))
  {
    my $ch = compare_feature_translations($feat1, $feat2, 
                                      $fact->codon_table, $fact->blast_factory);
    push @Changes, $ch if ($ch);
  }

  ##If check sequence is on, see if the feature sequence has changed
  if ($fact->checkseq || $fact->alignseq)
  {
    my $ch = compare_feature_sequences($feat1, $feat2, $fact->blast_factory);
    push @Changes, $ch if ($ch);
  }
  return @Changes;  
}

=head1 ACCESSOR METHODS

=head2 olddb

=cut

sub olddb
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
    unless ($value->isa("Bio::DB::SeqFeature::Store"))
    {
      $fact->throw("olddb value is not a Bio::DB::SeqFeature::Store object");
    }
	  $fact->{'olddb'} = $value;
  }
  return $fact->{'olddb'};
}

=head2 newdb

=cut

sub newdb
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
    unless ($value->isa("Bio::DB::SeqFeature::Store"))
    {
      $fact->throw("newdb value is not a Bio::DB::SeqFeature::Store object");
    }
	  $fact->{'newdb'} = $value;
  }
  return $fact->{'newdb'};
}

=head2 checktrx

=cut

sub checktrx
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'checktrx'} = $value;
  }
  return $fact->{'checktrx'};
}

=head2 aligntrx

=cut

sub aligntrx
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'aligntrx'} = $value;
  }
  return $fact->{'aligntrx'};
}

=head2 alignseq

=cut

sub alignseq
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'alignseq'} = $value;
  }
  return $fact->{'alignseq'};
}

=head2 checkseq

=cut

sub checkseq
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'checkseq'} = $value;
  }
  return $fact->{'checkseq'};
}


=head2 common_features

=cut
sub common_features
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'common_features'} = $value;
  }
  return $fact->{'common_features'};
}

=head2 only_old_features

=cut

sub only_old_features
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'only_old_features'} = $value;
  }
  return $fact->{'only_old_features'};
}

=head2 only_new_features

=cut

sub only_new_features
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'only_new_features'} = $value;
  }
  return $fact->{'only_new_features'};
}

=head2 codon_table

=cut

sub codon_table
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
	  $fact->{'codon_table'} = $value;
  }
  return $fact->{'codon_table'};
}

=head2 blast_factory

=cut

sub blast_factory
{
  my ($fact, $value) = @_;
  if (defined $value)
  {
    unless ($value->isa("Bio::Tools::Run::StandAloneBlastPlus"))
    {
      $fact->throw("blast_factory value is not a " .
                   "Bio::Tools::Run::StandAloneBlastPlus object");
    }
	  $fact->{'blast_factory'} = $value;
  }
  return $fact->{'blast_factory'};
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
