#
# BioStudio module for sequence segmentation
#
# POD documentation - main docs before the code

=head1 NAME

Bio::BioStudio::RestrictionEnzyme

=head1 VERSION

Version 1.04

=head1 DESCRIPTION

BioStudio object that represents a restriction enzyme - inherits from 
Bio::GeneDesign::RestrictionEnzyme and adds attributes for feature annotation
awareness

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::BioStudio::RestrictionEnzyme;

use Switch;

use strict;

use base qw(Bio::GeneDesign::RestrictionEnzyme);

my $VERSION = 1.04;

=head1 CONSTRUCTORS

=head2 new

 Title   : new
 Function:
 Returns : a new Bio::BioStudio::RestrictionEnzyme object
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  bless $self, $class;
  
  my ($name, $presence, $eligible, $end, $feature, $featureid, $overhangs, 
      $strand, $peptide, $offset, $dbid, $phang) =
     $self->_rearrange([qw(NAME
                           PRESENCE
                           IGNORE
                           END
                           FEATURE
                           FEATUREID
                           OVERHANGS
                           STRAND
                           PEPTIDE
                           OFFSET
                           DBID
                           PHANG)], @args);

  $self->throw("No name defined") unless ($name);
  $self->{'name'} = $name;

  $self->throw("No presence defined") unless ($presence);
  switch ($presence)
  {
    case "p"          {$self->{'presence'} = "potential"}
    case "potential"  {$self->{'presence'} = "potential"}
    case "i"          {$self->{'presence'} = "intergenic"}
    case "intergenic" {$self->{'presence'} = "intergenic"}
    case "e"          {$self->{'presence'} = "existing"}
    case "existing"   {$self->{'presence'} = "existing"}
    case "a"          {$self->{'presence'} = "appended"}
    case "appended"   {$self->{'presence'} = "appended"}
    else              {$self->throw("Cannot determine presence from $presence")} 
  };

  $eligible && $self->eligible($eligible);

  $end && $self->end($end);

  if ($feature)
  {
    $self->feature($feature);
    $self->featureid($feature->id);
  }
  elsif ($featureid)
  {
    $self->featureid($featureid);
  }
  $overhangs && $self->overhangs($overhangs);

  $strand && $self->strand($strand);
  
  $peptide && $self->peptide($peptide);
  
  $offset && $self->offset($offset);
  
  $dbid && $self->dbid($dbid);
  
  $phang && $self->phang($phang);

  return $self;
}

=head1 FUNCTIONS

=head2 dump

=cut

sub dump
{
  my ($self, $fieldterm, $lineterm) = @_;
  my $pep = $self->peptide ?  $self->peptide  : "NULL";
  my $line = $self->name . $fieldterm;
  $line .= $self->presence . $fieldterm;
  $line .= $self->start . $fieldterm;
  $line .= $self->end . $fieldterm;
  $line .= $self->id . $fieldterm;
  $line .= $self->feature->id . $fieldterm;
  $line .= $self->peptide . $fieldterm;
  $line .= join(",", keys %{$self->overhangs}) . $fieldterm;
  $line .= $self->strand . $fieldterm;
  $line .= $self->offset . $lineterm;
  return $line;
}

=head1 ACCESSORS

=head2 name

=cut

sub name
{
  my ($self) = @_;
  return $self->{'name'};
}

=head2 presence

=cut

sub presence
{
  my ($self) = @_;
  return $self->{'presence'};
}

=head2 end

=cut

sub end
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'end'} = $value;
  }
  return $self->{'end'};
}

=head2 eligible

=cut

sub eligible
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'eligible'} = $value;
  }
  return $self->{'eligible'};
}

=head2 feature

=cut

sub feature
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->throw("object of class " . ref($value) . " does not implement ".
		    "Bio::DB::SeqFeature.")
		  unless $value->isa("Bio::DB::SeqFeature");
    $self->{'feature'} = $value;
  }
  return $self->{'feature'};
}

=head2 featureid

=cut

sub featureid
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'featureid'} = $value;
  }
  return $self->{'featureid'};
}

=head2 overhangs

=cut

sub overhangs
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'overhangs'} = $value;
  }
  return $self->{'overhangs'};
}

=head2 strand

=cut

sub strand
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'strand'} = $value;
  }
  return $self->{'strand'};
}

=head2 peptide

=cut

sub peptide
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'peptide'} = $value;
  }
  return $self->{'peptide'};
}

=head2 offset

=cut

sub offset
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'offset'} = $value;
  }
  return $self->{'offset'};
}

=head2 movers

=cut

sub movers
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'movers'} = $value;
  }
  return $self->{'movers'};
}

=head2 creates

=cut

sub creates
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'creates'} = $value;
  }
  return $self->{'creates'};
}

=head2 dbid

=cut

sub dbid
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'dbid'} = $value;
  }
  return $self->{'dbid'};
}

=head2 phang

=cut

sub phang
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'phang'} = $value;
  }
  return $self->{'phang'};
}

1;

__END__

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, BioStudio developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Johns Hopkins nor the
      names of the developers may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
