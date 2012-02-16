#
# BioStudio module for sequence segmentation
#
# POD documentation - main docs before the code

=head1 NAME

Bio::BioStudio::Marker

=head1 VERSION

Version 1.04

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::BioStudio::Marker;

use strict;

use base qw(Bio::Root::Root);

my $VERSION = 1.04;

=head1 CONSTRUCTORS

=head2 new

 Title   : new
 Function:
 Returns : a new Bio::BioStudio::Marker object
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  bless $self, $class;
  
  my ($name, $db) =
     $self->_rearrange([qw(NAME
                           DB)], @args);

  $self->throw("No name defined") unless ($name);
  $self->{'name'} = $name;

  $self->throw("No db defined") unless ($db);
  $self->{'db'} = $db;
  
  my @regions = $db->get_features_by_type("region");
  my $region = $regions[0];
  if ($region->has_tag("color"))
  {
    $self->{'color'} = "#" . join("", $region->get_tag_values("color"));
  }
  $self->sequence($region->seq->seq);

  return $self;
}

=head1 FUNCTIONS

=head1 ACCESSORS

=head2 db

=cut

sub db
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->throw("object of class " . ref($value) . " does not implement ".
		    "Bio::DB::SeqFeature::Store.")
		  unless $value->isa("Bio::DB::SeqFeature::Store");
    $self->{'db'} = $value;
  }
  return $self->{'db'};
}

=head2 sequence

=cut

sub sequence
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'sequence'} = $value;
  }
  return $self->{'sequence'};
}

=head2 name

=cut

sub name
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'name'} = $value;
  }
  return $self->{'name'};
}

=head2 color

=cut

sub color
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'color'} = $value;
  }
  return $self->{'color'};
}

=head2 removeable_enzymes

=cut

sub removeable_enzymes
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'removeable_enzymes'} = $value;
  }
  return $self->{'removeable_enzymes'} ? $self->{'removeable_enzymes'}  : {};
}

=head2 static_enzymes

=cut

sub static_enzymes
{
  my ($self, $value) = @_;
  if (defined $value)
  {
    $self->{'static_enzymes'} = $value;
  }
  return $self->{'static_enzymes'}  ? $self->{'static_enzymes'} : {};
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
