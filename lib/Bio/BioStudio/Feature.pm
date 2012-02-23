#
# BioStudio module for sequence segmentation
#

=head1 NAME

Bio::BioStudio::Feature

=head1 VERSION

Version 1.05

=head1 DESCRIPTION

BioStudio object that represents a restriction enzyme - inherits from 
Bio::GeneDesign::Feature and adds attributes for feature annotation
awareness

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::BioStudio::Feature;

use Switch;

use strict;

use base qw(Bio::Root::Root);

my $VERSION = 1.05;

=head1 CONSTRUCTORS

=head2 new

 Title   : new
 Function:
 Returns : a new Bio::BioStudio::Feature object
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  bless $self, $class;
  
  my ($name, $type, $source, $sequence) =
     $self->_rearrange([qw(NAME
                           TYPE
                           SOURCE
                           SEQUENCE)], @args);

  $self->{'name'} = $name;

  $self->{'type'} = $type;

  $self->{'source'} = $source;

  $self->{'sequence'} = $sequence;    

  return $self;
}

=head1 ACCESSORS

=head2 name

=cut

sub name
{
  my ($self) = @_;
  return $self->{'name'};
}

=head2 type

=cut

sub type
{
  my ($self) = @_;
  return $self->{'type'};
}

=head2 source

=cut

sub source
{
  my ($self) = @_;
  return $self->{'source'};
}

=head2 sequence

=cut

sub sequence
{
  my ($self) = @_;
  return $self->{'sequence'};
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
