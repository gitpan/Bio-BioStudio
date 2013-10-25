#
# BioStudio object
#

=head1 NAME

Bio::BioStudio

=head1 VERSION

Version 2.00

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <smrichardson@lbl.gov>

=cut

package Bio::BioStudio;
use base qw(Bio::Root::Root);

use Bio::GeneDesign;
use Bio::BioStudio::ConfigData;
use Bio::BioStudio::Chromosome;
use Bio::BioStudio::Repository qw(:BS);
use Bio::BioStudio::DB qw(:BS);
use Bio::BioStudio::Mask qw(:BS);
use File::Path qw(make_path);
use YAML::Tiny;
use Carp;
use DBI;

use strict;
use warnings;

our $VERSION = 2.00;

=head1 CONSTRUCTORS

=head2 new

 Title   : new
 Function:
 Returns : a new Bio::BioStudio object
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
 
  my ($repo) = $self->_rearrange([qw(repo)], @args);
  bless $self, $class;
  
  $self->{bioperl_path} = Bio::BioStudio::ConfigData->config('bioperl_path');
  $self->{script_path} = Bio::BioStudio::ConfigData->config('script_path');
  $self->{tmp_path} = Bio::BioStudio::ConfigData->config('tmp_path');

  $self->{conf} = Bio::BioStudio::ConfigData->config('conf_path');

  $self->{path_to_features} = $self->{conf} . "features/";
  $self->{path_to_markers} = $self->{conf} . "markers/";

  $self->throw("$repo does not exist") if ($repo && ! -e $repo);
  $repo = $repo || $self->{conf} . 'genome_repository/';
  $self->{repo} = $repo;
  
  my $gbrowse = Bio::BioStudio::ConfigData->config('gbrowse_support');
  if ($gbrowse eq "Y")
  {
    $self->{gbrowse} = 1;
  }
  else
  {
    $self->{gbrowse} = 0;
  }
  
  my $cairo = Bio::BioStudio::ConfigData->config('cairo_support');
  if ($cairo && $cairo eq "Y")
  {
    $self->{cairo} = 1;
    $self->{cairo_conf_path} = $self->conf . "cairo/";
    $self->{cairo_colors_path} = $self->{cairo_conf_path} . "Cairo_colors.yaml";
  }
  else
  {
    $self->{cairo} = 0;
  }

  $self->{db_engine} = _engine();

  return $self;
}

=head1 FUNCTIONS

=head1 ACCESSORS

=cut

=head2 path_to_markers

=cut

sub path_to_markers
{
  my ($self) = @_;
  return $self->{path_to_markers};
}

=head2 path_to_repo

=cut

sub path_to_repo
{
  my ($self) = @_;
  return $self->{repo};
}

=head2 path_to_features

=cut

sub path_to_features
{
  my ($self) = @_;
  return $self->{path_to_features};
}

=head2 tmp_path

=cut

sub tmp_path
{
  my ($self) = @_;
  return $self->{tmp_path};
}

=head2 script_path

=cut

sub script_path
{
  my ($self) = @_;
  return $self->{script_path};
}

=head2 gbrowse

=cut

sub gbrowse
{
  my ($self) = @_;
  return $self->{gbrowse};
}

=head2 cairo

=cut

sub cairo
{
  my ($self) = @_;
  return $self->{cairo};
}

=head2 bioperl_path

=cut

sub bioperl_path
{
  my ($self) = @_;
  return $self->{bioperl_path};
}

=head2 conf

=cut

sub conf
{
  my ($self) = @_;
  return $self->{conf};
}

=head1 FUNCTIONS

=head2 set_chromosome

=cut

sub set_chromosome
{
  my ($self, @args) = @_;
  
  my ($chrname, $gbrowse) = $self->_rearrange([qw(chromosome gbrowse)], @args);
  
  $self->throw('No chromosome name provided') unless ($chrname);
  
  $gbrowse = $gbrowse || 0;
  
  my $chr = Bio::BioStudio::Chromosome->new(
    -name       => $chrname,
    -repo       => $self->{repo},
    -db_engine  => $self->{db_engine},
    -gbrowse    => $gbrowse
  );
  
  return $chr;
}

=head2 gather_latest_genome

=cut

sub gather_latest_genome
{
  my ($self, @args) = @_;
  
  my ($species) = $self->_rearrange([qw(species)], @args);
  
  $self->throw('No species provided') unless ($species);
  
  my @objset = ();
  
  my $latest_names = _gather_latest($species, $self->{repo});
  
  foreach my $chrname (@{$latest_names})
  {
    my $chr = Bio::BioStudio::Chromosome->new(-name => $chrname);
    push @objset, $chr->seqobj;
  }
  return \@objset;
}

=head2 fetch_custom_features()

Returns a hashref of custom features in the BioStudio configuration directory.
Each key is a feature name, each value is a Bio::SeqFeature object.

=cut

sub fetch_custom_features
{
  my ($self) = @_;
  opendir(my $FDIR, $self->{path_to_features});
  my @features = grep {$_ =~ m{\.yaml\z}msix} readdir($FDIR);
  closedir $FDIR;
  my %features;
  foreach my $feature (@features)
  {
    my $path = $self->path_to_features . $feature;
    
    $self->throw("$path does not exist!") unless (-e $path);
    my $yaml = YAML::Tiny->read($path);
    
    my $name     = $yaml->[0]->{name};
    my $type     = $yaml->[0]->{type};
    my $source   = $yaml->[0]->{source};
    my $sequence = $yaml->[0]->{sequence};
    my $length   = length $sequence;
    
    my $feat = Bio::SeqFeature::Generic->new(
      -start         => 1,
      -end           => $length,
      -primary_tag   => $type,
      -source_tag    => $source,
      -display_name  => $name,
    );
    $feat->attach_seq( Bio::Seq->new(-id => $name, -seq => $sequence) );
      
    $features{$name} = $feat;
  }
  return \%features;
}

=head2 fetch_custom_markers

=cut

sub fetch_custom_markers
{
  my ($self) = @_;
  opendir(my $FDIR, $self->{path_to_markers});
  my @markers = grep {$_ =~ m{\.gff\z}msix} readdir($FDIR);
  closedir $FDIR;
  my %markers;
  my $n = 1;
  foreach my $marker (@markers)
  {
    my $name = 'marker' . $n;
    $n++;
    $name = $1 if ($marker =~ m{([\w\d]+)\.gff\z}msix);
    my $path = $self->{path_to_markers} . $marker;
    my $db = Bio::DB::SeqFeature::Store->new(
      -adaptor => 'memory',
      -gff     => $path,
      -index_subfeatures => "true"
    );
    $markers{$name} = Bio::BioStudio::Marker->new(
      -name => $name,
      -db => $db
    );
  }
  return \%markers;
}

=head2 species_list

=cut

sub species_list
{
  my ($self) = @_;
  return _list_species($self->{repo});
}

=head2 source_list

=cut

sub source_list
{
  my ($self) = @_;
  return _list_repository($self->{repo});
}

=head2 db_list

=cut

sub db_list
{
  my ($self) = @_;
  return _list_databases($self->{db_engine});
}

=head2 gv_increment_warning

=cut

sub gv_increment_warning
{
  my ($self, $chromosome) = @_;

  my $species = $chromosome->species();
  my $seqid   = $chromosome->seq_id();
  my $cver    = $chromosome->chromosome_version();
  my $ngver   = $chromosome->genome_version() + 1;
  my $gscale  = $species . "_" . $seqid . "_" . $ngver . "_" . $cver;
  my $dblist  = $self->source_list();
  return $gscale if (exists $dblist->{$gscale});
  return 0;
}

=head2 cv_increment_warning

=cut

sub cv_increment_warning
{
  my ($self, $chromosome) = @_;

  my $species = $chromosome->species();
  my $seqid   = $chromosome->seq_id();
  my $gver    = $chromosome->genome_version();
  my $cver    = $chromosome->chromosome_version();
  my $ncver   = $chromosome->GD->pad($cver + 1, 2);
  my $cscale = $species . "_" . $seqid . "_" . $gver  . "_" . $ncver;
  my $dblist = $self->source_list();
  return $cscale if (exists $dblist->{$cscale});
  return 0;
}

=head2 prepare_repository

=cut

sub prepare_repository
{
  my ($self, @args) = @_;
  
  my ($species, $chrname) = $self->_rearrange([qw(species chrname)], @args);
  
  my $path = _prepare_repository($self->{repo}, $species, $chrname);
  return $path;
}

=head1 ACCESSORS

=head2 genome_repository

A path to a directory that can be used as a genome repository. Defaults to
the config dir as set at install slash genome_repository.

=cut

sub genome_repository
{
  my ($self, $value) = @_;
  if ($value)
  {
    $self->throw("$value does not exist") unless (-e $value);
    $value .= q{/} unless substr($value, -1, 1) eq q{/};
    $self->{repo} = $value;
  }
  return $self->{repo};
}

=head2 db_engine

A path to a directory that can be used as a genome repository. Defaults to
the config dir as set at install slash genome_repository.

=cut

sub db_engine
{
  my ($self, $value) = @_;
  if ($value)
  {
    $self->{db_engine} = $value;
  }
  return $self->{db_engine};
}

1;

__END__

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