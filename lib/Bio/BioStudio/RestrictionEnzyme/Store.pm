#
# BioStudio module for sequence segmentation
#
# POD documentation - main docs before the code

=head1 NAME

Bio::BioStudio::RestrictionEnzyme::Store

=head1 VERSION

Version 1.04

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>

=cut

package Bio::BioStudio::RestrictionEnzyme::Store;

use Switch;

use strict;

use base qw(Bio::Root::Root);

my $VERSION = 1.04;
my $tblname = "positions";

=head1 CONSTRUCTORS

=head2 new

 Title   : new
 Function:
 Returns :
 Args    :

=cut

sub new
{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($name, $user, $pass, $file, $create, $RES) =
     $self->_rearrange([qw(NAME
                           USER
                           PASS
                           FILE
                           CREATE
                           ENZYME_DEFINITIONS)], @args);

  $self->throw("No name defined") unless ($name);
  $self->{'name'} = $name;

  $self->throw("No enzymes defined") unless ($RES);
  $self->{'enzyme_definitions'} = $RES;
    
  $file && $self->dumpfile($file);
  
  $user && $self->user($user);
  $self->{'pass'} = $pass if ($pass);
  
  if ($create)
  {
    my $dbh = DBI->connect('dbi:mysql:mysql', $user, $pass,
  	                { RaiseError => 1, AutoCommit => 1});
  	$dbh->do("drop database IF EXISTS $name;");
    $dbh->do("create database $name;");
  	$dbh->disconnect();
  }
  
  my $dsn = "dbi:mysql:" . $self->name;
  my $dbh = DBI->connect($dsn, $self->user, $self->{'pass'}, 
          {RaiseError => 1, AutoCommit => 1, mysql_auto_reconnect => 1});
  $self->{'dbh'} = $dbh;
  
  $self->_initialize() if ($create);
  
  return $self;
}

=head1 FUNCTIONS

=head2 _initialize

=cut

sub _initialize
{
  my ($self) = @_;
  my $def = $self->_table_definition;
  my $command = "CREATE table IF NOT EXISTS $tblname $def->{$tblname}";
  $self->dbh->do($command) or die ($self->dbh->errstr."\n");
  return;
}

=head2 load

=cut

sub load
{
  my ($self) = @_;
  $self->throw("No dumpfile described!") unless ($self->dumpfile);
  $self->throw("Can't find dumpfile!") unless(-e $self->dumpfile);
  my $command = "LOAD DATA INFILE \"" . $self->dumpfile . "\" INTO TABLE ";
  $command .= $self->name . "." . $tblname;
  $command .= " FIELDS TERMINATED BY '.' LINES TERMINATED BY '\n' ";
  $command .= "(name, presence, start, end, enzyme, feature, peptide, ";
  $command .= "overhangs, strand, overhangoffset);";
  $self->dbh->do($command) or die ($self->dbh->errstr."\n");
  return;
}

=head2 search

=cut

sub search
{
  my ($self, @args) = @_;
  my ($name, $left, $right, $enzyme) = $self->_rearrange([qw(
       NAME   LEFT   RIGHT   ENZYME)], @args);
  my $command = "SELECT id, name, presence, eligible, start, end, enzyme, ";
  $command .= "feature, strand, overhangs, peptide, overhangoffset FROM ";
  $command .= $tblname;
  my @wheres = ();
  push @wheres, "name = '$name'" if ($name);
  push @wheres, "start >= $left" if ($left);
  push @wheres, "end <= $right" if ($right);
  push @wheres, "enzyme = '$enzyme'" if ($enzyme);
  $command .= " WHERE " . join(" AND ", @wheres) if (scalar(@wheres));
  $command .= ";";
  my $sth = $self->dbh->prepare($command) or die ($self->dbh->errstr."\n");
  $sth->execute or die ($self->dbh->errstr."\n");
  my @res = ();
  while (my $aref = $sth->fetchrow_arrayref)
  {
    my ($id, $name, $pre, $elg, $start, $end, $eid, $featid, $strand, 
        $overhangs, $peptide, $overhangoffset) = @$aref;
    my %overhangs = map {$_ => 1} split(",", $overhangs);
    my $enzyme = $self->enzyme_definitions->{$eid};
    $self->throw("$eid is undefined for this database") unless $enzyme;
    push @res, Bio::BioStudio::RestrictionEnzyme->new(
                -enzyme => $enzyme,
                -dbid => $id,
                -name => $name,
                -presence => $pre,
                -eligible => $elg,
                -start => $start,
                -end => $end,
                -featureid => $featid,
                -strand => $strand,
                -overhangs => \%overhangs,
                -peptide => $peptide,
                -offset => $overhangoffset
    );
  }
  $sth->finish; 
  return \@res;
}

=head2 cull

=cut

sub cull
{
  my ($self, $cullref) = @_;  
  my @list = @{$cullref};
  while (my @portion = splice(@list, 0, 500))
  {
    my $command = "DELETE from positions where id in (";
    $command .= join(",", @portion) . ");";
    
    my $sth = $self->dbh->prepare($command) or die ($self->dbh->errstr."\n");
    $sth->execute or die ($self->dbh->errstr."\n");    
  }
  return;
}

=head2 screen

=cut

sub screen
{
  my ($self, $screenref) = @_;
  my @list = @{$screenref};
  while (my @portion = splice(@list, 0, 500))
  {
    my $command = "UPDATE positions set `eligible` = \"no\" where id in (";
    $command .= join(",", @portion) . ");";
    
    my $sth = $self->dbh->prepare($command) or die ($self->dbh->errstr."\n");
    $sth->execute or die ($self->dbh->errstr."\n");    
  }
  return;
}

=head1 ACCESSORS

=head2 dbh

=cut

sub dbh
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'dbh'} = $value;
  }
  return $self->{'dbh'};
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

=head2 user

=cut

sub user
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'user'} = $value;
  }
  return $self->{'user'};
}

=head2 dumpfile

=cut

sub dumpfile
{
  my ($self, $value) = @_;
  if (defined $value)
  {
	  $self->{'dumpfile'} = $value;
  }
  return $self->{'dumpfile'};
}

=head2 enzyme_definitions

=cut

sub enzyme_definitions
{
  my ($self) = @_;
  return $self->{'enzyme_definitions'};
}

=head2 _table_definition

=cut

sub _table_definition
{
  my $hsh = {
    $tblname => <<END,
(
  id              int           not null auto_increment primary key,
  name            varchar(45)   not null, 
  presence        varchar(40)   not null, 
  eligible        varchar(3)    null,      
  start           int           not null,         
  end             int           not null,         
  enzyme          varchar(45)   not null, 
  feature         varchar(100)  not null,
  peptide         varchar(45)   null,     
  overhangs       text          null,            
  strand          varchar(3)    not null,
  overhangoffset  int           null,
  index ENZYME (enzyme ASC),
  index POSITION (start ASC, end ASC)
)
END

  };
  return $hsh;
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
