=head1 NAME

Bio::BioStudio::DB

=head1 VERSION

Version 2.00

=head1 DESCRIPTION

BioStudio functions for database interaction.

=head1 AUTHOR

Sarah Richardson <SMRichardson@lbl.gov>.

=cut

package Bio::BioStudio::DB;
require Exporter;

use Bio::BioStudio::ConfigData;
use Bio::DB::SeqFeature::Store;
use DBI;
use English qw(-no_match_vars);
use Carp;

use base qw(Exporter);

use strict;
use warnings;

our $VERSION = '2.00';

our @EXPORT_OK = qw(
  _engine
  _pass
  _user
  _dbh
  _list_databases
  _create_database
  _load_database
  _drop_database
  _fetch_database
  _db_execute
  _db_search
);
our %EXPORT_TAGS = (BS => \@EXPORT_OK);

=head1 DATABASE FUNCTIONS

=head2 _engine

=cut

sub _engine
{
  my $engine = lc Bio::BioStudio::ConfigData->config('db_engine');
  if ($engine =~ m{^y}msix)
  {
    if (Bio::BioStudio::ConfigData->config('mysql_support') =~ m{^y}msix)
    {
      return 'mysql';
    }
    elsif (Bio::BioStudio::ConfigData->config('pg_support') =~ m{^y}msix)
    {
      return 'pg';
    }
  }
  return 'memory';
}

=head2 _pass

=cut

sub _pass
{
  my ($engine) = @_;
  $engine = $engine || _engine();
  if ($engine ne 'memory')
  {
    return Bio::BioStudio::ConfigData->config($engine . "_pass");
  }
  return;
}

=head2 _user

=cut

sub _user
{
  my ($engine) = @_;
  $engine = $engine || _engine();
  if ($engine ne 'memory')
  {
    return Bio::BioStudio::ConfigData->config($engine . "_user");
  }
  return;
}

=head2 _dbh

=cut

sub _dbh
{
  my ($name, $engine) = @_;
  $engine = $engine || _engine();
  if ($engine ne 'memory')
  {
    $name = $name || $engine;
    my $dbh = DBI->connect
    (
      "dbi:$engine:$name", _user($engine), _pass($engine),
      {
        RaiseError => 1,
        AutoCommit => 1
      }
    );
    return $dbh;
  }
  return q{};
}

=head2 _fetch_database

Fetches a Bio::DB::SeqFeature::Store interface for a database containing
the annotations of the argument chromosome. An optional write flag sets whether
or not the interface will support adding, deleting, or modifying features.

  Returns: A L<Bio::DB::SeqFeature::Store> object.

=cut

sub _fetch_database
{
  my ($chromosome, $refresh) = @_;
  $refresh = $refresh || 0;
  my $name = $chromosome->name();
  my $engine = $chromosome->db_engine() || _engine();
  if ($engine ne 'memory')
  {
    my $dblist = _list_databases($engine);
    if (exists $dblist->{$name} && $refresh)
    {
      _drop_database($name, $engine);
      _create_database($name, $engine);
      _load_database($chromosome);
    }
    if (! exists $dblist->{$name})
    {
      _create_database($name, $engine);
      _load_database($chromosome)
    }
    my $db = Bio::DB::SeqFeature::Store->new
    (
      -adaptor  => "DBI::$engine",
      -dsn      => "dbi:$engine:$name",
      -user     => _user($engine),
      -pass     => _pass($engine),
      -write    => 1
    );
    return $db;
  }
  my $db = Bio::DB::SeqFeature::Store->new
  (
    -adaptor  => 'memory',
    -dsn     => $chromosome->path_to_GFF()
  );
  return $db;
}

=head2 _list_databases

This function returns a hash reference containing all BioStudio database names
as keys.

=cut

sub _list_databases
{
  my ($engine) = @_;
  my %dblist;
  $engine = $engine || _engine();
  my $dbh = _dbh(undef, $engine);
  if ($dbh)
  {
    my $sth = $dbh->prepare(q{SHOW DATABASES})
      or croak "Unable to prepare show databases: ". $dbh->errstr."\n";
    $sth->execute or croak "Unable to exec show databases: ". $dbh->errstr."\n";
    my $aref;
    while ($aref = $sth->fetchrow_arrayref)
    {
      $dblist{$aref->[0]}++;
    }
    $sth->finish;
    $dbh->disconnect();
  }
  return \%dblist;
}

=head2 _create_database

This function creates a database that is ready to be loaded with
chromosome data. It does NOT load that database.

=cut

sub _create_database
{
  my ($name, $engine) = @_;
  $engine = $engine || _engine();
  my $dbh = _dbh(undef, $engine);
  if ($dbh)
  {
    $dbh->do("create database $name;");
    $dbh->do("grant select on $name.* to nobody\@localhost;");
    $dbh->do("flush privileges;");
    $dbh->disconnect();
  }
  return;
}

=head2 _load_database

This function loads a database (which must have been previously created,
see create_database()) with a GFF file. The file is the one corresponding to the
first argument provided unless the alternate is defined, in which case the file
corresponding to the third argument is loaded into a database named after the
first argument.

  Arguments: Optionally, the name of an alternate chromosome to be loaded using
               the database name provided in the first argument

=cut

sub _load_database
{
  my ($chromosome) = @_;
  my @args = ("bp_seqfeature_load.pl", "--noverbose");
  push @args, "-f";
  push @args, "-c";
  push @args, "-d", $chromosome->name;
  push @args, $chromosome->path_to_GFF;
  push @args, "--user", _user($chromosome->db_engine);
  push @args, "-p", _pass($chromosome->db_engine);
  local $SIG{CHLD} = 'DEFAULT';
  system(@args) == 0 or croak "system @args failed: $OS_ERROR";
  return;
}

=head2 _drop_database

This function drops the database associated with a BioStudio chromosome.
 
=cut

sub _drop_database
{
  my ($name, $engine) = @_;
  $engine = $engine || _engine();
  my $dblist = _list_databases($engine);
  if (exists $dblist->{$name})
  {
    my $dbh = _dbh($name, $engine);
    $dbh->do("drop database $name;");
    $dbh->do("flush privileges;");
    $dbh->disconnect();
  }
  return;
}

=head2 db_execute

Execute an arbitrary command on an arbitrary database

=cut

sub _db_execute
{
  my ($dbname, $command, $engine) = @_;
  $engine = $engine || _engine();
  my $dbh = _dbh($dbname, $engine);
  if ($dbh)
  {
    $dbh->{'mysql_auto_reconnect'} = 1;
    my $sth = $dbh->prepare($command) or croak ($dbh->errstr . "\n");
    $sth->execute or croak ($dbh->errstr . "\n");
    $sth->finish;
    $dbh->disconnect();
  }
  return;
}

=head2 db_search

Execute a search command on an arbitrary database

=cut

sub _db_search
{
  my ($dbname, $command, $engine) = @_;
  $engine = $engine || _engine();
  my $dbh = _dbh($dbname, $engine);
  my @rowrefs = ();
  if ($dbh)
  {
    $dbh->{'mysql_auto_reconnect'} = 1;
    my $sth = $dbh->prepare($command) or croak ($dbh->errstr . "\n");
    $sth->execute or croak ($dbh->errstr . "\n");
    while (my $aref = $sth->fetchrow_arrayref)
    {
      my @instance = @{$aref};
      push @rowrefs, \@instance;
    }
    $sth->finish;
    $dbh->disconnect();
  }
  return \@rowrefs;
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
