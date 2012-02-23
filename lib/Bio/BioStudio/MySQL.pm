#
# BioStudio module for MySQL database interactions
#

=head1 NAME

Bio::BioStudio::MySQL

=head1 VERSION

Version 1.05

=head1 DESCRIPTION

BioStudio functions for MySQL interaction.

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::BioStudio::MySQL;

use Exporter;
use Bio::BioStudio::Basic qw($VERNAME);
use Bio::BioStudio::RestrictionEnzyme;
use DBI;
use Carp;
use Bio::DB::SeqFeature::Store;

use strict;
use warnings;

our $VERSION = '1.05';

our @ISA = qw(Exporter);
our @EXPORT_OK = qw( 
  drop_database
  list_databases
  create_database
  load_database
  fetch_database
  db_execute
);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

=head2 list_databases

This function returns a hash reference containing all BioStudio database names 
as keys.

  Arguments: The BioStudio configuration hashref

=cut

sub list_databases
{
	my ($BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', $BS->{mysql_user}, $BS->{mysql_pass},
	  { RaiseError => 1, AutoCommit => 1});
	my %dblist;
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
	return \%dblist;
}

=head2 create_database

This function creates a MySQL database that is ready to be loaded with 
chromosome data. It does NOT load that database.

  Arguments: A database name,
             the BioStudio configuration hashref

=cut

sub create_database
{
	my ($chrname, $BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', 
	  $BS->{mysql_user}, 
	  $BS->{mysql_pass},
	  { RaiseError => 1, AutoCommit => 1});
	$dbh->do("create database $chrname;");
	$dbh->do("grant select on $chrname.* to nobody\@localhost;");
	$dbh->do("flush privileges;");
	$dbh->disconnect();
	return;
}

=head2 load_database

This function loads a MySQL database (which must have been previously created, 
see create_database()) with a GFF file. The file is the one corresponding to the
first argument provided unless the alternate is defined, in which case the file
corresponding to the third argument is loaded into a database named after the 
first argument.

  Arguments: Database name as string
             BioStudio configuration hashref
             Optionally, the name of an alternate chromosome to be loaded using
               the database name provided in the first argument

=cut

sub load_database
{
	my ($chrname, $BS, $altchrname) = @_;
  my $gffget = $altchrname  ? $altchrname : $chrname;
	my ($species, $chromname) = ($1, $2) if ($gffget =~ $VERNAME);
	my $fileloc = "$BS->{genome_repository}/$species/chr$chromname/$gffget.gff";
	my @args = ($BS->{bioperl_bin} . "/bp_seqfeature_load.pl", "--noverbose");
	push @args, "-c", "-d", $chrname, $fileloc, "-f";
	push @args, "--user", $BS->{mysql_user}, "-p", "$BS->{mysql_pass}";
	$SIG{CHLD} = 'DEFAULT';
  system(@args) == 0 or croak "system @args failed: $!";
	return;
}

=head2 fetch_database

Fetches a Bio::DB::SeqFeature::Store interface for a MySQL database containing 
the annotations of the argument chromosome. An optional write flag sets whether 
or not the interface will support adding, deleting, or modifying features.

  Arguments: Database name 
             BioStudio configuration hashref
             Optionally, a write flag
             
  Returns: A L<Bio::DB::SeqFeature::Store> object.

=cut

sub fetch_database
{
  my ($chrname, $BS, $write) = @_;
  my $writeflag = $write  ? 1 : 0;
  my $db = Bio::DB::SeqFeature::Store->new(
        -adaptor => "DBI::mysql",
        -dsn => "dbi:mysql:$chrname",
        -user => $BS->{mysql_user},
        -pass => $BS->{mysql_pass},
        -write => $writeflag
  );
  return $db;
}

=head2 drop_database

This function deletes a MySQL database.
 
  Arguments: Database name
             BioStudio configuration hashref

=cut

sub drop_database
{
	my ($chrname, $BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', 
	  $BS->{mysql_user}, 
	  $BS->{mysql_pass},
	  { RaiseError => 1, AutoCommit => 1});
	$dbh->do("drop database $chrname;");
	$dbh->do("flush privileges;");
	$dbh->disconnect();
	return;
}

=head2 db_execute

Execute an arbitrary command on an arbitrary database
  
=cut

sub db_execute
{
  my ($dbname, $command, $BS) = @_;
  my $dbh = DBI->connect("dbi:mysql:$dbname", 
        $BS->{mysql_user},
        $BS->{mysql_pass}, 
        { RaiseError => 1, AutoCommit => 1});
  $dbh->{'mysql_auto_reconnect'} = 1;
  my $sth = $dbh->prepare($command) or die ($dbh->errstr . "\n");
  $sth->execute or die ($dbh->errstr . "\n");
  $sth->finish;
  $dbh->disconnect();
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
