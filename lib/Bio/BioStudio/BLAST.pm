#
# BioStudio BLAST interface
#

=head1 NAME

Bio::BioStudio::BLAST

=head1 VERSION

Version 1.05

=head1 DESCRIPTION

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=cut

package Bio::BioStudio::BLAST;

use Exporter;
use Bio::BioStudio::Basic qw($VERNAME);
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::GeneDesign::Codons qw(translate);
use Env;

use strict;
use warnings;

our $VERSION = '1.05';

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
  make_blast_factory
  make_BLAST_db
  make_megaBLAST_index
  bl2seq_blastn
  bl2seq_blastp
);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

=head1 FUNCTIONS

=head2 make_blast_factory

Returns a L<Bio::Tools::Run::StandAloneBlastPlus> object that can be used to 
run BLAST queries.

=cut

sub make_blast_factory
{
  my ($BS) = @_;
  $ENV{BLASTPLUSDIR} = $BS->{blast_executables};
  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                          -DB_DIR => $BS->{blast_directory});
  return $factory;
}

=head2 _make_FASTA

creates a FASTA file for BLAST database creation.

  Arguments: a hashref of chromosome, where the key is the name and the value is 
                the path, see L<Bio::BioStudio::Basic::gather_versions>
             the BioStudio configuration hashref
             a label

=cut

sub _make_FASTA
{
  my ($gff_hsh, $BS, $label) = @_;
  my $FASTAfile = $BS->{blast_directory} . "/" . $label . ".fa";
  print "seeking $FASTAfile\n";
  unless (-e $FASTAfile)
  {
    my $out = Bio::SeqIO->new(-file => ">$FASTAfile", -format => 'Fasta')
      || die "can't make fasta file $FASTAfile for output ($!)";
    foreach my $chr (keys %$gff_hsh)
    {
      print "Loading $chr...\n";
      my $db = Bio::DB::SeqFeature::Store->new(
          -adaptor => 'memory',
          -dir     => $gff_hsh->{$chr} );
      my $chrname = $2 if ($chr =~ $VERNAME);
      my $seqid = "chr$chrname";
      my $bases = $db->fetch_sequence($seqid);
      unless ($bases)
      {
        print "\t CAN'T FIND BASES for $seqid!!!\n";
        next;
      }
      my $seq = Bio::Seq->new( -id => $seqid, -seq => $bases);
      $out->write_seq($seq);
    }
  }
  return;
}

=head2 make_BLAST_db

creates a BLAST database.

  Arguments: a hashref of chromosomes where the key is the name and the value is
                the path, see L<Bio::BioStudio::Basic::gather_versions>
             the BioStudio configuration hashref
             a label

=cut

sub make_BLAST_db
{
  my ($gff_hsh, $BS, $label) = @_;
  my $BLASTdb = $BS->{blast_directory} . "/" . $label . ".nsq";
  my $FASTAfile = $BS->{blast_directory} . "/" . $label. ".fa";
  if (! -e $BLASTdb)
  {
    if (! -e $FASTAfile)
    {
      _make_FASTA($gff_hsh, $BS, $label);
      system("wait");
    }
    my @args = ("$BS->{makeblastdb} -input_type fasta -in $FASTAfile -dbtype nucl -parse_seqids -out $BS->{blast_directory}/$label");
    $SIG{CHLD} = 'DEFAULT';
    system(@args) == 0 or die "system @args failed: $!";
  }
  return $BS->{blast_directory} . "/" . $label;
}

=head2 make_megaBLAST_index
 
creates a megaBLAST index to speed BLASTing

  Arguments: the name of a BLAST database
             the BioStudio configuration hashref
             a label

=cut

sub make_megaBLAST_index
{
  my ($BLASTdb, $BS, $label) = @_;
  my $megaBLASTchk = $BS->{blast_directory} . "/mb" . $label . ".00.idx";
  my $megaBLASTidx = $BS->{blast_directory} . "/mb" . $label;
  if (! -e $megaBLASTchk)
  {
    my @args = ("$BS->{makembindex}", "-input", $BLASTdb);
    push @args, "-output", $megaBLASTidx, "-iformat", "blastdb";
    $SIG{CHLD} = 'DEFAULT';
    system(@args) == 0 or die "system @args failed: $!";
  }
  return $megaBLASTidx;
}

=head2 bl2seq_blastn

Runs a bl2seq on the nucleotide sequences of two features.

  Arguments: two L<Bio::DB::SeqFeature> objects
             a L<Bio::Tools::Run::StandAloneBlastPlus> object, probably from 
                the make_blast_factory function
           
=cut

sub bl2seq_blastn
{
  my ($feat1, $feat2, $blast) = @_;
  my $alnz = $blast->bl2seq(
                          -method  => 'blastn',
                          -query   => $feat2->seq,
                          -subject => $feat1->seq,
                          -method_args => [
                            gapopen => 11,
                            gapextend => 2]);
  $blast->cleanup;
  return $alnz;
}

=head2 bl2seq_blastp

Runs a bl2seq on the translated sequences of two features.

  Arguments: two L<Bio::DB::SeqFeature> objects
             a GeneDesign codon table hashref, from define_codon_table()
             a L<Bio::Tools::Run::StandAloneBlastPlus> object, probably from 
                the make_blast_factory function

=cut

sub bl2seq_blastp
{
  my ($feat1, $feat2, $CODON_TABLE, $blast) = @_;
  my $f1phase = $feat1->phase ? $feat1->phase : 0;
  my $f2phase = $feat2->phase ? $feat2->phase : 0;
  my $newfeat1 = Bio::Seq->new(
                -id => $feat1->id,
                -seq => translate($feat1->seq->seq, $f1phase + 1, $CODON_TABLE),
                -alphabet => 'protein');
  my $newfeat2 = Bio::Seq->new(
                -id => $feat2->id,
                -seq => translate($feat2->seq->seq, $f2phase + 1, $CODON_TABLE),
                -alphabet => 'protein');
  my $alnz = $blast->bl2seq(
                          -method  => 'blastp',
                          -query   => $newfeat2,
                          -subject => $newfeat1);
  $blast->cleanup;
  return $alnz;
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
