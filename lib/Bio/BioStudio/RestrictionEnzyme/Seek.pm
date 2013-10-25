#
# BioStudio restriction finding functions
#

=head1 NAME

Bio::BioStudio::RestrictionEnzyme::Seek

=head1 VERSION

Version 2.00

=head1 DESCRIPTION

BioStudio functions

=head1 AUTHOR

Sarah Richardson <smrichardson@lbl.gov>.

=cut

package Bio::BioStudio::RestrictionEnzyme::Seek;
require Exporter;

use Bio::BioStudio::RestrictionEnzyme;
use File::Find;
use Carp;
use English qw(-no_match_vars);

use base qw(Exporter);

use strict;
use warnings;

our $VERSION = '2.00';

our @EXPORT_OK = qw(
  find_enzymes_in_CDS
  find_enzymes_in_igen
);
our %EXPORT_TAGS = (BS => \@EXPORT_OK);

=head1 Functions

=cut

=head2 find_enzymes_in_CDS

=cut

sub find_enzymes_in_CDS
{
  my ($chr, $tree, $feat) = @_;
  my $GD = $chr->GD;
  my $phase = $feat->phase || 0;
  my $nucseq = $feat->seq->seq;
  my $aaseq = $GD->translate(-sequence => $nucseq, -frame => 1 + $phase);
  my $hits = $GD->search_suffix_tree(-tree => $tree, -sequence => $aaseq);
  my $results = [];
  foreach my $rog (@{$hits})
  {
    my $enz       = $rog->[0];
    my $nucstart  = $rog->[1] * 3;
    my $realpos   = $feat->start + $nucstart;
    my $peptide   = $rog->[2];
    my @notes     = @{$rog->[3]};
    my $enzyme    = $GD->enzyme_set->{$enz};
    my $flip = 0;
    my $recsite   = $enzyme->recseq();
    if ($notes[0] ne $recsite)
    {
      $flip++;
      $recsite   = $GD->complement(-sequence => $recsite, -reverse => 1);
    }
    my $presence  = "p";
    my $ohang     = {};
    my ($mattersbit, $ohangstart, $ohangend, $fabric) = (q{}, 0, 0, q{});
    my ($offset, $situ, $pproof, $sitelen) = (0, q{}, q{}, 0);

  ##Figure Possible overhangs
    if ($enzyme->type() eq "b" || $enzyme->class() eq "IIB")
    {
      $situ = substr($nucseq, $nucstart, length($peptide)*3);
      $pproof = $GD->translate(-sequence => $situ);
      $ohangstart = 0;
    }
    elsif ($enzyme->class() eq "IIP")
    {
      $sitelen = $enzyme->len;
      if ($sitelen + $realpos <= $feat->end)
      {
        my ($lef, $rig) = (length($1), length($2))
          if ($enzyme->cutseq() =~ $enzyme->classex());
        ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
        $ohangstart = $enzyme->len - $rig + 1;
        $ohangend = $enzyme->len - $lef;
        $situ = substr($nucseq, $nucstart, length($peptide)*3);
        ($fabric, $offset) = $GD->pattern_aligner(
          -sequence => $situ,
          -pattern => $recsite,
          -peptide => $peptide,
          -offset => 1
        );
        $mattersbit = substr($situ,
                             $ohangstart + $offset -1,
                             $ohangend - $ohangstart + 1);
        $pproof = $GD->translate(-sequence => $situ);
      }
    }
    elsif ($enzyme->class() eq "IIA")
    {
      my ($lef, $rig) = ($1, $2)  if ($enzyme->cutseq() =~ $enzyme->classex());
      ($rig, $lef) = ($lef, $rig) if ($rig < $lef);
      $sitelen = $rig >= 0 ? $enzyme->len + $rig  : $enzyme->len;
      if ($sitelen + $realpos <= $feat->end)
      {
        my $nuclen = length($peptide)*3;
        $nuclen++ while($nuclen % 3 != 0);
        $situ = substr($nucseq, $nucstart, $nuclen);
        ($fabric, $offset) = $GD->pattern_aligner(
          -sequence => $situ,
          -pattern => $recsite,
          -peptide => $peptide,
          -offset => 1
        );
        my $add;
        if ($flip == 1)
        {
          $ohangstart = $enzyme->len + $lef + 1;
          if ($rig > 0)
          {
            $add = $rig - (length($fabric) - ($offset + $enzyme->len));
            $add ++ while ($add % 3 != 0);
            $situ .= substr($nucseq, $nucstart + $nuclen, $add);
            $fabric .= "N" x $add;
          }

        }
        else
        {
          if ($rig > 0)
          {
            $add =  $rig - $offset;
            $add ++ while ($add % 3 != 0);
            $situ = substr($nucseq, $nucstart - $add, $add) . $situ;
            $fabric = "N" x $add . $fabric;
            $nucstart = $nucstart - $add;
            $ohangstart = $add - $rig + 1;
          }
          else
          {
            $ohangstart = $offset + abs($rig) + 1;
          }
        }
        $mattersbit = substr($nucseq, $nucstart + $ohangstart+1, $rig-$lef);
        $pproof = $GD->translate(-sequence => $situ);
      }
    }
    else
    {
      print "I don't recognize this type of enzyme: " . $enzyme->id;
    }
    if ($enzyme->type() ne "b" && $enzyme->class() ne "IIB"
     && $realpos + $sitelen <= $feat->end)
    {
      if ($fabric eq "0")
      {
        print "oh no bad fabric, $enz, $fabric, $peptide\n";
        next;
      }
      my $lenm = $mattersbit  ? length($mattersbit) : 0;
      my $matstart = $ohangstart + $offset - 1;
         $matstart-- while($matstart % 3 != 0);
      my $matend = $ohangstart + $offset + $lenm - 1;
         $matend++ while($matend % 3 != 2);
      my $matlen = $matend - $matstart + 1;
      my $peproof = substr($pproof, ($matstart / 3), $matlen / 3);
      my $what = substr($fabric, $matstart, $matlen);
#print "\n\n$flip $enz, $realpos, $peptide, fab: $fabric, situ: $situ ";
#print "pproof: $pproof, ohangstart: $ohangstart, ohangend: $ohangend, ";
#print "matters: $mattersbit, offset: $offset, matstart: $matstart, ";
#print "matlen: $matlen, peproof: $peproof, what: $what\n";
      my $transcs = $GD->ambiguous_transcription($what);
      foreach my $swapseq (@{$transcs})
      {
        next unless ($GD->translate(-sequence => $swapseq) eq $peproof);
        substr($fabric, $matstart, $matlen, $swapseq);
        my $tohang = substr($fabric, $ohangstart +  $offset - 1, $lenm);
#print "\t\t\tsubbing $swapseq into $fabric ($peproof) \tconsidering $tohang\n";
        $ohang->{$tohang}++ if ($tohang);#if ($tohang ne _complement($tohang,1));
      }
    }

  ##Determine Presence
    $presence = "e" if ($situ =~ $enzyme->regex()->[$flip - 1]);
#print "e $enz @ $nucstart $presence $peptide $fabric\n" if ($presence eq "e");
    unless ($enzyme->class() eq "IIB" && $presence eq "p")
    {
      next if ($presence eq 'p' && ! $pproof);
      my $ohangoffset = $ohangstart + $offset - 1;
      #push @{$results}, [$enzyme, $nucstart, $presence, $sitelen, $ohang,
      #                 $pproof, $ohangoffset, $mattersbit, $flip];
      #                 
      #my ($enz, $start, $pre, $len, $hangref, $pep, $offset, $hang, $flip)
      #    = @{$list};
      my $sitestart;
      $sitestart = $feat->start + $nucstart if ($feat->strand != -1);
      $sitestart = $feat->end - $nucstart + 1 - $enzyme->len
        if ($feat->strand == -1);
      my $name = $enzyme->id . "_" . $sitestart . "_" . $pproof;
      my $strand = $flip == 2 ? -1  : 1;
      if ($sitestart + $sitelen - 1 > $feat->end)
      {
        $ohang = $mattersbit  ? {$mattersbit => 1}  : {};
      }
      push @{$results}, Bio::BioStudio::RestrictionEnzyme->new
      (
        -enzyme => $enzyme,
        -name => $name,
        -presence => $presence,
        -start => $sitestart,
        -end => $sitestart + $enzyme->len - 1,
        -feature  => $feat,
        -overhangs  => $ohang,
        -strand => $strand,
        -peptide => $pproof,
        -offset => $ohangoffset
      );
     }
  }
  return $results;
}

=head2 find_enzymes_in_igen

=cut

sub find_enzymes_in_igen
{
  my ($chr, $feat, $pad) = @_;
  $pad = $pad || 0;
  my $chrseq = $chr->sequence();
  my $chrlen = length $chrseq;
  my $GD = $chr->GD;
  my $RES = $GD->enzyme_set();
  my $results = [];
  my $fstart = $feat->start - $pad;
  $fstart = 1 if ($fstart < 1);
  my $fend = $feat->end + $pad;
  $fstart = $chrlen if ($fend > $chrlen);
  my $seq = substr($chrseq, $fstart - 1, $fend - $fstart + 1);
  my $SITESTATUS = $GD->restriction_status(-sequence => $seq);
  foreach my $enzid ( grep {$SITESTATUS->{$_} >= 1} keys %{$SITESTATUS})
  {
    my $enz = $RES->{$enzid};
    my $enzlocs = $enz->positions($seq);
    foreach my $enzpos (keys %{$enzlocs})
    {
      my $siteseq = $enzlocs->{$enzpos};
      my $strand = ($siteseq =~ $enz->regex->[0]) ? 1  : -1;
      my $sitestart = $enzpos + $fstart + 1;
      my ($ohangoffset, $ohangseq) = (0, q{});
      if ($enz->class eq 'IIP')
      {
        ($ohangoffset, $ohangseq) = $enz->overhang($siteseq);
      }
      elsif ($enz->class eq 'IIA')
      {
        my ($lef, $rig) = ($1, $2) if ($enz->cutseq =~ $enz->class_regexes->{IIA});
        ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
        my $newseq;
        if ($strand == 1)
        {
          $newseq = substr($chrseq, $sitestart-2, $enz->len + $rig + 5);
        }
        else
        {
          $newseq = substr($chrseq, $sitestart - ($rig+7), $enz->len + $rig+5);
        }
        ($ohangoffset, $ohangseq) = $enz->overhang($siteseq, $newseq, $strand);
      }
      elsif ($enz->class eq 'IIB')
      {
        my ($rlef, $rrig) = ($3, $4) if ($enz->cutseq =~ $enz->class_regexes->{IIB});
        ($rlef, $rrig) = ($rrig, $rlef) if ($rrig < $rlef);
        ($ohangoffset, $ohangseq) = ($rlef, 'NULL');
      }
      my $ohang = $ohangseq ?  {$ohangseq => 1} : {};
      push @{$results}, Bio::BioStudio::RestrictionEnzyme->new(
        -enzyme => $enz,
        -name => $enz->id . q{_} . $sitestart,
        -presence => 'i',
        -start => $sitestart,
        -end => $sitestart + $enz->len - 1,
        -overhangs  => $ohang,
        -strand => $strand,
        -offset => $ohangoffset
      );
    }
  }
  return $results;
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
