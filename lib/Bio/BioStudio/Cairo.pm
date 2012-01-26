package Bio::BioStudio::Cairo;
require Exporter;

use Bio::BioStudio::Basic qw($VERNAME);
use DBI;
use Carp;
use Cairo;
use POSIX;
use Math::Trig ':pi';

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.03';

@ISA = qw(Exporter);
@EXPORT_OK = qw(  
  draw_scale
  draw_feature
  draw_RE
  draw_stop
  draw_centromere
  draw_ARS
  draw_SSTR
  draw_CDS
  draw_intron
  draw_amplicon
  draw_repeats
  draw_UTC
);
%EXPORT_TAGS = (
  all => [qw(draw_scale draw_RE draw_stop draw_ARS draw_SSTR draw_feature 
          draw_centromere draw_CDS draw_intron draw_amplicon draw_repeats
          draw_UTC)],
  basic => [qw(draw_scale draw_feature)]
  );

my %featsize = ("site_specific_recombination_target_region" => 50,
      "stop_retained_variant" => 20);
################################################################################
########################### Cairo Drawing  Functions ###########################
################################################################################

sub draw_scale
{
  my ($ctx, $pa) = @_;
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->set_line_width(6);
  my $bigScaleMark = int($pa->{U_BIG_SCALE_MARK} / $pa->{FACTOR});
  my $bigScaleMarkHeight = 100;
  my $bigScaleMarkText = $pa->{DATASTART};
  my $bigScaleMarkEnd = $pa->{LEFT_MARGIN} + $pa->{SCALE_WIDTH};
  my $totalLevel = ceil(($pa->{DATAEND} - $pa->{DATASTART}) / $pa->{LEVEL_WIDTH}) -1;
  for (my $i = 0; $i <= $totalLevel; $i++)
  {
    $ctx->move_to($pa->{LEFT_MARGIN}, $pa->{SCALE_HEIGHT} + $i * $pa->{LEVEL_HEIGHT});
    if ($i == $totalLevel)
    {
      $ctx->rel_line_to($pa->{LEFT_REP_WIDTH} + fmod(($pa->{DATAEND} - $pa->{DATASTART}), $pa->{LEVEL_WIDTH}) / $pa->{FACTOR} + $pa->{RIGHT_REP_WIDTH}, 0);
      $bigScaleMarkEnd = $pa->{LEFT_DATA_MARGIN} + ceil(fmod(($pa->{DATAEND} - $pa->{DATASTART}), $pa->{LEVEL_WIDTH}) / $pa->{FACTOR});
    }
    else
    {
      $ctx->rel_line_to($pa->{SCALE_WIDTH}, 0);
    }
    $ctx->stroke();
    for (my $j = $pa->{LEFT_DATA_MARGIN}; $j <= $bigScaleMarkEnd; $j += $bigScaleMark)
    {
      $ctx->move_to($j, $pa->{SCALE_HEIGHT} + $i * $pa->{LEVEL_HEIGHT});
      $ctx->rel_line_to(0, 0-$bigScaleMarkHeight);
      $ctx->stroke();

      my $thref = $ctx->text_extents(int($bigScaleMarkText / 1000) . " kb");
      $ctx->move_to($j - $thref->{width} / 2, $pa->{SCALE_HEIGHT} + $i * $pa->{LEVEL_HEIGHT} - ($bigScaleMarkHeight + 10));
      $ctx->show_text (int($bigScaleMarkText / 1000) . " kb");
      $bigScaleMarkText += $bigScaleMark * $pa->{FACTOR};

      #comment out after debugging
      for (my $k = $j; $k <= $j + $bigScaleMark; $k+= int(1000/$pa->{FACTOR}))
      {
        #$ctx->move_to($k, $pa->{SCALE_HEIGHT} + $i * $pa->{LEVEL_HEIGHT});
        #$ctx->rel_line_to(0, -50);
        #$ctx->stroke();
      }
    }
    $bigScaleMarkText -= $bigScaleMark * $pa->{FACTOR};
    if ($i == $totalLevel)
    {
      my $thref = $ctx->text_extents(int($pa->{DATAEND}/1000) . " kb");
      if ($pa->{DATAEND} > $bigScaleMarkText + $thref->{width} * $pa->{FACTOR})
      {
        $ctx->move_to($pa->{LEFT_DATA_MARGIN} + fmod(($pa->{DATAEND} - $pa->{DATASTART}), $pa->{LEVEL_WIDTH}) / $pa->{FACTOR}, $pa->{SCALE_HEIGHT} + $i*$pa->{LEVEL_HEIGHT});
        $ctx->rel_line_to(0, 0-$bigScaleMarkHeight);
        $ctx->stroke();

        $ctx->move_to($pa->{LEFT_DATA_MARGIN} + fmod(($pa->{DATAEND}-$pa->{DATASTART}), $pa->{LEVEL_WIDTH}) / $pa->{FACTOR} - $thref->{width} / 2, $pa->{SCALE_HEIGHT} + $i*$pa->{LEVEL_HEIGHT} - ($bigScaleMarkHeight + 10));
        $ctx->show_text (int($pa->{DATAEND}/1000) . " kb");
      }
    }
  }
  return $totalLevel;
}

sub draw_feature
{
  my ($ctx, $obj, $xbeg, $xend, $isdash, $pa) = @_;
  my $feat = $obj->{feat};
  $ctx->set_dash(15,7,7) if $isdash;
  $ctx->set_dash(0) unless $isdash;
  my $y = $obj->{LevelNum} * $pa->{LEVEL_HEIGHT};
  my $flag;
  my $clipEnd = $obj->{DrawEnd};
  if ($clipEnd >= $xend)
  {
    $clipEnd = $xend;
    $flag = 1;
  }
  my $clipStart = $obj->{DrawStart};  
  if ($clipStart <= $xbeg)
  {
     $clipStart = $xbeg;
     $flag = 1;
  }
  my $clipWidth = $clipEnd - $clipStart;
  if (exists $featsize{$feat->primary_tag})
  {
    $clipWidth = $featsize{$feat->primary_tag};
    $clipStart -= .5 *$clipWidth;
  }
  if ($flag)
  {
    $ctx->rectangle($clipStart, $y, $clipWidth, $pa->{LEVEL_HEIGHT});
    #print "\t drawing a $isdash $clipWidth wide rectangle at ($clipStart - $clipEnd), $y for ", 
     #         $feat, " from ", $obj->{DrawStart}, "..", $obj->{DrawEnd}, " ($xbeg..$xend)\n";
    $ctx->clip();    
  }
  
  #Check what feature is it and start drawing (Put it in alphabetical order) 
  if ($feat->primary_tag eq "CDS")
  {
    draw_CDS($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "PCR_product")
  {
    draw_amplicon($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "site_specific_recombination_target_region")
  {
    draw_SSTR($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "ARS")
  {
    draw_ARS($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "stop_retained_variant")
  {
    draw_stop($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "centromere")
  {
    draw_centromere($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "restriction_enzyme_recognition_site")
  {
    draw_RE($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "intron")
  {
    draw_intron($ctx, $obj, $pa);
  }
  elsif ($feat->primary_tag eq "repeat_family")
  {
    draw_repeats($ctx, $obj, $pa);
  }  
  elsif ($feat->primary_tag eq "universal_telomere_cap")
  {
    draw_UTC($ctx, $obj, $pa);
  }  
  #Reset all drawing setting
  $ctx->set_dash(0);
  $ctx->reset_clip();
}

sub draw_RE
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgba($colref->[0], $colref->[1], $colref->[2], .95);
  
  my $RELineHeight = 50.0;  
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $radius = ($end-$start)/2;
  $ctx->move_to($start + $radius, $pa->{SCALE_HEIGHT} + $obj->{LevelNum}*$pa->{LEVEL_HEIGHT});
  $ctx->rel_line_to(0, -$RELineHeight);
  my $thref = ($ctx->text_extents($feat->Tag_enzyme));
  my $tx = $start-$thref->{width}/2;
  my $ty = $pa->{SCALE_HEIGHT} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT} - ($RELineHeight+10);
  $ctx->move_to($tx, $ty);
  $ctx->show_text ($feat->Tag_enzyme);
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
  return; 
}

sub draw_centromere
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  #print "\t\tDrawing " . $feat->Tag_load_id, "\n";
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgb($colref->[0], $colref->[1], $colref->[2]);
  
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $radius = ($end-$start) / 2;
  my $ypos = $obj->{LevelNum} * $pa->{LEVEL_HEIGHT};
    #Centromere radius
  my $centromereRadius = 200/$pa->{FACTOR};
  $centromereRadius = $radius if ($radius > $centromereRadius);
  $centromereRadius = 7 if ($centromereRadius < 7);
  
  #Draw Circle Shape
  $ctx->arc($start + $radius, 
            $pa->{STRAND_Y_POS} + $pa->{STRAND_DISTANCE}/2 + $ypos, 
            $centromereRadius, 
            0, 
            2*pi);
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_stop
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  my $parent = $feat->Tag_parent_id;
  my $pobj = $pa->{FEATURES}->{$parent};
  my $pfeat = $pobj->{feat};
  #print "\t\tDrawing " . $feat->Tag_load_id, "\n";
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgb($colref->[0], $colref->[1], $colref->[2]);
  
  #Preset Size of the Diamond
  my $CodonSide = 80 / $pa->{FACTOR};
  $CodonSide = 4 if ($CodonSide < 4);
  my $smove = $CodonSide / sqrt(2);
  
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $radius = ($end-$start)/2;
  ($start,$end) = ($end, $start) if ($pfeat->strand == -1);
  my $ypos = $obj->{LevelNum} * $pa->{LEVEL_HEIGHT};
  my ($xmove, $ymove) = ($start + $radius, $pa->{STRAND_Y_POS} + $ypos);
  $ymove += $pa->{STRAND_DISTANCE} if ($pfeat->strand == -1);
  $ctx->move_to($xmove, $ymove);
  
  #Draw diamond octagon
  $ctx->rel_move_to(-$CodonSide / 2, -($CodonSide / 2+$smove));
  $ctx->rel_line_to($CodonSide, 0);
  $ctx->rel_line_to($smove, $smove);
  $ctx->rel_line_to(0, $CodonSide);
  $ctx->rel_line_to(-$smove, $smove);
  $ctx->rel_line_to(-$CodonSide, 0);
  $ctx->rel_line_to(-$smove, -$smove);
  $ctx->rel_line_to(0, -$CodonSide);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_ARS
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
#  print "\t\tDrawing " . $feat->Tag_load_id, "\n";

  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgba($colref->[0], $colref->[1], $colref->[2], .95);
  my $ARSHeight = 50;

  #Draw Block Shape
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $ypos = $obj->{LevelNum} * $pa->{LEVEL_HEIGHT};
  $ctx->move_to($start, $pa->{SCALE_HEIGHT}+$ypos);
  $ctx->rel_line_to(0, -$ARSHeight/2);
  $ctx->rel_line_to($end-$start, 0);
  $ctx->rel_line_to(0, $ARSHeight);
  $ctx->rel_line_to(-($end-$start), 0);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_SSTR
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  #print "\t\tDrawing " . $feat->Tag_load_id, "\n";
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgba($colref->[0], $colref->[1], $colref->[2], .95);
             
  #Draw diamond sharp
  my $SSRTsize = $pa->{U_BIG_SCALE_MARK} / ($pa->{FACTOR} * 15);
  $SSRTsize = 7 if ($SSRTsize < 7);
  $SSRTsize = 20 if ($SSRTsize > 20);
  
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $ypos = $obj->{LevelNum}*$pa->{LEVEL_HEIGHT};
  $ctx->move_to($start+($end-$start)/2, $pa->{SCALE_HEIGHT}+$ypos-$SSRTsize);
  $ctx->rel_line_to($SSRTsize, $SSRTsize);
  $ctx->rel_line_to(-$SSRTsize, $SSRTsize);
  $ctx->rel_line_to(-$SSRTsize, -$SSRTsize);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_CDS
{
  my ($ctx, $obj, $pa) = @_;
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $feat = $obj->{feat};
  my $triLen = 50;
  my $CDSHeight = 50;
  #$ctx->rectangle()
  #print "\t\tDrawing " . $feat->Tag_load_id, "\n";
  #Set Color
  my $key = $feat->Tag_essential_status eq "Essential"
          ? "Essential"
          : $feat->Tag_essential_status eq "fast_growth"
            ? "fast_growth"
            : $feat->Tag_orf_classification;
  my $colref = $pa->{FEAT_RGB}->{$key};
  $ctx->set_source_rgb($colref->[0], $colref->[1], $colref->[2]);
  
  #inital the drawing by identify which direction is processing and move to the start point
  my ($smove, $emove);
  if ($feat->strand == 1)
  {
    ($smove, $emove) = ($start, $pa->{STRAND_Y_POS} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT});
  }
  elsif ($feat->strand == -1)
  {
    ($start, $end) = ($end, $start);
    ($smove, $emove) = ($start, $pa->{STRAND_Y_POS} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT} + $pa->{STRAND_DISTANCE});
    my $unitVec = (($end-$start) / abs($end-$start));
    $triLen = $triLen * $unitVec;
  }  
  #Calculate x2 = rectangle's Len , x3 = triangle part's Len |----x2----->
  my ($x2, $x3);
  if (abs($end-$start) <= abs($triLen * 1.3))
  {   
    $x2 = $end-($end-$start)/2-$start;
    $x3 = ($end-$start)/2;
  } 
  else
  {
    $x2 = $end-$triLen-$start;
    $x3 = $triLen;
  }
  $pa->{CDSDATA}->{$feat->Tag_parent_id} = [$smove, $emove, $x2, $x3];
  $ctx->move_to($smove, $emove);
  
  $ctx->rel_line_to(0, 0-$CDSHeight/2);
  $ctx->rel_line_to($x2, 0);
  $ctx->rel_line_to($x3, $CDSHeight/2);
  $ctx->rel_line_to(0-$x3, $CDSHeight/2);
  $ctx->rel_line_to(0-$x2, 0);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_intron
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  my $parent = $feat->Tag_parent_id;
  my $pobj = $pa->{FEATURES}->{$parent};
  my $pfeat = $pobj->{feat};
  my ($triLen, $CDSHeight) = (50, 50);
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgb($colref->[0], $colref->[1], $colref->[2]);
    
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $radius = ($end-$start)/2;
  my $movey = $pa->{STRAND_Y_POS} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT};
  $movey+= $pa->{STRAND_DISTANCE} if ($pfeat->strand == -1);
  $ctx->move_to($start, $movey);

  #Draw two Line
  $ctx->rel_line_to($radius, -$CDSHeight/3.0);
  $ctx->rel_line_to($radius, $CDSHeight/3.0);
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_amplicon
{
  my ($ctx, $obj, $pa) = @_;
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $feat = $obj->{feat};
  my $gene = $feat->Tag_ingene;
  #print "\t\tDrawing " . $feat->Tag_load_id, "\n";
  
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgba($colref->[0], $colref->[1], $colref->[2], .5);
  
  #Use Part of cDS drawing part to get the exact shape
  my ($triLen, $CDSHeight) = (50, 50);
  
  my ($psmove, $pemove, $px2, $px3) = @{$pa->{CDSDATA}->{$feat->Tag_ingene}};
  my ($smove, $emove) = ($start, $pa->{STRAND_Y_POS} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT});
  $emove+= $pa->{STRAND_DISTANCE} if ($feat->strand == -1);
  my $mwidth = abs($end-$start);
  {
    $ctx->move_to($psmove, $pemove);
    $ctx->rel_line_to(0, 0-$CDSHeight/2);
    $ctx->rel_line_to($px2, 0);
    $ctx->rel_line_to($px3, $CDSHeight/2);
    $ctx->rel_line_to(0-$px3, $CDSHeight/2);
    $ctx->rel_line_to(0-$px2, 0);
    $ctx->clip();
  }
  $ctx->move_to($smove, $emove);
  $ctx->rel_line_to(0, $CDSHeight/2);
  $ctx->rel_line_to($mwidth, 0);
  $ctx->rel_line_to(0, -$CDSHeight);
  $ctx->rel_line_to(- $mwidth, 0); 
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_repeats
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  my $repeatfamilyHeight = 50;
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgba($colref->[0], $colref->[1], $colref->[2], .95);
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $width = $end - $start;    
  $ctx->move_to($start, $pa->{SCALE_HEIGHT} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT});
  $ctx->rel_line_to(0, -$repeatfamilyHeight/2);
  $ctx->rel_line_to($width, 0);
  $ctx->rel_line_to(0, $repeatfamilyHeight);
  $ctx->rel_line_to(-$width, 0);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

sub draw_UTC
{
  my ($ctx, $obj, $pa) = @_;
  my $feat = $obj->{feat};
  my $radius = 50;
  my $UTCHeight = 50;
  my $colref = $pa->{FEAT_RGB}->{$feat->primary_tag};
  $ctx->set_source_rgb($colref->[0], $colref->[1], $colref->[2]);
  my ($start, $end) = ($obj->{DrawStart}, $obj->{DrawEnd});
  my $width = $end - $start;
  #inital the drawing by identify which direction is processing and move to the start point
  if ($feat->strand == 1)
  {
    $ctx->move_to($start, $pa->{SCALE_HEIGHT} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT});
  }
  elsif ($feat->strand == -1)
  {
    ($start, $end) = ($end, $start);
    $ctx->move_to($start, $pa->{SCALE_HEIGHT} + $obj->{LevelNum} * $pa->{LEVEL_HEIGHT} + $pa->{STRAND_DISTANCE});
    my $unitVec = ($width/abs($width));
    $radius = $radius*$unitVec;
  }
  my ($x2, $x3);
  #Calculate x2 = rectangle's Len , x3 = triangle part's Len |----x2----->
  if (abs($width) <= abs($radius*1.3))
  {
    $x2 = $end-($width)*2/3-$start;
    $x3 = ($width)*2/3;
  }
  else
  {
    $x2 = $end - $radius - $start;
    $x3 = $radius;
  }

  $ctx->rel_line_to(0, -$UTCHeight/2);
  $ctx->rel_line_to($x2, 0);
  $ctx->rel_curve_to(0, 0, $x3, $UTCHeight/2, 0, $UTCHeight);
  $ctx->rel_line_to(-$x2, 0);
  $ctx->close_path();
  $ctx->fill_preserve();
  $ctx->set_source_rgb(0, 0, 0);
  $ctx->stroke();
}

1;
__END__

=head1 NAME

BioStudio::Cairo

=head1 VERSION

Version 1.03

=head1 DESCRIPTION

BioStudio functions for Cairo interaction.

=head1 FUNCTIONS

=head2 draw_scale()

=head2 draw_feature()

=head2 draw_RE()

=head2 draw_stop()

=head2 draw_centromere()

=head2 draw_ARS()

=head2 draw_SSTR()

=head2 draw_CDS()

=head2 draw_intron()

=head2 draw_amplicon()

=head2 draw_repeats()

=head2 draw_UTC()

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

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
