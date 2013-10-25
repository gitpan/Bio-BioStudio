use 5.006;

use Module::Build;
use Carp;
use File::Basename;

use strict;
use warnings;

print "\n";

my $llt = 0;
my ($cpath, $spath, $tpath) = (q{}, q{}, q{});
my ($dcpath, $dspath, $dtpath) = ('/etc/BioStudio', '/usr/local/bin/', '/tmp');
my ($gb, $gba, $bl, $c, $db) = (0, 0, 0, 0, 0);
my ($dgb, $dgba, $dbl, $dc, $ddb) = ('Y', '127.0.0.1', 'Y', 'Y', 'Y');
my ($ps, $pu, $pp, $ms, $mu, $mp) = (0, 0, 0, 0, 0, 0);
my ($dps, $dpu, $dpp) = ('N', 'bss', 'bsspass');
my ($dms, $dmu, $dmp) = ('Y', 'bss', 'bsspass');

my $check = eval
{
  require Bio::BioStudio::ConfigData;
};
if ($check)
{
  $dcpath = Bio::BioStudio::ConfigData->config('conf_path')        || $dcpath;
  $dtpath = Bio::BioStudio::ConfigData->config('tmp_path')         || $dtpath;
  $dspath = Bio::BioStudio::ConfigData->config('script_path')      || $dspath;
  $dgb    = Bio::BioStudio::ConfigData->config('gbrowse_support')  || $dgb;
  $dgba   = Bio::BioStudio::ConfigData->config('gbrowse_address')  || $dgba;
  $dbl    = Bio::BioStudio::ConfigData->config('blast_support')    || $dbl;
  $dc     = Bio::BioStudio::ConfigData->config('cairo_support')    || $dc;
  $ddb    = Bio::BioStudio::ConfigData->config('db_engine')        || $ddb;
  $dps    = Bio::BioStudio::ConfigData->config('pg_support')       || $dps;
  $dpu    = Bio::BioStudio::ConfigData->config('pg_user')          || $dpu;
  $dpp    = Bio::BioStudio::ConfigData->config('pg_pass')          || $dpp;
  $dms    = Bio::BioStudio::ConfigData->config('mysql_support')    || $dms;
  $dmp    = Bio::BioStudio::ConfigData->config('mysql_pass')       || $dmp;
  $dmu    = Bio::BioStudio::ConfigData->config('mysql_user')       || $dmu;
}

#BUILDER
my $BSB = Module::Build->new(
    module_name         => 'Bio::BioStudio',
    license             => 'bsd',
    dist_author         => q{Sarah Richardson <smrichardson@lbl.gov>},
    dist_version        => '2.00',
    dist_abstract       => 'perl for the design of synthetic chromosomes',
    add_to_cleanup      => [ 'Bio::BioStudio-*' ],
    create_makefile_pl  => 'traditional',
    build_requires      =>
    {
      'Test::More' => 0,
    },
    requires            =>
    {
      'perl'                    => 5.006,
      'Text::Diff'              => 1.3,
      'Time::Format'            => 1.11,
      'Digest::MD5'             => 2.51,
      'URI::Escape'             => 1.59,
      'YAML::Tiny'              => 1.50,
      'Bio::Root::Root'         => 0,
      'Bio::GeneDesign'         => 5.00,
    },
    configure_requires  =>
    {
      'Module::Build' => 0.38
    },
    auto_features       =>
    {
      gbrowse =>
      {
        description => "Use GBrowse2 as a BioStudio GUI",
        requires    =>
        {
          'Bio::Graphics::Browser2' => 2.42,
          'Cache::FileCache'        => 0
        }
      },
      mysql =>
      {
        description => "Use MySQL as a database engine",
        requires    =>
        {
          'DBI' => 1.616,
          'DBD::mysql' => 4.020
        }
      },
      pg =>
      {
        description => "Use Postgre as a database engine",
        requires    =>
        {
          'DBI' => 1.616,
          'DBD::Pg' => 2
        }
      },
      cairo =>
      {
        description => "Use Cairo to generate chromosome maps",
        requires    =>
        {
          'Cairo'           => 1.081,
          'Font::FreeType'  => 0.03,
        }
      },
      blast =>
      {
        description => "Use BLAST+ to design primers and make alignments",
        requires    =>
        {
          'Bio::Tools::Run::StandAloneBlastPlus' => 0
        }
      }
    },
    script_files        =>
    [
      'bin/BS_AddFromGenbank.pl',
      'bin/BS_AddFromGFF.pl',
      'bin/BS_Cartoonist.pl',
      'bin/BS_ChromosomeAnalyzer.pl',
      'bin/BS_ChromosomeCutter.pl',
      'bin/BS_ChromosomeDiff.pl',
      'bin/BS_ChromosomeSegmentationPlanner.pl',
      'bin/BS_ChromosomeSplicer.pl',
      'bin/BS_CodonJuggler.pl',
      'bin/BS_FirstGBrowse.pl',
      'bin/BS_GlobalREMarkup.pl',
      'bin/BS_GlobalREMarkup.pl',
      'bin/BS_PCRTagDumper.pl',
      'bin/BS_PCRTagger.pl',
      'bin/BS_RepeatSmasher.pl',
      'bin/BS_SeqChecker.pl',
    ],
    get_options =>
    {
      like_last_time =>
      {
        store => \$llt,
        default => 0,
      },
      conf_path =>
      {
        store => \$cpath,
        type => '=s',
      },
      script_path =>
      {
        store => \$spath,
        type => '=s',
      },
      tmp_path =>
      {
        store => \$tpath,
        type => '=s',
      },
      gbrowse_support =>
      {
        store => \$gb,
        type => '=s',
      },
      gbrowse_address =>
      {
        store => \$gba,
        type => '=s',
      },
      blast_support =>
      {
        store => \$bl,
        type => '=s',
      },
      cairo_support =>
      {
        store => \$c,
        type => '=s',
      },
      db_engine =>
      {
        store => \$db,
        type => '=s',
      },
      pg_support =>
      {
        store => \$ps,
        type => '=s',
      },
      pg_user =>
      {
        store => \$pu,
        type => '=s',
      },
      pg_pass =>
      {
        store => \$pp,
        type => '=s',
      },
      mysql_support =>
      {
        store => \$ms,
        type => '=s',
      },
      mysql_user =>
      {
        store => \$mu,
        type => '=s',
      },
      mysql_pass =>
      {
        store => \$mp,
        type => '=s',
      },
    }
);

if ($llt == 0)
{
  if (! $cpath)
  {
    $cpath = $BSB->prompt('Where should configuration files be installed?', $dcpath);
  }
  if (! $spath)
  {
    $spath = $BSB->prompt('Where should scripts be installed?', $dspath);
  }
  if (! $tpath)
  {
    $tpath = $BSB->prompt('Where should BioStudio write tmp files?', $dtpath);
  }
  if (! $gb && $BSB->feature('gbrowse'))
  {
    $gb = $BSB->y_n('Enable GBrowse interaction?', $dgb);
    $gba = $BSB->prompt(' Address of GBrowse server?', $dgba);
  }
  if (! $bl && $BSB->feature('blast'))
  {
    $bl = $BSB->y_n('Enable BLAST+ support?', $dbl);
  }
  if (! $c && $BSB->feature('cairo'))
  {
    $c = $BSB->y_n('Enable Cairo graphics support?', $dc);
  }
  if (! $db && ($BSB->feature('mysql') && $BSB->feature('pg')))
  {
    $db = $BSB->y_n('Use a database engine for annotation?', $ddb);
  }
  if ($db eq 'Y' && $BSB->feature('mysql'))
  {
    $ms = $BSB->y_n('Enable MySQL support?', $dms);
    $mu = $BSB->prompt(' Username for MySQL?', $dmu);
    $mp = $BSB->prompt(' Password for MySQL?', $dmp);
  }
  if ($db eq 'Y' && $ms eq 'N' && $BSB->feature('pg'))
  {
    $ps = $BSB->y_n('Enable Postgres support?', $dps);
    $pu = $BSB->prompt(' Username for Postgres?', $dpu);
    $pp = $BSB->prompt(' Password for Postgres?', $dpp);
  }
}
else
{
  $cpath = $cpath || $dcpath;
  $tpath = $tpath || $dtpath;
  $spath = $spath || $dspath;
  $gb  = $gb  || $dgb;
  $gba = $gba || $dgba;
  $bl  = $bl  || $dbl;
  $c   = $c   || $dc;
  $db  = $db  || $ddb;
  $ps  = $ps  || $dps;
  $pu  = $pu  || $dpu;
  $pp  = $pp  || $dpp;
  $ms  = $ms  || $dms;
  $mu  = $mu  || $dmu;
  $mp  = $mp  || $dmp;
}


my $ttp = $tpath;
$ttp .= q{/} unless substr($ttp, -1, 1) eq q{/};
$BSB->config_data(tmp_path => $ttp);

my $tsp = $spath;
$tsp .= q{/} unless substr($tsp, -1, 1) eq q{/};
$BSB->config_data(script_path => $tsp);

#Prepare configuration directory
my $tcp = $cpath;
$tcp .= q{/} unless substr($tcp, -1, 1) eq q{/};
$BSB->config_data(conf_path => $tcp);
my $confs =
[
  'cairo/Armata-Regular.otf',
  'cairo/Armata_license.txt',
  'cairo/Aurulent_license.txt',
  'cairo/AurulentSansMono-Regular.otf',
  'cairo/Bree_license.txt',
  'cairo/BreeSerif-Regular.otf',
  'cairo/Cairo_colors.yaml',
  'cairo/Inconsolata.otf',
  'cairo/Inconsolata_license.txt',
  'cairo/Pecita.otf',
  'cairo/Pecita_license.txt',
  'cairo/Verily_license.txt',
  'cairo/VerilySerifMono.otf',
  'features/loxPsym.yaml',
  'features/rox.yaml',
  'features/UTC_left.yaml',
  'features/UTC_right.yaml',
  'foswiki/WebAtom.txt',
  'foswiki/WebChanges.txt',
  'foswiki/WebCreateNewTopic.txt',
  'foswiki/WebHome.txt',
  'foswiki/WebIndex.txt',
  'foswiki/WebLeftBar.txt',
  'foswiki/WebNotify.txt',
  'foswiki/WebPreferences.txt',
  'foswiki/WebRss.txt',
  'foswiki/WebSearch.txt',
  'foswiki/WebSearchAdvanced.txt',
  'foswiki/WebStatistics.txt',
  'foswiki/WebTopicList.txt',
  'gbrowse/BS_GBrowse_chromosome.conf',
  'gbrowse/GBrowse_colors.yaml',
  'gbrowse/conf_repository/repository_README.txt',
  'genome_repository/repository_README.txt',
  'genome_repository/Escherichia_coli/chr01/Escherichia_coli_chr01_0_00.gff',
  'markers/LEU2.gff',
  'markers/URA3.gff'
];
$BSB->add_build_element('BioStudio');
process_conf_files($BSB, $confs);
$BSB->install_path(BioStudio => $BSB->config_data('conf_path'));
print "\tWill install configuration files in ";
print $BSB->config_data('conf_path') . "\n";

#Check for bioperl scripts
my $bp = `which bp_seqfeature_load.pl`;
die ("Either BioPerl scripts are not installed or they are not in my PATH")
  unless ($bp);
my ($x, $bppath) = fileparse($bp);
$BSB->config_data('bioperl_path' => $bppath);
print "\tWill use BioPerl scripts found in $bppath\n";

#Optionally configure GBrowse if GBrowse is installed
if ($BSB->feature('gbrowse'))
{
  $BSB->config_data(gbrowse_support => $gb);
  $BSB->config_data(gbrowse_address => $gba);
}
if ( $BSB->config_data('gbrowse_support') )
{
  my $plugins =
  [
    'gbrowse_plugins/BS_Cartoonist.pm',
    'gbrowse_plugins/BS_ChromosomeAnalyzer.pm',
    'gbrowse_plugins/BS_ChromosomeCutter.pm',
    'gbrowse_plugins/BS_ChromosomeDiff.pm',
    'gbrowse_plugins/BS_ChromosomeSegmentationPlanner.pm',
    'gbrowse_plugins/BS_ChromosomeSegmenter.pm',
    'gbrowse_plugins/BS_ChromosomeSplicer.pm',
    'gbrowse_plugins/BS_CodonJuggler.pm',
    'gbrowse_plugins/BS_PCRTagDumper.pm',
    'gbrowse_plugins/BS_PCRTagger.pm',
    'gbrowse_plugins/BS_RepeatSmasher.pm',
  ];
  require GBrowse::ConfigData;
  my $plugin_path = GBrowse::ConfigData->config('conf') . q{/} . 'plugins/';
  $BSB->config_data(gbrowse_plugin_path => $plugin_path);
  $BSB->add_build_element('gbrowse_plugins');
  $BSB->install_path(gbrowse_plugins => $plugin_path);
  print "\tWill install gbrowse plugins in $plugin_path\n";
  print "\tWill run with gbrowse installation in $gba\n";
  process_gbrowse_files($BSB, $plugins);
  chmod 0777, GBrowse::ConfigData->config('conf') . q{/} . 'GBrowse.conf';
}

#Optionally configure BLAST+ if Standalonetools is installed
$BSB->config_data(blast_support => $bl)   if ($BSB->feature('blast'));
if ( $BSB->config_data('blast_support') )
{
  my $blp = `which blastn`;
  die ("Either BLAST+ is not installed or it is not on my PATH") unless ($blp);
  my ($y, $blpath) = fileparse($blp);
  print "\tWill use BLAST+ executables found in $blpath\n";
}

#Optionally configure Cairo if Cairo is installed
$BSB->config_data(cairo_support => $c)   if ($BSB->feature('cairo'));

$BSB->config_data(db_engine => $db);
if ($BSB->feature('mysql') && $db eq 'Y')
{
  $BSB->config_data(mysql_support => $ms);
  $BSB->config_data(mysql_user => $mu);
  $BSB->config_data(mysql_pass => $mp);
  print "\tYou must ensure that BioStudio can interact with mysql:\n";
  print "\t\t\$ mysql -u root\n";
  print "\t\t" . '> CREATE USER "' . $mu . '"@"localhost" IDENTIFIED BY "' ;
  print $mp . q{";} . "\n";
  print "\t\t" . '> GRANT SELECT, INSERT, UPDATE, DELETE, CREATE, GRANT ';
  print 'OPTION, ALTER, DROP ON *.* TO "' . $mu . q{"@"localhost";} . "\n";
  print "\tYeah, that's not the most secure, sorry! Working on it!\n\n";
}
if ($BSB->feature('pg') && $db eq 'Y' && $ms eq 'N')
{
  $BSB->config_data(pg_support => $ps);
  $BSB->config_data(pg_user => $pu);
  $BSB->config_data(pg_pass => $pp);
}

$BSB->create_build_script();

print "\n\n Remember to make the configuration directory writeable so ";
print "BioStudio can edit its own repository:\n";
print "\n\t\t sudo chmod -R 777 " . $BSB->config_data('conf_path') . "\n\n\n";

sub process_conf_files
{
  my $BSB = shift;
  my $files = shift;
  return unless $files;

  my $conf_dir = File::Spec->catdir($BSB->blib, 'BioStudio');
  File::Path::mkpath( $conf_dir );

  foreach my $file (@{$files})
  {
    my $result = $BSB->copy_if_modified($file, $conf_dir) or next;
    $BSB->fix_shebang_line($result) unless $BSB->is_vmsish;
  }
}

sub process_gbrowse_files
{
  my $BSB = shift;
  my $files = shift;
  return unless $files;

  my $gbrowse_dir = File::Spec->catdir($BSB->blib, 'gbrowse_plugins');
  File::Path::mkpath( $gbrowse_dir );

  foreach my $file (@{$files})
  {
    my $result = $BSB->copy_if_modified($file, $gbrowse_dir, 'flatten') or next;
    $BSB->fix_shebang_line($result) unless $BSB->is_vmsish;
    $BSB->make_executable($result);
  }
}
