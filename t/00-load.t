#!perl -T

use Test::More tests => 8;

BEGIN {
    use_ok( 'Bio::BioStudio::Basic' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::GFF3' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::MySQL' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::GBrowse' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::Diff' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::Foswiki' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::BLAST' ) || print "Bail out!\n";
    use_ok( 'Bio::BioStudio::Cairo' ) || print "Bail out!\n";
}

diag( "Testing Bio::BioStudio $Bio::BioStudio::Basic::VERSION, Perl $], $^X" );
