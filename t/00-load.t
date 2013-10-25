#!perl -T

use Test::More tests => 21;

BEGIN {
    use_ok('Bio::BioStudio')                 || print "Can't use BioStudio!\n";
    use_ok('Bio::BioStudio::Cairo')          || print "Can't use Cairo!\n";
    use_ok('Bio::BioStudio::Chromosome')     || print "Can't use Chromosome!\n";
    use_ok('Bio::BioStudio::Chunk')          || print "Can't use Chunk!\n";
    use_ok('Bio::BioStudio::DB')             || print "Can't use DB!\n";
    use_ok('Bio::BioStudio::Diff')           || print "Can't use Diff!\n";
    use_ok('Bio::BioStudio::Exceptions')     || print "Can't use Exceptions!\n";
    use_ok('Bio::BioStudio::Foswiki')        || print "Can't use Foswiki!\n";
    use_ok('Bio::BioStudio::GBrowse')        || print "Can't use GBrowse!\n";
    use_ok('Bio::BioStudio::Git')            || print "Can't use Git!\n";
    use_ok('Bio::BioStudio::Marker')         || print "Can't use Marker!\n";
    use_ok('Bio::BioStudio::Mask')           || print "Can't use Mask!\n";
    use_ok('Bio::BioStudio::Megachunk')      || print "Can't use Megachunk!\n";
    use_ok('Bio::BioStudio::Repository')     || print "Can't use Repository!\n";

    use_ok('Bio::BioStudio::Diff::Difference')
      || print "Can't use Diff::Difference!\n";
      
    use_ok('Bio::BioStudio::RestrictionEnzyme')
      || print "Can't use RestrictionEnzyme!\n";
    
    use_ok('Bio::BioStudio::RestrictionEnzyme::Store')
      || print "Can't use RestrictionEnzyme::Store!\n";
    
    use_ok('Bio::BioStudio::RestrictionEnzyme::Seek')
      || print "Can't use RestrictionEnzyme::Seek!\n";

    use_ok( 'Bio::BioStudio::Analyze::ArbitraryFeatures' )
      || print "Can't use ArbitraryFeatures\n";

    use_ok( 'Bio::BioStudio::Analyze::RestrictionEnzymes' )
      || print "Can't use RestrictionEnzymes\n";

    use_ok( 'Bio::BioStudio::Analyze::ProteinCodingGenes' )
      || print "Can't use ProteinCodingGenes\n";

}
