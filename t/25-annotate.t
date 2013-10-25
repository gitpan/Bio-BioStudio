#! /usr/bin/perl

use Test::More tests => 2;
use Test::Deep;

use Bio::BioStudio;

use strict;
use warnings;

my $BS = Bio::BioStudio->new();
$BS->genome_repository();
$BS->db_engine();

my $achr = Bio::BioStudio::Chromosome->new(
  -name       => 'test_chr01_0_00',
  -repo       => 't/test_repo/',
  -db_engine  => 'memory',
  -gbrowse     => 0
);
my $bchr = $achr->iterate;


#Testing add_feature, no new seq
my $nname = 'test_tag';
my $nfeat = Bio::SeqFeature::Generic->new(
  -start         => 2680,
  -end           => 2700,
  -primary_tag   => $nname,
  -display_name  => $nname,
);
my $rchr = bless
{
    'chromosome_id' => '01',
    'sequence' =>
'TATGGGTACCACAAAGCGAGGTCGCTTTTGAAGAGCCCTCGGTAGCATAACATTTTTAATTATTACGACTGTTTTTTTTATTCATTATGTAGAGATAATTAAATGTTATAGATGCTCTATACTCAAACGGTGGAAGAAAAACAGCGAAAAAAAATAACCGATACCCCCTTTTCGAATACAAATGCTTGTATATTCAATTATGAATTATTTTTTTTTTTTTTCATTTCTTATATTATTTTTTGTTCGAGAATCACTTTTTCAAGATGGTAACAACATCTTCGTCTTCCAAAATGTGACTCAACCCCACGTATTGAGGTTGATGTTTGACACTGCTACCGTAAACCAGAGCATTTCTAAAGTCGTCCACTAAAGATTTATGAATTTGGTTACAAAAATCCTTGACACTGCAACGGTCTGATCTTAGCACCACAGGGTCGGTAAAATCTGGTATTTGGCCCTTTGGTTTAGTGTAAATACGGACTAGATTTAGTCTATCCCACATGACTTGCAACAGCTCGTCCAAGTTCCAATCTTGACCAGACGAAATAGGCACGGCATTAGGAATTCGGTAAAGTAATTCCAATTCCTCTATTGACAGAGAATCAATCTTGTTTAACACATAGATGGCAGGCATGTATCTTCTTGACGAAGCTTCCAAAACATCAATCAAATCATCCACAGTGGCATCACACCTGAAGGCAATCTCAGCGCTATTTATTCTGTACTCGCTCATAACGGCTCTGATTTCGTCATTCCCCAGATGGGTCAATGGGACTGTGTTTGTGATGGAAATACCACCTTTCTCTTTTTTTTTGATCAAGATATCTGGCGGAGTTTTATTCAGACGAATCCCCACACCTTCCAGTTCCTTCTCAATGATTTGCTTATGATGCAAGGGTTTGTTCACATCTAGGATGATAAATAACAGGTTACAGGTTCTTGCCACGGCAATAACTTGCTTACCTCTACCTCTACCATCCTTAGCACCATCGATAATACCAGGTAAATCCAACATTTGGATCTTGGCACCTTTATAACGAATGACACCGGGGACGGTAACCAGGGTGGTAAACTCGTACTCAGCTGCTTCAGACTCAGTACCAGTCAACTTGGACAGTAATGTAGATTTCCCCACCGACGGGAACCCGACAAACCCCACACTGGCCACACCAGTTCTAGCCACATCAAAACCAATACCAGCACCACCACCGCTGCCGGATGAAGCACTGGTCAACAATTCTCTTCTCAGTTTGGCCAGCTTGGCCTTCAGTTGACCCAAATGGAAAGATGTGGCCTTGTTCTTTTGGGTACGGGCCATTTCATCTTCGATAGCTTTGATTTTTTCAACTGTAGTAGACATTTTTGCTCAATCAACAACTCTACGCTTGCACCTACTGCATCTAGCTTCAAACACTTCCTATCATTGCGCCCTCATCACACCGTAATATCCCATCTTAAAAGTGGAAAACTCTTATAGCTCATCGATGAAAAAAACGGGCCCTCGTCGCTTGTGATGTGAAAAAATTTTTCAAGCTTTAAGCCCATTGAAAGCAAGAGATCTTGCACTAGAATAAGTGGCAAAGGTGAACTTTGAGGGGATAAGAAGGGCATCTCCTCCGGAGTCATTGCCATCTGCGTTGAGTACCAAAGCTTTGAGCCCGTCAGAATCCTTGGCCACCGGACATGCTTCACAGATATAGAACGTAGCATGGTCTGTGGGAGCTTCATTTCTATGTTTTACCTTCTCTTTTCGCTTTTATGGTTCTCAGTGACCAAATAAAGAAACTTATATATGTTCCGGAATGACGAATCAAAAAGAGAATAGCATCGTTAGCAGCAAACGAAAGTGGAAAGAGAATAATGTTCAAGAGAGCAATGAGCACAGATGGTCCCGTGGCACGTACCATCCTGAAGAGACTGGAATGCGGCTTTCCAGATTACAAGAACTTTGCGTTTGGCCTCTACAACGATTCTCACAAGCATAAGGGCCATGCTGGTGTACAGGGAAATGTCTCTGCTGAGACACATTTCCGGATTGAGATGGTCAGTAAAAAGTTCGAAGGCCTGAAACTTCCACAACGCCATCGTATGGTTTATTCCCTCTTGCAAGACGAGATGGCTCAGGCGAACGGTATCCATGCTTTACAATTGTCACTAAAGACCCCACAGGAGTATGAATCCAAAGCGAAATAGGAATGCATAAGCATAAGTGTACACGTTGAGTTTATTGTTTTATTTCCCCTACATATATATACATATATATGAAATTACTTTACGTACGTATAAGCTTTGTTCAGTCATCATGAACCAGTGTCTTTTCGTACTGTTCTAAGGACATTAGACCCTCGACCTGTTCCACATTAACGCCCTCACCAAGCTTCATTTTGACTAGCCAGCCGTCACCCATAGGATCTTCGTTCACCACACCTGGATTTTCCTCAAGTATTATATATTTTGCTTTATGATCCTCGGCTTGATGCTCGCCAACGTGAGATAGCTGGTCATCACAATAGATCAGCCGGGACGCTTTTCGATCACATCGAATCCCTTCGGGACGTTGCAACAATACGTGAAAAATGCCTCAAAAATAATAAATACAATGGTGAACAACGTTAAAAAAGCATAAAACAGCTGGCTATTTTGATCAGGATAACATCTATAAGTGCCATATTAAGGCAAGATATCAATTGACCATGCAAACACCTTCAGAAAATACCGACGTAAAGATGGATACTCTCGACGAACCCAGTGCACATTTAATCGAAGAGAATGTAGCTCTTCCCGAAGACACATTCAGTTCACATCTGAGTTATGTACTTTATGAAATTGCTCATTGTAAACCGATCATGTTTATGATCATCATAATCGTGAGTTTGATCTCATTGATTGTGCTTTTTCATGATAACGACGGGTGCACTGTGATCTTAGTGATGTCCCTTATAGTAGCCTCCATGGCTTTAATGGTGGTTGCAGCATTCACATTCGGGAAAGCGATCACTGAACAGGAGTTCATGATAAAGCTTTTAGTGGAGGTGATCGCACGCAAGCCTGCGGGGAAGGAATGGGGTACTGTCGCATATAATATGAACCAATATCTATTCATGAAGAGACTATGGTATACCCCGTACTATTTCTATAGCGGCAAGAAGTGCCATGAGTTCTTCACCACTCTTATCAAGGAAGTGAATTCTGGTTCGCACTCGGATTCCTCATCGAATAGTGCCGAGGATACACAATCACCTGTCTCAGCAGGGAAGACTTCAAATGGTCTAAACAACTTTTATAGTATTAGATCAGACCCTATTTTGATGGCATATGTTTTGAAGGCAACACAAATAGAAAAGGAGGCTCAAAGTGAATACTGGAGAAAGCAATATCCTGACGCTGATTTACCTTGAAAGCGGAAGCATTTTATTCACCAAGTATACTTACTTTTCTTTAAAACGAGAACAAGAATCGAATTCAAGAACATCTCGAAGCCAGAATTGAGCATCATATATTCGAGCTGTACAAACATCATGGCCTACAACTATCGTATTTGTAAGTTTTTTTAGAGGTTTTCATATTTGTTTAATAAGGGTTCTGTCAGTTTTTGTCACATTCTATTGTTGCGCTTCGCATAATGCAGCCAAGAAAATCCAAACAATACCTTTCTACATACTACTACATAATATATATATATAGTATAGAAATTGGTATATCACTACTTGTACAAATATCATATTGTACGATAATCGCGAAGAACGACGCACTGGTGGGAAGAAGTGGAAAACAGAAGCTTTAAGGTAGAAACAGAACAAGAATGTGGCTATGGTAGGATAGCAAAAGAGTACCATTGCTGTTATCATTTGTTGCCTAGCCCTATCAAGACCTGTCTGCTAATCCAACCCGAGAGATCATGGCGATCCAAACCCGTTTTGCCTCGGGCACATCTTTATCCGATTTGAAACCAAAACCAAGTGCAACTTCCATCTCCATACCCATGCAAAATGTCATGAACAAGCCTGTCACGGAACAGGACTCACTGTTCCATATATGCGCAAACATCCGGAAAAGACTGGAGGTGTTACCTCAACTCAAACCTTTTTTACAATTGGCCTACCAATCGAGCGAGGTTTTGAGTGAAAGGCAATCTCTTTTGCTATCCCAAAAGCAGCATCAGGAACTGCTCAAGTCCAATGGCGCTAACCGGGACAGTAGCGACTTGGCACCAACTTTAAGGTCTAGCTCTATCTCCACAGCTACCAGTCTCATGTCGATGGAAGGTATATCATACACGAATTCGAATCCCTCGGCCACCCCAAATATGGAGGACACTTTACTGACTTTTAGTATGGGTATTTTGCCCATTACCATGGATTGCGACCCTGTGACACAACTATCACAGCTGTTTCAACAAGGTGCGCCCCTCTGTATACTTTTCAACTCTGTGAAGCCGCAATTTAAATTACCGGTAATAGCATCTGACGATTTGAAAGTCTGTAAAAAATCCATTTATGACTTTATATTGGGCTGCAAGAAACACTTTGCATTTAACGATGAGGAGCTTTTCACTATATCCGACGTTTTTGCCAACTCTACTTCCCAGCTGGTCAAAGTGCTAGAAGTAGTAGAAACGCTAATGAATTCCAGCCCTACTATTTTCCCCTCTAAGAGTAAGACACAGCAAATCATGAACGCAGAAAACCAACACCGACATCAGCCTCAGCAGTCTTCGAAGAAGCATAACGAGTATGTTAAAATTATCAAGGAATTCGTTGCAACGGAAAGAAAATATGTTCACGATTTGGAAATTTTGGATAAATATAGACAGCAGTTATTAGACAGCAATCTAATAACGTCTGAAGAGTTGTACATGTTGTTCCCTAATTTGGGTGATGCTATAGATTTTCAAAGAAGATTTCTAATATCCTTGGAAATAAATGCTTTAGTAGAACCTTCCAAGCAAAGAATCGGGGCTCTTTTCATGCATTCCAAACATTTTTTTAAGTTGTATGAGCCTTGGTCTATTGGCCAAAATGCAGCCATCGAATTTCTCTCTTCAACTTTGCACAAGATGAGGGTTGATGAATCGCAGCGGTTCATAATTAACAATAAACTGGAATTGCAATCCTTCCTTTATAAACCCGTGCAAAGGCTTTGTAGATATCCCCTGTTGGTCAAAGAATTGCTTGCTGAATCGAGTGACGATAATAATACGAAAGAACTTGAAGCTGCTTTAGATATTTCTAAAAATATTGCGAGAAGTATCAACGAAAATCAAAGAAGAACAGAAAATCATCAAGTGGTGAAGAAACTTTATGGTAGAGTGGTCAACTGGAAGGGTTATAGAATTTCCAAGTTCGGTGAGTTATTATATTTCGATAAAGTGTTCATTTCAACAACAAATAGCTCCTCGGAACCTGAAAGAGAATTTGAGGTTTATCTTTTTGAAAAAATCATCATCCTTTTTTCAGAGGTAGTGACTAAGAAATCTGCATCATCACTAATCCTTAAGAAGAAATCCTCAACCTCAGCATCAATCTCCGCCTCGAACATAACGGACAACAATGGCAGCCCTCACCACAGTTACCATAAGAGGCATAGCAATAGTAGTAGCAGTAATAATATCCATTTATCTTCGTCTTCAGCAGCGGCGATAATACATTCCAGTACCAATAGTAGTGACAACAATTCCAACAATTCATCATCATCCTCATTATTCAAGCTGTCCGCTAACGAACCTAAGCTGGATCTAAGAGGTCGAATTATGATAATGAATCTGAATCAAATCATACCGCAAAACAACCGGTCATTAAATATAACATGGGAATCCATAAAAGAGCAAGGTAATTTCCTTTTGAAATTCAAAAATGAGGAAACAAGAGATAATTGGTCATCGTGTTTACAACAGTTGATTCATGATCTGAAAAATGAGCAGTTTAAGGCAAGACATCACTCTTCAACATCGACGACTTCATCGACAGCCAAATCATCTTCAATGATGTCACCCACCACAACTATGAATACACCGAATCATCACAACAGCCGCCAGACACACGATAGTATGGCTTCTTTCTCAAGTTCTCATATGAAAAGGGTTTCGGATGTCCTGCCTAAACGGAGGACCACTTCATCAAGTTTCGAAAGTGAAATTAAATCCATTTCAGAAAATTTCAAGAACTCTATTCCAGAATCTTCCATACTCTTCAGGATATCATATAATAACAACTCTAATAATACCTCTAGTAGCGAGATCTTCACACTTTTGGTAGAAAAAGTTTGGAATTTTGACGACTTGATAATGGCGATCAATTCTAAAATTTCGAATACACATAATAACAACATTTCACCAATCACCAAGATCAAATATCAGGACGAAGATGGGGATTTTGTTGTGTTAGGTAGCGATGAAGATTGGAATGTTGCTAAAGAAATGTTGGCGGAAAACAATGAGAAATTCTTGAACATTCGTCTGTATTGAATAAATAAAACTAGTATACAGCAAATACTAAATAATTCAAGAAAAAAACATTAGATAGAGAGGGGCAGATGTTCAAGCTATACCCATTATATTGATCCACACTTAGTATTAAGATACGTCTGTGAAGGATGAAAAAAAATGTATAATGTGACTAGAGGAAGTAAGGAGAAAAAACGATAGTAATCGTATTTTAGGTTGTGCGTTTTTATAATTTTTTTTTTTTTGTAATTCTATGCAAATGTAATATAAGGGTCAGTAAAAAGTTCGAAGGCCTGAAACTTCCACAACGCCATCGTATGGTTTATTCCCTCTTGCAAGACGAGATGGCTCAGGCGAACGGTATCCATGCTTTACAATTGTCACTAAAGACCCCACAGGAGTATGAATCCAAAGCGAAATAGAATGCATAAGCATAAGTGTACACGTTGAGTTTATTGTTTTATTTCCCCTACATATATATACATATATATGAAATTACTTTACGTACGTATAAGCTTTGTTCAGTCATCATGAACCAGTGTCTTTTCGTACTGTTCTAAGGACATTAGACCCTCGACCTGTTCCACATTAACGCCCTCACCAAGCTTCATTTTGACTAGCCAGCCGTCACCCATAGGATCTTCGTTCACCACACCTGGATTTTCCTCAAGATTAGTGTTAATTTCCTCTACGGTACCATCGGCAGGCTGGTAGATCTCGGAGGCTGACTTGACGGACTCAATGGACCCTAGCGACTCACCTTGGGAAATCTCAGTGCCCACTTCTGGCAACTCAACATAGGTAGCGTCCCCTAAGGAATCAGTGGCGTATTTTGTAATTCCGACAAAGGCAGTCTTGTCCTGATGCACAGCTATCCACTCATGTTGGGAAGTGTACCTCACGGCTTGAGGTCCTTGGGATGAGTACAAAAATGGTAGTTTATTCTTGTTTAGGGCATTGCCGGAGCTGTTTCTCAAAAACAATTTGCTCACAGCGGGCATGCGGGTGGTCCATAGTCTAGTAGTGCGTAACATTTGTCGATGTGGTATGCTTCATGTGGAGATTCCCTTTCCCATTAGATACTTGTTTGTTGGTCTGTATATATAGAAGAAAGAGTTAGCGAAAGTGACTCCGCCGCTGAATGACTCCTTACGGAAGTGTCAAAATTGCGAGGTCCCTATAGCACAGAATGATAGATAAAACATTGATTTGCAAGTTGAAGGAAGACCCTACACATGCGTATATATGATGTATGTAATGGTTGTGATCATTTTAGCCTGTCAA',
    'chromosome_version' => '01',
    'comments'           => [],
    'provisional'        => -1,
    'repo' => 't/test_repo/',
    'GD'                 => ignore(),
    'tag'            => undef,
    'mask'           => ignore(),
    'gbrowse'        => 0,
    'seq_id'         => 'chr01',
    'species'        => 'test',
    'db_engine'      => 'memory',
    'genome_version' => '0',
    'database'       => bless(
        {
            'fasta_file' => ignore(),
            '_index'     => {
                'attribute' => {
                    'essential_status' => {
                        'essential' => {
                            '5' => undef
                        },
                        'nonessential' => {
                            '6' => undef,
                            '4' => undef,
                            '1' => undef,
                            '2' => undef
                        }
                    },
                    'intro' => {
                        '0_01' => {
                            '3' => undef
                        }
                    },
                    'version' => {
                        '1' => {
                            '3' => undef
                        }
                    },
                    'wtpos' => {
                        '2680' => {
                            '3' => undef
                        }
                    },
                    'gene' => {
                        'cdc24' => {
                            '5' => undef
                        },
                        'gcv3' => {
                            '6' => undef
                        },
                        'rbg1' => {
                            '1' => undef
                        }
                    },
                    'load_id' => {
                        'yal036c' => {
                            '1' => undef
                        },
                        'test_tag' => {
                            '3' => undef
                        },
                        'yar028w' => {
                            '4' => undef
                        },
                        'yal044w-a' => {
                            '2' => undef
                        },
                        'yal044c' => {
                            '6' => undef
                        },
                        'yal041w' => {
                            '5' => undef
                        }
                    },
                    'dbxref' => {
                        'sgd:s000000034' => {
                            '1' => undef
                        },
                        'sgd:s000007586' => {
                            '2' => undef
                        },
                        'sgd:s000000039' => {
                            '5' => undef
                        },
                        'sgd:s000000042' => {
                            '6' => undef
                        },
                        'sgd:s000000076' => {
                            '4' => undef
                        }
                    },
                    'note' => {
'10-methylene-thf; expression is regulated by levels of levels of 5'
                          => {
                            '6' => undef
                          },
                        ' putative dna repair protein' => {
                            '2' => undef
                        },
'h subunit of the mitochondrial glycine decarboxylase complex'
                          => {
                            '6' => undef
                          },
' and mutants have morphological defects in bud formation and shmooing'
                          => {
                            '5' => undef
                          },
                        ' required for the catabolism of glycine to 5' => {
                            '6' => undef
                        },
'member of the drg family of gtp-binding proteins; interacts with translating ribosomes and with tma46p'
                          => {
                            '1' => undef
                          },
                        'putative integral membrane protein' => {
                            '4' => undef
                        },
                        'similar to pombe uvi31' => {
                            '2' => undef
                        },
                        ' member of dup240 gene family' => {
                            '4' => undef
                        },
'guanine nucleotide exchange factor (gef or gdp-release factor) for cdc42p; required for polarity establishment and maintenance'
                          => {
                            '5' => undef
                          },
                        '10-methylene-thf in the cytoplasm' => {
                            '6' => undef
                        }
                    },
                    'orf_classification' => {
                        'uncharacterized' => {
                            '4' => undef,
                            '2' => undef
                        },
                        'verified' => {
                            '6' => undef,
                            '1' => undef,
                            '5' => undef
                        }
                    },
                    'ontology_term' => {
                        'go:0006546' => {
                            '6' => undef
                        },
                        'go:0001403' => {
                            '5' => undef
                        },
                        'go:0030468' => {
                            '5' => undef
                        },
                        'go:0000723' => {
                            '6' => undef
                        },
                        'go:0004375' => {
                            '6' => undef
                        },
                        'go:0007119' => {
                            '5' => undef
                        },
                        'go:0000004' => {
                            '4' => undef,
                            '1' => undef,
                            '2' => undef
                        },
                        'go:0007264' => {
                            '5' => undef
                        },
                        'go:0005739' => {
                            '6' => undef
                        },
                        'go:0005960' => {
                            '6' => undef
                        },
                        'go:0005634' => {
                            '5' => undef
                        },
                        'go:0005737' => {
                            '1' => undef
                        },
                        'go:0008372' => {
                            '4' => undef,
                            '2' => undef
                        },
                        'go:0043332' => {
                            '5' => undef
                        },
                        'go:0007124' => {
                            '5' => undef
                        },
                        'go:0005525' => {
                            '1' => undef
                        },
                        'go:0000131' => {
                            '5' => undef
                        },
                        'go:0006033' => {
                            '5' => undef
                        },
                        'go:0006730' => {
                            '6' => undef
                        },
                        'go:0000750' => {
                            '5' => undef
                        },
                        'go:0004871' => {
                            '5' => undef
                        },
                        'go:0005089' => {
                            '5' => undef
                        },
                        'go:0005554' => {
                            '4' => undef,
                            '2' => undef
                        },
                        'go:0007096' => {
                            '5' => undef
                        },
                        'go:0000753' => {
                            '5' => undef
                        },
                        'go:0007118' => {
                            '5' => undef
                        },
                        'go:0005935' => {
                            '5' => undef
                        }
                    },
                    'alias' => {
                        'cdc24' => {
                            '5' => undef
                        },
                        'gcv3' => {
                            '6' => undef
                        },
                        'cls4' => {
                            '5' => undef
                        },
                        'fun11' => {
                            '1' => undef
                        },
                        'rbg1' => {
                            '1' => undef
                        }
                    }
                },
                'location' => {
                    'chr01' => {
                        '0' => {
                            '6' => undef,
                            '4' => undef,
                            '1' => undef,
                            '3' => undef,
                            '2' => undef,
                            '5' => undef
                        }
                    }
                },
                'ids' => {
                    '6' => undef,
                    '4' => undef,
                    '1' => undef,
                    '3' => undef,
                    '2' => undef,
                    '5' => undef
                },
                'name' => {
                    'yal036c' => {
                        '1' => 1
                    },
                    'test_tag' => {
                        '3' => 1
                    },
                    'cdc24' => {
                        '5' => 2
                    },
                    'gcv3' => {
                        '6' => 2
                    },
                    'yal044w-a' => {
                        '2' => 1
                    },
                    'rbg1' => {
                        '1' => 2
                    },
                    'yar028w' => {
                        '4' => 1
                    },
                    'cls4' => {
                        '5' => 2
                    },
                    'fun11' => {
                        '1' => 2
                    },
                    'yal044c' => {
                        '6' => 1
                    },
                    'yal041w' => {
                        '5' => 1
                    }
                },
                'type' => {
                    'test_tag' => {
                        'bio' => {
                            '3' => undef
                        }
                    },
                    'gene' => {
                        'sgd' => {
                            '6' => undef,
                            '4' => undef,
                            '1' => undef,
                            '2' => undef,
                            '5' => undef
                        }
                    }
                }
            },
            '_data' => {
                '6' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 6,
                        'store'       => ignore(),
                        'stop'        => 7476,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL044C',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => -1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Nonessential' ],
                            'load_id'          => [ 'YAL044C' ],
                            'gene'             => [ 'GCV3' ],
                            'Note'             => [
'H subunit of the mitochondrial glycine decarboxylase complex',
                                ' required for the catabolism of glycine to 5',
'10-methylene-THF; expression is regulated by levels of levels of 5',
                                '10-methylene-THF in the cytoplasm'
                            ],
                            'dbxref'             => [ 'SGD:S000000042' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'GCV3' ],
                            'Ontology_term'      => [
                                'GO:0006730', 'GO:0006546',
                                'GO:0005960', 'GO:0005739',
                                'GO:0004375', 'GO:0000723'
                            ]
                        },
                        'start' => 6964
                    },
                    'Bio::DB::SeqFeature'
                ),
                '4' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 4,
                        'store'       => ignore(),
                        'stop'        => 3398,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAR028W',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'orf_classification' => [ 'Uncharacterized' ],
                            'essential_status'   => [ 'Nonessential' ],
                            'load_id'            => [ 'YAR028W' ],
                            'Ontology_term' =>
                              [ 'GO:0008372', 'GO:0005554', 'GO:0000004' ],
                            'Note' => [
                                'Putative integral membrane protein',
                                ' member of DUP240 gene family'
                            ],
                            'dbxref' => [ 'SGD:S000000076' ]
                        },
                        'start' => 2694
                    },
                    'Bio::DB::SeqFeature'
                ),
                '1' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 1,
                        'store'       => ignore(),
                        'stop'        => 1360,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL036C',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => -1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Nonessential' ],
                            'load_id'          => [ 'YAL036C' ],
                            'gene'             => [ 'RBG1' ],
                            'Note'             => [
'Member of the DRG family of GTP-binding proteins; interacts with translating ribosomes and with Tma46p'
                            ],
                            'dbxref'             => [ 'SGD:S000000034' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'RBG1', 'FUN11' ],
                            'Ontology_term' =>
                              [ 'GO:0005737', 'GO:0005525', 'GO:0000004' ]
                        },
                        'start' => 251
                    },
                    'Bio::DB::SeqFeature'
                ),
                '3' => bless(
                    {
                        'source'      => 'BIO',
                        'primary_id'  => 3,
                        'store'       => ignore(),
                        'stop'        => 2700,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'test_tag',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 0,
                        'type'        => 'test_tag',
                        'attributes'  => {
                            'version' => [ '1' ],
                            'intro'   => [ '0_01' ],
                            'load_id' => [ 'test_tag' ],
                            'wtpos'   => [ '2680' ]
                        },
                        'start' => 2680
                    },
                    'Bio::DB::SeqFeature'
                ),
                '2' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 2,
                        'store'       => ignore(),
                        'stop'        => 2193,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL044W-A',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'orf_classification' => [ 'Uncharacterized' ],
                            'essential_status'   => [ 'Nonessential' ],
                            'load_id'            => [ 'YAL044W-A' ],
                            'Ontology_term' =>
                              [ 'GO:0008372', 'GO:0005554', 'GO:0000004' ],
                            'Note' => [
                                'Similar to pombe uvi31',
                                ' putative DNA repair protein'
                            ],
                            'dbxref' => [ 'SGD:S000007586' ]
                        },
                        'start' => 1861
                    },
                    'Bio::DB::SeqFeature'
                ),
                '5' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 5,
                        'store'       => ignore(),
                        'stop'        => 6463,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL041W',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Essential' ],
                            'load_id'          => [ 'YAL041W' ],
                            'gene'             => [ 'CDC24' ],
                            'Note'             => [
'Guanine nucleotide exchange factor (GEF or GDP-release factor) for Cdc42p; required for polarity establishment and maintenance',
' and mutants have morphological defects in bud formation and shmooing'
                            ],
                            'dbxref'             => [ 'SGD:S000000039' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'CDC24', 'CLS4' ],
                            'Ontology_term'      => [
                                'GO:0043332', 'GO:0030468',
                                'GO:0000753', 'GO:0000750',
                                'GO:0007264', 'GO:0007124',
                                'GO:0007119', 'GO:0007118',
                                'GO:0007096', 'GO:0006033',
                                'GO:0005935', 'GO:0005634',
                                'GO:0005089', 'GO:0004871',
                                'GO:0001403', 'GO:0000131'
                            ]
                        },
                        'start' => 3899
                    },
                    'Bio::DB::SeqFeature'
                )
            },
            '_children' => {},
            'setting'   => {
                'serializer'        => 'Storable',
                'compress'          => undef,
                'index_subfeatures' => 1
            },
            'fasta_fh'        => ignore(),
            'seqfeatureclass' => 'Bio::DB::SeqFeature',
            'fasta_db'        => ignore(),
        },
        'Bio::DB::SeqFeature::Store::memory'
    ),
    '_root_verbose' => 0
}, 'Bio::BioStudio::Chromosome';

$bchr->add_feature(
    -feature => $nfeat
);
$bchr->write_chromosome();
cmp_deeply( $bchr, $rchr, "add without sequence change" );

my $cchr = $bchr->iterate;


#Testing add_feature, new seq
my $yname = 'test_tag';
my $yfeat = Bio::SeqFeature::Generic->new(
  -start         => 2680,
  -end           => 2700,
  -primary_tag   => $yname,
  -display_name  => $yname,
  -tag           => {
    newseq => 'GATATCAATTGACCATGCATA'
  },
);
$rchr = bless
{
    'chromosome_id' => '01',
    'sequence' =>
'TATGGGTACCACAAAGCGAGGTCGCTTTTGAAGAGCCCTCGGTAGCATAACATTTTTAATTATTACGACTGTTTTTTTTATTCATTATGTAGAGATAATTAAATGTTATAGATGCTCTATACTCAAACGGTGGAAGAAAAACAGCGAAAAAAAATAACCGATACCCCCTTTTCGAATACAAATGCTTGTATATTCAATTATGAATTATTTTTTTTTTTTTTCATTTCTTATATTATTTTTTGTTCGAGAATCACTTTTTCAAGATGGTAACAACATCTTCGTCTTCCAAAATGTGACTCAACCCCACGTATTGAGGTTGATGTTTGACACTGCTACCGTAAACCAGAGCATTTCTAAAGTCGTCCACTAAAGATTTATGAATTTGGTTACAAAAATCCTTGACACTGCAACGGTCTGATCTTAGCACCACAGGGTCGGTAAAATCTGGTATTTGGCCCTTTGGTTTAGTGTAAATACGGACTAGATTTAGTCTATCCCACATGACTTGCAACAGCTCGTCCAAGTTCCAATCTTGACCAGACGAAATAGGCACGGCATTAGGAATTCGGTAAAGTAATTCCAATTCCTCTATTGACAGAGAATCAATCTTGTTTAACACATAGATGGCAGGCATGTATCTTCTTGACGAAGCTTCCAAAACATCAATCAAATCATCCACAGTGGCATCACACCTGAAGGCAATCTCAGCGCTATTTATTCTGTACTCGCTCATAACGGCTCTGATTTCGTCATTCCCCAGATGGGTCAATGGGACTGTGTTTGTGATGGAAATACCACCTTTCTCTTTTTTTTTGATCAAGATATCTGGCGGAGTTTTATTCAGACGAATCCCCACACCTTCCAGTTCCTTCTCAATGATTTGCTTATGATGCAAGGGTTTGTTCACATCTAGGATGATAAATAACAGGTTACAGGTTCTTGCCACGGCAATAACTTGCTTACCTCTACCTCTACCATCCTTAGCACCATCGATAATACCAGGTAAATCCAACATTTGGATCTTGGCACCTTTATAACGAATGACACCGGGGACGGTAACCAGGGTGGTAAACTCGTACTCAGCTGCTTCAGACTCAGTACCAGTCAACTTGGACAGTAATGTAGATTTCCCCACCGACGGGAACCCGACAAACCCCACACTGGCCACACCAGTTCTAGCCACATCAAAACCAATACCAGCACCACCACCGCTGCCGGATGAAGCACTGGTCAACAATTCTCTTCTCAGTTTGGCCAGCTTGGCCTTCAGTTGACCCAAATGGAAAGATGTGGCCTTGTTCTTTTGGGTACGGGCCATTTCATCTTCGATAGCTTTGATTTTTTCAACTGTAGTAGACATTTTTGCTCAATCAACAACTCTACGCTTGCACCTACTGCATCTAGCTTCAAACACTTCCTATCATTGCGCCCTCATCACACCGTAATATCCCATCTTAAAAGTGGAAAACTCTTATAGCTCATCGATGAAAAAAACGGGCCCTCGTCGCTTGTGATGTGAAAAAATTTTTCAAGCTTTAAGCCCATTGAAAGCAAGAGATCTTGCACTAGAATAAGTGGCAAAGGTGAACTTTGAGGGGATAAGAAGGGCATCTCCTCCGGAGTCATTGCCATCTGCGTTGAGTACCAAAGCTTTGAGCCCGTCAGAATCCTTGGCCACCGGACATGCTTCACAGATATAGAACGTAGCATGGTCTGTGGGAGCTTCATTTCTATGTTTTACCTTCTCTTTTCGCTTTTATGGTTCTCAGTGACCAAATAAAGAAACTTATATATGTTCCGGAATGACGAATCAAAAAGAGAATAGCATCGTTAGCAGCAAACGAAAGTGGAAAGAGAATAATGTTCAAGAGAGCAATGAGCACAGATGGTCCCGTGGCACGTACCATCCTGAAGAGACTGGAATGCGGCTTTCCAGATTACAAGAACTTTGCGTTTGGCCTCTACAACGATTCTCACAAGCATAAGGGCCATGCTGGTGTACAGGGAAATGTCTCTGCTGAGACACATTTCCGGATTGAGATGGTCAGTAAAAAGTTCGAAGGCCTGAAACTTCCACAACGCCATCGTATGGTTTATTCCCTCTTGCAAGACGAGATGGCTCAGGCGAACGGTATCCATGCTTTACAATTGTCACTAAAGACCCCACAGGAGTATGAATCCAAAGCGAAATAGGAATGCATAAGCATAAGTGTACACGTTGAGTTTATTGTTTTATTTCCCCTACATATATATACATATATATGAAATTACTTTACGTACGTATAAGCTTTGTTCAGTCATCATGAACCAGTGTCTTTTCGTACTGTTCTAAGGACATTAGACCCTCGACCTGTTCCACATTAACGCCCTCACCAAGCTTCATTTTGACTAGCCAGCCGTCACCCATAGGATCTTCGTTCACCACACCTGGATTTTCCTCAAGTATTATATATTTTGCTTTATGATCCTCGGCTTGATGCTCGCCAACGTGAGATAGCTGGTCATCACAATAGATCAGCCGGGACGCTTTTCGATCACATCGAATCCCTTCGGGACGTTGCAACAATACGTGAAAAATGCCTCAAAAATAATAAATACAATGGTGAACAACGTTAAAAAAGCATAAAACAGCTGGCTATTTTGATCAGGATAACATCTATAAGTGCCATATTAAGGCAAGATATCAATTGACCATGCATACACCTTCAGAAAATACCGACGTAAAGATGGATACTCTCGACGAACCCAGTGCACATTTAATCGAAGAGAATGTAGCTCTTCCCGAAGACACATTCAGTTCACATCTGAGTTATGTACTTTATGAAATTGCTCATTGTAAACCGATCATGTTTATGATCATCATAATCGTGAGTTTGATCTCATTGATTGTGCTTTTTCATGATAACGACGGGTGCACTGTGATCTTAGTGATGTCCCTTATAGTAGCCTCCATGGCTTTAATGGTGGTTGCAGCATTCACATTCGGGAAAGCGATCACTGAACAGGAGTTCATGATAAAGCTTTTAGTGGAGGTGATCGCACGCAAGCCTGCGGGGAAGGAATGGGGTACTGTCGCATATAATATGAACCAATATCTATTCATGAAGAGACTATGGTATACCCCGTACTATTTCTATAGCGGCAAGAAGTGCCATGAGTTCTTCACCACTCTTATCAAGGAAGTGAATTCTGGTTCGCACTCGGATTCCTCATCGAATAGTGCCGAGGATACACAATCACCTGTCTCAGCAGGGAAGACTTCAAATGGTCTAAACAACTTTTATAGTATTAGATCAGACCCTATTTTGATGGCATATGTTTTGAAGGCAACACAAATAGAAAAGGAGGCTCAAAGTGAATACTGGAGAAAGCAATATCCTGACGCTGATTTACCTTGAAAGCGGAAGCATTTTATTCACCAAGTATACTTACTTTTCTTTAAAACGAGAACAAGAATCGAATTCAAGAACATCTCGAAGCCAGAATTGAGCATCATATATTCGAGCTGTACAAACATCATGGCCTACAACTATCGTATTTGTAAGTTTTTTTAGAGGTTTTCATATTTGTTTAATAAGGGTTCTGTCAGTTTTTGTCACATTCTATTGTTGCGCTTCGCATAATGCAGCCAAGAAAATCCAAACAATACCTTTCTACATACTACTACATAATATATATATATAGTATAGAAATTGGTATATCACTACTTGTACAAATATCATATTGTACGATAATCGCGAAGAACGACGCACTGGTGGGAAGAAGTGGAAAACAGAAGCTTTAAGGTAGAAACAGAACAAGAATGTGGCTATGGTAGGATAGCAAAAGAGTACCATTGCTGTTATCATTTGTTGCCTAGCCCTATCAAGACCTGTCTGCTAATCCAACCCGAGAGATCATGGCGATCCAAACCCGTTTTGCCTCGGGCACATCTTTATCCGATTTGAAACCAAAACCAAGTGCAACTTCCATCTCCATACCCATGCAAAATGTCATGAACAAGCCTGTCACGGAACAGGACTCACTGTTCCATATATGCGCAAACATCCGGAAAAGACTGGAGGTGTTACCTCAACTCAAACCTTTTTTACAATTGGCCTACCAATCGAGCGAGGTTTTGAGTGAAAGGCAATCTCTTTTGCTATCCCAAAAGCAGCATCAGGAACTGCTCAAGTCCAATGGCGCTAACCGGGACAGTAGCGACTTGGCACCAACTTTAAGGTCTAGCTCTATCTCCACAGCTACCAGTCTCATGTCGATGGAAGGTATATCATACACGAATTCGAATCCCTCGGCCACCCCAAATATGGAGGACACTTTACTGACTTTTAGTATGGGTATTTTGCCCATTACCATGGATTGCGACCCTGTGACACAACTATCACAGCTGTTTCAACAAGGTGCGCCCCTCTGTATACTTTTCAACTCTGTGAAGCCGCAATTTAAATTACCGGTAATAGCATCTGACGATTTGAAAGTCTGTAAAAAATCCATTTATGACTTTATATTGGGCTGCAAGAAACACTTTGCATTTAACGATGAGGAGCTTTTCACTATATCCGACGTTTTTGCCAACTCTACTTCCCAGCTGGTCAAAGTGCTAGAAGTAGTAGAAACGCTAATGAATTCCAGCCCTACTATTTTCCCCTCTAAGAGTAAGACACAGCAAATCATGAACGCAGAAAACCAACACCGACATCAGCCTCAGCAGTCTTCGAAGAAGCATAACGAGTATGTTAAAATTATCAAGGAATTCGTTGCAACGGAAAGAAAATATGTTCACGATTTGGAAATTTTGGATAAATATAGACAGCAGTTATTAGACAGCAATCTAATAACGTCTGAAGAGTTGTACATGTTGTTCCCTAATTTGGGTGATGCTATAGATTTTCAAAGAAGATTTCTAATATCCTTGGAAATAAATGCTTTAGTAGAACCTTCCAAGCAAAGAATCGGGGCTCTTTTCATGCATTCCAAACATTTTTTTAAGTTGTATGAGCCTTGGTCTATTGGCCAAAATGCAGCCATCGAATTTCTCTCTTCAACTTTGCACAAGATGAGGGTTGATGAATCGCAGCGGTTCATAATTAACAATAAACTGGAATTGCAATCCTTCCTTTATAAACCCGTGCAAAGGCTTTGTAGATATCCCCTGTTGGTCAAAGAATTGCTTGCTGAATCGAGTGACGATAATAATACGAAAGAACTTGAAGCTGCTTTAGATATTTCTAAAAATATTGCGAGAAGTATCAACGAAAATCAAAGAAGAACAGAAAATCATCAAGTGGTGAAGAAACTTTATGGTAGAGTGGTCAACTGGAAGGGTTATAGAATTTCCAAGTTCGGTGAGTTATTATATTTCGATAAAGTGTTCATTTCAACAACAAATAGCTCCTCGGAACCTGAAAGAGAATTTGAGGTTTATCTTTTTGAAAAAATCATCATCCTTTTTTCAGAGGTAGTGACTAAGAAATCTGCATCATCACTAATCCTTAAGAAGAAATCCTCAACCTCAGCATCAATCTCCGCCTCGAACATAACGGACAACAATGGCAGCCCTCACCACAGTTACCATAAGAGGCATAGCAATAGTAGTAGCAGTAATAATATCCATTTATCTTCGTCTTCAGCAGCGGCGATAATACATTCCAGTACCAATAGTAGTGACAACAATTCCAACAATTCATCATCATCCTCATTATTCAAGCTGTCCGCTAACGAACCTAAGCTGGATCTAAGAGGTCGAATTATGATAATGAATCTGAATCAAATCATACCGCAAAACAACCGGTCATTAAATATAACATGGGAATCCATAAAAGAGCAAGGTAATTTCCTTTTGAAATTCAAAAATGAGGAAACAAGAGATAATTGGTCATCGTGTTTACAACAGTTGATTCATGATCTGAAAAATGAGCAGTTTAAGGCAAGACATCACTCTTCAACATCGACGACTTCATCGACAGCCAAATCATCTTCAATGATGTCACCCACCACAACTATGAATACACCGAATCATCACAACAGCCGCCAGACACACGATAGTATGGCTTCTTTCTCAAGTTCTCATATGAAAAGGGTTTCGGATGTCCTGCCTAAACGGAGGACCACTTCATCAAGTTTCGAAAGTGAAATTAAATCCATTTCAGAAAATTTCAAGAACTCTATTCCAGAATCTTCCATACTCTTCAGGATATCATATAATAACAACTCTAATAATACCTCTAGTAGCGAGATCTTCACACTTTTGGTAGAAAAAGTTTGGAATTTTGACGACTTGATAATGGCGATCAATTCTAAAATTTCGAATACACATAATAACAACATTTCACCAATCACCAAGATCAAATATCAGGACGAAGATGGGGATTTTGTTGTGTTAGGTAGCGATGAAGATTGGAATGTTGCTAAAGAAATGTTGGCGGAAAACAATGAGAAATTCTTGAACATTCGTCTGTATTGAATAAATAAAACTAGTATACAGCAAATACTAAATAATTCAAGAAAAAAACATTAGATAGAGAGGGGCAGATGTTCAAGCTATACCCATTATATTGATCCACACTTAGTATTAAGATACGTCTGTGAAGGATGAAAAAAAATGTATAATGTGACTAGAGGAAGTAAGGAGAAAAAACGATAGTAATCGTATTTTAGGTTGTGCGTTTTTATAATTTTTTTTTTTTTGTAATTCTATGCAAATGTAATATAAGGGTCAGTAAAAAGTTCGAAGGCCTGAAACTTCCACAACGCCATCGTATGGTTTATTCCCTCTTGCAAGACGAGATGGCTCAGGCGAACGGTATCCATGCTTTACAATTGTCACTAAAGACCCCACAGGAGTATGAATCCAAAGCGAAATAGAATGCATAAGCATAAGTGTACACGTTGAGTTTATTGTTTTATTTCCCCTACATATATATACATATATATGAAATTACTTTACGTACGTATAAGCTTTGTTCAGTCATCATGAACCAGTGTCTTTTCGTACTGTTCTAAGGACATTAGACCCTCGACCTGTTCCACATTAACGCCCTCACCAAGCTTCATTTTGACTAGCCAGCCGTCACCCATAGGATCTTCGTTCACCACACCTGGATTTTCCTCAAGATTAGTGTTAATTTCCTCTACGGTACCATCGGCAGGCTGGTAGATCTCGGAGGCTGACTTGACGGACTCAATGGACCCTAGCGACTCACCTTGGGAAATCTCAGTGCCCACTTCTGGCAACTCAACATAGGTAGCGTCCCCTAAGGAATCAGTGGCGTATTTTGTAATTCCGACAAAGGCAGTCTTGTCCTGATGCACAGCTATCCACTCATGTTGGGAAGTGTACCTCACGGCTTGAGGTCCTTGGGATGAGTACAAAAATGGTAGTTTATTCTTGTTTAGGGCATTGCCGGAGCTGTTTCTCAAAAACAATTTGCTCACAGCGGGCATGCGGGTGGTCCATAGTCTAGTAGTGCGTAACATTTGTCGATGTGGTATGCTTCATGTGGAGATTCCCTTTCCCATTAGATACTTGTTTGTTGGTCTGTATATATAGAAGAAAGAGTTAGCGAAAGTGACTCCGCCGCTGAATGACTCCTTACGGAAGTGTCAAAATTGCGAGGTCCCTATAGCACAGAATGATAGATAAAACATTGATTTGCAAGTTGAAGGAAGACCCTACACATGCGTATATATGATGTATGTAATGGTTGTGATCATTTTAGCCTGTCAA',
    'chromosome_version' => '02',
    'comments'           => [],
    'provisional'        => -1,
    'repo' => 't/test_repo/',
    'GD'                 => ignore(),
    'tag'            => undef,
    'mask'           => ignore(),
    'gbrowse'        => 0,
    'seq_id'         => 'chr01',
    'species'        => 'test',
    'db_engine'      => 'memory',
    'genome_version' => '0',
    'database'       => bless(
        {
            'fasta_file' => ignore(),
            '_index'     => {
                'attribute' => {
                    'essential_status' => {
                        'essential' => {
                            '6' => undef
                        },
                        'nonessential' => {
                            '1' => undef,
                            '7' => undef,
                            '2' => undef,
                            '5' => undef
                        }
                    },
                    'newseq' => {
                        'gatatcaattgaccatgcata' => {
                            '4' => undef
                        }
                    },
                    'intro' => {
                        '0_01' => {
                            '3' => undef
                        },
                        '0_02' => {
                            '4' => undef
                        }
                    },
                    'version' => {
                        '1' => {
                            '4' => undef,
                            '3' => undef
                        }
                    },
                    'wtpos' => {
                        '2680' => {
                            '4' => undef,
                            '3' => undef
                        }
                    },
                    'gene' => {
                        'cdc24' => {
                            '6' => undef
                        },
                        'gcv3' => {
                            '7' => undef
                        },
                        'rbg1' => {
                            '1' => undef
                        }
                    },
                    'load_id' => {
                        'yal036c' => {
                            '1' => undef
                        },
                        'test_tag' => {
                            '4' => undef,
                            '3' => undef
                        },
                        'yar028w' => {
                            '5' => undef
                        },
                        'yal044w-a' => {
                            '2' => undef
                        },
                        'yal044c' => {
                            '7' => undef
                        },
                        'yal041w' => {
                            '6' => undef
                        }
                    },
                    'dbxref' => {
                        'sgd:s000000034' => {
                            '1' => undef
                        },
                        'sgd:s000007586' => {
                            '2' => undef
                        },
                        'sgd:s000000039' => {
                            '6' => undef
                        },
                        'sgd:s000000042' => {
                            '7' => undef
                        },
                        'sgd:s000000076' => {
                            '5' => undef
                        }
                    },
                    'note' => {
'10-methylene-thf; expression is regulated by levels of levels of 5'
                          => {
                            '7' => undef
                          },
                        ' putative dna repair protein' => {
                            '2' => undef
                        },
'h subunit of the mitochondrial glycine decarboxylase complex'
                          => {
                            '7' => undef
                          },
' and mutants have morphological defects in bud formation and shmooing'
                          => {
                            '6' => undef
                          },
                        ' required for the catabolism of glycine to 5' => {
                            '7' => undef
                        },
'member of the drg family of gtp-binding proteins; interacts with translating ribosomes and with tma46p'
                          => {
                            '1' => undef
                          },
                        'putative integral membrane protein' => {
                            '5' => undef
                        },
                        'similar to pombe uvi31' => {
                            '2' => undef
                        },
                        ' member of dup240 gene family' => {
                            '5' => undef
                        },
'guanine nucleotide exchange factor (gef or gdp-release factor) for cdc42p; required for polarity establishment and maintenance'
                          => {
                            '6' => undef
                          },
                        '10-methylene-thf in the cytoplasm' => {
                            '7' => undef
                        }
                    },
                    'orf_classification' => {
                        'uncharacterized' => {
                            '2' => undef,
                            '5' => undef
                        },
                        'verified' => {
                            '6' => undef,
                            '1' => undef,
                            '7' => undef
                        }
                    },
                    'ontology_term' => {
                        'go:0006546' => {
                            '7' => undef
                        },
                        'go:0001403' => {
                            '6' => undef
                        },
                        'go:0030468' => {
                            '6' => undef
                        },
                        'go:0000723' => {
                            '7' => undef
                        },
                        'go:0004375' => {
                            '7' => undef
                        },
                        'go:0007119' => {
                            '6' => undef
                        },
                        'go:0000004' => {
                            '1' => undef,
                            '2' => undef,
                            '5' => undef
                        },
                        'go:0007264' => {
                            '6' => undef
                        },
                        'go:0005739' => {
                            '7' => undef
                        },
                        'go:0005960' => {
                            '7' => undef
                        },
                        'go:0005634' => {
                            '6' => undef
                        },
                        'go:0005737' => {
                            '1' => undef
                        },
                        'go:0008372' => {
                            '2' => undef,
                            '5' => undef
                        },
                        'go:0043332' => {
                            '6' => undef
                        },
                        'go:0007124' => {
                            '6' => undef
                        },
                        'go:0005525' => {
                            '1' => undef
                        },
                        'go:0000131' => {
                            '6' => undef
                        },
                        'go:0006033' => {
                            '6' => undef
                        },
                        'go:0006730' => {
                            '7' => undef
                        },
                        'go:0000750' => {
                            '6' => undef
                        },
                        'go:0004871' => {
                            '6' => undef
                        },
                        'go:0005089' => {
                            '6' => undef
                        },
                        'go:0005554' => {
                            '2' => undef,
                            '5' => undef
                        },
                        'go:0007096' => {
                            '6' => undef
                        },
                        'go:0000753' => {
                            '6' => undef
                        },
                        'go:0007118' => {
                            '6' => undef
                        },
                        'go:0005935' => {
                            '6' => undef
                        }
                    },
                    'alias' => {
                        'cdc24' => {
                            '6' => undef
                        },
                        'gcv3' => {
                            '7' => undef
                        },
                        'cls4' => {
                            '6' => undef
                        },
                        'fun11' => {
                            '1' => undef
                        },
                        'rbg1' => {
                            '1' => undef
                        }
                    }
                },
                'location' => {
                    'chr01' => {
                        '0' => {
                            '6' => undef,
                            '4' => undef,
                            '1' => undef,
                            '3' => undef,
                            '7' => undef,
                            '2' => undef,
                            '5' => undef
                        }
                    }
                },
                'ids' => {
                    '6' => undef,
                    '4' => undef,
                    '1' => undef,
                    '3' => undef,
                    '7' => undef,
                    '2' => undef,
                    '5' => undef
                },
                'name' => {
                    'yal036c' => {
                        '1' => 1
                    },
                    'test_tag' => {
                        '4' => 1,
                        '3' => 1
                    },
                    'cdc24' => {
                        '6' => 2
                    },
                    'gcv3' => {
                        '7' => 2
                    },
                    'yal044w-a' => {
                        '2' => 1
                    },
                    'rbg1' => {
                        '1' => 2
                    },
                    'yar028w' => {
                        '5' => 1
                    },
                    'cls4' => {
                        '6' => 2
                    },
                    'fun11' => {
                        '1' => 2
                    },
                    'yal044c' => {
                        '7' => 1
                    },
                    'yal041w' => {
                        '6' => 1
                    }
                },
                'type' => {
                    'test_tag' => {
                        'bio' => {
                            '4' => undef,
                            '3' => undef
                        }
                    },
                    'gene' => {
                        'sgd' => {
                            '6' => undef,
                            '1' => undef,
                            '7' => undef,
                            '2' => undef,
                            '5' => undef
                        }
                    }
                }
            },
            '_data' => {
                '6' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 6,
                        'store'       => ignore(),
                        'stop'        => 6463,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL041W',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Essential' ],
                            'load_id'          => [ 'YAL041W' ],
                            'gene'             => [ 'CDC24' ],
                            'Note'             => [
'Guanine nucleotide exchange factor (GEF or GDP-release factor) for Cdc42p; required for polarity establishment and maintenance',
' and mutants have morphological defects in bud formation and shmooing'
                            ],
                            'dbxref'             => [ 'SGD:S000000039' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'CDC24', 'CLS4' ],
                            'Ontology_term'      => [
                                'GO:0043332', 'GO:0030468',
                                'GO:0000753', 'GO:0000750',
                                'GO:0007264', 'GO:0007124',
                                'GO:0007119', 'GO:0007118',
                                'GO:0007096', 'GO:0006033',
                                'GO:0005935', 'GO:0005634',
                                'GO:0005089', 'GO:0004871',
                                'GO:0001403', 'GO:0000131'
                            ]
                        },
                        'start' => 3899
                    },
                    'Bio::DB::SeqFeature'
                ),
                '4' => bless(
                    {
                        'source'      => 'BIO',
                        'primary_id'  => 4,
                        'store'       => ignore(),
                        'stop'        => 2700,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'test_tag',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 0,
                        'type'        => 'test_tag',
                        'attributes'  => {
                            'newseq'  => [ 'GATATCAATTGACCATGCATA' ],
                            'version' => [ '1' ],
                            'intro'   => [ '0_02' ],
                            'load_id' => [ 'test_tag' ],
                            'wtpos'   => [ '2680' ]
                        },
                        'start' => 2680
                    },
                    'Bio::DB::SeqFeature'
                ),
                '1' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 1,
                        'store'       => ignore(),
                        'stop'        => 1360,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL036C',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => -1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Nonessential' ],
                            'load_id'          => [ 'YAL036C' ],
                            'gene'             => [ 'RBG1' ],
                            'Note'             => [
'Member of the DRG family of GTP-binding proteins; interacts with translating ribosomes and with Tma46p'
                            ],
                            'dbxref'             => [ 'SGD:S000000034' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'RBG1', 'FUN11' ],
                            'Ontology_term' =>
                              [ 'GO:0005737', 'GO:0005525', 'GO:0000004' ]
                        },
                        'start' => 251
                    },
                    'Bio::DB::SeqFeature'
                ),
                '3' => bless(
                    {
                        'source'      => 'BIO',
                        'primary_id'  => 3,
                        'store'       => ignore(),
                        'stop'        => 2700,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'test_tag',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 0,
                        'type'        => 'test_tag',
                        'attributes'  => {
                            'version' => [ '1' ],
                            'intro'   => [ '0_01' ],
                            'load_id' => [ 'test_tag' ],
                            'wtpos'   => [ '2680' ]
                        },
                        'start' => 2680
                    },
                    'Bio::DB::SeqFeature'
                ),
                '7' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 7,
                        'store'       => ignore(),
                        'stop'        => 7476,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL044C',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => -1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'essential_status' => [ 'Nonessential' ],
                            'load_id'          => [ 'YAL044C' ],
                            'gene'             => [ 'GCV3' ],
                            'Note'             => [
'H subunit of the mitochondrial glycine decarboxylase complex',
                                ' required for the catabolism of glycine to 5',
'10-methylene-THF; expression is regulated by levels of levels of 5',
                                '10-methylene-THF in the cytoplasm'
                            ],
                            'dbxref'             => [ 'SGD:S000000042' ],
                            'orf_classification' => [ 'Verified' ],
                            'Alias'              => [ 'GCV3' ],
                            'Ontology_term'      => [
                                'GO:0006730', 'GO:0006546',
                                'GO:0005960', 'GO:0005739',
                                'GO:0004375', 'GO:0000723'
                            ]
                        },
                        'start' => 6964
                    },
                    'Bio::DB::SeqFeature'
                ),
                '2' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 2,
                        'store'       => ignore(),
                        'stop'        => 2193,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAL044W-A',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'orf_classification' => [ 'Uncharacterized' ],
                            'essential_status'   => [ 'Nonessential' ],
                            'load_id'            => [ 'YAL044W-A' ],
                            'Ontology_term' =>
                              [ 'GO:0008372', 'GO:0005554', 'GO:0000004' ],
                            'Note' => [
                                'Similar to pombe uvi31',
                                ' putative DNA repair protein'
                            ],
                            'dbxref' => [ 'SGD:S000007586' ]
                        },
                        'start' => 1861
                    },
                    'Bio::DB::SeqFeature'
                ),
                '5' => bless(
                    {
                        'source'      => 'SGD',
                        'primary_id'  => 5,
                        'store'       => ignore(),
                        'stop'        => 3398,
                        'ref'         => 'chr01',
                        'is_circular' => 0,
                        'name'        => 'YAR028W',
                        'score'       => undef,
                        'phase'       => undef,
                        'strand'      => 1,
                        'type'        => 'gene',
                        'attributes'  => {
                            'orf_classification' => [ 'Uncharacterized' ],
                            'essential_status'   => [ 'Nonessential' ],
                            'load_id'            => [ 'YAR028W' ],
                            'Ontology_term' =>
                              [ 'GO:0008372', 'GO:0005554', 'GO:0000004' ],
                            'Note' => [
                                'Putative integral membrane protein',
                                ' member of DUP240 gene family'
                            ],
                            'dbxref' => [ 'SGD:S000000076' ]
                        },
                        'start' => 2694
                    },
                    'Bio::DB::SeqFeature'
                )
            },
            '_children' => {},
            'setting'   => {
                'serializer'        => 'Storable',
                'compress'          => undef,
                'index_subfeatures' => 1
            },
            'fasta_fh'        => ignore(),
            'seqfeatureclass' => 'Bio::DB::SeqFeature',
            'fasta_db'        => ignore(),
        },
        'Bio::DB::SeqFeature::Store::memory'
    ),
    '_root_verbose' => 0
}, 'Bio::BioStudio::Chromosome';

$cchr->add_feature(
    -feature => $yfeat
);
$cchr->write_chromosome();
cmp_deeply( $cchr, $rchr, "add with sequence change" );

