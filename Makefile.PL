use ExtUtils::MakeMaker;
use strict;
use warnings;

my @scripts = qw (Split_senstive_read_counter.pl);

WriteMakefile(
	'NAME'	=>	'SSRC',
	'AUTHOR'	=>	q{Konstantin Riege},
	'LICENSE'	=>	'perl',
	'PREREQ_PM'	=>	{
		'Getopt::Long'	=>	'2.0',
		'Bio::Root::Version'	=>	'1.00690001',
		'Bio::DB::Sam'	=>	'1.39',
		'Bio::DB::SeqFeature::Store'	=>	'0',
		'Try::Tiny'	=>	'0',
		'Pod::Usage'	=>	'0',
		'IO::Select'	=>	'0',
		'IO::Pipe'	=>	'0',
		'POSIX'	=>	'0'
	},
	'MIN_PERL_VERSION'	=>	'5.010',
	'EXE_FILES'	=>	[ map { "$_" } @scripts ],
	'CONFIGURE_REQUIRES'	=>	{
		'ExtUtils::MakeMaker'	=>	0
	},
	'BUILD_REQUIRES'	=>	{
		'Test::More' => 0
	},
	dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', }
);
