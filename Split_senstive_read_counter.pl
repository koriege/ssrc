#! /usr/bin/perl

use strict;
use warnings;
use IO::Select;
use IO::Pipe;
use Bio::DB::Sam;
use Bio::DB::SeqFeature::Store;
use v5.10;
use Try::Tiny;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(:sys_wait_h);

(Getopt::Long::Parser->new)->getoptions (
	'gff|gff=s' => \my $gfffile,
	'gtf|gtf=s' => \my $gtffile,
	'gix|gtfididx=i' => \my $gix,
	'exidsep|exonidseperator=s' => \my $exidsep,
	'bam|bam=s' => \my $bamfile,
	'gene|gene=s' => \my $genetag,	
	'exon|exon=s' => \my $exontag,
	'gc|genecount' => \my $genecount,
	't|threads=i' => \my $threads,
	'o|out=s' => \my $outfile,
	'ol|overlap=i' => \my $overlap,
	'h|help' => \my $help
) or pod2usage(-exitval => 2, -verbose => 1);
pod2usage(-exitval => 2, -verbose => 2) if $help;
pod2usage(-exitval => 2, -verbose => 1) unless $gfffile || $gtffile || $bamfile || $outfile || $genetag || $exontag;
pod2usage(-exitval => 2, -verbose => 1) if $gtffile && ! $gix;
pod2usage(-exitval => 2, -verbose => 1) if $gfffile && ! -e $gfffile;
pod2usage(-exitval => 2, -verbose => 1) if $gtffile && ! -e $gtffile;
pod2usage(-exitval => 2, -verbose => 1) unless -e $bamfile;
$threads = 1 unless $threads;
$overlap = 20 unless $overlap;

sub add_gff_gene(){
	my ($db,$line) = @_;

	chomp $line;
	return if ! $line || $line=~/^#/ || $line=~/^\s*$/;	
	my @line = split /\s+/ , $line;	
	pod2usage(-exitval => 1, -verbose => 3) unless $line[-1]=~/id=/i;

	return unless $line[2] eq $genetag;
	
	my %attributes = map {
		my @k_v = split /=/ , $_; 
		if($k_v[0]=~/^id$/i){
			ID => $k_v[1];
		} else {
			$k_v[0] => $k_v[1];
		}
	} split /;/ , $line[-1];

	$db->new_feature(
		-start => $line[3],
		-stop => $line[4],
		-seq_id => $line[0],
		-strand => $line[6],
		-primary_tag => $line[2],
		-source => $line[1],
		-score => $line[5],
		-index => 1,
		-attributes => \%attributes
	);	
}

sub add_gff_exon(){
	my ($db,$line) = @_;

	chomp $line;
	return if ! $line || $line=~/^#/ || $line=~/^\s*$/;	
	my @line = split /\s+/ , $line;	
	pod2usage(-exitval => 1, -verbose => 3) unless $line[-1]=~/id=/i;

	return unless $line[2] eq $exontag;

	my %attributes = map {
		my @k_v = split /=/ , $_; 
		if($k_v[0]=~/^id$/i){
			ID => $exidsep ? (split /$exidsep/ , $k_v[1])[0] : $k_v[1];
		} else {
			$k_v[0] => $k_v[1];
		}
	} split /;/ , $line[-1];

	my ($parent) = $db->get_features_by_attribute(ID => $attributes{ID});
	my $child = $db->new_feature(
		-start => $line[3],
		-stop => $line[4],
		-seq_id => $line[0],
		-strand => $line[6],
		-primary_tag => $line[2],
		-source => $line[1],
		-score => $line[5],
		-index => 0
	);
	$db->add_SeqFeature($parent,($child));
}

sub add_gtf_gene(){
	my ($db,$line) = @_;

	chomp $line;

	return if ! $line || $line=~/^#/ || $line=~/^\s*$/;	
	my @line  = split /\s+/ , $line;

	return unless $line[2] eq $genetag;

	my $id = substr($line[$gix],1,-2);	
	my %attributes = $gix == 9 ? (ID => $id) : ($line[8] => substr($line[9],1,-2), ID => $id);
	my @attributes = split /;/ , $line;
	for (1..$#attributes){		
		next if $_ == $gix - 11;
		my @k_v = split /\s+/ , $attributes[$_];		
		$attributes{$k_v[1]}=substr($k_v[2],1,-1);
	}
	
	$db->new_feature(
		-start => $line[3],
		-stop => $line[4],
		-seq_id => $line[0],
		-strand => $line[6],
		-primary_tag => $line[2],
		-source => $line[1],
		-score => $line[5],
		-index => 1,
		-attributes => \%attributes
	);	
}

sub add_gtf_exon(){
	my ($db,$line) = @_;

	chomp $line;

	return if ! $line || $line=~/^#/ || $line=~/^\s*$/;	
	my @line  = split /\s+/ , $line;

	return unless $line[2] eq $exontag;

	my $id = $exidsep ? (split /$exidsep/ , substr($line[$gix],1,-2))[0] : substr($line[$gix],1,-2);
	
	my ($parent) = $db->get_features_by_attribute(ID => $id);
	my $child = $db->new_feature(
		-start => $line[3],
		-stop => $line[4],
		-seq_id => $line[0],
		-strand => $line[6],
		-primary_tag => $line[2],
		-source => $line[1],
		-score => $line[5],
		-index => 0		
	);
	$db->add_SeqFeature($parent,($child));	
}

my $gff = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
$gff->init_database([1]);
say "reading reference annotation";
if ($gfffile){	
	open GFF , '<'.$gfffile or die $!;
	&add_gff_gene($gff,$_) while(<GFF>);
	unless ($genecount) {
		seek GFF, 0, 0;
		&add_gff_exon($gff,$_) while(<GFF>);
	}
	close GFF;
} else {
	$gix--;
	$gff = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
	$gff->init_database([1]);
	open GTF , '<'.$gtffile or die $!;
	&add_gtf_gene($gff,$_) while(<GTF>);
	unless ($genecount) {
		seek GTF, 0, 0;
		&add_gtf_exon($gff,$_) while(<GTF>);
	}
	close GTF;
}

my $libsize = 0;
my @features = $gff->features;
my @data;
my @pids;
my $chunks = ($#features+1)/($threads-1);
push @data, [ splice @features, 0, $chunks ] while @features;
my $select = IO::Select->new();


my $bam = Bio::DB::Sam->new(-bam => $bamfile, -autindex => 1, -verbose => -1);
my $next;		
my %chr;
for (`samtools idxstats $bamfile`){ #todo counts all fragments
	next if $_=~/^\*/;
	next if $_=~/fail/;			
	my ($header , $size, $mapped, $unmapped) = split /\s+/ , $_;
	$chr{$header}=1;
	$libsize += $mapped;	
}		
while (my $features = shift @data){
	my $pipe = IO::Pipe->new;
	my $pid = fork();			
	unless ($pid) {
		$pipe->writer();
		$pipe->autoflush(1);
		$bam->clone;

		for my $f (@$features){					
			next unless exists $chr{$f->seq_id};

			my ($count,$tpm,$strand,$length) = (0,0,1,0);
			
			my @reads;
			for my $r ($bam->get_features_by_location(-type => 'read_pair', -seq_id => $f->seq_id, -start => $f->start, -end => $f->stop)){
				next if $r->stop - $f->start < 20 || $f->stop - $r->start < 20;
				for ($r->segments){
					$strand = -1 if $_->get_tag_values('FLAGS')=~/REVERSED/;
				}
				push @reads , $r if $strand == $f->strand;				
			}
			my ($id) = $f->get_tag_values('ID');
			if ($genecount){
				$count += $#reads + 1;
				$length = $f->stop - $f->start;
			} else {
				my $readdb = Bio::DB::SeqFeature::Store->new( -adaptor => 'memory', -verbose => -1 );
				$readdb->add_features(\@reads);

				for my $s ( sort {$a->start <=> $b->start} $f->segments){
					$length += $s->stop - $s->start;
					my @o = $readdb->features(-range_type => 'overlaps', -start => $s->start, -stop => $s->stop); 					
					for my $r (@o){						
						my ($start, $stop, $strand) = $s->intersection($r);
						if ($start && $stop-$start >= $overlap){
							$count++;
							$readdb->delete(($r));
						}	
					}
				}				
			}
			say $pipe "$id $length $count";
		}
		exit;
	} else {				
		push @pids, $pid;
		$pipe->reader();
		$select->add($pipe);
	}
}

waitpid($_,0) for @pids;

my %counts;
my %lengths;
while( my @responses = $select->can_read(0) ){
	for my $pipe (@responses){			
		while(<$pipe>){	
			my ($id, $length, $count)  = split /\s+/ , $_;
			$counts{$id} += $count;
			$lengths{$id} = $length;
		}			
		$select->remove($pipe->fileno);
	}
}

open OUT , '>'.$outfile or die $!;
for my $f ($gff->features){
	my ($id) = $f->get_tag_values('ID');
	my ($count, $tpm ) = (0,0);
	my @out = split /\s+/ , $f->gff_string;
	$#out=8;
	$out[8]='';
	my @tags = $f->get_all_tags;
	for (sort @tags){
		my ($v) = $f->get_tag_values($_);
		$out[8].=$_.'='.$v.';';
	}
	try {		
		if (exists $counts{$id}){
			$count = $counts{$id};
			$tpm = ($count / ($lengths{$id}/10^3)) / ($libsize/10^6);
		}	
	} catch {

	};
	$out[8].="reads=$count;TPM=$tpm";
	say OUT join "\t" , @out;
}
close OUT;

__END__

=head1 NAME

SSRC - counts mapped reads for given annotation file

counts splitted reads only once
- counts strand specific
- can count on whole genes (default: exons only)
- calculates TPM values
- input can be gff or gtf
- output is gene based gff
- is multithreaded

=head1 DEPENDENCIES

Samtools, Bio::Perl, Bio::DB::Sam
 

=head1 SYNOPSIS

bamsplitpaircount.pl [-h]
  
example: bamsplitpaircount.pl -bam /my/bam -gff /my/gff -gene CDS -exon exon -o /my/out

=head1 DESCRIPTION

B<-h>, B<--help>	
	this (help) message

B<-o>, B<--out>=F<FILE>

	(requierd)
	output gff file with read counts and tpm values

B<-bam>, B<--bam>=F<FILE>

	(requierd)
	single bam file
	
B<-gff>, B<--gff>=F<FILE>

	(requierd)
	input is gff file - needs case insensitive id=uniqGeneID001 tag

B<-gtf>, B<--gtf>=F<FILE>

	(requierd)
	input is gtf file - requiers option -gix

B<-gix>, B<--gtfididx>=I<INT>
	
	(requierd if input is gtf)
	column of "uniqGeneID001";

B<-gene>, B<--gene>=F<STRING>

	(requierd)
	for feature type gene

B<-exon>, B<--exon>=F<STRING>

	(requierd)
	for feature type exon

B<-ol>, B<--overlap>=I<INT>
	
	(optional, default: 20)	
	minimum nucleotide overlap of read with gene/exon to count

B<-t>, B<--threads>=I<INT>

	(optional)
	number of threads to use

B<-gc>, B<--genecount>

	(optional)
	count reads mapped on whole gene instead of exons only


=head1 AUTHOR

Konstantin Riege, E<lt>konstantin.riege@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2016 by Konstantin Riege

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

The full description of all in-/outputs and parameters is maintained as PDF manual. 
You can access it on L<www.rna.uni-jena.de/software.php>.

=cut

