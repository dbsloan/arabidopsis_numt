#!/usr/bin/perl

use strict;
use warnings;
use sloan;
use List::Util qw(min max);

my $usage  = "\nUSAGE: perl $0 fasta_alignment\n\n";

my $file = shift or die ($usage);

my %fasta = fasta2hash($file);

my @headers = sort keys (%fasta);

my %cur_pos; 

foreach (@headers){
	my @sh = split (/\_/, $_);
	my @sr = split (/\-/, $sh[-1]);
	$cur_pos{$_} = $sr[0] - 1; 
}



print "SNP1\tSNP2\n";

for (my $i = 0; $i < length ($fasta{$headers[0]}); ++$i){

	my %allele_counts;
	foreach my $name (@headers){
		my $base = substr ($fasta{$name}, $i, 1);		
		$base eq '-' or ++$cur_pos{$name};
		++$allele_counts{$base};
	}
	
	my @alleles = sort keys (%allele_counts);
	
	if (scalar (@alleles) == 2 and $allele_counts{$alleles[0]} == 2){
		my $pair_found = 0;
		if (substr ($fasta{$headers[1]}, $i, 1) eq substr ($fasta{$headers[2]}, $i, 1)){
			print min ($cur_pos{$headers[1]}, $cur_pos{$headers[2]}), "\t", max ($cur_pos{$headers[1]}, $cur_pos{$headers[2]}), "\n";
			++$pair_found;
		}
		if (substr ($fasta{$headers[1]}, $i, 1) eq substr ($fasta{$headers[3]}, $i, 1)){
			print min ($cur_pos{$headers[1]}, $cur_pos{$headers[3]}), "\t", max ($cur_pos{$headers[1]}, $cur_pos{$headers[3]}), "\n";
			++$pair_found;
		}
		if (substr ($fasta{$headers[2]}, $i, 1) eq substr ($fasta{$headers[3]}, $i, 1)){
			print min ($cur_pos{$headers[2]}, $cur_pos{$headers[3]}), "\t", max ($cur_pos{$headers[2]}, $cur_pos{$headers[3]}), "\n";
			++$pair_found;
		}
		
		$pair_found == 0 and die ("\nERROR: No PI pairing found at $i\n");
		$pair_found > 1 and die ("\nERROR: More than 1 PI pairing found at $i\n");
		
	}
}