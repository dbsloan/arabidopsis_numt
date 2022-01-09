#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_fasta_pairwise_alignment\n\n";

my $file = shift or die ($usage);

my ($header_ref, $seqs_ref) = get_fasta_names_and_seqs($file);

my @headers = @{$header_ref};
my @seqs = @{$seqs_ref};

scalar (@headers) == 2 or die ("\nERROR: fasta file should include 2 and only 2 sequences.\n\n");

print "Position\t$headers[0]\t$headers[1]\n";

my $length = length ($seqs[0]);

$length == length ($seqs[1]) or die ("\nERROR: the two sequences should be of the same length.\n\n");
my $pos = 0;

for (my $i = 0; $i < $length; ++$i){
	
	substr ($seqs[1], $i, 1) eq '-' or ++$pos;
	
	unless (substr ($seqs[0], $i, 1) eq substr ($seqs[1], $i, 1)){
		print $pos, "\t", substr ($seqs[0], $i, 1), "\t", substr ($seqs[1], $i, 1), "\n";
	}
}