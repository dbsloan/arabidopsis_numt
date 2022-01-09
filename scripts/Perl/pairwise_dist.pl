#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 fasta_file (pairwise alignment)\n\n";

my $file = shift or die ($usage);

my ($header_ref, $seq_ref) = get_fasta_names_and_seqs ($file);

my @headers = @{$header_ref};
my @seqs = @{$seq_ref};

my $alignable = 0;
my $vars = 0;

print "Alignment Position\tNuc1\tNuc2\n";

for (my $i = 0; $i < length ($seqs[0]); ++$i){

	if (substr ($seqs[0], $i, 1) eq '-' or substr ($seqs[1], $i, 1) eq '-'){
		next;
	}
	
	++$alignable;
	
	unless (substr ($seqs[0], $i, 1) eq substr ($seqs[1], $i, 1)){
		++$vars;
		print $i + 1, "\t", substr ($seqs[0], $i, 1), "\t", substr ($seqs[1], $i, 1), "\n";
	} 

}

print "\n\nVariants: $vars\n$alignable Positions: $alignable\nPairwise Identity: ", 1 - ($vars/$alignable), "\n";