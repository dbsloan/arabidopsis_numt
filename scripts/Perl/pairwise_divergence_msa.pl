#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 fasta_input\n\n";

my $file = shift or die ($usage);

my %fasta = fasta2hash($file);

my @headers = sort keys (%fasta);

print "Sequence1\tSequence2\tBase1\tBase2\tAlignmentPosition\n";

my $summary_string = "\n\n\nSequence1\tSequence2\tMismatches\tAlignedPositions\n";

for (my $i = 0; $i < scalar(@headers) - 1; ++$i){
	for (my $j = $i+1; $j < scalar(@headers); ++$j){
		
		my $align_count = 0;
		my $mismatches = 0;
		
		my $seq1 = $fasta{$headers[$i]};
		my $seq2 = $fasta{$headers[$j]};
		
		my $length = length ($seq1);
		$length == length ($seq2) or die ("\nERROR: sequences in the alignment must be the same length: see $headers[$i] and $headers[$j]\n\n");
		
		for (my $k = 0; $k < $length; ++$k){
			
			my $base1 = uc (substr ($seq1, $k, 1));
			my $base2 = uc (substr ($seq2, $k, 1));
			
			$base1 eq '-' and next;
			$base2 eq '-' and next;
			
			++$align_count;
			
			unless ($base1 eq $base2){
				print "$headers[$i]\t$headers[$j]\t$base1\t$base2\t", $k+1, "\n";
				++$mismatches;
			}
			
		}
		
		$summary_string .= "$headers[$i]\t$headers[$j]\t$mismatches\t$align_count\n";
		
		
	}
}

print $summary_string;