# arabidopsis_numt
Analysis of *Arabidopsis thaliana* Chromosome 2 numt. This repository contains scripts, sequence data, and sequence alignments associated with the de novo assembly and analysis of the large numt on chromosome 2 of the *A. thaliana* Col-0 genome.

## Scripts

### R

#### numt_map.R

- This R script takes the input files in the `numt_map` subdirectory and generates a figure summarizing the structure of the numt. It requires the `diagram` package.

### Perl

All Perl scripts require the [sloan.pm](https://github.com/dbsloan/perl_modules) module.

#### pairwise_divergence_msa.pl

- This scripts calculates divergence between pairs of sequences (repeat copies) in a multi-fasta alignment. It was called as follows. Output was then manually updated to exclude MNVs (dinucleotide subs and microinversions) and unalignable regions including SNVs right at indel boundaries.

- `perl  pairwise_divergence_msa.pl  repeat1.aligned.fas`

- `perl  pairwise_divergence_msa.pl  repeat2.aligned.fas`

#### extract_PI_positions.pl

- This script extracts pair of sites for which two repeat copies within the numt share an allele that is different than the mitogenome and the other repeat copy. It was called as follows.

- `perl  extract_PI_positions.pl  repeat1.aligned.fas`

- `perl  extract_PI_positions.pl  repeat2.aligned.fas`

#### pairwise_dist.pl

- This script calculates pairwise nucleotide sequence identity between two sequences in a fasta alignment. It ignores indel gaps. It was run on a pairwise alignment between the numt and a manually constructed concatenation of corresponding regions in the Col-0 mitogenome after removal of structural variants and multinucleotide variants.

- `perl  pairwise_dist.pl  alignment_emboss_stretcher_noStructOrMNVs.fas`

#### summarize_pairwise_variants.pl

- This script provides a list of all variants between two sequences in a pairwise fasta alignment. It was called on multiple different EMBOSS Stretcher alignment files as follows. Where necessary for counting purposes, variants spanning multiple alignments positions (MNVs, indels, and complex structural variants) were then manually collapsed.

- `perl  summarize_pairwise_variants.pl  input_fasta_file`

## Sequence data and alignments

The `sequences_and_alignments` directory contains data and input files used with the above scripts

### Repeat alignments

#### repeat1.aligned.fas and repeat2.aligned.fas

- These are alignments in fasta format of regions present in three copies in the numt along with the corresponding sequence from the Col-0 mitogenome.

### Pairwise EMBOSS Stretcher alignments

Pairwise alignments in fasta format generated with EMBOSS Stretcher

#### Col-CEN-denovo_vs_Col-XJTU-denovo.fas

- An alignment of the numt assembled from Col-CEN HiFi reads and the numt assembled from Col-XJTU HiFi reads

#### Col-CEN-denovo_vs_Col-XJTU-published.fas

- An alignment of the numt assembled from Col-CEN HiFi reads and the numt from the Col-XJTU assembly published by Wang et al. 2021.

#### Col-CEN-published_vs_Col-XJTU-published.fas

- An alignment of the numt from the Col-CEN assembly published by Naish et al. 2021 and the numt from the Col-XJTU assembly published by Wang et al. 2021.

#### Col-CEN-denovo_vs_mitogenome.fas

- An alignment of the numt assembled from Col-CEN HiFi reads and a manually constructed concatenation of corresponding regions in the Col-0 mitogenome.

#### alignment_emboss_stretcher_noStructOrMNVs.fas

- A modified version of Col-CEN-denovo_vs_mitogenome.fas that removes structural variants and multinucleotide variants.

### numt sequences

The numt sequences extracted from our de novo assemblies.

#### numt_Col-CEN.fas and numt_Col-CEN_withFlanking.fas

- Each file contains the 641kb numt sequence obtained from a de novo assembly of the HiFi reads from the Col-CEN project. The "withFlanking" version contains 10kb if flanking sequence on either side of the numt.

#### numt_Col-XJTU.fas and numt_Col-XJTU_withFlanking.fas

- Each file contains the 641kb numt sequence obtained from a de novo assembly of the HiFi reads from the Col-CEN project. The "withFlanking" version contains 10kb if flanking sequence on either side of the numt.

