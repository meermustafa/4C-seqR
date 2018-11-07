# 4C-seqR
Processing &amp; analysing 4C-seq to get comparisons between experiments & call significantly interacting sites

## Processing
FASTA:
Digesting your full genome into a reduced FASTA file. 
Subsetting the enzyme I viewpoints to only the 5' and 3' edges of the viewpoint.

FASTQ:
Trim the reads to only the interaction sequence.
Trim the reads further to have uniform length (e.g. 50 bp).

## Alignment
Create bowtie index of the FASTA. Align the FASTQ reads to the digested genome.


## Checks
Does 50% of reads map to the same chromosome as the viewpoint? (cis interactions)

Optional: What percentage of reads map within the topologically associated domain (TAD) as the viewpoint? (cis interactions)


## Analysis
Normalization

Comparisons across conditions (e.g. across cell types, experiments, treatments)

Calling significantly interacting regions with associated p-value
