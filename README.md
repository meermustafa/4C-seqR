# 4C-seqR
Repo for processing &amp; analysing 4C-seq

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


## Analysis
