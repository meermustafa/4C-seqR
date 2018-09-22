#!/bin/bash
#
# name it
#$ -N reducdGen
#
# execute the job from current directory
#$ -cwd
#
# notification to email address
#$ -m abe
#$ -M mm8805@nyu.edu
#
# preserve the environment
#$ -V
#
# Modify these depending on your job
#
# specify number of CPUs on a single host
#$ -pe smp 8
#
# specify the maximum amount of memory
#$ -l h_vmem=32G
#
# specify the time needed
# -l h_rt=4:00:00
#
# put your script below this line


# genome digestion as done by 4C-seq 2-step digestion & ligation

# this script takes in 3 inputs:
# 1) enzyme 1 motif FASTA file
# 2) enzyme 2 motif FASTA file
# 3) genome FASTA file (e.g. hg19.fa)

# this script has 2 outputs:
# 1) n bp read length FASTA file of possible DNA interaction sequences
# 2) n bp read length BED file with coordinates of possible DNA interaction sequences


module load kentutils/302.1
module load bedtools/2.26.0
module load R/3.5.1
#cd firstPassMM


# modify fragment length, enzyme for naming, genome for naming

# fragment length
fl=50
# primary enzyme for file name
enzyme=csp6i
# genome for file name
genome=hg19
# enzyme 2 for file name
enzyme2=dpnii

# digest the genome and produce a bed file of nonblind VP intervals
Rscript 4C_createFragmentLibrary.R ${enzyme}.fa ${enzyme2}.fa 50 ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints.bed 

# if numbers are round (e.g. 7000000), they will be concatenated into 7e6 and bedtools cannot recognize those numbers for parsing
Rscript disableScientificNotationR.R ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints.bed ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints.bed

# Fragment library adjusted to reflect 0-based BED coordinates
# remove first row of headers, and keep only column 1,2,3, subtract 1 from column 2 to put the 1-based R indices onto the 0-based BED indices
awk -v shift=1 '{print $1"\t"$2-shift"\t"$3}' ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints.bed | awk '{if (NR!=1) {print}}' > ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints_adjusted.bed 

# for new files containing oligoMatch coord & VP information
# get coordinates of downstream fragments
awk -v fl=$fl '{print $1"\t"$2"\t"$2+fl}' ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints_adjusted.bed > Prime5_viewpoint_intervals.txt

# get coordinates of upstream fragments
awk -v fl=$fl '{print $1"\t"$3-fl"\t"$3}' ${genome}_${enzyme}_restriction_sites_nonblind_viewpoints_adjusted.bed > Prime3_viewpoint_intervals.txt

# combine up and downstream fragments
cat Prime5_viewpoint_intervals.txt Prime3_viewpoint_intervals.txt > ${genome}_${enzyme}_flanking_sites_${fl}_2.bed

# remove any fragments with negative coordinates, unwanted unknown chromosomes, sort
awk '{if($2 >= 0 && $3 >=0) print $0}' ${genome}_${enzyme}_flanking_sites_${fl}_2.bed | grep -v -E 'random|JH|GL' - | sort -k1,1 -k2,2n | uniq  > ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed

# get the sequence of unique flanking coordinates
fastaFromBed -fi ${genome}.fa -bed ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed -fo ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa

# TODO
# is this 'get only sequences from FASTA file' only for bowtie2 and not for bowtie
# get only sequences from FASTA file
grep -v '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | sort | uniq -i -u | grep -xF -f - -B 1 ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | grep -v '^--' > ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa



# remove unwanted intermediate files
# rm Prime5_viewpoint_intervals.txt
# rm Prime3_viewpoint_intervals.txt


# make a BED file of the unique FASTA sequences
grep '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa > ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
# remove the '>'
sed -i 's/>//g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
# convert to bed file of unique possible interactions - this accompanies the FASTA file of n bp seqeunce
sed -i 's/:\|-/\t/g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed

exit 0;

