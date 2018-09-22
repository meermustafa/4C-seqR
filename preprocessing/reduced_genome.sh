#!/bin/bash
#
# name it
#$ -N reducdGen
#
# Join stdout and stderr into one file
#
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
# -l h_rt=2:00:00
#
# put your script below this line



module load kentutils/302.1
module load bedtools/2.26.0
module load R/3.3.1
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


# get coordinates for all RE sites in the genome
# DON'T FORGET THAT BED FILES ARE 0-BASED!
oligoMatch ${enzyme}.fa ${genome}.fa ${genome}_${enzyme}_restriction_sites_oligomatch.bed
oligoMatch ${enzyme2}.fa ${genome}.fa ${genome}_${enzyme2}_restriction_sites_oligomatch.bed

# TODO Meer
# if you want only non-blind fragments (fragments where two consecutive enzyme 1 sites have an enzyme 2 site in between), we need to find oligo matches of enzyme two in the genome
# then we take each consecutive pair of enzyme 1 sites (with if statement making sure the chr's are the same because bed file is continuous), and ask perform an overlap in R or ask shell if there is a row in the enzyme 2 oligo match file that is positioned in between these two pairs
# if there is a enzyme 2 site positioned between an enzyme 2, save the distance from the closest enzyme 1 to enzyme 2
Rscript 4C_processing-MM.R ${genome}_${enzyme}_restriction_sites_oligomatch.bed ${genome}_${enzyme2}_restriction_sites_oligomatch.bed ${enzyme} ${enzyme2} OnlyNonBlind


# for new files containing oligoMatch coord & VP information
# get coordinates of upstream fragments
awk -v fl=$fl '{print $7"\t"$8-fl"\t"$8}' filtered_${genome}_${enzyme}_restriction_sites_oligomatch.bed > up.txt

# get coordinates of downstream fragments
awk -v fl=$fl '{print $7"\t"$9"\t"$9+fl}' filtered_${genome}_${enzyme}_restriction_sites_oligomatch.bed > down.txt

# for only oligoMatch files
# # get coordinates of upstream fragments
# awk -v fl=$fl '{print $1"\t"$2-fl"\t"$2}' ${genome}_${enzyme}_restriction_sites_oligomatch.bed > up.txt

# # get coordinates of downstream fragments
# awk -v fl=$fl '{print $1"\t"$3"\t"$3+fl}' ${genome}_${enzyme}_restriction_sites_oligomatch.bed > down.txt


# combine up and downstream fragments
cat up.txt down.txt > ${genome}_${enzyme}_flanking_sites_${fl}_2.bed


# remove any fragments with negative coordinates, unwanted unknown chromosomes, sort
awk '{if($2 >= 0 && $3 >=0) print $0}' ${genome}_${enzyme}_flanking_sites_${fl}_2.bed | grep -v -E 'random|JH|GL' - | sort -k1,1 -k2,2n | uniq  > ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed

# TODO
# inspect this .bed in IGV


# use this file which contains coordinates for both sides of the 


# get the sequence of unique flanking coordinates
fastaFromBed -fi ${genome}.fa -bed ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed -fo ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa

# TODO
# is this 'get only sequences from FASTA file' only for bowtie2 and not for bowtie
# get only sequences from FASTA file
grep -v '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | sort | uniq -i -u | grep -xF -f - -B 1 ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | grep -v '^--' > ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa



# remove unwanted intermediate files
# rm up.txt
# rm down.txt


# make a BED file of unique sequences
grep '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa > ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
sed -i 's/>//g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
sed -i 's/:\|-/\t/g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed

exit 0;
