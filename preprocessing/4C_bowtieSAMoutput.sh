#!/bin/bash

# Meer Mustafa
# Sept 12, 2018

module purge
module load bowtie/1.1.2


# input files as a list. must iterate over these
INPUT=($(ls *trim))

# specify # of cores used
THREADS=6


BOWTIEFASTA='hg19_csp6i_flanking_sequences_50_unique.fa'  #ls *unique.fa
BOWTIEIDX='hg19_csp6i_flanking_sequences_50_unique' # ls *unique.fa stripped of fa name

# build index
bowtie2-build ${BOWTIEFASTA} ${BOWTIEIDX}


# use the fq env. variable to name all other files
for fq in ${INPUT[@]}; do
    echo "About to work on fastq file ${fq}"

    echo "bowtie -v 1 -m 1 ${BOWTIEIDX} -q ${fq} -S ${fq}.btie.aligned.sam " | qsub -V -cwd -m abe -M mmustafa@nygenome.org -N ${fq}_job -pe smp ${THREADS}

    echo "Submitted!"
done
