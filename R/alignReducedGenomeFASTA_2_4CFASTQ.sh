#!/bin/bash

# Meer Mustafa
# Sept 12, 2018

module purge
module load bowtie/1.1.2

# input files as a list. must iterate over these
INPUT=($(ls *trim))

# specify # of cores used
THREADS=6


# use the fq env. variable to name all other files
for fq in ${INPUT[@]}; do
    echo "about to work on fastq file ${fq}"

    echo "bowtie -v 1 -m 1 hg19_fragment_50_reducedGenome -q ${fq} ${fq}.btie.fq; cat ${fq}.btie.fq | cut -f 3 | sort | uniq -c | sort -rn > ${fq}.trim.bowtie.results " | qsub -V -cwd -m abe -M mmustafa@nygenome.org -N ${fq}_job -pe smp ${THREADS}

    echo "all done!"
done

