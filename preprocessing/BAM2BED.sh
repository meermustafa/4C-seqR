#!/bin/bash

module load bedtools/2.25.0

# save the 1st argument on the command line as an env. variable
#directory=$1
# use it to grep the ls
#INPUT=($(ls *${directory}))

INPUT=($(ls *.sorted.bam))

THREADS=6


for sam in ${INPUT[@]}; do
    echo "About to work on SAM file ${sam}"
    echo "bedtools bamtobed -i ${sam} > ${sam}.bed " | qsub -V -cwd -m abe -M mustafa@nygenome.org -N ${sam}_job -pe smp ${THREADS}
done