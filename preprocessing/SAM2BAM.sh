#!/bin/bash

INPUT=($(ls *.sam))

THREADS=6

module purge
module load samtools/1.7

# # convert SAM -> BAM
# samtools view -S -b sample.sam > sample.bam

# # convert BAM -> sorted BAM
# samtools sort sample.bam -o sample.sorted.bam

# # index sorted BAM
# samtools index sample.sorted.bam


# use the fq env. variable to name all other files
for sam in ${INPUT[@]}; do
    echo "About to work on SAM file ${sam}"

    echo "
    # convert SAM -> BAM
	samtools view -S -b ${sam} > ${sam}.bam \

	# convert BAM -> sorted BAM
	samtools sort ${sam}.bam -o ${sam}.sorted.bam \

	# index sorted BAM
	samtools index ${sam}.sorted.bam " | qsub -V -cwd -m abe -M mmustafa@nygenome.org -N ${sam}_job -pe smp ${THREADS}

    echo "Submitted!"
done
