#!/bin/bash

# Meer Mustafa
# Sept 12, 2018

module purge
module load bowtie2/2.2.8


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


    echo " bowtie2 -p ${THREADS} -N 0 \
	--un ${fq}.sam \
	-x ${BOWTIEIDX} \
	-U ${fq} \
	-S ${fq}.sam " | qsub -V -cwd -m abe -M mmustafa@nygenome.org -N ${fq}_job -pe smp ${THREADS}

    echo "Submitted!"
done
	



bowtie2 -p ${THREADS} -N 0 \
	--un ${fq}.sam \
	-x ${BOWTIEIDX} \
	-U {fq} \
	-S ${fq}.sam



# after alignment with bowtie, perform wait command and then convert the final aligned.sam files into bedgraphs
wait
bash sam_bedGraph.sh input_dir_forSAMfiles/





# WAIT SHELL CMD
command1 & #run command1 in background
PID=$! #catch the last PID, here from command1
command2 #run command2 while command1 is running in background
wait $PID #wait for command1, in background, to end
command3 #execute once command1 ended





# bowtie2
# -5 option can be used to trim the first ‘x’ bps that contain the barcode (if present)
# -un write out unaligned reads to a file
# -U is the input FASTQ file name
# -S is the output SAM file name
# -x is prefix database name