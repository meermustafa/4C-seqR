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




# use the fq env. variable to name all other files
for fq in ${INPUT[@]}; do
    echo "about to work on fastq file ${fq}"


    echo " ; \
    bowtie -v 1 -m 1 hg19_fragment_50_reducedGenome -q ${fq} ${fq}.btie.fq; cat ${fq}.btie.fq | cut -f 3 | sort | uniq -c | sort -rn > ${fq}.trim.bowtie.results " | qsub -V -cwd -m abe -M mmustafa@nygenome.org -N ${fq}_job -pe smp ${THREADS}

    echo "all done!"
done
	


bowtie2-build mm10_hindiii_flanking_sequences_25_unique.fa mm10_hindiii_flanking_sequences_25_unique


bowtie2 -p 12 -N 0 \
	--un sample_unaligned.sam \
	-x mm10_hindiii_flanking_sequences_25_unique \
	-U sample.fastq \
	-S sample_aligned.sam



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