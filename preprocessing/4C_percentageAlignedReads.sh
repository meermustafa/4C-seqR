#!/bin/bash

module load bedtools/2.25.0


# loop through the sorted.bam.bed files containing all aligned reads
INPUT=($(ls *.sorted.bam.bed))



########## % reads aligning within viewpoints fragments ##########

for bed in ${INPUT[@]}; do
    echo "About to work on BED file ${bed}"

	totalreads=$(wc -l ${bed}| awk '{print $1}')

	# how many total reads overlap the 5 and 3 prime sides of the VPs
	prime5=$(bedtools intersect -a ${bed} -b Prime5_viewpoint_intervals.txt | wc -l)

	prime3=$(bedtools intersect -a ${bed} -b Prime3_viewpoint_intervals.txt | wc -l)

	# add the two overlapped read counts
	sum=$(echo | awk "{print $prime5 + $prime3}")

	# percentage reads inside the viewpoint edges
	percentageReadsInVPedges=$( echo | awk "{print $sum / $totalreads * 100}" )

	# append this BED file name & % read count
	echo ${bed} >> BEDfiles.txt
	echo $percentageReadsInVPedges >> percentageReadsInVPedges.txt

done

# cbind the %s and the 
awk '{print $1}'  BEDfiles.txt | paste -d' ' percentageReadsInVPedges.txt - > BEDpercentagesWithinVPedges.txt



########## proxmity of closest read to any enzyme 1 site ########## ----- 

# option 1 --
bedtools closest -a ${bed} -b /gpfs/commons/groups/sanjana_lab/MYCScreen/4C-seq/MM/4Cker/hg19_csp6i_restriction_sites_oligomatch.bed  -d > closest.bed


for bed in ${INPUT[@]}; do
    echo "About to work on BED file ${bed}"

    # find the closest enzyme 1 edge (5' or 3') to the read (5' or 3')
	bedtools closest -a ${bed} -b /gpfs/commons/groups/sanjana_lab/MYCScreen/4C-seq/MM/4Cker/hg19_csp6i_restriction_sites_oligomatch.bed  -d > ${bed}.closest.bed
	echo "Created ${bed}.closest.bed..."

	# run R script to load in the file and produce a histogram of distance from read to closest enzyme 1 site
	Rscript getHistograms.R *closest.bed 13

done

# cbind the %s and the 
awk '{print $1}'  BEDfiles.txt | paste -d' ' percentageReadsInVPedges.txt - > BEDpercentagesWithinVPedges.txt









# option 2 --
# read in BED file of enzyme 1 cut sites
# read in BED file of aligned reads
# for these reads:
# if the read is on the + strand, the 3' coordinate on this read should be closest to the 5' coordinate of the enzyme 1 site
# subtract column 3 of reads BED from column 1 of enzyme 1 sites
# sort by increasing values



















BED='VL_4C-a549-rep1-distal1_S5_L001_R1_001.fastq.gz.trim_temp.fq_50bp_trim.btie.fq.aligned.sam.sorted.bam.bed'

totalreads=$(wc -l VL_4C-a549-rep1-distal1_S5_L001_R1_001.fastq.gz.trim_temp.fq_50bp_trim.btie.fq.aligned.sam.sorted.bam.bed | awk '{print $1}')

# how many total reads overlap the 5 and 3 prime sides of the VPs
prime5=$(bedtools intersect -a VL_4C-a549-rep1-distal1_S5_L001_R1_001.fastq.gz.trim_temp.fq_50bp_trim.btie.fq.aligned.sam.sorted.bam.bed -b Prime5_viewpoint_intervals.txt | wc -l)

prime3=$(bedtools intersect -a VL_4C-a549-rep1-distal1_S5_L001_R1_001.fastq.gz.trim_temp.fq_50bp_trim.btie.fq.aligned.sam.sorted.bam.bed -b Prime3_viewpoint_intervals.txt | wc -l)

# add the two overlapped read counts
sum=$(echo | awk "{print $prime5 + $prime3}")

# percentage reads inside the viewpoint edges
percentageReadsInVPedges=$( echo | awk "{print $sum / $totalreads * 100}" )

# append this BED file name & % read count
echo $BED >> BEDfiles.txt
echo $percentageReadsInVPedges >> percentageReadsInVPedges.txt

# cbind the %s and the 
awk '{print $1}'  BEDfiles.txt | paste -d' ' percentageReadsInVPedges.txt - > BEDpercentagesWithinVPedges