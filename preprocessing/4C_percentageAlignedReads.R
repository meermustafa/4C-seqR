# Meer Mustafa
# Sep 26, 2018

# utility script to get diagnostic measurements about alignments to various sites

########## % reads aligning within viewpoints fragments ########## ----- 

list.files(pattern = '*.sorted.bam')
percentageReadsWithinViewpoints = function(readsBEDfile, Prime5UpstreamViewpointBED, Prime3DownstreamViewpointBED) {
  
  # read in files
  readsBEDfile = read.table(readsBEDfile, header = F)
  Prime5UpstreamViewpointBED = read.table(Prime5UpstreamViewpointBED, header = F)
  Prime3DownstreamViewpointBED = read.table(Prime3DownstreamViewpointBED, header = F)
  
  
  
  
}


########## % reads aligning within viewpoint fragment edges ########## ----- 

########## proxmity of closest read to any enzyme 1 site ########## ----- 
# read in BED file of aligned reads,


# function that reads in bed file