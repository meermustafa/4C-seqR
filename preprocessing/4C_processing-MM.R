# Meer Mustafa
# September 17, 2018

# processing genome FASTA file to convert to reduced genome

# input:
# 1. Kentutils oligoMatch output of enzyme 1 on full genome (.bed)
# 2. Kentutils oligoMatch output of enzyme 2 on full genome (.bed)
# 3. name of enzyme 1 e.g. csp6i
# 4. name of enzyme 2 e.g. dpnii
# 5. selection of fragments: 'all' to include blind and nonblind, or 'OnlyNonBlind' to keep only nonblind fragments


# output:
# 1. Kentutils oligoMatch output of enzyme 1 on full genome (.bed) with 3 fields added: 
# 1) NonblindStatus on 5' side of enzyme 1, 2) and 3' side of enzyme1 3) length VP (distance between the consecutive enzyme 1 sites)



#########################
######## Outline ########
#########################

# 1. find pairs of enzyme 1 sites that have enzyme 2 sites in between them

# first take 2 two consecutive rows of the enz1 bed file
# check if they have the same chr
# if yes, continue
# set the interval to seek overlaps with to be the 3' end of the 1st enzyme 1 site and the 5' end of the 2nd enzyme 1 enzyme
# if not, next

# lessen the search space for overlaps
# subset the enz2 bed file to be




#setwd("~/Dropbox (Sanjana Lab)/SLab histone/4C-seq/Analysis/MM/")
#list.files()


# enzyme1sites = read.table('test_hg19_csp6i_restriction_sites_oligomatch.bed');head(enzyme1sites)
# # add index
# enzyme1sites$index = seq_len(nrow(enzyme1sites))
# enzyme2sites = read.table('test_hg19_dpnii_restriction_sites_oligomatch.bed');head(enzyme2sites)

### passing in command line arguments ----
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ptm = proc.time() # calc time it takes to run on system

cat('Reading in data...\n')
enzyme1sites = read.table(args[1])#;head(enzyme1sites)
# add index
enzyme1sites$index = seq_len(nrow(enzyme1sites))
enzyme2sites = read.table(args[2])#;head(enzyme2sites)


library(data.table);library(GenomicRanges)

rowIndexList = list()
ViewPointDNAfragmentLength = list()
ViewPointDNAfragment = list()
NonBlindStatus = list()
alloverlaps = list()

cat('Finding enzyme II sites within paired enzyme I sites (viewpoints)...\n')
# iterate from row 1 to the second to last row
for (rowIndex in 1:(nrow(enzyme1sites)-1)) {
  
  # progress diagnostic at every 10th percentile
  if (rowIndex %% round((nrow(enzyme1sites)*0.1), 0) == 0) {
    cat((rowIndex/nrow(enzyme1sites))*100,'% done...\n')
  } 
  
  # check if enz1 sites on same chr
  if( enzyme1sites[rowIndex, 1] == enzyme1sites[(rowIndex+1), 1] ) {
    
    rowIndexList[[rowIndex]] = rowIndex
    
    site1 = enzyme1sites[rowIndex, ]
    site2 = enzyme1sites[(rowIndex+1), ]
    
    Prime5_viewpoint = site1[, 2]
    Prime3_viewpoint = site2[, 3]
    ViewPointDNAfragmentLength[[rowIndex]] = Prime3_viewpoint - Prime5_viewpoint
    
    # save the chr, start, end of viewpoint
    ViewPointDNAfragment[[rowIndex]] = data.table(VPchr = enzyme1sites[rowIndex, 1],
                                                  VPstart = Prime5_viewpoint,
                                                  VPend = Prime3_viewpoint)
    
    
    # create table with keys
    enzyme1Intervals = data.table(V1 = enzyme1sites[rowIndex, 1],
                                  V2 = Prime5_viewpoint,
                                  V3 = Prime3_viewpoint
                                  ,key = c('V1','V2','V3')
                                  )
    
    # narrow down enzyme 2 table to those in the near vicinity of this enzyme 1 interval
    enzyme2PossibleMatches = data.table(V1 = enzyme1Intervals$V1,
                                        V2 = enzyme2sites[ ( enzyme2sites$V1 == enzyme1Intervals$V1 
                                                             & enzyme2sites$V2 > (enzyme1Intervals$V2 - 10000) 
                                                             & enzyme2sites$V2 < (enzyme1Intervals$V2 + 10000) ), 
                                                           2],
                                        V3 = enzyme2sites[ ( enzyme2sites$V1 == enzyme1Intervals$V1 
                                                             & enzyme2sites$V2 > (enzyme1Intervals$V2 - 10000) 
                                                             & enzyme2sites$V2 < (enzyme1Intervals$V2 + 10000) ),
                                                           3]
                                        , key = c('V1','V2','V3')
                                        )
    
    # overlap the enzyme 1 interval (viewpoint) with the possible enzyme 2 matches
    overlaps  = foverlaps(
              enzyme2PossibleMatches, 
      enzyme1Intervals, 
              type = 'within', 
              nomatch = 0L)
    
    alloverlaps[[rowIndex]] = overlaps
    
    if (nrow(overlaps) > 0) {
      NonBlindStatus[[rowIndex]] = 'TRUE'
    } else {
      NonBlindStatus[[rowIndex]] = 'FALSE'
    }
    
    
  }

  
}



cat('Unlisting all lists & assembling together...\n')
ViewPointIndices = unlist(rowIndexList)
ViewPointDNAfragmentLength = unlist(ViewPointDNAfragmentLength)
ViewPointDNAfragment = rbindlist(ViewPointDNAfragment)
NonBlindStatus = unlist(NonBlindStatus)
NonBlindStatus_VPindex = data.frame(ViewPointDNAfragment,
                                    NonBlindStatus,
                                    ViewPointIndices, 
                                    ViewPointDNAfragmentLength)
#NonBlindStatus_VPindex



# merge back to original enzyme 1 file
cat('Merging original enzyme I file with new viewpoint info...\n')
NonBlindStatus_VPindex = merge(enzyme1sites,NonBlindStatus_VPindex,by.x='index',by.y='ViewPointIndices', all = F)
#NonBlindStatus_VPindex

# save merged DF as RDS
cat('Saving DF as RDS...\n')
save(list = 'NonBlindStatus_VPindex',file = paste0(args[3],args[4],NonBlindStatus_VPindex,'DF',collapse = "_"))


# save merged DF as txt 
cat('Writing DF out...\n')
if(args[5] == 'OnlyNonBlind') {
  # name file like ${enzyme1}_${enzyme2}_
  write.table(NonBlindStatus_VPindex[NonBlindStatus_VPindex$NonBlindStatus == 'TRUE', ], 
              paste0('filtered',args[1], collapse = "_"), quote = F, sep = '\t', col.names = F, row.names = F)
}
if(args[5] == 'all') {
  # name file like ${enzyme1}_${enzyme2}_
  write.table(NonBlindStatus_VPindex, 
              paste0('filtered',args[1], collapse = "_"), quote = F, sep = '\t', col.names = F, row.names = F)
}

# how long did this script take?
(proc.time() - ptm)[1]/60/60



