# Meer Mustafa
# October 8, 2018

# correlations between signal vectors



# 1st run 4C-analysis to make the bedgraphs of smoothed normalized 4C signal


ls(pattern = 'smoothed')

# concatenate all subset libraries together
allMergedScreensOnUID = Reduce(function(x,y) merge(x = x, y = y, by = "UniqueID"), list(
  `A549 Rep1 Distal1_normalized_bdg4C`[,'NormReadCounts'], 
  `A549 Rep1 Distal2_normalized_bdg4C`[,'NormReadCounts'], 
  `A549 Rep1 Myc1_normalized_bdg4C`[,'NormReadCounts'],
  `A549 Rep2 Myc1_normalized_bdg4C`[,'NormReadCounts'],
  `A549 Rep1 Myc2_normalized_bdg4C`[,'NormReadCounts'],
  `A549 Rep1 Pvt1_normalized_bdg4C`[,'NormReadCounts'],
  `A549 Rep1 Pvt2_normalized_bdg4C`[,'NormReadCounts']
  #, 
))

allMergedScreensOnUID = data.frame(`A549 Rep1 Distal1_smoothed_normalized_bdg4C`[,2], 
                                     `A549 Rep1 Distal2_smoothed_normalized_bdg4C`[,2], 
                                     `A549 Rep1 Myc1_smoothed_normalized_bdg4C`[,2],
                                     `A549 Rep2 Myc1_smoothed_normalized_bdg4C`[,2],
                                     `A549 Rep1 Myc2_smoothed_normalized_bdg4C`[,2],
                                     `A549 Rep1 Pvt1_smoothed_normalized_bdg4C`[,2],
                                     `A549 Rep1 Pvt2_smoothed_normalized_bdg4C`[,2]
                                     )

head(allMergedScreensOnUID)
nrow(allMergedScreensOnUID)



# subset to only the mean 20 guide CS
# allMergedScreensOnUID_OnlyCS = allMergedScreensOnUID[, grep('CS',colnames(allMergedScreensOnUID)) ]
# head(allMergedScreensOnUID_OnlyCS)

# round to 2 digits
allMergedScreensOnUID_OnlyCS = round(allMergedScreensOnUID, 5)
head(allMergedScreensOnUID_OnlyCS)

# compute correlation
allMergedScreensOnUID_OnlyCS_Corr = cor(allMergedScreensOnUID_OnlyCS, use = 'pairwise.complete.obs')
allMergedScreensOnUID_OnlyCS_Corr

# change the rownames and colnames of the correlation matrix
colnames(allMergedScreensOnUID_OnlyCS_Corr) = c('Rep1 Distal 1',
                                                'Rep1 Distal 2',
                                                'Rep1 MYC 1',
                                                'Rep2 MYC 1',
                                                'Rep1 MYC 2',
                                                'Rep1 PVT1 1',
                                                'Rep1 PVT1 2'
)
rownames(allMergedScreensOnUID_OnlyCS_Corr) = c('Rep1 Distal 1',
                                                'Rep1 Distal 2',
                                                'Rep1 MYC 1',
                                                'Rep2 MYC 1',
                                                'Rep1 MYC 2',
                                                'Rep1 PVT1 1',
                                                'Rep1 PVT1 2'
)
allMergedScreensOnUID_OnlyCS_Corr

library(pheatmap)
breaksList = seq(0, 1, 0.1)

library(RColorBrewer)
# set colour
color = colorRampPalette(rev (brewer.pal(n = 9, name = "Reds"))) (length(breaksList))

pdf(useDingbats = F, width = 4.5, height = 4, file = '4C normlized correlation heat map w noms.pdf')

pheatmap(allMergedScreensOnUID_OnlyCS_Corr,
         color = rev(color), 
         #breaks = seq(0,1,0.05),
         border_color = 'black',
         fontsize = 12,
         fontsize_row = 16,
         fontsize_col = 16, 
         display_numbers = T,
         number_format = "%.2f", 
         fontsize_number = 12,
         number_color = 'grey8',
         show_rownames = T,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F)

dev.off()


pdf(useDingbats = T, onefile = T, width = 4.5, height = 4, file = 'subset ess. gene controls myc screens 200 bp from TSS corr.pdf')
pheatmap(allMergedScreensOnUID_OnlyCS_Corr,
         color = rev(color), 
         #breaks = seq(0,1,0.05),
         border_color = 'black',
         fontsize = 12,
         fontsize_row = 16,
         fontsize_col = 16, 
         display_numbers = T,
         number_format = "%.2f", 
         fontsize_number = 12,
         number_color = 'black',
         show_rownames = T,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F)

dev.off()
