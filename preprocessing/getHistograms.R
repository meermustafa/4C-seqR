# Meer Mustafa
# Sep 26, 2018

# function that reads in bed file, file pattern to match in directory,t akes column index of which file to make histogram of


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pattern = paste0('*',args[1],collapse = "")

# black white pie chart
getBWpie = function(pieDF, title, sizeOfPie) {
  
  ### takes in DF with V1 = vector of labels, V2 = counts
  pie(table( c(rep(pieDF[,1], times = pieDF[,2]) ) ),
      labels = pieDF[,1],
      main = title, # title
      cex.main = 1.5, # change title size
      cex = 1.5,
      radius = sizeOfPie, # change size of pie
      col = c('black','white')  # change colour of pie
  ) 
  
}


getHistograms = function (filePattern, columnIndexToMakeHistogram) {
  
  listOfFilesAndValues = list()
  
  
  pdf(useDingbats = F, file = paste0(
    #CellLineAndRepAndViewpointName
    '4C-seq'
    , ' aligned read distance to nearest enzyme1 site.pdf')
      , width = 10, onefile = T)
  
  # plot multiple plots
  par(mfrow = c( length(list.files(pattern = filePattern))/2, 
                 (length(list.files(pattern = filePattern)))/2 + 1)
  )
  
  # loop through files that have the pattern
  for (fileIndex in 1:length(list.files(pattern = filePattern))) {
    
    # get file name
    filename = list.files(pattern = filePattern)[fileIndex]
    cat('Getting values for',filename,'...\n')
    
    # names of samples for plots
    CellLineandRep = paste0(sapply(strsplit(filename, split='-',fixed=T),FUN = function(x) x[c(2,3)] ), collapse = " ")
    temp = paste0(sapply(strsplit(filename, split='-',fixed=T),FUN = function(x) x[4] ))
    CellLineAndRepAndViewpointName = paste0( c(CellLineandRep, sapply(strsplit(temp, split='_',fixed=T),FUN = function(x) x[1] )),collapse = ' ')
    
    # read in files
    file = read.table(list.files(pattern = filePattern)[fileIndex], header = F)
    
    currentDF = data.frame(Values = file[, columnIndexToMakeHistogram],
               Name = CellLineAndRepAndViewpointName)
    
    # plots -----
    # pdf(useDingbats = F, file = paste0(CellLineAndRepAndViewpointName, ' aligned read distance to nearest enzyme1 site.pdf')
    #     , width = 10)
    # par(mar = c(5,10,0.1,2)
    #     ,mgp=c(0.2,0.2,0)
    #     )
    
    hist(currentDF$Values,
         breaks = c(50),
         las = 0,
         main = paste0(CellLineAndRepAndViewpointName),
         xlab = "Distance from aligned reads to closest enzyme 1 site"
         #cex.lab = 1.6, cex.axis = 1.5, cex.main = 1.8
         )
    #dev.off()
    
    # pie chart
    percentage1bp = round ( nrow(currentDF[currentDF$Values <= 1,]) / nrow(currentDF) * 100 , 0)
    cat("1 bp %:",percentage1bp,'\n')
    percentageGreater1bp = round ( nrow(currentDF[currentDF$Values > 1,]) / nrow(currentDF) * 100 , 0)
    cat("> 1 bp %:",percentageGreater1bp,'\n')
    
    pieDF = data.frame(Labels = c('1 bp', '> 1 bp'),
                       Percentages = c(percentage1bp, percentageGreater1bp))
    # pdf(useDingbats = F, file = paste0(CellLineAndRepAndViewpointName, ' pie aligned read distance to nearest enzyme1 site.pdf')
    #     , width = 10)
    getBWpie(pieDF = pieDF, title = paste0(CellLineAndRepAndViewpointName), sizeOfPie = 0.5)
    #dev.off()
    
    
    # save DF in list
    listOfFilesAndValues[[fileIndex]] = currentDF
    
  }
  
  dev.off()
  
  # save to env
  assign(x = 'listOfFilesAndValues', value = listOfFilesAndValues, envir = globalenv())
  
}

# call function
getHistograms(filePattern = , 
              columnIndexToMakeHistogram = as.numeric(args[2])
              )


# unbind the list into a DT
library(data.table)
FilesAndValues = rbindlist(listOfFilesAndValues)








# ggplot
library(ggplot2)
ggplot(data =  FilesAndValues, 
       aes(x = Values
           ,colour = Name
           #, fill = allThresholdsAllDistancesToNearestSigGuide$ThresholdName
       )) +
  
  # geom_histogram(binwidth=0.01
  #                #, position="dodge"
  #                , alpha = 0.4
  #                ) +
  
  geom_density(na.rm = T)+
  
  # geom_freqpoly(binwidth = 10) + 
  
  scale_color_brewer(palette="Set1") + # Accent
  
  # scale_color_manual(values = c(
  #   "#0000FF",
  #   '#FF0000',
  #   '#800080',
  #   "#008000", 
  #   '#000000'
  # )) +
  
  
  
  labs (
    #title = paste0(screen), 
        x = "Distance from aligned reads to \nclosest enzyme 1 site"
        #,y = "Frequency"
        , color = ''
        # , fill = 'Library Group'
  ) +
  
  geom_vline(xintercept = 0, linetype = 'dashed') +
  
  #guides(fill=guide_legend("my awesome title")) +
  
  theme_classic() +
  
  theme(title = element_text(size=22), 
        legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size=18),
        # legend.position="bottom",
        axis.text=element_text(size=18), 
        axis.title=element_text(size=20, face="bold"))
