# Meer Mustafa
# October 9, 2018

# 4C - differentially interacting DNA regions

# take in smoothed 4C signals and compare contact to sites in various conditions



# create list of 4C bedgraph signal bedgraphs
createList = function(patternInEnvToCombineToList) {
  
  newList = list()
  
  # find objects in env through listing
  # search for items inside list using grep
  namesToRetrieve = names(as.list(environment())) [ grep (patternInEnvToCombineToList, names(as.list(environment())) ) ]
  
  # retrieves the R objects
  # TODO
  #newList = sapply(namesToRetrieve, get, envir = sys.frame(sys.parent(0)), simplify = F) # sys.frame(sys.parent(0))

  newList = sapply(ls(pattern = patternInEnvToCombineToList), get, simplify = F)
  
  return(newList) 
}

a = createList(patternInEnvToCombineToList = 'smoothed_normalized')
a



# function to take in n 4C bedgraph samples, VP boundaries, the contact that is in question's boundaries
# output: DF (that can be boxplotted) that contains the contact boundary's signal pileup distribution, the sample name

# full left peak colon enhancer
contactBoundary = data.frame(chr = 'chr8', start = 128225000, end = 128240000)


find4CContactFrequency = function (ListOf4CBedGraphs, contactBoundary) {
  # ListOf4CBedGraphs is a list of data frames containing 4C signal bedgraphs
  # contactBoundary is a data frame containing chr, start and end of the contact boundary in question
  
  # empty list for saving the contacts of each sample inside of 
  ListOf4CContactFrequencies = list()
  
  # iterate through the 
  for (listIndex in 1:length(ListOf4CBedGraphs)) {
    
    # get name of sample
    sampleName = sapply(strsplit(names(ListOf4CBedGraphs)[listIndex], split = '_', fixed = T), function(x)x[1])
    
    current4CBedGraph = ListOf4CBedGraphs[[listIndex]]
    
    # find all signal in the contact boundaries
    contactsInCurrent4CBedGraph = current4CBedGraph [current4CBedGraph[,'chr'] == contactBoundary[,1] & 
                                                       current4CBedGraph[,'start'] > contactBoundary[,2] & 
                                                       current4CBedGraph[,'end'] < contactBoundary[,3] 
                                                     , ]
    # add column with sample name
    contactsInCurrent4CBedGraph$sampleName = sampleName
    
    # save to list
    ListOf4CContactFrequencies[[listIndex]] = contactsInCurrent4CBedGraph
  }
  
  names(ListOf4CContactFrequencies) = names(ListOf4CBedGraphs)
  
  return(ListOf4CContactFrequencies)
  
}


contactfreqs = find4CContactFrequency(ListOf4CBedGraphs = newList, contactBoundary = contactBoundary)
contactfreqs
rbindlist(contactfreqs)
rbind.data.frame(contactfreqs, make.row.names = T)

# normalize contact frequency
normalize4CContactFrequency = function (ListOf4CContactFrequencies, sampleWithMeanToNormalizeAllSamplesTo, columnWithSignal) {
  
  # ListOf4CContactFrequencies is result of find4CContactFrequency function
  
  # unbind the list into a DF
  DF = as.data.frame(rbindlist(ListOf4CContactFrequencies))
  print(head(DF))
  
  # subset all samples to the sample in question
  DFsubset = DF[DF$sampleName == sampleWithMeanToNormalizeAllSamplesTo, ]
  
  # get mean
  meanOfSample = mean(DFsubset[ , columnWithSignal], na.rm = T)
  
  # divide all values by mean to get signal normalized to one sample
  DF[,columnWithSignal] = DF[,columnWithSignal] / meanOfSample
  
  return(DF)
  
}

normalizedcontactfreqDF = normalize4CContactFrequency(ListOf4CContactFrequencies = contactfreqs, 
                            sampleWithMeanToNormalizeAllSamplesTo = 'A549 Rep1 Myc1', 
                            columnWithSignal = 2
                            )

# plot box plots
with(normalizedcontactfreqDF, boxplot(`libOverlappedWithWindows$NormReadCounts` ~ sampleName, 
                                      col = palette()[1:length(unique(sampleName))]
                                      )
     )


     