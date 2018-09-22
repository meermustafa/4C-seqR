# Meer Mustafa
# September 19, 2018

library(Basic4Cseq)
library(BSgenome.Hsapiens.UCSC.hg19)

enzymeImotif = 'GTAC'
enzymeIImotif = 'GATC'
readLengthOfFragment = 50

createVirtualFragmentLibrary(chosenGenome = BSgenome.Hsapiens.UCSC.hg19 ,
                             firstCutter = enzymeImotif,
                             secondCutter = enzymeIImotif,
                             readLength = readLength,
                             libraryName = paste0("hg19_fragment_",readLengthOfFragment,".csv")
                             )

