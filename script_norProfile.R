############################ Get BlastP data ###########################
setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/')

# get all human protein names
library('Biostrings')

# get annotation
humanPro <- readBStringSet('hsa.fasta')
humanNames <- names(humanPro)
humanNamesSpli <- strsplit(humanNames, split = ' ', fixed = TRUE)
humanProNames <- sapply(humanNamesSpli, function(x) x[1])
humanAnno <- sapply(humanNamesSpli, function(x){paste(x[-1], collapse = ' ')})
humanProAnno <- cbind(humanProNames, humanAnno)

# process phyloRes data
phyloResFile <- dir('blastRes', full.names = TRUE)
phyloResSpe <- sapply(strsplit(phyloResFile, split = '_', fixed = TRUE), function(x) x[2])
phyloResSpe <- sapply(strsplit(phyloResSpe, split = '.', fixed = TRUE), function(x) x[1])

## make NA matrix
dfhit <- 1
wholePhyloData <- matrix(dfhit, ncol = length(phyloResFile), nrow = nrow(humanProAnno))
rownames(wholePhyloData) <- humanProAnno[, 1]
colnames(wholePhyloData) <- phyloResSpe

for (i in 1:length(phyloResFile)){

  print(paste('It is running ', phyloResSpe[i], '.', sep = ''))
  phyloResData <- read.table(phyloResFile[i])
  phyloResData <- aggregate(phyloResData[, 12], list(phyloResData[, 1]), max)
  gc()

  wholePhyloData[match(phyloResData[, 1], humanProAnno[, 1]), i] <- phyloResData[, 2]
}

save(humanProAnno, wholePhyloData, file = 'wholePhyloDataHit.RData')
#####################################################################

#############################prepare NRR profile######################
library('PhyloProfile') ## version 0.3.13

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholePhyloDataHit.RData')

norProfile <- NPPNor(wholePhyloData, bitCutoff = 70, bitReset = 1)
save(norProfile, file = 'NPP70_profile.RData')
#####################################################################

##########################SVD preprocess############################
library('PhyloProfile') ## version 0.3.13

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholePhyloDataHit.RData')

norProfile100 <- SVDPhy(wholePhyloData, bitCutoff = 60, bitReset = 1, trimming = 1, minConserve = -0.1)
norProfile15 <- SVDPhy(wholePhyloData, bitCutoff = 60, bitReset = 0.15, trimming = 1, minConserve = -0.1)
norProfile15 <- SVDPhy(wholePhyloData, bitCutoff = 60, bitReset = 0.3, trimming = 1, minConserve = -0.1)
norProfile15 <- SVDPhy(wholePhyloData, bitCutoff = 60, bitReset = 0.45, trimming = 1, minConserve = -0.1)

save(norProfile100, norProfile45, norProfile30, norProfile15, file = 'SVD_profile.RData')
####################################################################


