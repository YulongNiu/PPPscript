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
setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholePhyloDataHit.RData')

## ref: http://www.nature.com/nature/journal/v493/n7434/extref/nature11779-s1.pdf
## step1: set hit < 50 to 1
norProfile <- apply(wholePhyloData, 1:2, function(x){
  x <- ifelse(x < 50, 1, x)
  return(x)
})

## step2: log2(x[i, j]/max(x[i, ]))
norProfile <- apply(norProfile, 1, function(x) {
  x <- log2(x/max(x))
  return(x)
})
norProfile <- t(norProfile)

## step3: z-score for each column
norProfile <- apply(norProfile, 2, function(x) {
  x <- scale(x)
  return(x)
})
rownames(norProfile) <- rownames(wholePhyloData)

save(norProfile, file = 'NPP_profile.RData')
#####################################################################

###############################NPP_cor##############################
setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

library('PhyloProfile') ## version 0.3.10
library('pROC')

load('complexAll/allRS_cutInf_seed123.RData')
load('NPP_profile.RData')

profile <- t(norProfile)

## correlation
simcor <- SimDistBatch(allRS, profile, SimCor, n = 4)
corMatNPP <- data.frame(simcor = simcor, status = allRS[, 3])
corRocNPP <- roc(status ~ simcor, corMatNPP, levels = c('TP', 'TN'))

save(corMatNPP, corRocNPP, file = 'complexAll/NPPCorROC_cutInf_seed123.RData')
####################################################################


