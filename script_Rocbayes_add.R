#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')
##################################LR matrix###########################
library('PhyloProfile') ## version 0.3.7
library('pROC')
library('ape')

load('complexAll/allRS_cut40_seed123.RData')
allRS40 <- allRS
load('complexAll/allRS_cutInf_seed123.RData')
allRSInf <- allRS
load('wholePhyloData.RData')

allRS40Vec <- apply(allRS40, 1, paste, collapse = '|')
allRSInfVec <- apply(allRSInf, 1, paste, collapse = '|')

allRS <- allRSInf[!(allRSInfVec %in% allRS40Vec), ]

profile <- t(wholePhyloDataNet)
tree <- read.nexus('lessSpeciesLR/RAxML_bestTree.nexus')
profile <- profile[match(tree$tip, rownames(profile)), ]

simLR <- BayesTraitsBatch(ftMat = allRS,
                          profileMat = profile,
                          n = 8,
                          BayesTraitsPath = '/home/Yulong/Biotools/BayesTraitsV2/BayesTraitsV2',
                          treeFilePath = 'lessSpeciesLR/RAxML_bestTree.nexus')
simLR <- sapply(simLR, '[[', 3)

save(simLR, allRS, file = 'complexAll/addRS402Inf.RData')
#################################################################

## ###############################add LR#############################
## library(pROC)

## load('complexAll/LRROC_cut30_seed123.RData')
## load('complexAll/allRS_cut30_seed123.RData')
## allRS30 <- allRS
## load('complexAll/allRS_cut40_seed123.RData')
## allRS40 <- allRS
## load('complexAll/addRS30240.RData')


## allRS30Vec <- apply(allRS30, 1, paste, collapse = '|')
## allRS40Vec <- apply(allRS40, 1, paste, collapse = '|')
## addRSVec <- apply(allRS, 1, paste, collapse = '|')
## addLRVec <- c(as.numeric(LRMat[, 1]), simLR)
## names(addLRVec) <- c(allRS30Vec, addRSVec)

## LRMat <- data.frame(simLR = addLRVec[match(allRS40Vec, names(addLRVec))],
##                     status = allRS40[, 3])
## LRRoc <- roc(status ~ simLR, LRMat, levels = c('TP', 'TN'))

## save(LRMat, LRRoc, file = 'complexAll/LRROC_cut40_seed123.RData')
## ##################################################################

