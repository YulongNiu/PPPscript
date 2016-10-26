#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')
##################################LR matrix###########################
library('PhyloProfile') ## version 0.3.12
library('pROC')
library('ape')

load('complexAll/allRS_cutInf_seed123.RData')
allRS123 <- allRS
load('complexAll/allRS_Bioinfo.RData')
allRSBio <- allRS
load('wholePhyloData.RData')

allRS123Vec <- apply(allRS123, 1, paste, collapse = '|')
allRSBioVec <- apply(allRSBio, 1, paste, collapse = '|')

allRS <- allRSBio[!(allRSBioVec %in% allRS123Vec), ]

profile <- t(wholePhyloDataNet)
tree <- read.nexus('lessSpeciesLR/RAxML_bestTree.nexus')
profile <- profile[match(tree$tip, rownames(profile)), ]

simLR <- BayesTraitsBatch(ftMat = allRS,
                          profileMat = profile,
                          n = 8,
                          BayesTraitsPath = '/home/Yulong/Biotools/BayesTraitsV2/BayesTraitsV2',
                          treeFilePath = 'lessSpeciesLR/RAxML_bestTree.nexus')
simLR <- sapply(simLR, '[[', 3)

save(simLR, allRS, file = 'complexAll/addRS1232Bio.RData')
#################################################################

###############################add LR#############################
library(pROC)

load('complexAll/LRROC_cutInf_seed123.RData')
load('complexAll/allRS_cutInf_seed123.RData')
allRS123 <- allRS
load('complexAll/allRS_Bioinfo.RData')
allRSBio <- allRS
load('complexAll/addRS1232Bio.RData')

allRS123Vec <- apply(allRS123, 1, paste, collapse = '|')
allRSBioVec <- apply(allRSBio, 1, paste, collapse = '|')
addRSVec <- apply(allRS, 1, paste, collapse = '|')
addLRVec <- c(as.numeric(LRMat[, 1]), simLR)
names(addLRVec) <- c(allRS123Vec, addRSVec)

LRMat <- data.frame(simLR = addLRVec[match(allRSBioVec, names(addLRVec))],
                    status = allRSBio[, 3])
LRRoc <- roc(status ~ simLR, LRMat, levels = c('TP', 'TN'))

save(LRMat, LRRoc, file = 'complexAll/LRROC_Bioinfo.RData')
##################################################################

