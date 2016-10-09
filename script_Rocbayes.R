#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

library('PhyloProfile') ## version 0.3.7
library('ape')
library('pROC')

load('complexAll/allRS_cut40_seed123.RData')
load('wholePhyloData.RData')

profile <- t(wholePhyloDataNet)
tree <- read.nexus('lessSpeciesLR/RAxML_bestTree.nexus')
profile <- profile[match(tree$tip, rownames(profile)), ]

simLR <- BayesTraitsBatch(ftMat = allRS,
                          profileMat = profile,
                          n = 8,
                          BayesTraitsPath = '/home/Yulong/Biotools/BayesTraitsV2/BayesTraitsV2',
                          treeFilePath = 'lessSpeciesLR/RAxML_bestTree.nexus')
simLR <- sapply(simLR, '[[', 3)

LRMat <- data.frame(simLR = simLR, status = allRS[, 3])
LRRoc <- roc(status ~ simLR, LRMat, levels = c('TP', 'TN'))

save(LRMat, LRRoc, file = 'complexAll/LRROC_cut30_seed123.RData')
