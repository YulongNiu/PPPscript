#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

library('PhyloProfile') ## version 0.3.8
library('pROC')
library('ape')

load('complexAll/allRS_cutInf_seed123.RData')
load('wholePhyloData.RData')

profile <- t(wholePhyloDataNet)
tree <- read.nexus('lessSpeciesLR/RAxML_bestTree.nexus')
profile <- profile[match(tree$tip, rownames(profile)), ]
pathList <- nodepath(tree)

distDollo <- DolloDistBatch(ftMat = allRS,
                            profileMat = profile,
                            edgeMat = tree$edge,
                            tipPath = pathList,
                            n = 8)

DolloMat <- data.frame(distDollo = distDollo, status = allRS[, 3])
DolloRoc <- roc(status ~ distDollo, DolloMat, levels = c('TP', 'TN'))

save(DolloMat, DolloRoc, file = 'complexAll/DolloROC_cutInf_seed123.RData')
