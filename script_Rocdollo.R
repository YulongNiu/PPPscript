#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

library('PhyloProfile') ## version 0.3.12
library('pROC')
library('ape')

load('complexAll/allRS_Bioinfo.RData')
load('wholePhyloData.RData')

profile <- t(wholePhyloDataNet)
tree <- read.nexus('lessSpeciesLR/RAxML_bestTree.nexus')
profile <- profile[match(tree$tip, rownames(profile)), ]
pathList <- nodepath(tree)

distDollo <- DolloDistBatch(ftMat = allRS,
                            profileMat = profile,
                            edgeMat = tree$edge,
                            tipPath = pathList,
                            n = 4)

DolloMat <- data.frame(distDollo = distDollo, status = allRS[, 3])
DolloRoc <- roc(status ~ distDollo, DolloMat, levels = c('TP', 'TN'))

save(DolloMat, DolloRoc, file = 'complexAll/DolloROC_Bioinfo.RData')

