#! /usr/bin/Rscript --vanilla

## setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

## #############################construct TP and TN###################
## load('complexAll/interMat_cutInf.RData')
## load('complexAll/randomMat.RData')

## set.seed(123)
## PRS <- interMat
## NRS <- randomMat[sample(1:nrow(randomMat), 10 * nrow(PRS)), ]
## NRS <- apply(NRS, 1:2, as.character)
## allRS <- rbind(PRS, NRS)
## allRS <- cbind(allRS, rep(c('TP', 'TN'), c(nrow(PRS), nrow(NRS))))
## rownames(allRS) <- NULL
## colnames(allRS) <- c('From', 'To', 'Status')
## save(allRS, file = 'complexAll/allRS_cutInf_seed123.RData')

## set.seed(456)
## PRS <- interMat
## NRS <- randomMat[sample(1:nrow(randomMat), 10 * nrow(PRS)), ]
## NRS <- apply(NRS, 1:2, as.character)
## allRS <- rbind(PRS, NRS)
## allRS <- cbind(allRS, rep(c('TP', 'TN'), c(nrow(PRS), nrow(NRS))))
## rownames(allRS) <- NULL
## colnames(allRS) <- c('From', 'To', 'Status')
## save(allRS, file = 'complexAll/allRS_cutInf_seed456.RData')
## ###################################################################

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

###########################simdist################################
library('PhyloProfile') ## version 0.3.11
library('pROC')

load('complexAll/allRS_Bioinfo.RData')
load('NPP70_profile.RData')

profile <- t(norProfile)

## correlation
simcor <- SimDistBatch(allRS, profile, SimCor, n = 8)
corMat <- data.frame(simcor = simcor, status = allRS[, 3])
corRoc <- roc(status ~ simcor, corMat, levels = c('TP', 'TN'))

## ## jaccard
## simjac <- SimDistBatch(allRS, profile, SimJaccard, n = 8)
## jacMat <- data.frame(simjac = simjac, status = allRS[, 3])
## jacRoc <- roc(status ~ simjac, jacMat, levels = c('TP', 'TN'))

## MI
simMI <- SimDistBatch(allRS, profile, SimMIConti, n = 8)
MIMat <- data.frame(simMI = simMI, status = allRS[, 3])
MIRoc <- roc(status ~ simMI, MIMat, levels = c('TP', 'TN'))

## ## hamming
## distham <- SimDistBatch(allRS, profile, DistHamming, n = 8)
## hamMat <- data.frame(distham = distham, status = allRS[, 3])
## hamRoc <- roc(status ~ distham, hamMat, levels = c('TP', 'TN'))

## euclidean
disteu <- SimDistBatch(allRS, profile, DistEuclidean, n = 8)
euMat <- data.frame(distham = disteu, status = allRS[, 3])
euRoc <- roc(status ~ disteu, euMat, levels = c('TP', 'TN'))

## save(corMat, corRoc, jacMat, jacRoc, MIMat, MIRoc, hamMat, hamRoc, euMat, euRoc, file = 'complexAll/simdistROC_Bioinfo.RData')
save(corMat, corRoc, MIMat, MIRoc, euMat, euRoc, file = 'complexAll/simdistROCNPP70_Bioinfo.RData')
##################################################################
