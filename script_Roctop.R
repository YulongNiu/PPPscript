#! /usr/bin/Rscript --vanilla

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

library('PhyloProfile') ## version 0.3.8
library('pROC')

load('complexAll/allRS_cut40_seed123.RData')
load('complexAll/simdistROC_cut40_seed123.RData')

GetPosJac <- function(jacMat, corMat, linkVec, corSet){
  # USE: choose the top position in Jaccard mat
  # INPUT: 'jacMat' Jaccard matrix. 'corMat' correlation matrix. 'jacMat' and 'corMat' should have the same rownames and colnames. 'linkVec' linkages (from, to, jaccard, cor), for example c("hsa:55851" "hsa:10569" "0.389313"  "0.5293421"). 'corSet' is the correlation cutoff.
  # OUTPUT: the minimum top number in 'jacMat'
  # EXEAMPLE: topInter(jaccardSim, wholeCor, c('hsa:1004', 'hsa:34235'), 0)

  geneNames <- rownames(jacMat)
  linkFromNum <- match(linkVec[1], geneNames)
  linkToNum <- match(linkVec[2], geneNames)
  jacValue <- as.numeric(linkVec[3])
  corValue <- as.numeric(linkVec[4])

  ## select with corSet
  if (corValue <= corSet) {
    ## the last one
    minTop <- length(geneNames)
  } else {
    jacFrom <- jacMat[, linkFromNum]
    corFrom <- corMat[, linkFromNum]
    jacFrom <- jacFrom[corFrom > corSet]

    jacTo <- jacMat[, linkToNum]
    corTo <- corMat[, linkToNum]
    jacTo <- jacTo[corFrom > corSet]
    
    fromTop <- sum(jacFrom >= jacValue)
    toTop <- sum(jacTo >= jacValue)
    minTop <- min(fromTop, toTop)
  }

   return(minTop)
}


GetPosJacBatch <- function(jacMatDesc, corMatDesc, linkMat, corSet, n = 1) {

  require('bigmemory')
  require('doParallel')
  require('foreach')
  
  jacMat <- attach.big.matrix(jacMatDesc)
  corMat <- attach.big.matrix(corMatDesc)

  ## register multiple core
  registerDoParallel(cores = n)
  linkNum <- nrow(linkMat)

  batchVec <- foreach(i = 1:linkNum, .combine = c) %dopar% {
    print(paste0('It is running ', i, ' in a total of ', linkNum, '.'))
    eachSD <- GetPosJac(jacMat, corMat, linkMat[i, ], corSet)
    return(eachSD)
  }

  ## stop multiple core
  stopImplicitCluster()

  return(batchVec)
  
}
linkMat <- cbind(allRS[, 1:2], jacMat[, 1], corMat[, 1])
TPNum <- sum(allRS[, 3] == 'TP')
distTop <- GetPosJacBatch('jaccardSim.desc', 'wholeCor.desc', linkMat, 0, n = 4)
N <- sum(allRS[, 3] == 'TP')
topMat <- data.frame(distTop = distTop, status = allRS[, 3])
topRoc <- roc(status ~ distTop, topMat, levels = c('TP', 'TN'))

save(topMat, topRoc, file = 'complexAll/topROC_cut40_seed123.RData')








