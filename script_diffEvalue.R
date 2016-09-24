###################different blast threshold for E-value####################
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholePhyloDataBlast.RData')

thresVec <- c(1e-3, 1e-4, 5e-4)

for (i in thresVec) {
  wholeProfile <- apply(wholePhyloData, 1:2, function(x) {
    if (x < i) {
      return(1)
    } else {
      return(0)
    }
  })

  hsaNum <- ncol(wholeProfile) + 1
  wholeProfile <- cbind(wholeProfile, rep(1, nrow(wholeProfile)))
  colnames(wholeProfile)[hsaNum] <- 'hsa'

  fileName <- paste0('diffEvalue/wholeProfile', i, '.RData')
  save(wholeProfile, file = fileName)
}
############################################################################


#########################generate cor and jaccard matrix##################
library('bigmemory')
library('vegan')
library('compiler')
enableJIT(3)

## file name
thresVec <- c('1e-03', '1e-04', '5e-04')
phyloFileVec <- paste0('diffEvalue/wholeProfile', thresVec, '.RData')

for (i in 1:3) {

  ## load profile data
  load(phyloFileVec[i])
  
  ## jaccard similariy
  jaccardSim <- 1 - vegdist(wholeProfile[1:10, ], method = 'jaccard', diag = TRUE, upper = TRUE)
  jaccardSim <- as.matrix(jaccardSim)
  jaccardSim <- diag(1, ncol(jaccardSim), nrow(jaccardSim)) + jaccardSim
  
  simFilePrefix <- paste0('jaccardSim', thresVec[i])
  as.big.matrix(jaccardSim, backingfile = paste0(simFilePrefix, '.bin'), descriptorfile=paste0(simFilePrefix, '.desc'), backingpath = 'diffEvalue')

  ## correlation coefficient
  wholeCor <- cor(t(wholeProfile))

  corFilePrefix <- paste0('wholeCor', thresVec[i])
  as.big.matrix(wholeCor, backingfile = paste0(corFilePrefix, '.bin'), descriptorfile=paste0(corFilePrefix, '.desc'), backingpath = 'diffEvalue')
}

enableJIT(0)
########################################################################

#######################select Jaccard top 500########################
library('bigmemory')
library('foreach')
library('doMC')
registerDoMC(8)

## file name
thresVec <- c('1e-03', '1e-04', '5e-04')
filePath <- 'diffEvalue'

for (thresi in 1:3) {
  jaccardDescPath <- paste0(filePath, '/jaccardSim', thresVec[thresi], '.desc')
  wholeCorDescPath <- paste0(filePath, '/wholeCor', thresVec[thresi], '.desc')

  jaccardSim <- attach.big.matrix(jaccardDescPath)
  wholeCor <- attach.big.matrix(wholeCorDescPath)

  interList <- foreach(i = 1:ncol(jaccardSim)) %dopar% {
    corVec <- wholeCor[, i]
    jacSimVec <- jaccardSim[, i]
    corSetNum <- which(corVec >= 0)
    jacSimVec <- jacSimVec[corSetNum]

    ## delete self nodes
    jacSimVec <- jacSimVec[names(jacSimVec) != rownames(jaccardSim)[i]]
    jacSimVec <- sort(jacSimVec, decreasing = TRUE)

    eachTop500 <- jacSimVec[1:500]
  }

  interNames <- sapply(interList, names)
  colnames(interNames) <- rownames(jaccardSim)
  interJac <- sapply(interList, unname)
  colnames(interJac) <- rownames(jaccardSim)

  ## get cor
  interCor <- matrix(ncol = ncol(interNames), nrow = nrow(interNames))
  for (i in 1:ncol(interNames)) {
    corIdx <- match(interNames[, i], rownames(wholeCor))
    interCor[, i] <- wholeCor[corIdx, i]
  }
  colnames(interCor) <- rownames(jaccardSim)

  ## round jaccard similarity and cor
  interJac <- round(interJac, 2)
  interCor <- round(interCor, 2)

  top500List <- list(interNames = interNames,
                     interJac = interJac,
                     interCor = interCor)

  save(top500List, file = paste0(filePath, '/top500List', thresVec[thresi], '.RData'))
}
#####################################################################
