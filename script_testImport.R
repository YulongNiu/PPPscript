##########################test the top 400 results for one single gene###############
library('bigmemory')
library('foreach')
library('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

testGene <- 'hsa:4929'
testGene <- 'hsa:5309'

interList <- foreach(i = 1:ncol(jaccardSim)) %dopar% {
  corVec <- wholeCor[, i]
  jacSimVec <- jaccardSim[, i]
  corSetNum <- which(corVec >= 0)
  jacSimVec <- jacSimVec[corSetNum]

  ## delete self nodes
  jacSimVec <- jacSimVec[names(jacSimVec) != rownames(jaccardSim)[i]]
  jacSimVec <- sort(jacSimVec, decreasing = TRUE)
  
  if (rownames(jaccardSim)[i] == testGene) {
    ## choose self inter genes
    interName <- names(jacSimVec[1:400])
  } else {
    interIdx <- match(testGene, names(jacSimVec))
    if (!is.na(interIdx)) {
      ## cross cor selectioin
      if (interIdx <= 400) {
        interName <- rownames(jaccardSim)[i]
      } else {
        interName <- NULL
      }
    } else {
      interName <- NULL
    }
  }
  return(interName)
}

interVec <- unique(unlist(interList))
#################################################################



################################select Jaccard top 500###############
library('bigmemory')
library('foreach')
library('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

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
save(top500List, file = 'top500List.RData')
######################################################################
