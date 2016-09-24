##########################List --> matrix interaction########################
# interList --> interMat is unique
# interMat --> interList is unique
# some interList is not undirected which need to be transfered at first, like
# mat2listInter(list2matInter(tmp1))

list2matInter <- function(interList){
  # USE: transfer list format to undirected matrix
  # INPUT: 'interList' is the list format interaction should have names. The vector in each list should be gene entrez, like 'hsa:1000'. No self-self interaction.
  # OUTPUT: undirected matrix.

  # list 2 matrix
  fromGenes <- rep(names(interList), sapply(interList, length))
  toGenes <- unlist(interList)
  interMat <- cbind(fromGenes, toGenes)

  # delete repeat
  fromNum <- sapply(strsplit(fromGenes, split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(toGenes, split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  interMat <- rbind(interMat[isSmall, ], interMat[!isSmall, 2:1])
  interMat <- interMat[!duplicated(interMat, MARGIN = 1), ]
  
  return(interMat)
}

list2matInterWithValue <- function(interList){
  # USE: transfer list format to undirected matrix with value
  # INPUT: 'interList' is the list format interaction should have names. The vector in each list should be gene entrez, like 'hsa:1000'. No self-self interaction.
  # OUTPUT: undirected matrix.

  # list 2 matrix
  fromGenes <- rep(names(interList), sapply(interList, ncol))
  toMat <- lapply(interList, function(x) {
    return(t(x))
  })
  toMat <- do.call(rbind, toMat)
  toGenes <- toMat[, 1]
  interMat <- cbind(fromGenes, toMat)

  ## # delete repeat
  fromNum <- sapply(strsplit(fromGenes, split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(toGenes, split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  interMat <- rbind(interMat[isSmall, ], interMat[!isSmall, c(2, 1, 3, 4)])
  interMat <- interMat[!duplicated(interMat, MARGIN = 1), ]

  return(interMat)
}


list2vecInter <- function(interList, removeDup = FALSE) {
  # USE: transfer list format to vector
  # INPUT: 'interList' is the list format interaction should have names. The vector in each list should be gene entrez, like 'hsa:1000'. No self-self interactions. 'removeDup' is whether to remove the duplicated elements, set 'FALSE' to increase speed.
  # OUTPUT: interaction vector.

  # list 2 matrix
  fromGenes <- rep(names(interList), sapply(interList, length))
  toGenes <- unlist(interList)

  fromNum <- sapply(strsplit(fromGenes, split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(toGenes, split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum

  interVec1 <- paste(fromGenes[isSmall], toGenes[isSmall], sep = '|')
  interVec2 <- paste(toGenes[!isSmall], fromGenes[!isSmall], sep = '|')

  interVec <- c(interVec1, interVec2)

  if (removeDup) {
    interVec <- unique(interVec)
  } else {}

  return(interVec)
}
  
tmp1 <- list()
tmp1[[1]] <- paste('hsa', 1:3, sep = ':')
tmp1[[2]] <- paste('hsa', 4:6, sep = ':')
names(tmp1) <- c('hsa:5', 'hsa:1')

mat2listInter <- function(interMat){
  # USE: transfer undirected matrix to list format interaction.
  # INPUT: 'interMat' is the matrix format interaction should have names.
  # OUTPUT: list format interaction.

  allGenes <- unique(c(interMat[, 1], interMat[, 2]))
  interList <- sapply(allGenes, function(x) {
    toInter <- interMat[interMat[, 1] %in% x, 2]
    toInter <- as.character(toInter)
    fromInter <- interMat[interMat[, 2] %in% x, 1]
    fromInter <- as.character(fromInter)

    allInter <- c(toInter, fromInter)
  })

  return(interList)
}

list2allmatInter <- function(interList, fileName = 'interAllmat') {
  # USE: transfer list format to undirected all matrix, in which row number is equal to column number
  # INPUT: 'interList' is the list format interaction should have names. The vector in each list should be gene entrez, like 'hsa:1000'. No self-self interactions. 'fileName' is the backup file prefix.
  # OUTPUT: undirected all matrix (file-backed big matrix object).

  require('bigmemory')
  
  allNodes <- c(names(interList), unlist(interList))
  allNodes <- unique(allNodes)

  interAllmat <- filebacked.big.matrix(nrow = length(allNodes),
                                       ncol = length(allNodes),
                                       dimnames = list(allNodes, allNodes),
                                       init = 0,
                                       type = 'integer',
                                       backingfile = paste(fileName, 'bin', sep = '.'),
                                       descriptorfile = paste(fileName, 'desc', sep = '.'))

  for (i in 1:length(interList)) {
    fromName <- names(interList)[i]
    fromIndex <- match(fromName, allNodes)
    toIndex <- match(interList[[i]], allNodes)

    interAllmat[fromIndex, toIndex] <- as.integer(interAllmat[fromIndex, toIndex] | rep(1, length(toIndex)))
    interAllmat[toIndex, fromIndex] <- as.integer(interAllmat[toIndex, fromIndex] | rep(1, length(toIndex)))
  }

  return(interAllmat)
}

mat2allMatInter <- function(interMat, fileName = 'interAllmat'){
  # USE: transfer undirected matrix to allMat interaction format.
  # INPUT: 'interMat' is the matrix format interaction should have names. 'fileName' is the backup file prefix.
  # OUTPUT: allMat interaction format.

  require('bigmemory')
  
  allNodes <- unique(c(interMat[, 1], interMat[, 2]))

  interAllmat <- filebacked.big.matrix(nrow = length(allNodes),
                                       ncol = length(allNodes),
                                       dimnames = list(allNodes, allNodes),
                                       init = 0,
                                       type = 'integer',
                                       backingfile = paste(fileName, 'bin', sep = '.'),
                                       descriptorfile = paste(fileName, 'desc', sep = '.'))

  for (i in 1:nrow(interMat)) {
    fromIndex <- match(interMat[i, 1], allNodes)
    toIndex <- match(interMat[i, 2], allNodes)
    interAllmat[fromIndex, toIndex] <- as.integer(interAllmat[fromIndex, toIndex] | c(1, 1))
    interAllmat[toIndex, fromIndex] <- as.integer(interAllmat[toIndex, fromIndex] | c(1, 1))
  }

  return(interAllmat)
}
#############################################################################


#######################test accuracy of percentage method################'
# run on the server
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')
load('complex/complexInterList.RData')
## load('jaccardSim.RData')
## load('wholeCor.RData')
## load('complexInterList.RData')

topInter <- function(jacMat, corMat, linkVec, corSet, withValue = FALSE) {
  # USE: choose the top interactions based on Jaccard with cor cutoff.
  # INPUT: 'jacMat' choosed Jaccard matrix. 'corMat' choosed correlation matrix. 'linkVec' is the topN --> topM, for example, c(1001, 2000). 'corSet' is the correlation cutoff. 'withValue' whether to return Jaccard and correlation values, and the default is FALSE (not return).
  # OUTPUT: interaction list. If no interaction, the fromNode will be removed from the final result. ATTENTION: the output interaction list may be not a undirected. 
  # EXEAMPLE: topInter(jaccardSim, wholeCor, c(1, 100), 0)
    
  jacTopInterList <- foreach (i = 1:nrow(jacMat)) %dopar% {
    ## print(paste('It is running ', i, '.', sep = ''))
    ## set corSet
    corVec <- corMat[i, ]
    jacSimVec <- jacMat[i, ]
    corSetNum <- which(corVec >= corSet)
    jacSimVec <- jacSimVec[corSetNum]

    ## delete self nodes
    jacSimVec <- jacSimVec[names(jacSimVec) != rownames(jacMat)[i]]
    jacSimVecLen <- length(jacSimVec)

    if (jacSimVecLen < linkVec[1]) {
      jacSimInter <- NA
    }
    else if (jacSimVecLen >= linkVec[1] & jacSimVecLen < linkVec[2]) {
      ## sort jacSimVec
      jacSimVec <- sort(jacSimVec, decreasing = TRUE)
      jacSimInter <- jacSimVec[linkVec[1] : jacSimVecLen]
    }
    else if (jacSimVecLen >= linkVec[2]) {
      ## sort jacSimVec
      jacSimVec <- sort(jacSimVec, decreasing = TRUE)
      jacSimInter <- jacSimVec[linkVec[1]:linkVec[2]]
    }

    # if jacSimInter == NA, will return NULL
    if (withValue) {
      if (is.na(jacSimInter)) {
        topRes <- NULL
      } else {
        topResCor <- corVec[match(names(jacSimInter), colnames(jacMat))]
        topRes <- rbind(jacSimInter, topResCor)
        topRes <- rbind(names(jacSimInter), topRes)
      }
    } else {
      return(names(jacSimInter))
    }
  }
  
  names(jacTopInterList) <- rownames(jacMat)

  # remove NULL
  isNull <- sapply(jacTopInterList, is.null)
  jacTopInterList <- jacTopInterList[which(!isNull)]

  
  return(jacTopInterList)
}


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
    jacTop <- jacMat[, linkToNum]
    fromTop <- sum(jacFrom >= jacValue)
    toTop <- sum(jacTop >= jacValue)
    minTop <- min(fromTop, toTop)
  }

  ## return(minTop)
  return(minTop)
}


## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## complexNum <- match(names(complexInterList), rownames(jaccardSim))
## jacSelectMat <- jaccardSim[complexNum, ]
## corSelectMat <- wholeCor[complexNum, ]
## top400 <- topInter(jaccardSim, wholeCor, c(1, 400), 0)
## save(top400, file = 'top400.RData')

## top400 <- topInter(jaccardSim, wholeCor, 400, 0)
## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test ATP synthase
# run on server
atpSub <- read.table('F1F0Gene.txt')
atpSub <- match(as.character(atpSub[, 2]), rownames(jaccardSim))
topATP400 <- topInter(jaccardSim, wholeCor, 400, 0, atpSub)
topATP400Mat <- list2matInter(topATP400)


# run on local
require('bigmemory')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')
load('wholePhyloData.RData')

entrez <- rownames(jaccardSim)
NRSindex <- apply(topATP400Mat, 1:2, function(x) {
  return(match(x, entrez))
})

jacsim <- apply(NRSindex, 1, function(x){
  return(jaccardSim[x[1], x[2]])
})

corValue <- apply(NRSindex, 1, function(x){
  return(wholeCor[x[1], x[2]])
})

fromVec <- annoFirst[match(topATP400Mat[, 1], names(annoFirst))]
toVec <- annoFirst[match(topATP400Mat[, 2], names(annoFirst))]

topATP400MatAnno <- data.frame(from = fromVec, to = toVec, jaccard = jacsim, cor = corValue)
rownames(topATP400MatAnno) <- NULL
write.csv(topATP400MatAnno, 'topATP400MatAnno.csv')

atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
topATP400MatAnno[topATP400MatAnno[, 2] %in% atpSub[, 1]]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('KEGGBioCycAPI')
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
## load('complex/our_complex.RData')
## load('complexAll/our_complex_567.RData')
## load('complex/Bioinfor_complex.RData')
load('complexAll/our_complex.RData')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')


RearrangeMat <- function(interMat){
  # USE: to rearrange a matrix that the entrez of first colum is smaller the second one.
  # INPUT: 'interMat' should be an undirected matrix with entrez like 'hsa:1000'.
  # OUTPUT: a rearranged matrix

  fromNum <- sapply(strsplit(interMat[, 1], split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(interMat[, 2], split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  arrangeMat <- rbind(interMat[isSmall, ], interMat[!isSmall, 2:1])

  return(arrangeMat)
  
}

# transfer TP and TN to inter list
NRSMat <- apply(NRSJacCor[, 1:2], 1:2, as.character)
NRSMat <- RearrangeMat(NRSMat)
PRSMat <- apply(PRSJacCor[, 1:2], 1:2, as.character)
PRSMat <- RearrangeMat(PRSMat)
PRSVec <- apply(PRSMat, 1, paste, collapse = '|')
NRSVec <- apply(NRSMat, 1, paste, collapse = '|')
rm(NRSJacCor, PRSJacCor)
gc()

allGenes <- unique(c(PRSMat[, 1], PRSMat[, 2], NRSMat[, 1], NRSMat[, 2]))
entrez <- rownames(jaccardSim)
complexNum <- match(allGenes, entrez)
jacSelectMat <- jaccardSim[complexNum, ]
corSelectMat <- wholeCor[complexNum, ]
rm(jaccardSim, wholeCor, NRSMat, PRSMat, entrez, allGenes, complexNum)
gc()

## tow part
## first part
topNumMat <- matrix(c(rep(0, 11), c(50, seq(100, 1000, 100))), nrow = 2, byrow = TRUE)
## second part
topNumMat <- cbind(CutSeqEqu(4000, 500) + 1000,
                   CutSeqEqu(10000, 1000) + 5000)

ROCMat <- apply(topNumMat, 2, function(x) {
  
  topInterList <- topInter(jacSelectMat, corSelectMat, x, 0)
  ## to undirected mat
  topInterVec <- list2vecInter(topInterList)
  rm(topInterList)
  gc()

  TP <- sum(PRSVec %in% topInterVec)
  ## FN <- length(PRSVec) - TP
  FP <- sum(NRSVec %in% topInterVec)
  ## TN <- length(NRSVec) - FP
  rm(topInterVec)
  gc()
  print(x)
  
  return(c(TP, FP))
})

ROCMat <- t(ROCMat)
rownames(ROCMat) <- topNumMat[2, ]
colnames(ROCMat) <- c('TP', 'FP')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test ROC top~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
## load('complex/our_complex.RData')
## load('complexAll/our_complex_567.RData')
## load('complex/Bioinfor_complex.RData')
load('complexAll/our_complex.RData')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')
NRSJacCor <- apply(NRSJacCor, 1:2, as.character)
PRSJacCor <- apply(PRSJacCor, 1:2, as.character)

PRSVec <- foreach(i = 1:nrow(PRSJacCor), .combine = c) %dopar% {
  print(paste('It is running ', i, sep = ''))
  GetPosJac(jaccardSim, wholeCor, PRSJacCor[i, ], 0)
}

NRSVec <- foreach(i = 1:nrow(NRSJacCor), .combine = c) %dopar% {
  print(paste('It is running ', i, sep = ''))
  GetPosJac(jaccardSim, wholeCor, NRSJacCor[i, ], 0)
}

PRStop <- cbind(PRSJacCor[, 1:2], PRSVec)
NRStop <- cbind(NRSJacCor[, 1:2], NRSVec)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hist plot~~~~~~~~~~~~~~~~~~~~~~~
topVec <- c(PRStop[, 3], NRStop[, 3])
topMat <- data.frame(Top = as.numeric(topVec), Status = rep(c('TP', 'TN'), c(nrow(PRStop), nrow(NRStop))))

p <- ggplot(data = topMat, aes(x = Top))
p +
  geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = factor(Status))) +
  stat_density(geom = 'line', position = 'identity', aes(colour = factor(Status)))
save(PRStop, NRStop, file = 'complexAll/our_complex_top.RData')
## save(PRStop, NRStop, file = 'complexAll/our_complex_top_567.RData')
## save(PRStop, NRStop, file = 'complex/Bioinfor_complex_top.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test top ROC~~~~~~~~~~~~~~~~~~~~
load('complexAll/our_complex_top.RData')
## load('complexAll/our_complex_top_567.RData')
## load('complex/Bioinfor_complex_top.RData')

# verse top
topVec <- 20127 - as.numeric(c(PRStop[, 3], NRStop[, 3]))
topMat <- data.frame(Top = as.numeric(topVec), Status = rep(c('pos', 'neg'), c(nrow(PRStop), nrow(NRStop))))

cutVec <- seq(min(topVec), max(topVec), 10)
topROCMat <- SingROCMat(cutVec, topMat)
save(topROCMat, file = 'complexAll/topROCMat.RData')
## save(topROCMat, file = 'complexAll/topROCMat_567.RData')
## save(topROCMat, file = 'complex/Bioinfor_topROCMat.RData')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~get whole genome interaction~~~~~~~~~~~~~~~
require('KEGGBioCycAPI')
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

sepMat <- CutSeqEqu(nrow(jaccardSim), 2000)

for (j in 1:ncol(sepMat)) {

  print(paste('it is running ', j, '.', sep = ''))
  jacSelectMat <- jaccardSim[sepMat[1, j]:sepMat[2, j], ]
  corSelectMat <- wholeCor[sepMat[1, j]:sepMat[2, j], ]

  topInterList <- topInter(jacSelectMat, corSelectMat, c(1, 400), 0, withValue = TRUE)
  rm(jacSelectMat, corSelectMat)
  gc()
  
  fileName <- paste('topInterList', sepMat[1, j], '_', sepMat[2, j], '.RData', sep = '')
  save(topInterList, file = fileName)
}

# to mat inter
topInterFiles <- dir(pattern = '^topInterList\\d')

for (i in 1:length(topInterFiles)) {
  load(topInterFiles[i])
  topInterMat <- list2matInterWithValue(topInterList)

  fileName <- paste(substr(topInterFiles[i], 1, 8), 'Mat', substring(topInterFiles[i], 13), sep = '')
  save(topInterMat, file = fileName)
}


# merged interMat and remove duplicated rows
topInterMatFiles <- dir(pattern = '^topInterMat\\d')
top400 <- foreach(i = 1:length(topInterMatFiles), .combine = rbind) %dopar% {
  load(topInterMatFiles[i])
  return(topInterMat)
}
top400 <- top400[!duplicated(top400, MARGIN = 1), ]
colnames(top400) <- c('From', 'To', 'Jaccard', 'Cor')
rownames(top400) <- NULL
save(top400, file = 'top400.RData')

# annotation
load('top400.RData')
load('wholePhyloData.RData')
fromAnno <- annoFirst[match(top400[, 1], names(annoFirst))]
toAnno <- annoFirst[match(top400[, 2], names(annoFirst))]
top400Anno <- cbind(fromAnno, toAnno)
top400Anno <- cbind(top400Anno, top400[, 3:4])
rownames(top400Anno) <- NULL
save(top400Anno, file = 'top400Anno.RData')
write.csv(top400Anno, 'top400Anno.csv')
write.table(top400Anno, 'top400Anno.txt', row.names = FALSE, quote = FALSE, sep = '\t')

atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
top400Anno[(top400Anno[, 1] %in% atpSub[, 1]) & (top400Anno[, 2] %in% atpSub[, 1]), ]
top400Anno[(top400Anno[, 1] %in% atpSub[, 1]) | (top400Anno[, 2] %in% atpSub[, 1]), ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################################################################


