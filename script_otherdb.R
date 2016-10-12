############# transfer Entrez to Uniprot##############################
require('KEGGBioCycAPI')
require('foreach')
require('doMC')
registerDoMC(4)
load('wholePhyloData.RData')

# transfer entrez to uniprot
uniProt <- KEGGConv('uniprot', 'hsa')

# select phylogenetic profiles
entrz <- names(annoFirst)
uniProt <- uniProt[uniProt[, 1] %in% entrz, ]

save(uniProt, file = 'uniProt.RData')
#######################################################################

#########################preprocess TP/TN interaction##########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~CCSB TP~~~~~~~~~~~~~~~~~~~~~~~~
# CCSB database
CCSBTP <- read.csv('Interact/HI2_2011.tsv', sep = '\t')
CCSBTP <- CCSBTP[, c(1, 3)]
CCSBTP <- apply(CCSBTP, 1:2, function(x){paste('hsa:', x, sep = '')})

load('wholePhyloData.RData')
# without anno
fromNA <- is.na(match(CCSBTP[, 1], names(annoFirst)))
toNA <- is.na(match(CCSBTP[, 2], names(annoFirst)))
totalNA <- fromNA | toNA
CCSBTP <- CCSBTP[!totalNA, ]


load('jaccardSim.RData')
require('foreach')
require('doMC')
registerDoMC(4)

CCSBTPJacSim <- foreach(i = 1:nrow(CCSBTP), .combine = rbind) %dopar% {
  rowNum <- match(CCSBTP[i, 1], rownames(jaccardSim))
  colNum <- match(CCSBTP[i, 2], rownames(jaccardSim))
  jacsim <- jaccardSim[rowNum, colNum]
  jacsimVec <- c(CCSBTP[i, ], jacsim)
  return(jacsimVec)
}

colnames(CCSBTPJacSim)[3] <- 'JaccardSim'

load('wholeCor.RData')

CCSBTPJacCor <- foreach(i = 1:nrow(CCSBTPJacSim), .combine = rbind) %dopar% {
  rowNum <- match(CCSBTP[i, 1], rownames(wholeCor))
  colNum <- match(CCSBTP[i, 2], rownames(wholeCor))
  cor <- wholeCor[rowNum, colNum]
  corVec <- c(CCSBTPJacSim[i, ], cor)
  return(corVec)
}

colnames(CCSBTPJacCor)[4] <- 'Cor'
write.csv(CCSBTPJacCor, 'CCSBTPJacCor.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~Negatome 2.0~~~~~~~~~~~~~~~~~~~~~~~~~~
require('foreach')
require('doMC')
registerDoMC(4)
load('uniProt.RData')
negTN <- read.csv('Interact/manualTN2.txt', sep = '\t', header = FALSE)
## negTN <- read.csv('Interact/manualTN1.txt', sep = '\t')
## negTN <- read.csv('Interact/manualTNPDB1.txt', sep = '\t')
negTN <- negTN[, 1:2]

uniProtIDs <- sapply(strsplit(uniProt[, 2], split = ':', fixed = TRUE), '[[', 2)
uniProt[, 2] <- uniProtIDs
col1Trans <- match(negTN[, 1], uniProt[, 2])
col2Trans <- match(negTN[, 2], uniProt[, 2])
negTN <- cbind(col1Trans, col2Trans)

# without anno
fromNA <- is.na(col1Trans)
toNA <- is.na(col2Trans)
totalNA <- fromNA | toNA
negTN <- negTN[!totalNA, ]

# without selfinter
notSelf <- apply(negTN, 1, function(x) {
  if (x[1] == x[2]) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
negTN <- negTN[notSelf, ]


load('jaccardSim.RData')
require('foreach')
require('doMC')
registerDoMC(4)

negTNJacSim <- foreach(i = 1:nrow(negTN), .combine = rbind) %dopar% {
  rowNum <- match(negTN[i, 1], rownames(jaccardSim))
  colNum <- match(negTN[i, 2], rownames(jaccardSim))
  jacsim <- jaccardSim[rowNum, colNum]
  jacsimVec <- c(negTN[i, ], jacsim)
  return(jacsimVec)
}

colnames(negTNJacSim)[3] <- 'JaccardSim'

rm(jaccardSim)
gc()

load('wholeCor.RData')

negTNJacCor <- foreach(i = 1:nrow(negTNJacSim), .combine = rbind) %dopar% {
  rowNum <- match(negTN[i, 1], rownames(wholeCor))
  colNum <- match(negTN[i, 2], rownames(wholeCor))
  cor <- wholeCor[rowNum, colNum]
  corVec <- c(negTNJacSim[i, ], cor)
  return(corVec)
}

colnames(negTNJacCor)[4] <- 'Cor'
write.csv(negTNJacCor, 'negTNJacCor.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('ggplot2')
posInt <- read.csv('Interact/CCSBTPJacCor.csv', row.names = 1)
posInt <- cbind(posInt, rep('pos', nrow(posInt)))
set.seed(123)
posInt <- posInt[sample(1:nrow(posInt), 1000), ]
colnames(posInt) <- c('geneA', 'geneB', 'JaccardSim', 'Cor', 'Fact')


negInt <- read.csv('Interact/negTNJacCor2.csv', row.names = 1)
negInt <- cbind(negInt, rep('neg', nrow(negInt)))
set.seed(123)
negInt <- negInt[sample(1:nrow(negInt), 1000), ]
colnames(negInt) <- c('geneA', 'geneB', 'JaccardSim', 'Cor', 'Fact')

cutVec <- seq(0, 1, 0.01)
# jaccard similarity
jacsimVec <- c(posInt[, 'JaccardSim'], negInt[, 'JaccardSim'])
jacsim <- data.frame(Jaccard = as.numeric(jacsimVec), Status = rep(c('pos', 'neg'), c(nrow(posInt), nrow(negInt))))
jacsimROCMat <- SingROCMat(cutVec, jacsim)
# correlation coefficient
corVec <- c(posInt[, 'Cor'], negInt[, 'Cor'])
corMat <- data.frame(Cor = corVec, Status = rep(c('pos', 'neg'), c(nrow(posInt), nrow(negInt))))
corROCMat <- SingROCMat(cutVec, corMat)

mergedROCMat <- rbind(jacsimROCMat, corROCMat)
mergedROCMat <- cbind(mergedROCMat, Cuttype = rep(c('jaccard', 'cor'), c(nrow(jacsimROCMat), nrow(corROCMat))))
ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Cuttype)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################################################################


#########################validate bioinfor 2011###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~deal with TP and TN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
load('uniProt.RData')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

# true positive interaction
PRS <- read.csv('complex/data_files/PRS_scores_human.tab', sep = '\t')
PRS <- PRS[, 1:2]

# whole genome entrz
entrez <- rownames(jaccardSim)

uniProtIDs <- sapply(strsplit(uniProt[, 2], split = ':', fixed = TRUE), '[[', 2)
uniProt[, 2] <- uniProtIDs
col1Trans <- match(PRS[, 1], uniProt[, 2])
col2Trans <- match(PRS[, 2], uniProt[, 2])
PRS <- cbind(col1Trans, col2Trans)

# without anno
fromNA <- is.na(col1Trans)
toNA <- is.na(col2Trans)
totalNA <- fromNA | toNA
PRS <- PRS[!totalNA, ]

# without selfinter
notSelf <- apply(PRS, 1, function(x) {
  if (x[1] == x[2]) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
PRS <- PRS[notSelf, ]
PRS <- apply(PRS, 1:2, function(x) {
  return(uniProt[x, 1])
})


# jaccard
PRSindex <- apply(PRS, 1:2, function(x) {
  return(match(x, entrez))
})

jacsim <- apply(PRSindex, 1, function(x){
  return(jaccardSim[x[1], x[2]])
})

corValue <- apply(PRSindex, 1, function(x){
  return(wholeCor[x[1], x[2]])
})

PRSJacCor <- cbind(PRS, jacsim, corValue)
colnames(PRSJacCor) <- c('proteinA', 'proteinB', 'JaccardSim', 'Cor')



# true negative interaction
NRS <- read.csv('complex/data_files/NRS_scores_human.tab', sep = '\t')
NRS <- NRS[, 1:2]

col1Trans <- match(NRS[, 1], uniProt[, 2])
col2Trans <- match(NRS[, 2], uniProt[, 2])
NRS <- cbind(col1Trans, col2Trans)

# without anno
fromNA <- is.na(col1Trans)
toNA <- is.na(col2Trans)
totalNA <- fromNA | toNA
NRS <- NRS[!totalNA, ]

# without selfinter
notSelf <- apply(NRS, 1, function(x) {
  if (x[1] == x[2]) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
NRS <- NRS[notSelf, ]
NRS <- apply(NRS, 1:2, function(x) {
  return(uniProt[x, 1])
})

# jaccard
NRSindex <- apply(NRS, 1:2, function(x) {
  return(match(x, entrez))
})

jacsim <- apply(NRSindex, 1, function(x){
  return(jaccardSim[x[1], x[2]])
})

corValue <- apply(NRSindex, 1, function(x){
  return(wholeCor[x[1], x[2]])
})

NRSJacCor <- cbind(NRS, jacsim, corValue)
colnames(NRSJacCor) <- c('proteinA', 'proteinB', 'JaccardSim', 'Cor')
save(PRSJacCor, NRSJacCor, file = 'Bioinfor_complex.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROC curve~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('ggplot2')
load('complex/Bioinfor_complex.RData')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~plot histogram~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# jaccard plot
jacsimVec <- c(PRSJacCor[, 'JaccardSim'], NRSJacCor[, 'JaccardSim'])
jacsim <- data.frame(Jaccard = as.numeric(jacsimVec), Status = rep(c('TP', 'TN'), c(nrow(PRSJacCor), nrow(NRSJacCor))))

p <- ggplot(data = jacsim, aes(x = Jaccard))
p +
  geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = factor(Status))) +
  stat_density(geom = 'line', position = 'identity', aes(colour = factor(Status)))


corVec <- c(PRSJacCor[, 'Cor'], NRSJacCor[, 'Cor'])
corVec <- sapply(as.numeric(corVec), function(x){
  if (x < 0) {
    x <- 0
  } else {
    x <- x
  }
  return(x)
})
corMat <- data.frame(Cor = as.numeric(corVec), Status = rep(c('TP', 'TN'), c(nrow(PRSJacCor), nrow(NRSJacCor))))

p <- ggplot(data = corMat, aes(x = Cor))
p +
  geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = factor(Status))) +
  stat_density(geom = 'line', position = 'identity', aes(colour = factor(Status)))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


SingROCMat <- function(ctVec, tptnMat) {
  # USE: use single threshold to generate ROC matrix.
  # INPUT: 'ctVec' is the vector of cutoff vector. 'tptnMat' is the true positive and true negative matrix.
  # OUTPUT: return the 'TPR' and 'FPR' matrix
  
  ROCMat <- sapply(ctVec, function(x){
    predMat <- table(tptnMat[, 1] >= x, tptnMat[, 2])
    
    # deal with 1 row
    if (nrow(predMat) == 1) {
      if (rownames(predMat) == 'TRUE') {
        predMat <- rbind(predMat, c(0, 0))
        rownames(predMat)[2] <- 'FALSE'
      }
      else if (rownames(predMat) == 'FALSE') {
        predMat <- rbind(predMat, c(0, 0))
        rownames(predMat)[2] <- 'TRUE'
      }
    } else {}
    
    TP <- predMat['TRUE', 'pos']
    TN <- predMat['FALSE', 'neg']
    FP <- predMat['TRUE', 'neg']
    FN <- predMat['FALSE', 'pos']
    TPR <- TP/(TP+FN)
    FPR <- FP/(FP+TN)
    return(c(TPR, FPR))
  })

  ROCMat <- t(ROCMat)
  colnames(ROCMat) <- c('TPR', 'FPR')
  rownames(ROCMat) <- ctVec
  ROCMat <- as.data.frame(ROCMat)
  return(ROCMat)
}

DoubleROCMat <- function(jcct, ctVec, tptnMat) {
  # USE: use single threshold to generate ROC matrix.
  # INPUT: 'jcct' is the cutoff of jaccard similarity. 'ctVec' is the vector of cutoff vector. 'tptnMat' is the true positive and true negative matrix.
  # OUTPUT: return the 'TPR' and 'FPR' matrix
  
  ROCMat <- sapply(ctVec, function(x){
    predMat <- table(tptnMat[, 1] > jcct & tptnMat[, 2] > x, tptnMat[, 3])
    
    # deal with 1 row
    if (nrow(predMat) == 1) {
      if (rownames(predMat) == 'TRUE') {
        predMat <- rbind(predMat, c(0, 0))
        rownames(predMat)[2] <- 'FALSE'
      }
      else if (rownames(predMat) == 'FALSE') {
        predMat <- rbind(predMat, c(0, 0))
        rownames(predMat)[2] <- 'TRUE'
      }
    } else {}
    
    TP <- predMat['TRUE', 'pos']
    TN <- predMat['FALSE', 'neg']
    FP <- predMat['TRUE', 'neg']
    FN <- predMat['FALSE', 'pos']
    TPR <- TP/(TP+FN)
    FPR <- FP/(FP+TN)
    return(c(TPR, FPR))
  })

  ROCMat <- t(ROCMat)
  colnames(ROCMat) <- c('TPR', 'FPR')
  rownames(ROCMat) <- ctVec
  ROCMat <- as.data.frame(ROCMat)
}

cutVec <- seq(0, 1, 0.01)
# jaccard similarity
jacsimVec <- c(PRSJacCor[, 'JaccardSim'], NRSJacCor[, 'JaccardSim'])
jacsim <- data.frame(Jaccard = as.numeric(jacsimVec), Status = rep(c('pos', 'neg'), c(nrow(PRSJacCor), nrow(NRSJacCor))))
jacsimROCMat <- SingROCMat(cutVec, jacsim)
# correlation coefficient
corVec <- c(PRSJacCor[, 'Cor'], NRSJacCor[, 'Cor'])
corVec <- sapply(as.numeric(corVec), function(x){
  if (x < 0) {
    x <- 0
  } else {
    x <- x
  }
  return(x)
})
corMat <- data.frame(Cor = as.numeric(corVec), Status = rep(c('pos', 'neg'), c(nrow(PRSJacCor), nrow(NRSJacCor))))
corROCMat <- SingROCMat(cutVec, corMat)

## # combine two method
## twoMat <- data.frame(Jaccard = as.numeric(jacsimVec), Cor = as.numeric(corVec), Status = rep(c('pos', 'neg'), c(nrow(PRSJacCor), nrow(NRSJacCor))))
## twoROCMat <- lapply(seq(0, 1, 0.1), function(x){
##   return(DoubleROCMat(x, cutVec, twoMat))
## })
## twoROCMat <- do.call(rbind, twoROCMat)
## ggplot(data = twoROCMat, mapping = aes(x = FPR, y = TPR)) +
##   geom_line() +
##   geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed")


mergedROCMat <- rbind(jacsimROCMat, corROCMat)
mergedROCMat <- cbind(mergedROCMat, Cuttype = rep(c('jaccard', 'cor'), c(nrow(jacsimROCMat), nrow(corROCMat))))
ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Cuttype)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed")

## ggplot(data = mergedROCMat, mapping = aes(x = 1-FPR, y = TPR, colour = Cuttype)) +
##   geom_line() +
##   geom_abline(intercept = 1, slope = -1, colour="grey", linetype = "dashed")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################################


##############################human GO################################
# human GO construct version:10/29/2014
# http://geneontology.org/page/download-annotations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~preprocess~~~~~~~~~~~~~~~~~~~~~~~~~~
GORaw <- read.csv('GO/gene_association.goa_human', skip = 11, header = FALSE, sep = '\t')
GORaw <- GORaw[, c(2, 5, 9)]
load('uniProt.RData')

uniProtIDs <- sapply(strsplit(uniProt[, 2], split = ':', fixed = TRUE), '[[', 2)
uniprotIndex <- match(GORaw[, 1], uniProtIDs)
GORaw <- cbind(GORaw, uniprotIndex)

# remove NA
GOMat <- GORaw[!is.na(uniprotIndex), ]
entrezVec <- uniProt[GOMat[, 4], 1]
GOMat <- cbind(GOMat, entrezVec)

# rearrage GOMat
GOMat <- GOMat[, c(5, 2, 3)]
GOMat <- GOMat[!duplicated(GOMat, MARGIN = 1), ]
GOMat <- apply(GOMat, 1:2, as.character)
GOMat <- unname(GOMat)


# BP
humanBP <- GOMat[GOMat[, 3] == 'P', ]
humanBPList <- split(humanBP[, 1], factor(humanBP[, 2]))
humanBPListRev <- split(humanBP[, 2], factor(humanBP[, 1]))
# CC
humanCC <- GOMat[GOMat[, 3] == 'C', ]
humanCCList <- split(humanCC[, 1], factor(humanCC[, 2]))
humanCCListRev <- split(humanCC[, 2], factor(humanCC[, 1]))
# MF
humanMF <- GOMat[GOMat[, 3] == 'F', ]
humanMFList <- split(humanMF[, 1], factor(humanMF[, 2]))
humanMFListRev <- split(humanMF[, 2], factor(humanMF[, 1]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################
