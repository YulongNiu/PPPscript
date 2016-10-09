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


##################################Complex################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~complex database~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read MIPS complex database
MIPSCom <- read.csv('complex/allComplexesCore.csv', sep = ';')
# select human complex
MIPSHum <- MIPSCom[MIPSCom[, 4] == 'Human', ]

MIPSHumList <- list()
for (i in 1:nrow(MIPSHum)) {
  # get entrez number
  # delet '(' and ')'
  mergedEntr <- MIPSHum[i, 6]
  mergedEntr <- sub('\\(', '', mergedEntr)
  mergedEntr <- sub('\\)', '', mergedEntr)
  entr <- unlist(strsplit(as.character(mergedEntr), split = ',', fixed = TRUE))
  entr <- paste('hsa:', entr, sep = '')
  entr <- unique(entr)
  # get uniprot
  mergedUniprot <- MIPSHum[i, 5]
  mergedUniprot <- sub('\\(', '', mergedUniprot)
  mergedUniprot <- sub('\\)', '', mergedUniprot)
  uniprot <- unlist(strsplit(as.character(mergedUniprot), split = ',', fixed = TRUE))
  uniprot <- unique(uniprot)
  # complex annotation
  anno <- unlist(strsplit(as.character(MIPSHum[i, 9]), split = ',', fixed = TRUE))
  # complex name
  nameCom <- as.character(MIPSHum[i, 2])
  MIPSHumList[[i]] <- list(Uniprot = uniprot, Entrez = entr, Annotation = anno, Name = nameCom)
}

# delete with one unit
len <- sapply(MIPSHumList, function(x) {
  length(x[[2]])
})
MIPSHumList <- MIPSHumList[which(len != 1)]

names(MIPSHumList) <- paste('Complex', 1:length(MIPSHumList), sep = '')

sink('MIPSHumList.txt')
MIPSHumList
sink()

save(MIPSHumList, file = 'MIPSHumList.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~add JacardSim and Cor~~~~~~~~~~~~~~~~~~~~~
load('complex/MIPSHumList.RData')
load('wholePhyloData.RData')

# without anno
for (i in 1:length(MIPSHumList)) {
  entr <- match(MIPSHumList[[i]][[2]], names(annoFirst))
  entr <- entr[!is.na(entr)]
  MIPSHumList[[i]][[2]] <- names(annoFirst)[entr]
}

# delete with one unit
len <- sapply(MIPSHumList, function(x) {
  length(x[[2]])
})
MIPSHumList <- MIPSHumList[which(len != 1)]
save(MIPSHumList, file = 'MIPSHumListNoRed.RData')

#~~~~~~~~~~~~tranfer complex to 'gene-list' interaction format~~~~~~~~~~~~~~~
require('foreach')
require('doMC')
registerDoMC(4)

# translate to 'gene-list' interaction formate
complexList <- sapply(MIPSHumList, '[[', 2)
interVec <- unique(unlist(complexList))

complexInterList <- foreach (i = 1:length(interVec)) %dopar% {
  hasGene <- sapply(complexList, function(x) {
    return(interVec[i] %in% x)
  })
  hasGeneList <- complexList[which(hasGene)]
  hasGeneVec <- unique(unlist(hasGeneList))
  # remove self
  hasGeneVec <- hasGeneVec[hasGeneVec != interVec[i]]
  return(hasGeneVec)
}
names(complexInterList) <- interVec
save(complexInterList, file = 'complexInterList.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~test Jac/Cor/decomposed~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

# all possible combination
for (i in 1:length(MIPSHumList)) {
  possCom <- t(combn(MIPSHumList[[i]][[2]], 2))
  MIPSHumList[[i]][[2]] <- possCom
}

# jaccard similarity
for (i in 1:length(MIPSHumList)) {
  print(paste('Now it is running ', i, ' in a total of ', length(MIPSHumList), '.', sep = ''))
  
  comInter <- MIPSHumList[[i]][[2]]
  
  comInterJacCor <- foreach(i = 1:nrow(comInter), .combine = rbind) %dopar% {
    rowNum <- match(comInter[i, 1], rownames(jaccardSim))
    colNum <- match(comInter[i, 2], rownames(jaccardSim))
    jacsim <- jaccardSim[rowNum, colNum]
    corValue <- wholeCor[rowNum, colNum]
    interMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% comInter[i, 1:2], ])
    cateVec <- unlist(CatePhylo(interMat, domain))
    jacCorVec <- c(comInter[i, ], jacsim, corValue, cateVec)
    return(jacCorVec)
  }

  # deal with row number is 1
  if (!is.matrix(comInterJacCor)) {
    comInterJacCor <- matrix(comInterJacCor, ncol = 16)
  } else {}
  
  colnames(comInterJacCor) <- c('proteinA', 'proteinB', 'JaccardSim', 'Cor', paste('Bacteria', c('11', '10', '01', '00'), sep = ''), paste('Eukaryotes', c('11', '10', '01', '00'), sep = ''), paste('Archaea', c('11', '10', '01', '00'), sep = ''))
  MIPSHumList[[i]][[2]] <- comInterJacCor
}

sink('comInterJacCor.txt')
MIPSHumList
sink()

save(MIPSHumList, file = 'comInterJacCor.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########################################################################



###############################Random select###############################
load('complex/MIPSHumListNoRed.RData')

getAllComb <- function(index1, index2, complexMat){
  # USE: all the combination
  vec1 <- complexMat[[index1]][[2]]
  vec2 <- complexMat[[index2]][[2]]

  combMat <- expand.grid(vec1, vec2)
  return(combMat)
}

getPetComb <- function(firstPet, secondPet, complexMat) {
  patIndex <- expand.grid(firstPet, secondPet)
  randomMat <- apply(patIndex, 1, function(x){
    return(getAllComb(x[1], x[2], complexMat))
  })
  randomMat <- do.call(rbind, randomMat)
}

# get complex annotation
comAnno <- sapply(MIPSHumList, '[[', 3)

has70List <- lapply(comAnno, function(x){
  has70 <- x[grepl('^70', x)]
  return(has70)
})
## table(unlist(has70List))


# cell membrane
pat1 <- c('70.02', '70.02.01')
pat1 <- sapply(has70List, function(x){sum(pat1 %in% x)})
pat1 <- which(pat1 > 0)
# cell junctions
pat2 <- c('70.06')
pat2 <- sapply(has70List, function(x){sum(pat2 %in% x)})
pat2 <- which(pat2 > 0)
# mito
pat3 <- c('70.16', '70.16.01', '70.16.03', '70.16.05', '70.16.07')
pat3 <- sapply(has70List, function(x){sum(pat3 %in% x)})
pat3 <- which(pat3 > 0)
# nuclear
pat4 <- c('70.10', '70.10.03', '70.10.04', '70.10.05', '70.10.06', '70.10.07', '70.10.09')
pat4 <- sapply(has70List, function(x){sum(pat4 %in% x)})
pat4 <- which(pat4 > 0)
# extracellular matrix component
pat5 <- c('70.27.01')
pat5 <- sapply(has70List, function(x){sum(pat5 %in% x)})
pat5 <- which(pat5 > 0)

patMat <- matrix(c('pat2', 'pat3',
                   'pat5', 'pat3',
                   'pat2', 'pat4',
                   'pat1', 'pat4',
                   'pat3', 'pat5',
                   'pat3', 'pat4'), ncol = 2, byrow = TRUE)

randomMat <- apply(patMat, 1, function(x){
  selectOrder <- paste('getPetComb(', x[1], ',', x[2], ',MIPSHumList)', sep = '')
  selectMat <- eval(parse(text = selectOrder))
  return(selectMat)
})

randomMat <- do.call(rbind, randomMat)

notSelf <- apply(randomMat, 1, function(x) {
  if (x[1] == x[2]) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
randomMat <- randomMat[notSelf, ]
save(randomMat, file = 'randomMat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')
load('complex/randomMat.RData')
load('complex/comInterJacCor.RData')

# TP
PRSJacCor <- sapply(MIPSHumList, '[[', 2)
PRSJacCor <- do.call(rbind, PRSJacCor)
PRSJacCor <- PRSJacCor[, 1:4]

# TN
set.seed(123)
set.seed(567)
NRS <- randomMat[sample(1:nrow(randomMat), nrow(PRSJacCor)), ]
# whole genome entrz
entrez <- rownames(jaccardSim)
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
save(PRSJacCor, NRSJacCor, file = 'our_complex.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('ggplot2')
load('complex/our_complex.RData')
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

mergedROCMat <- rbind(jacsimROCMat, corROCMat)
mergedROCMat <- cbind(mergedROCMat, Cuttype = rep(c('jaccard', 'cor'), c(nrow(jacsimROCMat), nrow(corROCMat))))

# add top TP/NP
## mergedROCMat <- rbind(jacsimROCMat, corROCMat, ROCMat1)
## mergedROCMat <- cbind(mergedROCMat, Cuttype = rep(c('jaccard', 'cor', 'top'), c(nrow(jacsimROCMat), nrow(corROCMat), nrow(ROCMat1))))

ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Cuttype)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################test accuaacy###############################
require('ggplot2')
load('complexAll/our_complex.RData')

cutVec <- seq(0, 1, 0.01)
# jaccard similarity
jacsimVec <- c(PRSJacCor[, 'JaccardSim'], NRSJacCor[, 'JaccardSim'])
jacsim <- data.frame(Jaccard = as.numeric(jacsimVec), Status = rep(c('pos', 'neg'), c(nrow(PRSJacCor), nrow(NRSJacCor))))

jacsimAcc <- sapply(cutVec, function(x) {
  statusVec <- jacsim[jacsim[, 'Jaccard'] >= x, 2]
  acc <- sum(statusVec == 'pos')/length(statusVec)

  return(cbind(x, acc))
})
jacsimAcc <- t(jacsimAcc)

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

corAcc <- sapply(cutVec, function(x) {
  statusVec <- corMat[corMat[, 'Cor'] >= x, 2]
  acc <- sum(statusVec == 'pos')/length(statusVec)

  return(cbind(x, acc))
})
corAcc <- t(corAcc)

mergedAcc <- rbind(jacsimAcc, corAcc)
mergedAcc <- data.frame(Threshold = mergedAcc[, 1],
                        Accuracy = mergedAcc[, 2] * 100,
                        Cuttype = rep(c('jaccard', 'cor'), c(nrow(jacsimAcc), nrow(corAcc))))

ggplot(data = mergedAcc, mapping = aes(x = Threshold, y = Accuracy, colour = Cuttype)) +
  geom_line() +
  scale_y_continuous(limits=c(0, 100))

#########################################################################

#########################################################################
#######################test accuracy of percentage method################
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
load('complex/complexInterList.RData')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
complexNum <- match(names(complexInterList), rownames(jaccardSim))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cutoff set
linkNum <- 400
corSet <- 0

## jacTopInterList <- foreach(i = 1:length(complexNum)) %dopar% {
## ## jacTopInterList <- foreach(i = 1:20) %dopar% {
##   print(paste('It is running ', i, '.', sep = ''))
##   # set corSet
##   rowNum <- complexNum[i]
##   corSetNum <- which(wholeCor[rowNum, ] > corSet)
##   jacSimVec <- jaccardSim[rowNum, ]
##   jacSimVec <- jacSimVec[corSetNum]

##   # delete self nodes
##   jacSimVec <- jacSimVec[names(jacSimVec) != rownames(jaccardSim)[rowNum]]
  
##   if (length(jacSimVec) < linkNum && length(jacSimVec) > 0) {
##     jacSimInter <- jacSimVec
##   }
##   else if (length(jacSimVec) == 0) {
##     jacSimInter <- NA
##   }
##   else if (length(jacSimVec) >= linkNum) {
##     ## sort jacSimVec
##     jacSimVec <- sort(jacSimVec, decreasing = TRUE)
##     jacSimInter <- jacSimVec[1:linkNum]
##   }
  
##   return(jacSimInter)
## }


jacTopInterList <- list()
for(i in 1:length(complexNum)) {
  print(paste('It is running ', i, '.', sep = ''))
  # set corSet
  rowNum <- complexNum[i]
  corSetNum <- which(wholeCor[rowNum, ] > corSet)
  jacSimVec <- jaccardSim[rowNum, ]
  jacSimVec <- jacSimVec[corSetNum]

  # delete self nodes
  jacSimVec <- jacSimVec[names(jacSimVec) != rownames(jaccardSim)[rowNum]]
  
  if (length(jacSimVec) < linkNum && length(jacSimVec) > 0) {
    jacSimInter <- jacSimVec
  }
  else if (length(jacSimVec) == 0) {
    jacSimInter <- NA
  }
  else if (length(jacSimVec) >= linkNum) {
    ## sort jacSimVec
    jacSimVec <- sort(jacSimVec, decreasing = TRUE)
    jacSimInter <- jacSimVec[1:linkNum]
  }

  jacTopInterList[[i]] <- jacSimInter
  
}

names(jacTopInterList) <- rownames(jaccardSim)[complexNum]


accVec <- rep(NA, length(jacTopInterList))
for(i in 1:length(jacTopInterList)) {
  accNum <- sum(names(jacTopInterList[[i]]) %in% complexInterList[[i]])
  accVec[i] <- accNum
}

sum(accVec) / sum(sapply(complexInterList, length))
#########################################################################





























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
