setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

##############################complex database#########################
## read MIPS complex database
MIPSCom <- read.csv('complexAll/allComplexes.csv', sep = ';')
## select human complex
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

save(MIPSHumList, file = 'complexAll/MIPSHumList.RData')
#######################################################################

##############################remove redondant subunits################
load('complexAll/MIPSHumList.RData')
load('wholePhyloData.RData')

# without anno
for (i in 1:length(MIPSHumList)) {
  entr <- MIPSHumList[[i]][[2]]
  entr <- entr[entr %in% names(annoFirst)]
  MIPSHumList[[i]][[2]] <- entr
}

# delete with one unit
len <- sapply(MIPSHumList, function(x) {
  length(x[[2]])
})
MIPSHumList <- MIPSHumList[which(len != 1)]

## sink('complexAll/MIPSHumListNoRed.txt')
## MIPSHumList
## sink()

save(MIPSHumList, file = 'complexAll/MIPSHumListNoRed.RData')
#######################################################################

###############################Random select###########################
load('complexAll/MIPSHumListNoRed.RData')

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
save(randomMat, file = 'complexAll/randomMat.RData')
######################################################################

########################Remove big complexes#########################
load('complexAll/MIPSHumListNoRed.RData')

lenAll <- sapply(MIPSHumList, function(x) length(x[[2]]))

## all possible complex interactions
complexInter <- lapply(MIPSHumList, function(x) {
  possCom <- t(combn(x[[2]], 2))
  return(possCom)
})

## remove big complexes
cutoffNum <- 30
complexInterCut <- complexInter[lenAll <= cutoffNum]
interMat <- do.call(rbind, complexInterCut)

## ## remove duplicate inters
## interMat <- t(apply(interMat, 1, sort))
## interMat <- interMat[!duplicated(interMat), ]

save(interMat, file = 'complexAll/interMat_30.RData')
######################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~test Jac/Cor/decomposed~~~~~~~~~~~~~~~~~~~~~~~
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

## sink('comInterJacCor.txt')
## MIPSHumList
## sink()

save(MIPSHumList, file = 'comInterJacCor.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################



#~~~~~~~~~~~~~~~~~~~~~~~~~~ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')
load('complexAll/randomMat.RData')
load('complexAll/comInterJacCor.RData')

# TP
PRSJacCor <- sapply(MIPSHumList, '[[', 2)
PRSJacCor <- do.call(rbind, PRSJacCor)
PRSJacCor <- PRSJacCor[, 1:4]

# TN
set.seed(123)
# second random
# set.seed(567)
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
save(PRSJacCor, NRSJacCor, file = 'complexAll/our_complex.RData')
# save(PRSJacCor, NRSJacCor, file = 'complexAll/our_complex_567.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('ggplot2')
load('complexAll/our_complex.RData')
## load('complexAll/our_complex_567.RData')
## load('complex/Bioinfor_complex.RData')
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

save(jacsimROCMat, file = 'complexAll/jacsimROCMat.RData')
save(corROCMat, file = 'complexAll/corROCMat.RData')
## save(jacsimROCMat, file = 'complexAll/jacsimROCMat_567.RData')
## save(corROCMat, file = 'complexAll/corROCMat_567.RData')
## save(jacsimROCMat, file = 'complex/Bioinfor_jacsimROCMat.RData')
## save(corROCMat, file = 'complex/Bioinfor_corROCMat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


