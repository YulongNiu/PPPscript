#############################Useful functions#################################
InterTest <- function(targetInter, wholeInterList, ncor = 4) {
  # USE: Fisher exact test for selected interactions
  # INPUT: 'targetInter' is the target interactions with KEGG IDs, and it is a two-column matrix. 'wholeInterList' is a list of interactions, and each element is the two-column matrix. 'wholeInterList' should have names. 'targetInter' and 'wholeInterList' matrix should be orderd with the same rule. 'ncor' is the threads number.
  # OUTPUT: matrix of Fisher exact test
  # REF: http://david.abcc.ncifcrf.gov/content.jsp?file=functional_annotation.html

  require('foreach')
  require('doMC')
  registerDoMC(ncor)
  
  wholeVecList <- lapply(wholeInterList, function(x){
    eachVec <- apply(x, 1, paste, collapse = '|')
    return(eachVec)
  })
  wholeVec <- unique(unlist(wholeVecList))

  tarVec <- apply(targetInter, 1, paste, collapse = '|')
  tarVec <- tarVec[tarVec %in% wholeVec]

  # user's geneList in pathways
  userIn <- foreach (i = 1:length(wholeInterList), .combine = c) %dopar% {
    print(paste('It is running ', i, ' in total of ', length(wholeInterList), '.', sep = ''))
    numInEach <- sum(wholeVecList[[i]] %in% tarVec)
    return(numInEach)
  }
  
  # user's geneList not in pathways
  userOut <- length(tarVec) - userIn

  # genome in pathways
  genomeIn <- sapply(wholeVecList, length)

  # genome not in pathways
  genomeOut <- length(wholeVec) - genomeIn

  # Fisher Test
  FisherMat <- cbind(userIn, userOut, genomeIn, genomeOut)
  FisherVal <- apply(FisherMat, 1, function(x){
    fMat <- matrix(x, ncol = 2)
    fVal <- fisher.test(fMat)
    return(fVal$p.value)
  })
  adFisherVal <- p.adjust(FisherVal)

  # FisherRes
  FisherRes <- cbind(adFisherVal, FisherVal, userIn)

  # sort FisherRes
  FisherRes <- FisherRes[order(FisherRes[, 3], decreasing = TRUE), ]

  return(FisherRes)

}

OrderHumMat <- function(humMat) {
  # USE: order the human matrix (pathways, genesets) by the first and second rows
  # INPUT: 'humMat' is the human matrix. The first and second rows should be KEGGIDs, like 'hsa:1000'.
  # OUTPUT: ordered matrix

  fromNum <- sapply(strsplit(humMat[, 1], split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(humMat[, 2], split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  humMat <- rbind(humMat[isSmall, ,drop = FALSE], humMat[!isSmall, 2:1, drop = FALSE])

  return(humMat)
}
##############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~insert the KEGG pathway test~~~~~~~~~~~~~~
require('ggplot2')
load('pathway/pathFilter.RData')
load('complexAll/comInterJacCor.RData')
load('top400.RData')
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
# get FATP complex
top400ATP <- top400[(top400[, 1] %in% atpSub[, 2]) | (top400[, 2] %in% atpSub[, 2]), ]

# deal with complex
humComp <- lapply(MIPSHumList, function(x) {
  eachHumMat <- x[[2]][, 1:2, drop = FALSE]
  # order each mat
  eachHumMat <- OrderHumMat(eachHumMat)
  return(eachHumMat)
})

# pathway test for KEGG/Biocarta/NCI/reactome/complex databases
pathTestList <- list()
pathTestList$KEGG <- InterTest(top400ATP[, 1:2], keggPathFilter)
pathTestList$Biocarta <- InterTest(top400ATP[, 1:2], biocartaPathFilter)
pathTestList$NCI <- InterTest(top400ATP[, 1:2], nciPathFilter)
pathTestList$Reactome <- InterTest(top400ATP[, 1:2], reactomePathFilter)
pathTestList$Complex <- InterTest(top400ATP[, 1:2], humComp)

# select 
pathTestMat <- lapply(pathTestList, function(x){
  eachMat <- x[x[, 3] >= 5, , drop = FALSE]
  return(eachMat)
})
pathTestMat <- do.call(rbind, pathTestMat)
pathTestMat[, 1] <- -log10(pathTestMat[, 1])
rownames(pathTestMat)[2] <- c('Respiratory electron transport')
rownames(pathTestMat)[3] <- c('The citric acid cycle')
rownames(pathTestMat)[8] <- c('F1F0-ATP synthase')

# plot significant 
sigPath <- pathTestMat[pathTestMat[, 1] > -log10(0.05), c(1, 3)]
sigPathMat <- data.frame(pathName = rownames(sigPath), logFDR = sigPath[, 1], InterNum = sigPath[, 2])

pdf('ATPsub_inter.pdf', width = 8, height = 2.5)
ggplot(sigPathMat, aes(x = pathName, y = logFDR, label = InterNum)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  geom_text(hjust = -0.2) +
  scale_y_continuous(limits=c(0, 76)) +
  xlab('') + ylab('-logFDR')
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot sigPath~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# significant pathways are chose from 5 databases with
# 1. predicted number >= 15
# 2. predicted percentage >= 50%
# 3. ordered by percentage

require('ggplot2')
require('scales')
topSigPath <- read.csv('topSigPath.csv')
topSigMat <- data.frame(pathName = topSigPath[, 1],
                       pred = paste(topSigPath[,2], topSigPath[, 3], sep = '/'),
                       percentage = topSigPath[, 4],
                       database = factor(topSigPath[, 5]))


pdf('top_sig_pathway.pdf', width = 9, height = 10)
ggplot(topSigMat[1:60, ], aes(x = pathName, y = percentage, label = pred)) +
  geom_point(aes(color = database), size = 3) +
  coord_flip() +
  geom_text(y = 1.07, size = 3.5) +
  scale_y_continuous(limits = c(0.6, 1.1), breaks = seq(0.6, 1, 0.1), labels = percent_format()) +
  xlab('') +
  ylab('Predicted percentage')
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################################################################

########################add NPP test###############################
##~~~~~~~~~~~~~~~~~~~~~~~~load complex~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('complexAll/comInterJacCor.RData')
complexPathFilter <- lapply(MIPSHumList, function(x) {
  x <- x[[2]][, 1:2, drop = FALSE]
  colnames(x) <- c('from', 'to')
  return(x)
})
names(complexPathFilter) <- sapply(MIPSHumList, '[[', 4)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~NPP pathways~~~~~~~~~~~~~~~~~~
library('foreach')
library('doMC')
registerDoMC(4)

load('pathway/pathFilter.RData')
load('NPP70_profile.RData')
topSigPath <- read.csv('topSigPath.csv', stringsAsFactors = FALSE, row.names = 1)
topSigPath <- topSigPath[1:60, ]
topSigPath <- topSigPath[order(topSigPath[, 4]), ]
allNames <- rownames(topSigPath)
geneNames <- rownames(norProfile)

## BioCarta
pathNames1 <- allNames[topSigPath[, 4] == 'BioCarta']
eachPath1 <- biocartaPathFilter[pathNames1]
pathCor1 <- lapply(eachPath1, function(x) {
  eachPathCor <- foreach(i = seq_len(nrow(x)), .combine = c) %dopar% {
    eachCor <- cor(t(norProfile[geneNames %in% x[i, 1:2], ]))[1, 2]
    return(eachCor)
  }
  return(eachPathCor)
})

## Complex
pathNames1 <- allNames[topSigPath[, 4] == 'Complex']
eachPath1 <- complexPathFilter[pathNames1]
pathCor2 <- lapply(eachPath1, function(x) {
  eachPathCor <- foreach(i = seq_len(nrow(x)), .combine = c) %dopar% {
    eachCor <- cor(t(norProfile[geneNames %in% x[i, 1:2], ]))[1, 2]
    return(eachCor)
  }
  return(eachPathCor)
})

## KEGG
pathNames1 <- allNames[topSigPath[, 4] == 'KEGG']
eachPath1 <- keggPathFilter[pathNames1]
pathCor3 <- lapply(eachPath1, function(x) {
  eachPathCor <- foreach(i = seq_len(nrow(x)), .combine = c) %dopar% {
    eachCor <- cor(t(norProfile[geneNames %in% x[i, 1:2], ]))[1, 2]
    return(eachCor)
  }
  return(eachPathCor)
})

## Reactome
pathNames1 <- allNames[topSigPath[, 4] == 'Reactome']
eachPath1 <- reactomePathFilter[pathNames1]
pathCor4 <- lapply(eachPath1, function(x) {
  eachPathCor <- foreach(i = seq_len(nrow(x)), .combine = c) %dopar% {
    eachCor <- cor(t(norProfile[geneNames %in% x[i, 1:2], ]))[1, 2]
    return(eachCor)
  }
  return(eachPathCor)
})

pathCor <- c(pathCor1, pathCor2, pathCor3, pathCor4)

## threshold
## top 400 sensitivity 0.9670519
## NPP 0.7274617 sensitivity 0.9670519
thres <- 0.73
pathCorNum <- sapply(pathCor, function(x){return(sum(x > thres))})
write.csv(cbind(topSigPath, pathCorNum), file = 'test.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################################################################
