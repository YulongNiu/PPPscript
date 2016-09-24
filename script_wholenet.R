#!/usr/bin/Rscript --vanilla

############# Preprocess the phylogenetic profiles for network ################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ process for network ~~~~~~~~~~~~~~~~~~~~~~~
# annotation 
load('wholePhyloData0001.RData')
wholePhyloDataAnnoRaw <- wholePhyloData0001
wholePhyloDataAnno <- as.character(wholePhyloDataAnnoRaw[, 2])
names(wholePhyloDataAnno) <- as.character(wholePhyloDataAnnoRaw[1:nrow(wholePhyloDataAnnoRaw), 1])

# phylo 0-1 data
wholePhyloDataNet <- wholePhyloDataAnnoRaw[, -1:-2]
wholePhyloDataNet <- apply(wholePhyloDataNet, 1:2, as.numeric)
# add hsa data
hsaNum <- ncol(wholePhyloDataNet) + 1
wholePhyloDataNet <- cbind(wholePhyloDataNet, rep(1, nrow(wholePhyloDataNet)))
colnames(wholePhyloDataNet)[hsaNum] <- 'hsa'


# If the phyloData is all 0, remove this row
preSum <- apply(wholePhyloDataNet, 1, sum)
hasZero <- preSum == 0
wholePhyloDataNet <- wholePhyloDataNet[!hasZero, ]
wholePhyloDataAnno <- wholePhyloDataAnno[!hasZero]

# get first anno name
getcontent <- function(s,g) {
    substring(s,g,g+attr(g,'match.length')-1)
  }
annoGrep <- regexpr('.+?(,|;|$)', wholePhyloDataAnno)
annoFirstRaw <- getcontent(wholePhyloDataAnno, annoGrep)
library('bigmemory')
library('foreach')
library('doMC')
registerDoMC(4)
annoFirst <- foreach(i = 1:length(annoFirstRaw), .combine = c) %dopar% {
  if (grepl('(,|;)', annoFirstRaw[i])) {
    annoClean <- substr(annoFirstRaw[i], 1, nchar(annoFirstRaw[i]) - 1)
    names(annoClean) <- names(annoFirstRaw[i])
  } else {
    annoClean <- annoFirstRaw[i]
    names(annoClean) <- names(annoFirstRaw[i])
  }
  return(annoClean)
}


repeatGenes <- sort(annoFirst[annoFirst %in% annoFirst[duplicated(annoFirst)]])
repeatGenesTab <- table(repeatGenes)
repeatGenesNum <- foreach(i = 1:length(repeatGenesTab), .combine = c) %dopar% {
  return(1:repeatGenesTab[i])
}
repeatGenesUni <- paste(repeatGenesNum, repeatGenes, sep = '_')
names(repeatGenesUni) <- names(repeatGenes)
annoFirst[match(names(repeatGenesUni), names(annoFirst))] <- repeatGenesUni

save(wholePhyloDataNet, wholePhyloDataAnno, annoFirst, file = 'wholePhyloData.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################################


######################### Construct Correlation network ###################
## adj2ft <- function(adMat){

##   # USE: this function is used to convert adjacency matrix to from-to matrix
##   # INPUT: 'adMat' is the input adjacency matrix
##   # OUTPUT: from-to matrix to denote the undirected network

##   require(graph)
##   adMatNEL <- as(adMat, 'graphNEL')
##   adBAM <- as(adMatNEL, 'graphBAM')
##   adft <- extractFromTo(adBAM)
##   return(adft)
## }

## adj2ftBig <- function(adMat, adAllRowCol){
  
##   # a big matrix version
##   # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also a bigmatrix.
  
##   library(bigmemory)

##   rowNames <- rownames(adMat)
##   colNames <- colnames(adMat)

##   big.matrix(nrow = nrow(adAllRowCol), ncol = 3, type = 'char', backingfile = 'adft.bin', descriptorfile = 'adft.desc')
##   adft <- attach.big.matrix('adft.desc')
  
##   for (i in 1:nrow(adAllRowCol)){
##     print(paste('It is running ', i, ' in total of ', nrow(adAllRowCol), '.', sep = ''))
##     linkData <- c(rowNames[adAllRowCol[i, 1]], colNames[adAllRowCol[i, 2]], adMat[adAllRowCol[i, 1], adAllRowCol[i, 2]])
##     adft[i, ] <- linkData
##   }
    
##   return(adft)
## }

## adj2ftBig2 <- function(adMat, adAllRowCol, n = 4){
  
##   # a big matrix version
##   # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also a bigmatrix.
  
##   library(bigmemory)
##   library(foreach)
##   library(doMC)
##   registerDoMC(n)

##   rowNames <- rownames(adMat)
##   colNames <- colnames(adMat)

##   adft <- foreach (i = 1:nrow(adAllRowCol), .combine = rbind, .inorder=TRUE) %dopar% {
##     print(paste('It is running ', i, ' in total of ', nrow(adAllRowCol), '.', sep = ''))
##     linkData <- c(rowNames[adAllRowCol[i, 1]], colNames[adAllRowCol[i, 2]], adMat[adAllRowCol[i, 1], adAllRowCol[i, 2]])
##     return(linkData)
##   }
    
##   return(adft)
## }

## adj2ftBig3 <- function(adMat, adAllRowCol, n = 2){
  
##   # a big matrix version
##   # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also a bigmatrix.
  
##   library(bigmemory)
##   library(parallel)
##   cl <- makeCluster(n, type = "MPI")
  
##   rowNames <- rownames(adMat)
##   colNames <- colnames(adMat)

##   ## adftBig <- big.matrix(nrow = nrow(adAllRowCol), ncol = 3, backingfile = 'adftBig.bin', descriptorfile = 'adftBig.desc')
##   ## adftBig <- attach.big.matrix('adftBig.desc')
##   ## adftBig[, 1] <- adAllRowCol[, 1]
##   ## adftBig[, 2] <- adAllRowCol[, 2]
    
##   adft <- parRapply(cl = cl, adAllRowCol, function(x) {
##     linkData <- c(rowNames[x[1]], colNames[x[2]], adMat[x[1], x[2]])
##     return(linkData)
##   })
##   stopCluster(cl)

##   return(adft)
## }

## adj2ftBig4 <- function(adMat, adAllRowCol, n = 2){
  
##   # a big matrix version
##   # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also a bigmatrix.
  
##   require(bigmemory)
##   require(parallel)
##   cl <- makeCluster(n, type = "MPI")

##   adMatDescFile <- describe(adMat)
##   adAllRowColDescFile <- describe(adAllRowCol)
  
##   rowNames <- rownames(adMat)
##   colNames <- colnames(adMat)

##   ignore <- clusterEvalQ(cl, {library(bigmemory); NULL})

##   Getft <- function(i, adAllRowColDesc, adMatDesc){
##     ## require(bigmemory)
##     adAllRowColData <- attach.big.matrix(adAllRowColDesc)
##     adMatData <- attach.big.matrix(adMatDesc)
##     rowIndex <- adAllRowColData[i, 1]
##     colIndex <- adAllRowColData[i, 2]
##     linkData <- c(rowNames[rowIndex], rowNames[colIndex], adMatData[rowIndex, colIndex])
##     return(linkData)
##   }

##   adft <- parSapply(cl, 1:nrow(adAllRowCol), Getft, adAllRowColDescFile, adMatDescFile)

##   stopCluster(cl)

##   return(adft)
## }

## # !!!!!!!!!!!!!!!!!!!!! test code!!!!!!!!!!!!!!!!!!!!!
## library(bigmemory)
## library(parallel)

## tmp1 <- attach.big.matrix('tmp1.desc')
## tmp1Grid <- attach.big.matrix('tmp1Grid.desc')
## adj2ftBig4(tmp1, tmp1Grid, n = 2)
## # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#~~~~~~~~~~~~~~~~~~ Pearson's correlation coefficient~~~~~~~~~~~~~~~
load('wholePhyloData.RData')
wholeCor <- cor(t(wholePhyloDataNet))
save(wholeCor, file = 'wholeCor.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ##################### Use the bigmemory and foreach packages #####################
## ## library(bigmemory)
## ## library(foreach)
## ## library(doMC)
## ## registerDoMC(8)
## ##
## ## wholeCor <- attach.big.matrix('wholeCor.desc')
## ## wholeCor[1:10, 1:10]
## ## wholeCorPos <- foreach(i = 1:nrow(wholeCor), .combine = rbind) %:%
## ##   foreach(j  = 1:ncol(wholeCor), .combine = c) %dopar% {
## ##     if (wholeCor[i, j] < 0) {
## ##       return(0)
## ##     } else {
## ##       return(wholeCor[i, j])
## ##     }
## ##   }
## ###############################################################################

######################## Use the parallel package ############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ parallel and openMPI ~~~~~~~~~~~~~~~~~~~~~
load('wholeCor.RData')
require('parallel')
require('KEGGBioCycAPI')
cl <- makeCluster(7, type = "MPI")

# cut mat
cutMat <- CutSeqEqu(nrow(wholeCor), 100)

for (i in 1:ncol(cutMat)) {
  print(paste('It is running ', i, ' in totally of ', ncol(cutMat), '.', sep = ''))
  wholeCorPost <- parApply(cl, wholeCor[cutMat[1, i]:cutMat[2, i], ], 1:2, function(x){if (x < 0) x = 0 else x = x})
  save(wholeCorPost, file = paste('wholeCorPos', i, '.RData', sep = ''))
  rm(wholeCorPost)
  gc()
}

stopCluster(cl)

# combine data
library('foreach')
library('doMC')
registerDoMC(8)

getcontent <- function(s,g) {
    substring(s,g,g+attr(g,'match.length')-1)
}

files <- dir()
fileList <- files[grep('wholeCorPos\\d+.RData', files)]
fileListNumGre <- gregexpr('\\d+', fileList)
fileListNum <- foreach(i = 1:length(fileList), .combine = c) %dopar% {
  getcontent(fileList[i], fileListNumGre[[i]])
}
fileList <- fileList[order(as.numeric(fileListNum))]

wholeCorPostAll <- foreach(i = 1:length(fileList), .combine = rbind) %dopar% {
  load(fileList[i])
  return(wholeCorPost)
}

library('bigmemory')
wholeCorPostAll <- as.big.matrix(wholeCorPostAll, backingfile = 'wholeCorPostAll.bin', descriptorfile = 'wholeCorPostAll.desc')
file.remove(fileList)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################################

#~~~~~~~~~~~~~~~~~~~~ Get Network suffix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('bigmemory')
library('KEGGBioCycAPI')
wholeCorPostAll <- attach.big.matrix('wholeCorPostAll.desc')
## ## # all row-column combination, and transfer it to bigmatrix
## ## allRowCol <- expand.grid(1:nrow(wholeCorPostAll), 1:ncol(wholeCorPostAll))
## ## # allRowCol <- apply(allRowCol, 1:2, as.numeric)
## ## allRowCol <- as.big.matrix(allRowCol, backingfile = 'allRowCol.bin', descriptorfile='allRowCol.desc')
## ## rm(allRowCol)
## ## gc()
## ## allRowCol <- attach.big.matrix('allRowCol.desc')

rowNum <- nrow(wholeCorPostAll)
big.matrix(nrow = rowNum^2, ncol = 2, type = 'integer', backingfile = 'allRowCol.bin', descriptorfile = 'allRowCol.desc')
cutMat <- CutSeqEqu(rowNum^2, rowNum)
cutMat <- as.big.matrix(cutMat, backingfile = 'cutMat.bin', descriptorfile = 'cutMat.desc')
allRowCol <- attach.big.matrix('allRowCol.desc')
gc()
for(i in 1:ncol(cutMat)) {
  print(paste('It is running ', i, ' in tatol of ', ncol(cutMat), '.', sep = ''))
  allRowCol[cutMat[1, i] : cutMat[2, i], 1] <- as.integer(1:rowNum)
  allRowCol[cutMat[1, i] : cutMat[2, i], 2] <- as.integer(rep(cutMat[2, i] %/% rowNum, rowNum))
}
allRowColLogic <- allRowCol[, 1] < allRowCol[, 2]
rm(cutMat, i, rowNum, wholeCorPostAll)
gc()
allRowColLogicData <- allRowCol[allRowColLogic, ]
rm(allRowCol, allRowColLogic)
gc()
allRowColLogicData <- as.big.matrix(allRowColLogicData, backingfile = 'allRowColLogicData.bin', descriptorfile = 'allRowColLogicData.desc')
gc()
## ## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ## wholeCorPostAll <- attach.big.matrix('wholeCorPostAll001.desc')
## ## allRowColLogicData <- attach.big.matrix('allRowColLogicData001.desc')
## ## wholeCorft <- adj2ftBig4(wholeCorPostAll, allRowColLogicData, n = 5)
## ## save(wholeCorft, file = 'wholeCorft001.RData')
## ## #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('bigmemory')
wholeCorPostAll <- attach.big.matrix('wholeCorPostAll.desc')
wholeCorPostAll <- wholeCorPostAll[, ]
upperSeq <- wholeCorPostAll[upper.tri(wholeCorPostAll)]
save(upperSeq, file = 'upperSeq.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('bigmemory')
allRowColLogicData <- attach.big.matrix('allRowColLogicData.desc')
wholeCorPostAll <- attach.big.matrix('wholeCorPostAll.desc')
rowNames <- rownames(wholeCorPostAll)
colNames <- colnames(wholeCorPostAll)
rowNamesUpper <- rowNames[allRowColLogicData[, 1]]
colNamesUpper <- colNames[allRowColLogicData[, 2]]
load('upperSeq.RData')
rm(allRowColLogicData)
rm(wholeCorPostAll)
rm(colNames, rowNames)
gc()
wholeCorft <- cbind(rowNamesUpper, colNamesUpper, upperSeq)
save(wholeCorft, file = 'wholeCorft.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################

####################### Get annotated from-to matrix ######################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ threshold value ~~~~~~~~~~~~~~~~~~~~~~
load('wholeCorft.RData')

cutoff04 <- which(as.numeric(wholeCorft[, 3]) >= 0.4)
cutoff04Mat <- wholeCorft[cutoff04, ]
save(cutoff04Mat, file = 'cutoff04Mat.RData')
rm(wholeCorft, cutoff04)
gc()

cutoff05 <- which(as.numeric(cutoff04Mat[, 3]) >= 0.5)
cutoff05Mat <- cutoff04Mat[cutoff05, ]
save(cutoff05Mat, file = 'cutoff05Mat.RData')
rm(cutoff04Mat, cutoff05)
gc()

cutoff06 <- which(as.numeric(cutoff05Mat[, 3]) >= 0.6)
cutoff06Mat <- cutoff05Mat[cutoff06, ]
save(cutoff06Mat, file = 'cutoff06Mat.RData')
rm(cutoff05Mat, cutoff06)
gc()

cutoff07 <- which(as.numeric(cutoff06Mat[, 3]) >= 0.7)
cutoff07Mat <- cutoff06Mat[cutoff07, ]
save(cutoff07Mat, file = 'cutoff07Mat.RData')
rm(cutoff06Mat, cutoff07)
gc()

cutoff08 <- which(as.numeric(cutoff07Mat[, 3]) >= 0.8)
cutoff08Mat <- cutoff07Mat[cutoff08, ]
save(cutoff08Mat, file = 'cutoff08Mat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~ Annotation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('cutoff08Mat.RData')
load('wholePhyloData.RData')
fromNames <- annoFirst[match(cutoff08Mat[, 1], names(annoFirst))]
toNames <- annoFirst[match(cutoff08Mat[, 2], names(annoFirst))]
cutoff08MatAnno <- cbind(fromNames, toNames, cutoff08Mat[, 3])
rownames(cutoff08MatAnno) <- NULL
save(cutoff08MatAnno, file = 'cutoff08MatAnno.RData')
## write.csv(cutoff08MatAnno, file = 'cutoff08MatAnno.csv')

load('cutoff07Mat.RData')
load('wholePhyloData.RData')
fromNames <- annoFirst[match(cutoff07Mat[, 1], names(annoFirst))]
toNames <- annoFirst[match(cutoff07Mat[, 2], names(annoFirst))]
cutoff07MatAnno <- cbind(fromNames, toNames, cutoff07Mat[, 3])
rownames(cutoff07MatAnno) <- NULL
save(cutoff07MatAnno, file = 'cutoff07MatAnno.RData')
## write.csv(cutoff07MatAnno, file = 'cutoff07MatAnno.csv')

load('cutoff06Mat.RData')
load('wholePhyloData.RData')
fromNames <- annoFirst[match(cutoff06Mat[, 1], names(annoFirst))]
toNames <- annoFirst[match(cutoff06Mat[, 2], names(annoFirst))]
cutoff06MatAnno <- cbind(fromNames, toNames, cutoff06Mat[, 3])
rownames(cutoff06MatAnno) <- NULL
save(cutoff06MatAnno, file = 'cutoff06MatAnno.RData')
## write.csv(cutoff06MatAnno, file = 'cutoff06MatAnno.csv')

load('cutoff05Mat.RData')
load('wholePhyloData.RData')
fromNames <- annoFirst[match(cutoff05Mat[, 1], names(annoFirst))]
toNames <- annoFirst[match(cutoff05Mat[, 2], names(annoFirst))]
cutoff05MatAnno <- cbind(fromNames, toNames, cutoff05Mat[, 3])
rownames(cutoff05MatAnno) <- NULL
save(cutoff05MatAnno, file = 'cutoff05MatAnno.RData')
## write.csv(cutoff05MatAnno, file = 'cutoff05MatAnno.csv')

load('cutoff04Mat.RData')
load('wholePhyloData.RData')
fromNames <- annoFirst[match(cutoff04Mat[, 1], names(annoFirst))]
toNames <- annoFirst[match(cutoff04Mat[, 2], names(annoFirst))]
cutoff04MatAnno <- cbind(fromNames, toNames, cutoff04Mat[, 3])
rownames(cutoff04MatAnno) <- NULL
save(cutoff04MatAnno, file = 'cutoff04MatAnno.RData')
## write.csv(cutoff04MatAnno, file = 'cutoff04MatAnno.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~ ATP synthase interaction~~~~~~~~~~~~~~~~~~~~~
load('cutoff04MatAnno.RData')
load('mitoPhyloData.RData')
fromList <- cutoff04MatAnno[, 1] %in% ATPgeneList[, 1]
toList <- cutoff04MatAnno[, 2] %in% ATPgeneList[, 1]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff04ATP <- cutoff04MatAnno[fullList, ]
save(cutoff04ATP, file = 'cutoff04ATPAnno.RData')
## write.csv(cutoff04ATP, file = 'cutoff04ATP.csv')

load('cutoff05MatAnno.RData')
load('mitoPhyloData.RData')
fromList <- cutoff05MatAnno[, 1] %in% ATPgeneList[, 1]
toList <- cutoff05MatAnno[, 2] %in% ATPgeneList[, 1]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff05ATP <- cutoff05MatAnno[fullList, ]
save(cutoff05ATP, file = 'cutoff05ATPAnno.RData')
## write.csv(cutoff05ATP, file = 'cutoff05ATP.csv')

load('cutoff06MatAnno.RData')
load('mitoPhyloData.RData')
fromList <- cutoff06MatAnno[, 1] %in% ATPgeneList[, 1]
toList <- cutoff06MatAnno[, 2] %in% ATPgeneList[, 1]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff06ATP <- cutoff06MatAnno[fullList, ]
save(cutoff06ATP, file = 'cutoff06ATPAnno.RData')
## write.csv(cutoff06ATP, file = 'cutoff06ATP001.csv')

load('cutoff07MatAnno.RData')
load('mitoPhyloData.RData')
fromList <- cutoff07MatAnno[, 1] %in% ATPgeneList[, 1]
toList <- cutoff07MatAnno[, 2] %in% ATPgeneList[, 1]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff07ATP <- cutoff07MatAnno[fullList, ]
save(cutoff07ATP, file = 'cutoff07ATPAnno.RData')
## write.csv(cutoff07ATP, file = 'cutoff07ATP.csv')

load('cutoff08MatAnno.RData')
load('mitoPhyloData.RData')
fromList <- cutoff08MatAnno[, 1] %in% ATPgeneList[, 1]
toList <- cutoff08MatAnno[, 2] %in% ATPgeneList[, 1]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff08ATP <- cutoff08MatAnno[fullList, ]
save(cutoff08ATP, file = 'cutoff08ATPAnno.RData')
## write.csv(cutoff08ATP, file = 'cutoff08ATP.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ATP with Entrez names~~~~~~~~~~~~~~~~~~~~~
load('cutoff08Mat.RData')
load('mitoPhyloData.RData')
fromList <- cutoff08Mat[, 1] %in% ATPgeneList[, 2]
toList <- cutoff08Mat[, 2] %in% ATPgeneList[, 2]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff08ATP <- cutoff08Mat[fullList, ]
save(cutoff08ATP, file = 'cutoff08ATP.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~ One protein to multiple ATP synthase subunits~~~~~~~~~~~~
library(foreach)
library(doMC)
registerDoMC(4)
load(file = 'cutoff08ATP.RData')
load('mitoPhyloData.RData')
MSVal <- as.character(unlist(read.table('MSVal.txt')))

cutoff08ATPUpdata <- foreach(i = 1:nrow(cutoff08ATP), .combine = rbind) %dopar% {
  
  if (cutoff08ATP[i, 1] %in% ATPgeneList[, 1]) {
    cutoff08ATP[i, 1:2] <- cutoff08ATP[i, 2:1]
  } else {}

  return(cutoff08ATP[i, ])
}

cutoff08Cprs <- aggregate(cutoff08ATPUpdata[, 2:3],
                  by = list(cutoff08ATPUpdata[, 1]),
                  paste, collapse = '|')
MShas <- rep(0, nrow(cutoff08Cprs))
MShas[match(MSVal, cutoff08Cprs[, 1])] <- 1
cutoff08Cprs <- cbind(cutoff08Cprs, MShas)
colnames(cutoff08Cprs) <- c('Symbol', 'Inter', 'Cor', 'MS')
write.csv(cutoff08Cprs, file = 'cutoff08Cprs.csv')


library('foreach')
library('doMC')
registerDoMC(4)
cutoff07ATP <- read.csv('cutoff07ATP.csv', row.names = 1)
load('mitoPhyloData.RData')
MSVal <- as.character(unlist(read.table('MSVal.txt')))

cutoff07ATPUpdata <- foreach(i = 1:nrow(cutoff07ATP), .combine = rbind) %dopar% {
  
  if (cutoff07ATP[i, 1] %in% ATPgeneList[, 1]) {
    cutoff07ATP[i, 1:2] <- cutoff07ATP[i, 2:1]
  } else {}

  return(cutoff07ATP[i, ])
}

cutoff07Cprs <- aggregate(cutoff07ATPUpdata[, 2:3],
                  by = list(cutoff07ATPUpdata[, 1]),
                  paste, collapse = '|')
MShas <- rep(0, nrow(cutoff07Cprs))
MShas[match(MSVal, cutoff07Cprs[, 1])] <- 1
cutoff07Cprs <- cbind(cutoff07Cprs, MShas)
colnames(cutoff07Cprs) <- c('Symbol', 'Inter', 'Cor', 'MS')
write.csv(cutoff07Cprs, file = 'cutoff07Cprs.csv')


library(foreach)
library(doMC)
registerDoMC(4)
load(file = 'cutoff06ATP.RData')
load('mitoPhyloData.RData')
MSVal <- as.character(unlist(read.table('MSVal.txt')))

cutoff06ATPUpdata <- foreach(i = 1:nrow(cutoff06ATP), .combine = rbind) %dopar% {
  
  if (cutoff06ATP[i, 1] %in% ATPgeneList[, 1]) {
    cutoff06ATP[i, 1:2] <- cutoff06ATP[i, 2:1]
  } else {}

  return(cutoff06ATP[i, ])
}

cutoff06Cprs <- aggregate(cutoff06ATPUpdata[, 2:3],
                  by = list(cutoff06ATPUpdata[, 1]),
                  paste, collapse = '|')
MShas <- rep(0, nrow(cutoff06Cprs))
MShas[match(MSVal, cutoff06Cprs[, 1])] <- 1
cutoff06Cprs <- cbind(cutoff06Cprs, MShas)
colnames(cutoff06Cprs) <- c('Symbol', 'Inter', 'Cor', 'MS')
write.csv(cutoff06Cprs, file = 'cutoff06Cprs.csv')

library(foreach)
library(doMC)
registerDoMC(4)
load(file = 'cutoff05ATP.RData')
load('mitoPhyloData.RData')
MSVal <- as.character(unlist(read.table('MSVal.txt')))

cutoff05ATPUpdata <- foreach(i = 1:nrow(cutoff05ATP), .combine = rbind) %dopar% {
  
  if (cutoff05ATP[i, 1] %in% ATPgeneList[, 1]) {
    cutoff05ATP[i, 1:2] <- cutoff05ATP[i, 2:1]
  } else {}

  return(cutoff05ATP[i, ])
}

cutoff05Cprs <- aggregate(cutoff05ATPUpdata[, 2:3],
                  by = list(cutoff05ATPUpdata[, 1]),
                  paste, collapse = '|')
MShas <- rep(0, nrow(cutoff05Cprs))
MShas[match(MSVal, cutoff05Cprs[, 1])] <- 1
cutoff05Cprs <- cbind(cutoff05Cprs, MShas)
colnames(cutoff05Cprs) <- c('Symbol', 'Inter', 'Cor', 'MS')
write.csv(cutoff05Cprs, file = 'cutoff05Cprs.csv')

library(foreach)
library(doMC)
registerDoMC(4)
load(file = 'cutoff04ATP.RData')
load('mitoPhyloData.RData')
MSVal <- as.character(unlist(read.table('MSVal.txt')))

cutoff04ATPUpdata <- foreach(i = 1:nrow(cutoff04ATP), .combine = rbind) %dopar% {
  
  if (cutoff04ATP[i, 1] %in% ATPgeneList[, 1]) {
    cutoff04ATP[i, 1:2] <- cutoff04ATP[i, 2:1]
  } else {}

  return(cutoff04ATP[i, ])
}

cutoff04Cprs <- aggregate(cutoff04ATPUpdata[, 2:3],
                  by = list(cutoff04ATPUpdata[, 1]),
                  paste, collapse = '|')
MShas <- rep(0, nrow(cutoff04Cprs))
MShas[match(MSVal, cutoff04Cprs[, 1])] <- 1
cutoff04Cprs <- cbind(cutoff04Cprs, MShas)
colnames(cutoff04Cprs) <- c('Symbol', 'Inter', 'Cor', 'MS')
write.csv(cutoff04Cprs, file = 'cutoff04Cprs.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##############################################################################

#############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Output cytoscape data type~~~~~~~~~~~~~~~~~~~
cutoff07ATP001 <- read.csv('cutoff07ATP001.csv', row.names = 1)
fromHas <- sapply(gregexpr(' ', as.character(cutoff07ATP001[, 1])), '[[', 1)
fromHas <- fromHas != -1
toHas <- sapply(gregexpr(' ', as.character(cutoff07ATP001[, 2])), '[[', 1)
toHas <- toHas != -1
wholeHas <- fromHas | toHas
cutoff07ATP001 <- cutoff07ATP001[!wholeHas, ]
write.table(cutoff07ATP001, 'cutoff07ATP001.txt', col.names = FALSE, row.names = FALSE, quote = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################################################################
#~~~~~~~~~~~~~~~~~~~~Test correlation coefficient~~~~~~~~~~~~~~~~~~~~~~~~~
load('corValues.RData')
load('cutoff06ATP.RData')
library(foreach)
library(doMC)
registerDoMC(8)

# set the bootstrap number
bootNum <- 10e6
corP <- foreach(i = 1:nrow(cutoff06ATP), .combine = c) %dopar% {
  print(paste('It is running ', i, ' in total of ', nrow(cutoff06ATP), '.', sep = ''))
  smpCor <- sample(corValues, bootNum, replace = TRUE)
  pCor <- sum(smpCor > as.numeric(cutoff06ATP[i, 3])) / bootNum

  return(pCor)
}
write.csv(cbind(cutoff06ATP, corP), 'cutoff06ATPcorP.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~Jaccard Similarity Coefficient~~~~~~~~~~~~~~~~~~~
library(foreach)
library(doMC)
registerDoMC(4)
load('jaccardSim.RData')
load('wholePhyloData.RData')
cutoff06ATPcorP <- read.csv('cutoff06ATPcorP.csv', row.names = 1)
colnames(cutoff06ATPcorP)[3] <- 'cor'

corSim <- foreach(i = 1:nrow(cutoff06ATPcorP), .combine = rbind) %dopar% {
  print(paste('It is running ', i, ' in total of ', nrow(cutoff06ATPcorP), sep = ''))
  corInputIndex1 <- which(annoFirst %in% cutoff06ATPcorP[i, 1])
  corInputIndex2 <- which(annoFirst %in% cutoff06ATPcorP[i, 2])

  corValue <- cor(t(jaccardSim[c(corInputIndex1, corInputIndex2), ]))
  corValue <- corValue[1, 2]

  jacSim <- jaccardSim[corInputIndex1, corInputIndex2]

  return(c(jacSim, corValue))

}

colnames(corSim) <- c('JaccardSim', 'SimCor')

cutoff06ATPcorSim <- cbind(cutoff06ATPcorP, corSim)
write.csv(cutoff06ATPcorSim, 'cutoff06ATPcorSim.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~Cluster whole phylogenetic profiles~~~~~~~~~~~~~~~~~~
load('wholePhyloData.RData')
library(compiler)
enableJIT(3)

print('It is running Manhattan distance.')
sampleDist <- dist(wholePhyloDataNet, method = 'manhattan')
save(sampleDist, file = 'manhattanDist.RData')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~cor/bootstrap/Jac/JacCor of cluster~~~~~~~~~~~~~~~
load('manhattanDist.RData')
load('mitoPhyloData.RData')
load('wholePhyloData.RData')

cutoff06ATPcorSim <- read.csv('cutoff06ATPcorSim.csv', row.names = 1)

sampleTree <- hclust(sampleDist, method = 'average')
memb <- cutree(sampleTree, k = 10)
## memb[names(memb) %in% ATPgeneList[, 2]]

# sort memb
memb <- memb[order(names(memb))]
memb <- memb[rank(names(annoFirst))]

# choose F1 subcomplex (Euclidean method includes alpha/beta/gamma/OSCP and a subunits)
F1ClNames <- annoFirst[memb == 1]
F1fromHas <- cutoff06ATPcorSim[, 'fromNames'] %in% F1ClNames
F1toHas <- cutoff06ATPcorSim[, 'toNames'] %in% F1ClNames
cutoff06ATPcorSimEu <- cutoff06ATPcorSim[F1fromHas | F1toHas, ]
write.csv(cutoff06ATPcorSimEu, 'cutoff06ATPcorSimEug1.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#~~~~~~~~~~~~~~~~~~~~~~~~~ circos ~~~~~~~~~~~~~~~~~~~~~~~~
library(biomaRt)
load('mitoPhyloData.RData')
load('wholePhyloData.RData')
cutoff06ATPcorSimEu <- read.csv('cutoff06ATPcorSimEug3.csv', row.names = 1)
## tmp1 <- read.csv('cutoff06ATPcorSimEug1.csv', row.names = 1)
## tmp2 <- read.csv('cutoff06ATPcorSimEug2.csv', row.names = 1)
## cutoff06ATPcorSimEu <- rbind(cutoff06ATPcorSimEu, tmp1)
## cutoff06ATPcorSimEu <- rbind(cutoff06ATPcorSimEu, tmp2)

getProID <- function(KEGGspec){
  # USE: get the whole KEGG IDs of certain species.
  # INPUT: 'KEGGspec' is the KEGG species, for example, 'hsa'.
  # OUTPU: matrix of KEGG ID

  require('RCurl')

  # get KEGG ID annotation list
  linkKEGGID <- paste('http://rest.kegg.jp/list/', KEGGspec, sep = '')
  webIDList <-getURL(linkKEGGID)

  # transfer webpage into a matrix
  IDAnno <- unlist(strsplit(webIDList, split = '\n', fixed = TRUE))
  IDAnno <- sapply(IDAnno, strsplit, split = '\t', fixed = TRUE)
  IDAnno <- matrix(unlist(IDAnno), ncol = 2, byrow = TRUE)

  return(IDAnno)
}

ft2circos <- function(ft, locaAnno, nodeName = 'all', showEdge = TRUE, thick = 5, showHeatmap = TRUE){
  # USE: convert the 'from-to matrix' to circos link file
  # INPUT: 'ft' if from-to matrix, here we used the undirected network; 'locaAnno' is the gene location annotation information; 'nodeName' is the node name list, users may change; 'showEdge' whether to show the thickness, which is corresponding to edge weight; 'thick' parameter that transfer weight to thickness.
  # ft file
##          from   to    weight
## 423837   UQCRQ ZCD1 0.2466398
## 423838   USMG5 ZCD1 0.2488576
## 423840   VAMP1 ZCD1 0.3520448
## 423841   VAMP8 ZCD1 0.3520448
## 423846 WBSCR16 ZCD1 0.2861196
## 423848   WDR22 ZCD1 0.2861196
  # locaAnno
##     SYM chromosome_name start_position end_position
## 1  ABAT              16        8768422      8878432
## 2 ACAA1               3       38144620     38178733
## 3 ACACA              17       35441923     35766909
## 4 ACACB              12      109554400    109706031
## 5 ACADL               2      211052663    211090215
## 6 ACADS              12      121163538    121177811
  # OUTPUT: list
  # first element is circos links file
## hs2 71175494 71176689 hs3 107187735 107188930 thickness=5p
## hs2 85623609 85625914 hs3 113896945 113899242 thickness=5p
  # second element is circos label file
## hs4 186064395 186068434 SLC25A4
## hsX 118602363 118605282 SLC25A5
## hs3 49396578 49450431 RHOA

  if (!identical(nodeName, 'all')){
    ftIn <- ft[ft[, 1] %in% nodeName, , drop = FALSE]
    ftOut <- ft[ft[, 2] %in% nodeName, , drop = FALSE]
    ftOut[, 1:2] <- ftOut[, 2:1]
    ft <- rbind(ftIn, ftOut)

  } else {}

  locaAnno[, 2] <- paste('hs', locaAnno[, 2], sep = '')

  nodeIn <- merge(ft, locaAnno, by.x = colnames(ft)[1], by.y = colnames(locaAnno)[1], sort = FALSE)

  # since some of nodes in 'ft' may not have local annotation, so we just use the out nodes in 'nodeIn'

  nodeOut <- merge(nodeIn[, 1:3], locaAnno, by.x = colnames(nodeIn)[2], by.y = colnames(locaAnno)[1], sort = FALSE)

  circosLink <- merge(nodeIn, nodeOut, by.x = colnames(nodeIn)[2], by.y = colnames(nodeOut)[1], sort = FALSE)

  # get circos label
  circosLabelName <- unique(c(as.character(circosLink[, 1]), as.character(circosLink[, 2])))
  circosLabel <- locaAnno[locaAnno[, 1] %in% circosLabelName, ]
  circosLabel <- circosLabel[, c(2:4, 1)]

  # get circos columns
  weight = as.numeric(as.character(circosLink[, 3]))
  circosLink <- circosLink[, c(4:6, 9:11)]

  # show edge
  # in case of no nodes and edges
  if (dim(circosLink)[1] != 0){
    if (showEdge){
      weight <- weight * thick
      weight <- paste('thickness=', weight, 'p', sep = '')
      circosLink <- cbind(circosLink, weight)
    } else {
      circosLink <- cbind(circosLink, paste('thickness=', thick, 'p', sep = ''))
    }
  } else {}

  # show heatmap
  if (dim(circosLink)[1] != 0){
    if (showHeatmap){
      heatVal <- cbind(circosLink[, 4:6], weight)
      heatVal <- apply(heatVal, 1:2, as.character)
      heatVal <- rbind(heatVal, c(as.character(circosLink[1, 1:3]), 1))
    } else {}
  } else {
    heatVal <- circosLink[, 1:4, drop = FALSE]
  }

  return(list(circosLink = circosLink, circosLabel = circosLabel, circosHeatmap = heatVal))

}


# annotation file
hsaWholeEntrez <- getProID('hsa')
hsaEntrez <- sapply(strsplit(hsaWholeEntrez[, 1], split = ':', fixed = TRUE), function(x) x[2])
hsaSymbol <- sapply(strsplit(hsaWholeEntrez[, 2], split = ',', fixed = TRUE), function(x) x[1])
hsaSymbol <- sapply(strsplit(hsaSymbol, split = ';', fixed = TRUE), function(x) x[1])
hsaEntrezSymbol <- cbind(hsaWholeEntrez, hsaEntrez, hsaSymbol)
colnames(hsaEntrezSymbol)[1:2] <- c('KEGGID', 'KEGGAnno')
hsaEntrezSymbol <- as.data.frame(hsaEntrezSymbol)

# biomaRt annotation to get gene position 'from' and 'to'
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
hsaAnno <- getBM(attributes = c('entrezgene', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'entrezgene', value = hsaEntrez, mart = ensembl)
# remove ones without chromesome annotation
chromName <- c(1:22, 'X', 'Y', 'MT')
hsaAnno <- hsaAnno[hsaAnno[, 3] %in% chromName, ]
hsaAnno <- aggregate(hsaAnno, by = list(hsaAnno[, 1]), function(x) x[1])
hsaAnno <- hsaAnno[, -1]

wholeAnno <- merge(hsaAnno, hsaEntrezSymbol, by.x = 'entrezgene', by.y = 'hsaEntrez')
wholeAnno[, 'hsaSymbol'] <- wholeAnno[, 'hgnc_symbol']
wholeAnno[which(wholeAnno[, 'hsaSymbol'] == ''), 'hsaSymbol'] <- paste('gene', 1:2207, sep = '')
wholeAnno[wholeAnno[, 'hsaSymbol'] == 'MT-ATP6', 'hsaSymbol'] <- 'ATP6'
wholeAnno[wholeAnno[, 'hsaSymbol'] == 'MT-ATP8', 'hsaSymbol'] <- 'ATP8'


# transfer the from-to file
fromNodes <- cutoff06ATPcorSimEu[, 'fromNames']
toNodes <- cutoff06ATPcorSimEu[, 'toNames']
fromNodesKEGG <- names(annoFirst[match(fromNodes, annoFirst)])
toNodesKEGG <- names(annoFirst[match(toNodes, annoFirst)])
fromNodesSymbol <- as.character(wholeAnno[match(fromNodesKEGG, wholeAnno[, 'KEGGID']), 'hsaSymbol'])
toNodesSymbol <- as.character(wholeAnno[match(toNodesKEGG, wholeAnno[, 'KEGGID']), 'hsaSymbol'])
cutoff06ATPcorSimEu[, 1:2] <- cbind(fromNodesSymbol, toNodesSymbol)
hasNA <- is.na(fromNodesSymbol) | is.na(toNodesSymbol)
cutoff06ATPcorSimEu <- cutoff06ATPcorSimEu[!hasNA, ]

labelGene <- ft2circos(cutoff06ATPcorSimEu[, 1:3], wholeAnno[, c(8, 3:5)], ATPgeneList[1:6, 1])[[2]]
write.table(labelGene, 'labelGene.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

for (i in 1:17){
  corPhyloCir <- ft2circos(cutoff06ATPcorSimEu[, 1:3], wholeAnno[, c(8, 3:5)], ATPgeneList[i, 1], showEdge = FALSE, thick = 3)
  write.table(corPhyloCir[[1]], paste(ATPgeneList[i, 1], 'txt', sep = '.'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(corPhyloCir[[3]], paste(ATPgeneList[i, 1], 'heatmap', '.txt', sep = ''), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# genome frequency
phyloSpe <- read.csv('wholeListFile.csv', row.names = 1)
phyloCode <- as.character(phyloSpe[, 4])
phyloCode <- sapply(strsplit(phyloCode, split = ';', fixed = TRUE), function(x) x[1:2])
phyloCode <- t(phyloCode)
phyloCode <- cbind(as.character(phyloSpe[, 2]), phyloCode)

arcSpe <- phyloCode[phyloCode[, 3] %in% 'Archaea', 1]
bacSpe <- phyloCode[phyloCode[, 3] %in% 'Bacteria', 1]
euSpe <- phyloCode[!(phyloCode[, 3] %in% c('Archaea', 'Bacteria')), 1]

arc <- wholePhyloDataNet[, colnames(wholePhyloDataNet) %in% arcSpe]
bac <- wholePhyloDataNet[, colnames(wholePhyloDataNet) %in% bacSpe]
eu <- wholePhyloDataNet[, colnames(wholePhyloDataNet) %in% euSpe]

arcFreq <- apply(arc, 1, function(x) {sum(x)/length(x)})
bacFreq <- apply(bac, 1, function(x) {sum(x)/length(x)})
euFreq <- apply(eu, 1, function(x) {sum(x)/length(x)})

mitoVecFreq <- cbind(arcFreq, bacFreq, euFreq)
mitoVecFeqRN <- cbind(as.character(wholeAnno[match(rownames(mitoVecFreq), wholeAnno[, 'KEGGID']), 'hsaSymbol']), mitoVecFreq)
colnames(mitoVecFeqRN)[1] <- 'symbol'

freqGene <- merge(labelGene, mitoVecFeqRN, by.x = 'hsaSymbol', by.y = 'symbol', sort = FALSE)
freqGene <- freqGene[, -1]

write.table(freqGene[, c(1:3, 4)], 'freqArc.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(freqGene[, c(1:3, 5)], 'freqBac.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(freqGene[, c(1:3, 6)], 'freqEu.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


######################  similarity method #################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(compiler)
library(foreach)
library(doMC)
registerDoMC(4)
load('wholePhyloData.RData')
cutoff06ATPcorSim <- read.csv('cutoff06ATPcorSim.csv', row.names = 1)

YuleqDiff <- function(seq1, seq2) {
  # USE: the yuleq difference 2bc/(ad+bc)
  # INPUT: seq1 and seq2 are 0-1 data with the same length
  # OUTPUT: yuleq difference

  sumSeq <- seq1 + seq2
  difSeq <- seq1 - seq2
  
  a <- sum(sumSeq == 2)

  b <- sum(difSeq == 1)

  c <- sum(difSeq == -1)

  d <- sum(sumSeq == 0)

  diff <- 2*b*c/(a*d + b*c)

  return(diff)
}

YuleqDiffMat1 <- function(matInput, n = 4) {

  # USE: generate the yuleq difference by a matrix
  # INPUT: 'matInput' is the input matrix. 'n' is the threads number
  # OUTPUT: a yuleq matrix

  library(foreach)
  library(doMC)
  registerDoMC(n)
  
  rowNum <- nrow(matInput)
  ## yuleqMat <- matrix(ncol = rowNum, nrow = rowNum)

  yuleqMat <- foreach(i = 1:rowNum, .combine = cbind) %:%
    foreach(j = 1:rowNum, .combine = c) %dopar% {
      YuleqDiff(matInput[i, ], matInput[j, ])
    }

  rownames(yuleqMat) <- rownames(matInput)
  colnames(yuleqMat) <- rownames(matInput)

  return(yuleqMat)
  
}

YuleqDiffMat2 <- function(matInput, n = 3) {

  library('Rmpi')
  library('parallel')
  library('bigmemory')
  cl <- makeCluster(n, type = "MPI")

  YuleqDiff <- function(seq1, seq2) {
  # USE: the yuleq difference 2bc/(ad+bc)
  # INPUT: seq1 and seq2 are 0-1 data with the same length
  # OUTPUT: yuleq difference

    sumSeq <- seq1 + seq2
    difSeq <- seq1 - seq2
    
    a <- sum(sumSeq == 2)

    b <- sum(difSeq == 1)

    c <- sum(difSeq == -1)

    d <- sum(sumSeq == 0)

    diff <- 2*b*c/(a*d + b*c)

    return(diff)
  }

  rowNum <- nrow(matInput)
  big.matrix(ncol = rowNum, nrow = rowNum, backingfile = 'yuleqMat.bin', descriptorfile = 'yuleqMat.desc', dimnames = list(rownames(matInput), rownames(matInput)))
  yuleqMat <- attach.big.matrix('yuleqMat.desc')
  
  for (j in 1:rowNum) {

    print(paste('It is running ', j, ' in a total of ', rowNum, sep = ''))
    diffEachRow <- parRapply(cl, matInput, function(i) {
      diffValue <- YuleqDiff(matInput[j, ], i)
      return(diffValue)
    })

    yuleqMat[j, ] <- diffEachRow
  }

  stopCluster(cl)
}

YuleqDiffMatCom2 <- cmpfun(YuleqDiffMat2)
YuleqDiffMatCom2(wholePhyloDataNet, n = 7)


## diffYuleq <- foreach(i = 1:nrow(cutoff06ATPcorSim), .combine = c) %dopar% {
##   print(paste('It is running ', i, ' in total of ', nrow(cutoff06ATPcorSim), sep = ''))
##   corInputIndex1 <- which(annoFirst %in% cutoff06ATPcorSim[i, 1])
##   corInputIndex2 <- which(annoFirst %in% cutoff06ATPcorSim[i, 2])

##   diffYuleqEach <- YuleqDiff(wholePhyloDataNet[corInputIndex1, ], wholePhyloDataNet[corInputIndex2, ])
  
##   return(diffYuleqEach)
## }

## cutoff06ATPYuleq <- cbind(cutoff06ATPcorSim, diffYuleq)
## colnames(cutoff06ATPYuleq)[7] <- 'YuleqDiff'
## write.csv(cutoff06ATPYuleq, 'cutoff06ATPYuleq.csv')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################################################################


#######################explore phyloData########################################
load('wholePhyloData.RData')
phyloSpe <- read.csv('wholeListFile.csv', row.names = 1)
cutoff07Cprs <- read.csv('cutoff07Cprs.csv', row.names = 1)

# genome frequency
phyloCode <- as.character(phyloSpe[, 4])
phyloCode <- sapply(strsplit(phyloCode, split = ';', fixed = TRUE), function(x) x[1:2])
phyloCode <- t(phyloCode)
phyloCode <- cbind(as.character(phyloSpe[, 2]), phyloCode)

arcSpe <- phyloCode[phyloCode[, 3] %in% 'Archaea', 1]
bacSpe <- phyloCode[phyloCode[, 3] %in% 'Bacteria', 1]
euSpe <- phyloCode[!(phyloCode[, 3] %in% c('Archaea', 'Bacteria')), 1]

## # explore the interaction proteins
## as.matrix(sort(table(cutoff07Cprs[, 2]), decreasing = TRUE))

PseudoSpeFreq <- function(geneName = 'ATP5E') {
  ATP5EInter <- as.character(cutoff07Cprs[cutoff07Cprs[, 2] == geneName, 1])
  ATP5EInter <- names(annoFirst[annoFirst %in% ATP5EInter])
  ATP5EInterPhyl <- wholePhyloDataNet[rownames(wholePhyloDataNet) %in% ATP5EInter, ]
  speNames <- as.character(phyloSpe[-1, 2])
  ATP5EInterPhyl <- ATP5EInterPhyl[ ,order(colnames(ATP5EInterPhyl))]
  ATP5EInterPhyl <- ATP5EInterPhyl[, rank(speNames)]

  arc <- ATP5EInterPhyl[, colnames(ATP5EInterPhyl) %in% arcSpe]
  bac <- ATP5EInterPhyl[, colnames(ATP5EInterPhyl) %in% bacSpe]
  eu <- ATP5EInterPhyl[, colnames(ATP5EInterPhyl) %in% euSpe]

  arcFreq <- apply(arc, 1, function(x) {sum(x)/length(x)})
  bacFreq <- apply(bac, 1, function(x) {sum(x)/length(x)})
  euFreq <- apply(eu, 1, function(x) {sum(x)/length(x)})

  mitoVecFreq <- cbind(euFreq, bacFreq, arcFreq)
  matplot(mitoVecFreq, type = 'l')
}



geneNameMat <- as.matrix(sort(table(cutoff07Cprs[, 2]), decreasing = TRUE))
geneNameList <- names(geneNameMat[1:12, 1])
par(mfrow=c(3,4))
for(i in geneNameList) {
  PseudoSpeFreq(i)
}
################################################################################

