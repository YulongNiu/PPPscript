###########################select linkages######################
library('stringr')

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/')

## whole 6422 linkages
MMMLinksWhole <- read.csv('Sup1_MMM12+.csv', stringsAsFactor = FALSE)
MMMLinksWhole <- MMMLinksWhole[, c(3:7)]

## whole annotation
MMMAnno <- read.csv('Sup2_OMAprotein_attributes.csv', stringsAsFactor = FALSE)
MMMAnno <- MMMAnno[, c(1, 3)]

## select linkages Anno
linkedGenes <- unique(c(MMMLinksWhole[, 1], MMMLinksWhole[, 2]))
MMMAnno <- MMMAnno[MMMAnno[, 1] %in% linkedGenes, ]
rownames(MMMAnno) <- NULL

## HNGC
HGNCAnno <- read.csv('HGNC_anno.txt', sep = '\t', stringsAsFactor = FALSE)
cHGNC <- HGNCAnno[, 2]
pHGNC <- HGNCAnno[, 5]
mAnno <- MMMAnno[, 2]
geneIdx <- numeric(nrow(MMMAnno))
for(i in seq_along(geneIdx)) {
  eachGene <- mAnno[i]
  eachGene <- sub('ORF', 'orf', eachGene)
  currentIdx <- match(eachGene, cHGNC)

  if (is.na(currentIdx)) {

    for (j in seq_along(pHGNC)) {
      eachP <- unlist(strsplit(pHGNC[j], split = ',', fixed = TRUE))
      eachP <- str_trim(eachP)
      if (sum(eachP %in% eachGene) > 0) {
        currentIdx <- j
        break
      } else {}
    }

  } else {}
  geneIdx[i] <- currentIdx
}

MMMAnnoEntrez <- cbind(MMMAnno, entrez = HGNCAnno[geneIdx, 11])
which(is.na(MMMAnnoEntrez[, 3]))

## manually check
MMMAnnoEntrez[70, 3] <- 23189
MMMAnnoEntrez[189, 3] <- 26354
MMMAnnoEntrez[264, 3] <- 6133
MMMAnnoEntrez[284, 3] <- 5783
MMMAnnoEntrez[289, 3] <- 84247
MMMAnnoEntrez[329, 3] <- NA
MMMAnnoEntrez[360, 3] <- 2966
MMMAnnoEntrez[489, 3] <- 4233
MMMAnnoEntrez[503, 3] <- 285966
MMMAnnoEntrez[528, 3] <- 5395
MMMAnnoEntrez[533, 3] <- NA
MMMAnnoEntrez[569, 3] <- NA
MMMAnnoEntrez[706, 3] <- 9738
MMMAnnoEntrez[748, 3] <- 54768
MMMAnnoEntrez[796, 3] <- 6143
MMMAnnoEntrez[920, 3] <- 56342
MMMAnnoEntrez[992, 3] <- 6122
MMMAnnoEntrez[1034, 3] <- 22925
MMMAnnoEntrez[1059, 3] <- 3329
MMMAnnoEntrez[1136, 3] <- 64395
MMMAnnoEntrez[1330, 3] <- 8078
MMMAnnoEntrez[1331, 3] <- 7167
MMMAnnoEntrez[1428, 3] <- 1650
MMMAnnoEntrez[1520, 3] <- 6202

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## update
noGene <- MMMAnnoEntrez[is.na(MMMAnnoEntrez[, 3]), 1]
## from 1608 --> 1605
MMMAnnoEntrez <- MMMAnnoEntrez[!is.na(MMMAnnoEntrez[, 3]), ]
noLinkLog <- (MMMLinksWhole[, 1] %in% noGene) | (MMMLinksWhole[, 2] %in% noGene)
## from 6422 --> 6386
MMMLinksWhole <- MMMLinksWhole[!noLinkLog, ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(MMMLinksWhole, MMMAnnoEntrez, file = 'MMM12.RData')
#################################################################

################select all present#####################
library('doParallel')
library('foreach')

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/wholePhyloData.RData')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/MMM12.RData')

## set threshold
thres <- 0
phyloP <- wholePhyloDataNet[rowSums(wholePhyloDataNet) >= floor(ncol(wholePhyloDataNet) * thres), , drop = FALSE]

## check MMM12
MMMAnnoEntrez[, 3] <- paste0('hsa:', MMMAnnoEntrez[, 3])
entrezLinks <- cbind(entrez1 = MMMAnnoEntrez[match(MMMLinksWhole[, 1], MMMAnnoEntrez[, 1]), 3],
                     entrez2 = MMMAnnoEntrez[match(MMMLinksWhole[, 2], MMMAnnoEntrez[, 1]), 3])
MMMLinksWhole <- cbind(MMMLinksWhole[, 1:4], entrezLinks, mmm = MMMLinksWhole[, 5])

## select thres links
thresLinks <- MMMLinksWhole[(MMMLinksWhole[, 5] %in% rownames(phyloP)) & (MMMLinksWhole[, 6] %in% rownames(phyloP)), , drop = FALSE]
thresPer <- cbind(p1 = rowSums(phyloP[match(thresLinks[, 5], rownames(phyloP)), ])/ncol(wholePhyloDataNet),
                  p2 = rowSums(phyloP[match(thresLinks[, 6], rownames(phyloP)), ])/ncol(wholePhyloDataNet))
rownames(thresPer) <- NULL
thresPer <- apply(thresPer, 1:2, function(x) {
  return(paste0(round(x*100, 1), '%'))
})
thresLinks <- cbind(thresLinks[, 1:6],
                    thresPer,
                    mmm = thresLinks[, 7])

##~~~~~~~~~~~~~~~~~~~~~~~~top number~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

SimJaccardR <- function(pairProfile) {
  f <- pairProfile[, 1]
  t <- pairProfile[, 2]

  A <- sum((f + 2*t) == 3)

  jac <- A / (sum(f) + sum(t) - A)
  return(jac)
}

linkMat <- thresLinks[, 5:6]
linkMat <- apply(linkMat, 1:2, as.character)
linkjc <- foreach(i = seq_len(nrow(linkMat)), .combine = rbind) %do% {
  eachP <- t(phyloP[match(linkMat[i, ], rownames(phyloP)), ])
  return(c(SimJaccardR(eachP),
           cor(eachP)[1, 2]))
}
linkMat <- cbind(linkMat, linkjc)

distTop <- GetPosJacBatch('../jaccardSim.desc', '../wholeCor.desc', linkMat, 0, n = 4)
thresLinksRes <- cbind(thresLinks[, 1:8],
                    Jaccard = round(as.numeric(linkMat[, 3]), 2),
                    Cor = round(as.numeric(linkMat[, 4]), 2),
                    mmm = thresLinks[, 9],
                    Top = distTop)

save(thresLinksRes, file = 'thresLinksRes.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################

#############################plot thresLink########################
library('ggplot2')

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/thresLinksRes.RData')

p1 <- as.numeric(sub('%', '', as.character(thresLinksRes[, 7])))/100
p2 <- as.numeric(sub('%', '', as.character(thresLinksRes[, 8])))/100

thresL <- data.frame(p1 = p1,
                     p2 = p2,
                     Top = thresLinksRes[, 12])

PresentAcc <- function(thresMat, pThres, tThres = 50) {
  passLogic <- (thresMat[, 1] > pThres) & (thresMat[, 2] > pThres)
  thresM <- thresMat[passLogic, ]
  topNum <- thresM[, 3]

  tVec <- seq(min(topNum), max(topNum), tThres)
  tAcc <- numeric(length(tVec))

  for (i in seq_along(tAcc)) {
    tAcc[i] <- sum(topNum <= tVec[i])
  }

  hitM <- data.frame(hitRate = tAcc/length(topNum),
                     top = tVec)

  return(hitM)
}

tmpM <- PresentAcc(thresL, 0.95)
linkM <- rbind(PresentAcc(thresL, 0.95),
               PresentAcc(thresL, 0.85),
               PresentAcc(thresL, 0.75),
               PresentAcc(thresL, 0.65),
               PresentAcc(thresL, 0.55))
pdf('MMM_compare.pdf')
linkM <- cbind(linkM, thres = rep(c(0.95, 0.85, 0.75, 0.65, 0.55), each = nrow(tmpM)))
ggplot(linkM, aes(top, hitRate, colour = factor(thres))) +
  geom_line() +
  xlim(0, 3000) +
  xlab('Top rank') +
  ylab('Hit rate')
dev.off()
###################################################################
