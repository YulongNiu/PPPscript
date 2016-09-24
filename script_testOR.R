###########################test OR###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~seperate species by taxa~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~seperate species by domain~~~~~~~~~~~~~~~~~~~~
library('foreach')
library('doMC')
registerDoMC(4)
load('wholePhyloData.RData')
KEGGSpeMat <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)

# separate species by domain
KEGGphylo <- as.character(KEGGSpeMat$Phylo)
domain <- foreach(i = 1:length(KEGGphylo), .combine = c) %dopar% {
  firstPhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[1]
  if (firstPhylo == 'Eukaryotes') {
    return(firstPhylo)
  } else {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
    return(sePhylo)
  }
}
KEGGID <- as.character(KEGGSpeMat$KEGGID)
names(KEGGID) <- domain

# sort KEGGID
KEGGID <- sort(KEGGID)
KEGGID <- KEGGID[rank(colnames(wholePhyloDataNet))]
domain <- names(KEGGID)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~seperate species by kindom~~~~~~~~~~~~~~~~~~~~~~~~
# separate species by kindom
KEGGphylo <- as.character(KEGGSpeMat$Phylo)
kingdom <- foreach(i = 1:length(KEGGphylo), .combine = c) %dopar% {
  firstPhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[1]
  if (firstPhylo == 'Eukaryotes') {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
  } else {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
  }
  return(sePhylo)
}
KEGGID <- as.character(KEGGSpeMat$KEGGID)
names(KEGGID) <- kingdom

# sort KEGGID
KEGGID <- sort(KEGGID)
KEGGID <- KEGGID[rank(colnames(wholePhyloDataNet))]
kingdom <- names(KEGGID)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~cate function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CatePhylo <- function(phyloMat, cateData){
  # USE: category the phylogenetic profiles data
  # INPUT: 'phyloMat' is a pair of phylogenetic profile data binded by column. 'catData' is vector in the same order of 'phyloMat'.
  # OUTPUT: list of binary matrix

  Bin2Mat <- function(binMat) {
    
    ## USE: transfer the binary 0-1 data to matrix
    ## INPUT: 'binMat' is the binary 0-1 data with column binded.
    ## OUTPUT: The binary matrix, which could be used for fisher.test()
    ## EXAMPLE (Barker 2005 PLoS Comp):
    ## CIN4 <- c(rep(0, 8), rep(1, 7))
    ## ORC3 <- c(rep(0, 7), rep(1, 8))
    ## Bin2Mat(cbind(CIN4, ORC3))
    ## L9A <- c(rep(0, 4), 1, rep(0, 2), 1, rep(0, 5), rep(1, 2))
    ## L42B <- c(rep(1, 8), rep(0, 5), rep(1, 2))
    ## Bin2Mat(cbind(L9A, L42B))

    seq1 <- binMat[, 1]
    seq2 <- binMat[, 2]
    
    sumSeq <- seq1 + seq2
    difSeq <- seq1 - seq2
    
    a <- sum(sumSeq == 2)

    b <- sum(difSeq == 1)

    c <- sum(difSeq == -1)

    d <- sum(sumSeq == 0)

    Mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    rownames(Mat) <- c('1', '0')
    colnames(Mat) <- c('1', '0')
    
    return(Mat)
  }

  cateNames <- unique(cateData)
  binList <- list()
  
  for (i in 1:length(cateNames)) {
    phyloMatCat <- phyloMat[cateData %in% cateNames[i], ,drop = FALSE]
    binList[[i]] <- Bin2Mat(phyloMatCat)
  }

  names(binList) <- cateNames
  return(binList)
}


# alpha and beta subunits
atpGenes <- c('hsa:498', 'hsa:506')

# c1c2
atpGenes <- c('hsa:516', 'hsa:517')

# neg 
atpGenes <- c('hsa:4536', 'hsa:1390')

# pos
atpGenes <- c('hsa:596', 'hsa:581')

atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% atpGenes, ])
CatePhylo(atpMat, kingdom)
sapply(CatePhylo(atpMat, kingdom), function(x) {fisher.test(x)$p.value})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test TP TN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
require('foreach')
require('doMC')
registerDoMC(4)
jaccardSim <- attach.big.matrix('jaccardSim.desc')
wholeCor <- attach.big.matrix('wholeCor.desc')

TP <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phylotree/true_positive.csv', header = FALSE)
TP <- TP[, -1:-3]
TP <- apply(TP, 1:2, as.character)
TN <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phylotree/true_negative.csv', header = FALSE)
TN <- TN[, -1:-3]
TN <- apply(TN, 1:2, as.character)
TN <- TN[1:13, ]

# merged jaccard/cor/kingdom information
comInterJacCor <- foreach(i = 1:nrow(TP), .combine = rbind) %dopar% {
  rowNum <- match(TP[i, 1], rownames(jaccardSim))
  colNum <- match(TP[i, 2], rownames(jaccardSim))
  jacsim <- jaccardSim[rowNum, colNum]
  corValue <- wholeCor[rowNum, colNum]
  atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% TP[i, 1:2], ])
  cateVec <- unlist(CatePhylo(atpMat, kingdom))
  jacCorVec <- c(TP[i, ], jacsim, corValue, cateVec)
  return(jacCorVec)
}

## colnames(comInterJacCor) <- c('proteinA', 'proteinB', 'JaccardSim', 'Cor', paste('Bacteria', c('11', '10', '01', '00'), sep = ''), paste('Eukaryotes', c('11', '10', '01', '00'), sep = ''), paste('Archaea', c('11', '10', '01', '00'), sep = ''))
colnames(comInterJacCor)[3:4] <- c('JaccardSim', 'Cor')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~Jaccard Similarity Coefficient~~~~~~~~~~~~~~~~~~~
load('cutoff08jaccardATP.RData')

jacCor <- foreach(i = 1:nrow(cutoff08ATP), .combine = rbind) %dopar% {
  rowNum <- match(cutoff08ATP[i, 1], rownames(jaccardSim))
  colNum <- match(cutoff08ATP[i, 2], rownames(jaccardSim))
  jacsim <- jaccardSim[rowNum, colNum]
  corValue <- wholeCor[rowNum, colNum]
  atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% cutoff08ATP[i, 1:2], ])
  cateVec <- unlist(CatePhylo(atpMat, kingdom))
  jacCorVec <- c(cutoff08ATP[i, 1:2], jacsim, corValue, cateVec)
  return(jacCorVec)
}

colnames(jacCor)[3:4] <- c('JaccardSim', 'Cor')

tmp1 <- foreach(i = 1:nrow(jacCor), .combine = c) %dopar% {
  if (jacCor[i, 1] %in% atpSub[, 1]) {
    return(jacCor[i, 1])
  } else {return(jacCor[i, 2])}
}


jacCor[, 1] <- annoFirst[match(jacCor[, 1], names(annoFirst))]
jacCor[, 2] <- annoFirst[match(jacCor[, 2], names(annoFirst))]
write.csv(jacCor, 'cutoff08ATPSimCate.csv')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~test ATP synthase~~~~~~~~~~~~~~~~~~~~~~~~~
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
atpSubMat <- t(combn(as.character(atpSub[, 2]), 2))
atpSubMatAnno <- t(combn(as.character(atpSub[, 1]), 2))

ATPCate <- foreach (i = 1:nrow(atpSubMat), .combine = rbind) %dopar% {
  atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% atpSubMat[i, 1:2], ])
  cateVec <- unlist(CatePhylo(atpMat, kindom))
  jaccardSim <- as.character(1 - vegdist(t(atpMat), method = 'jaccard'))
  corvalue <- cor(atpMat)
  corvalue <- corvalue[1, 2]
  return(c(cateVec, jaccardSim, corvalue))
}
colnames(ATPCate)[13:14] <- c('Jaccard', 'Cor')
ATPCate <- cbind(atpSubMatAnno, ATPCate)
write.csv(ATPCate, 'ATPCate.csv')


atpSubMat <- read.csv('cutoff08ATPcorSim.csv', row.names = 1)
atpSubMat <- apply(atpSubMat, 1:2, as.character)
ATPCate <- foreach (i = 1:nrow(atpSubMat), .combine = rbind) %dopar% {
  atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% atpSubMat[i, 1:2], ])
  cateVec <- unlist(CatePhylo(atpMat, kindom))
  return(cateVec)
}

atpSubMat <- cbind(atpSubMat, ATPCate)
atpSubMat[, 1] <- annoFirst[match(atpSubMat[, 1], names(annoFirst))]
atpSubMat[, 2] <- annoFirst[match(atpSubMat[, 2], names(annoFirst))]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~test ATP at Cor>0.8~~~~~~~~~~~~~~~~~~~~~~~~
require('foreach')
require('doMC')
require('vegan')
registerDoMC(4)
load('wholePhyloData.RData')
load('cutoff08ATPAnno.RData')

## test1 <- ('hsa:652991', 'hsa:522')
## test1 <- c('hsa:127069', 'hsa:9551')
## test1 <- c('hsa:81566', 'hsa:9551')

cutoff08ATPCate <- foreach (i = 1:nrow(cutoff08ATP), .combine = rbind) %dopar% {
  print(paste('It is running ', i, ' in a total of ', nrow(cutoff08ATP), '.', sep = ''))
  atpMat <- t(wholePhyloDataNet[rownames(wholePhyloDataNet) %in% cutoff08ATP[i, 1:2], ])
  cateVec <- unlist(CatePhylo(atpMat, kindom))
  jaccardSim <- as.character(1 - vegdist(t(atpMat), method = 'jaccard'))
  return(c(cateVec, jaccardSim))
}

cutoff08ATPCate <- cbind(cutoff08ATPAnno, cutoff08ATPCate)
colnames(cutoff08ATPCate)[3] <- 'Cor'
colnames(cutoff08ATPCate)[16] <- 'Jaccard'
write.csv(cutoff08ATPCate, 'cutoff08ATPCate.csv')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#################################################################


