#########################circos functions#################################
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

  # get nodes genome location information
  weight = as.numeric(as.character(ft[, 3]))
  nodeIn <- locaAnno[match(ft[, 1], locaAnno[, 1]), 1:4]
  nodeOut <- locaAnno[match(ft[, 2], locaAnno[, 1]), 1:4]

  # get gene lables data
  circosLabel <- rbind(nodeIn[, c(2:4, 1)], nodeOut[, c(2:4, 1)])
  circosLabel <- circosLabel[!duplicated(circosLabel, MARGIN = 1), ]
  rownames(circosLabel) <- NULL

  # get link data
  circosLink <- cbind(nodeIn[, 2:4], nodeOut[, 2:4])
  rownames(circosLink) <- NULL

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
      ## heatVal <- apply(heatVal, 1:2, as.character)
      ## heatVal <- rbind(heatVal, c(as.character(circosLink[1, 1:3]), 1))
    } else {}
  } else {
    heatVal <- circosLink[, 1:4, drop = FALSE]
  }

  return(list(circosLink = circosLink, circosLabel = circosLabel, circosHeatmap = heatVal))

}

SpeFreq <- function(KEGGPhylo, selectPhyloData, splitEu = FALSE) {
  # USE: calculate gene frequence of three domains of life.
  # INPUT: 'KEGGPhylo' is a 2-column matrix, the first column is the KEGGID, and the second column is species phylo. 'selectPhyloData' is a part/whole phylogeneti 0-1matrix, in which the colnames is KEGGIDs and the row names is the gene IDs. 'splitEu' whether to show the eukaryotic split.
  # OUTPUT:
  # INPUT EXAMPLE
  # KEGGPhylo file
  ##     KEGGID                                  Phylo
  ## 1    hsa Eukaryotes;Animals;Vertebrates;Mammals
  ## 2    ptr Eukaryotes;Animals;Vertebrates;Mammals
  ## 3    pps Eukaryotes;Animals;Vertebrates;Mammals
  ## 4    ggo Eukaryotes;Animals;Vertebrates;Mammals
  ## 5    pon Eukaryotes;Animals;Vertebrates;Mammals
  ## 6    mcc Eukaryotes;Animals;Vertebrates;Mammals
  # selectPhyloData file
  ##              aac aae aag aao aar aas aau aba abe abi
  ## hsa:60509    0   0   1   0   0   0   0   0   0   0
  ## hsa:57223    0   0   1   0   0   0   0   0   1   0
  ## hsa:56941    0   0   1   0   0   0   1   1   1   0
  ## hsa:66035    1   0   1   0   0   0   1   1   1   0
  ## hsa:55635    0   0   0   0   0   0   0   0   0   0
  ## hsa:66000    0   0   0   0   0   0   0   0   0   0
  ## hsa:171568   0   0   1   0   0   0   0   0   1   1
  ## hsa:65059    0   0   1   0   0   0   0   0   0   0
  ## hsa:9699     0   0   1   0   0   0   0   0   1   0
  ## hsa:317649   0   0   1   0   0   0   0   0   1   0

  phyloCode <- as.character(KEGGPhylo[, 2])
  phyloCode <- sapply(strsplit(phyloCode, split = ';', fixed = TRUE), function(x) x[1:2])
  phyloCode <- t(phyloCode)
  phyloCode <- cbind(as.character(KEGGPhylo[, 1]), phyloCode)

  phyloList <- list()
  phyloList$arcSpe <- phyloCode[phyloCode[, 3] %in% 'Archaea', 1]
  phyloList$bacSpe <- phyloCode[phyloCode[, 3] %in% 'Bacteria', 1]
  if (!splitEu) {
    phyloList$euSpe <- phyloCode[!(phyloCode[, 3] %in% c('Archaea', 'Bacteria')), 1]
  } else {
    phyloList$Animals <- phyloCode[phyloCode[, 3] %in% 'Animals', 1]
    phyloList$Plants <- phyloCode[phyloCode[, 3] %in% 'Plants', 1]
    phyloList$Fungi <- phyloCode[phyloCode[, 3] %in% 'Fungi', 1]
    phyloList$Protists <- phyloCode[phyloCode[, 3] %in% 'Protists', 1]
  }

  phyloFreqList <- lapply(phyloList, function(x) {
    phyloEachData <- selectPhyloData[, colnames(selectPhyloData) %in% x]
    phyloEachFreq <- apply(phyloEachData, 1, function(x) {sum(x)/length(x)})
    return(phyloEachFreq)
  })

  phyloFreqMat <- do.call(cbind, phyloFreqList)

  # asign names
  if (!splitEu) {
    colnames(phyloFreqMat) <- c('freqArc', 'freqBac', 'freqEu')
  } else {
    colnames(phyloFreqMat) <- c('freqArc', 'freqBac', 'freqAni', 'freqPlant', 'freqFun', 'freqProt')
  }
  
  return(phyloFreqMat)
}

CombinePredTrue <- function(predMat, trueMat, trueWeight = 0.4) {
  # USE: combine prediction interaction and true interaction
  # INPUT: 'predMat' is the prection interaction with weight. 'trueMat' is the true interaction with weight such as 1. Both 'predMat' and 'trueMat' should have the same nodes order. 'trueWeight' is FN weight
  # OUTOUT: combined mat, 'TP' use predition weight, 'FP' use prediction weight, 'FN' use 'trueWeight'.

  trueMat[, 3] <- rep(as.character(trueWeight), nrow(trueMat))
    
  predVec <- predMat[, 3]
  names(predVec) <- paste(predMat[, 1], predMat[, 2], sep = '|')

  trueVec <- trueMat[, 3]
  names(trueVec) <- paste(trueMat[, 1], trueMat[, 2], sep = '|')

  TPVec <- predVec[names(predVec) %in% names(trueVec)]
  FPVec <- predVec[!(names(predVec) %in% names(trueVec))]
  FNVec <- trueVec[!(names(trueVec) %in% names(predVec))]

  combineVec <- c(TPVec, FPVec, FNVec)
  combineMat <- cbind(combineVec,
                      rep(c('TP', 'FP', 'FN'), c(length(TPVec), length(FPVec), length(FNVec))))

  nodeList <- strsplit(names(combineVec), split = '|', fixed = TRUE)
  nodeMat <- do.call(rbind, nodeList)

  combineMat <- cbind(nodeMat, combineMat)

  rownames(combineMat) <- NULL
  colnames(combineMat) <- c('from', 'to', 'weight', 'pred')

  return(combineMat)
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
######################################################################



####################################plot circos#########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ annotation of gene postition~~~~~~~~~~~~~
require('biomaRt')
load('wholePhyloData.RData')

# annotation file
hsaWholeEntrez <- names(wholePhyloDataAnno)
hsaEntrez <- sapply(strsplit(hsaWholeEntrez, split = ':', fixed = TRUE), '[[', 2)
hsaSymbol <- annoFirst
hsaEntrezSymbol <- cbind(hsaWholeEntrez, hsaEntrez, hsaSymbol)
colnames(hsaEntrezSymbol)[1:2] <- c('KEGGID', 'KEGGAnno')
hsaEntrezSymbol <- as.data.frame(hsaEntrezSymbol)

# biomaRt annotation to get gene position 'from' and 'to'
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
hsaAnno <- getBM(attributes = c('entrezgene', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'entrezgene', value = hsaEntrez, mart = ensembl)
# remove ones without chromesome annotation
chromName <- c(1:22, 'X', 'Y', 'MT')
hsaAnno <- hsaAnno[hsaAnno[, 3] %in% chromName, ]
# get the first one
hsaAnno <- aggregate(hsaAnno, by = list(hsaAnno[, 1]), function(x) x[1])
hsaAnno <- hsaAnno[, -1]

# if no gene symbol, using 'hsa:***' instead.
wholeAnno <- merge(hsaAnno, hsaEntrezSymbol, by.x = 'entrezgene', by.y = 'KEGGAnno')
wholeAnno[, 'hsaSymbol'] <- wholeAnno[, 'hgnc_symbol']
withoutAnno <- as.character(wholeAnno[wholeAnno[, 'hsaSymbol'] == '', 6])
wholeAnno[which(wholeAnno[, 'hsaSymbol'] == ''), 'hsaSymbol'] <- withoutAnno

## transfer the MT-*** symbols to ***
mtGenes <- wholeAnno[wholeAnno[, 3] == 'MT', 7]
mtGenes <- strsplit(mtGenes, split = '-', fixed = TRUE)
mtGenes <- sapply(mtGenes, function(x){return(x[x != 'MT'])})
wholeAnno[wholeAnno[, 3] == 'MT', 7] <- mtGenes

## remove duplicate rows
wholeAnno <- wholeAnno[!duplicated(wholeAnno[, 7]), ]

save(wholeAnno, file = 'wholeAnno.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot ATP interacted~~~~~~~~~~~~~~~~~~~~~~
load('wholeAnno.RData')
load('top400.RData')
load('wholePhyloData.RData')
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
phyloSpe <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)
phyloSpe <- phyloSpe[, c(2, 4)]

# get FATP complex
top400ATP <- top400[(top400[, 1] %in% atpSub[, 2]) | (top400[, 2] %in% atpSub[, 2]), ]

# KEGG gene IDs to annotated symbols
annoNodesMat <- apply(top400ATP[, 1:2], 1:2, function(x){
  eachSymbol <- as.character(wholeAnno[match(x, wholeAnno[, 6]), 7])
  return(eachSymbol)
})
top400ATP[, 1:2] <- annoNodesMat
hasNA <- is.na(annoNodesMat[, 1]) | is.na(annoNodesMat[, 2])
top400ATP <- top400ATP[!hasNA, ]

labelGene <- ft2circos(top400ATP[, 1:3], wholeAnno[, c(7, 3:5)], atpSub[1:6, 1])[[2]]
write.table(labelGene, 'labelGene.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

for (i in 1:17){
  corPhyloCir <- ft2circos(top400ATP[, 1:3], wholeAnno[, c(7, 3:5)], atpSub[i, 1], showEdge = FALSE, thick = 3)
  write.table(corPhyloCir[[1]], paste(atpSub[i, 1], 'txt', sep = '.'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(corPhyloCir[[3]], paste(atpSub[i, 1], 'heatmap', '.txt', sep = ''), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
}


freqGene <- SpeFreq(phyloSpe, wholePhyloDataNet)
labelKEGGID <- as.character(wholeAnno[match(labelGene[, 4], wholeAnno[, 7]), 6])
for (i in 1:ncol(freqGene)) {
  freqVec <- freqGene[match(labelKEGGID, rownames(freqGene)), i]
  freqMat <- cbind(labelGene[, 1:3], freqVec)
  write.table(freqMat, paste(colnames(freqGene)[i], 'txt', sep = '.'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################


###################################plot complex/pathway#######################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ plot F-ATP complex~~~~~~~~~~~~~~~~~~~~~~~~~
load('wholeAnno.RData')
load('top400.RData')
load('wholePhyloData.RData')
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
phyloSpe <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)
phyloSpe <- phyloSpe[, c(2, 4)]


# get FATP complex
FATPcomplexPred <- top400[(top400[, 1] %in% atpSub[, 2]) & (top400[, 2] %in% atpSub[, 2]), ]
FATPcomplexPred[(FATPcomplexPred[, 1] %in% 'hsa:6648') & (FATPcomplexPred[, 2] %in% 'hsa:6648'), ]


# combine all the true interaction
FATPcomplexTrue <- read.csv('complex/FATPcomplex_true.csv')
FATPcomplexTrue <- apply(FATPcomplexTrue, 1:2, as.character)
FATPcomplexTrue <- cbind(FATPcomplexTrue, rep(1, nrow(FATPcomplexTrue)))

# combine prediction data and true data
FATPcomplex <- CombinePredTrue(FATPcomplexPred, FATPcomplexTrue, trueWeight = 0.5)

# KEGG gene IDs to annotated symbols
annoNodesMat <- apply(FATPcomplex[, 1:2], 1:2, function(x){
  eachSymbol <- as.character(wholeAnno[match(x, wholeAnno[, 6]), 7])
  return(eachSymbol)
})
FATPcomplex[, 1:2] <- annoNodesMat
hasNA <- is.na(annoNodesMat[, 1]) | is.na(annoNodesMat[, 2])
FATPcomplex <- FATPcomplex[!hasNA, ]

# gene labels
labelGene <- ft2circos(FATPcomplex[, 1:3], wholeAnno[, c(7, 3:5)])[[2]]
write.table(labelGene, 'labelGene.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

# links with weight
corPhyloCir <- ft2circos(FATPcomplex[, 1:3], wholeAnno[, c(7, 3:5)], showEdge = TRUE, thick = 4)
corPhyloCir <- apply(corPhyloCir[[1]], 1:2, as.character)
for (i in 1:nrow(FATPcomplex)) {
  if (FATPcomplex[i, 4] == 'TP') {
    corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=dgreen_a1', sep = ',')
  }
  else if (FATPcomplex[i, 4] == 'FP') {
    corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=grey_a1', sep = ',')
  }
  else if (FATPcomplex[i, 4] == 'FN') {
    corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=dpurple_a1', sep = ',')
  }
}
write.table(corPhyloCir, 'links.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

# gene freqency
freqGene <- SpeFreq(phyloSpe, wholePhyloDataNet, splitEu = TRUE)
labelKEGGID <- as.character(wholeAnno[match(labelGene[, 4], wholeAnno[, 7]), 6])
for (i in 1:ncol(freqGene)) {
  freqVec <- freqGene[match(labelKEGGID, rownames(freqGene)), i]
  freqMat <- cbind(labelGene[, 1:3], freqVec)
  write.table(freqMat, paste(colnames(freqGene)[i], 'txt', sep = '.'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot MAPK/TCA/GABA~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('pathway/pathFilter.RData')
load('complexAll/comInterJacCor.RData')
load('wholeAnno.RData')
load('top400.RData')
load('wholePhyloData.RData')
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
phyloSpe <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)
phyloSpe <- phyloSpe[, c(2, 4)]

# deal with complex
humComp <- lapply(MIPSHumList, function(x) {
  eachHumMat <- x[[2]][, 1:2, drop = FALSE]
  # order each mat
  eachHumMat <- OrderHumMat(eachHumMat)
  return(eachHumMat)
})

## allInter <- apply(top400[, 1:2], 1, paste, collapse = '|')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get true pathway
pathwayTrue <- biocartaPathFilter$`mapkinase signaling pathway`
pathwayTrue <- keggPathFilter$`Citrate cycle (TCA cycle)`
pathwayTrue <- reactomePathFilter$`GABA A receptor activation`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathwayTrue <- cbind(pathwayTrue, rep(1, nrow(pathwayTrue)))
pathwayNodes <- unique(c(pathwayTrue[, 1], pathwayTrue[, 2]))
pathwayPred <- top400[(top400[, 1] %in% pathwayNodes) & (top400[, 2] %in% pathwayNodes), ]
## pathwayPred <- top400[allInter %in% pathwayTrueVec, ]
pathwayMat <- CombinePredTrue(pathwayPred, pathwayTrue, trueWeight = 0.5)

# KEGG gene IDs to annotated symbols
annoNodesMat <- apply(pathwayMat[, 1:2], 1:2, function(x){
  eachSymbol <- as.character(wholeAnno[match(x, wholeAnno[, 6]), 7])
  return(eachSymbol)
})
pathwayMat[, 1:2] <- annoNodesMat
hasNA <- is.na(annoNodesMat[, 1]) | is.na(annoNodesMat[, 2])
pathwayMat <- pathwayMat[!hasNA, ]

# gene labels
labelGene <- ft2circos(pathwayMat[, 1:3], wholeAnno[, c(7, 3:5)])[[2]]
write.table(labelGene, 'labelGene.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

# links with weight
corPhyloCir <- ft2circos(pathwayMat[, 1:3], wholeAnno[, c(7, 3:5)], showEdge = TRUE, thick = 4)
corPhyloCir <- apply(corPhyloCir[[1]], 1:2, as.character)
for (i in 1:nrow(pathwayMat)) {
  if (pathwayMat[i, 4] == 'TP') {
    corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=dgreen_a1', sep = ',')
  }
  ## else if (pathwayMat[i, 4] == 'FP') {
  ##   corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=grey_a1', sep = ',')
  ## }
  else if (pathwayMat[i, 4] == 'FN') {
    corPhyloCir[i, 7] <- paste(corPhyloCir[i, 7], 'color=dpurple_a1', sep = ',')
  }
}
write.table(corPhyloCir, 'links.txt', sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)

# gene freqency
freqGene <- SpeFreq(phyloSpe, wholePhyloDataNet, splitEu = TRUE)
labelKEGGID <- as.character(wholeAnno[match(labelGene[, 4], wholeAnno[, 7]), 6])
for (i in 1:ncol(freqGene)) {
  freqVec <- freqGene[match(labelKEGGID, rownames(freqGene)), i]
  freqMat <- cbind(labelGene[, 1:3], freqVec)
  write.table(freqMat, paste(colnames(freqGene)[i], 'txt', sep = '.'), sep = ' ', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################################




#################### generate combn################

allOrderdComb <- function(nodeVec) {
  # USE: generate all the combination of nodes and order the interaction (order character, like '3412' and '345251')
  # INPUT: 'nodeVec' is the KEGG ID vectors, for example, c('hsa:2345', 'hsa:3435')

  allComb <- t(combn(nodeVec, 2))
  
  fromNum <- sapply(strsplit(allComb[, 1], split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(allComb[, 2], split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  for (i in 1:length(isSmall)) {
    if (!isSmall[i]) {
      allComb[i, 1:2] <- allComb[i, 2:1]
    } else {}
  }
  ## interMat <- rbind(allComb[isSmall, ], allComb[!isSmall, 2:1])

  return(allComb)
}

# FATP synthase combination
allOrderdComb(as.character(atpSub[, 2]))
###################################################




######################generate data for web-based software#########
library(stringr)
geneAnno <- wholeAnno[, c(6:7, 3:5)]
geneAnno <- apply(geneAnno, 1:2, as.character)
geneAnno <- apply(geneAnno, 1:2, str_trim)
####################################################################
