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
setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/wholePhyloData.RData')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/MMM12.RData')

## set threshold 95%
thres <- 0.95
phyloP <- wholePhyloDataNet[rowSums(wholePhyloDataNet) >= floor(ncol(wholePhyloDataNet) * thres), , drop = FALSE]

## check MMM12
MMMAnnoEntrez[, 3] <- paste0('hsa:', MMMAnnoEntrez[, 3])
sum(rownames(phyloP) %in% MMMAnnoEntrez[, 3])
#######################################################
