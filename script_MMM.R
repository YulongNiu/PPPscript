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
  currentIdx <- match(mAnno[i], cHGNC)
  if (is.na(currentIdx)) {

    for (j in seq_along(pHGNC)) {
      eachP <- unlist(strsplit(pHGNC[j], split = ',', fixed = TRUE))
      eachP <- str_trim(eachP)
      if (sum(eachP %in% mAnno[i]) > 0) {
        currentIdx <- j
        break
      } else {}
    }

  } else {}
  geneIdx[i] <- currentIdx
}

MMMAnnoEntrez <- cbind(MMMAnno, entrez = HGNCAnno[geneIdx, 11])

## manually check

save(MMMLinksWhole, MMMAnno, file = 'MMM12.RData')


#################################################################

################select all present#####################
library('biomaRt')

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/MMM/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/wholePhyloData.RData')

## set threshold 95%
thres <- 0.95
phyloP <- wholePhyloDataNet[rowSums(wholePhyloDataNet) >= floor(ncol(wholePhyloDataNet) * thres), , drop = FALSE]

## annotation select genes
geneEntrez <- strsplit(rownames(phyloP),
                       split = ':',
                       fixed = TRUE)
geneEntrez <- sapply(geneEntrez, '[[', 2)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
hsaAnno <- getBM(attributes = c('entrezgene', 'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'entrezgene', value = geneEntrez, mart = ensembl)

## HGNC annotation
HGNCAnno <- read.csv('HGNC_anno.txt', sep = '\t', stringsAsFactor = FALSE)

#######################################################
