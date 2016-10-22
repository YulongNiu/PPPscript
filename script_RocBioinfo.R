setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

######################get TP/NP###########################
library('bigmemory')

load('uniProt.RData')
load('wholePhyloData.RData')

uniProtIDs <- sapply(strsplit(uniProt[, 2], split = ':', fixed = TRUE), '[[', 2)
uniProt[, 2] <- uniProtIDs

entrez <- rownames(wholePhyloDataNet)

# true positive interaction
PRS <- read.csv('complex/data_files/PRS_scores_human.tab', sep = '\t')
PRS <- PRS[, 1:2]
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

allRS <- rbind(PRS, NRS)
allRS <- cbind(allRS, rep(c('TP', 'TN'), c(nrow(PRS), nrow(NRS))))
colnames(allRS) <- c('From', 'To', 'Status')

save(allRS, file = 'complexAll/allRS_Bioinfo.RData')
##########################################################
