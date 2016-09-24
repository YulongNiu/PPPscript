##########################top400 for other genes#######################
#~~~~~~~~~~~~~~~~~~TSPO interaction~~~~~~~~~~~~~~~~~~~~~~~
load('top400Anno.RData')

TSPOInter <- top400Anno[top400Anno[, 1] == 'TSPO' | top400Anno[, 2] == 'TSPO', ]
write.csv(TSPOInter, 'TSPOInter.csv')
write.table(TSPOInter, file = 'TSPOInter.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

## select other interaction
interVec <- unique(c(TSPOInter[, 1], TSPOInter[, 2]))
interVec <- interVec[interVec != 'TSPO']
otherInterLogic <- (top400Anno[, 1] %in% interVec) & (top400Anno[, 2] %in% interVec)
otherInter <- top400Anno[otherInterLogic, ]
wholeInter <- rbind(TSPOInter, otherInter)
write.table(wholeInter, file = 'TSPOwholeInter.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

## for cytocape color
mitoColor <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/mito/cytonetwork/top400mitoNode.txt', header = TRUE)
mitoColor <- apply(mitoColor, 1:2, as.character)

TSPOGenes <- c(interVec, 'TSPO')
mitoGenesComplex <- rep('other', length(TSPOGenes))
mitoGenesMS <- rep('N', length(TSPOGenes))
for(i in 1:length(TSPOGenes)) {
  ind <- match(TSPOGenes[i], mitoColor[, 1])
  if (!is.na(ind)) {
    mitoGenesMS[i] <- 'Y'
   if (mitoColor[ind, 2] != 'other') {
     mitoGenesComplex[i] <- mitoColor[ind, 2]
   } else {}
  } else {}
}
TSPOColorMat <- cbind(TSPOGenes, mitoGenesComplex, mitoGenesMS)
colnames(TSPOColorMat)[1] <- 'mitoGenesAnno'
write.table(TSPOColorMat, file = 'TSPOColorMat.txt', row.names = FALSE, quote = FALSE, sep = '\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################################################################
