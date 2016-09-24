#~~~~~~~~~~~~~~~~~~~~~~~~~~Jaccard Similarity~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('vegan')
require('compiler')
load('wholePhyloData.RData')
enableJIT(3)

jaccardSim <- 1 - vegdist(wholePhyloDataNet[1:10, ], method = 'jaccard', diag = TRUE, upper = TRUE)
jaccardSim <- as.matrix(jaccardSim)
jaccardSim <- diag(1, ncol(jaccardSim), nrow(jaccardSim)) + jaccardSim
save(jaccardSim, file = 'jaccardSim.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ parallel and openMPI ~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
load('jaccardSim.RData')
upperSeq <- jaccardSim[upper.tri(jaccardSim)]
save(upperSeq, file = 'upperjaccardSeq.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('bigmemory')
load('jaccardSim.RData')
load('upperjaccardSeq.RData')
allRowColLogicData <- attach.big.matrix('allRowColLogicData.desc')
rowNames <- rownames(jaccardSim)
colNames <- colnames(jaccardSim)
rowNamesUpper <- rowNames[allRowColLogicData[, 1]]
colNamesUpper <- colNames[allRowColLogicData[, 2]]
rm(allRowColLogicData)
rm(jaccardSim)
rm(colNames, rowNames)
gc()
jaccardft <- cbind(rowNamesUpper, colNamesUpper, upperSeq)
save(jaccardft, file = 'jaccardft.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('jaccardft.RData')
cutoff04 <- which(as.numeric(jaccardft[, 3]) >= 0.4)
cutoff04jaccardMat <- jaccardft[cutoff04, ]
save(cutoff04jaccardMat, file = 'cutoff04jaccardMat.RData')
rm(jaccardft, cutoff04)
gc()

cutoff05 <- which(as.numeric(cutoff04jaccardMat[, 3]) >= 0.5)
cutoff05jaccardMat <- cutoff04jaccardMat[cutoff05, ]
save(cutoff05jaccardMat, file = 'cutoff05jaccardMat.RData')
rm(cutoff04jaccardMat, cutoff05)
gc()

cutoff06 <- which(as.numeric(cutoff05jaccardMat[, 3]) >= 0.6)
cutoff06jaccardMat <- cutoff05jaccardMat[cutoff06, ]
save(cutoff06jaccardMat, file = 'cutoff06jaccardMat.RData')
rm(cutoff05jaccardMat, cutoff06)
gc()

cutoff07 <- which(as.numeric(cutoff06jaccardMat[, 3]) >= 0.7)
cutoff07jaccardMat <- cutoff06jaccardMat[cutoff07, ]
save(cutoff07jaccardMat, file = 'cutoff07jaccardMat.RData')
rm(cutoff06jaccardMat, cutoff07)
gc()

cutoff08 <- which(as.numeric(cutoff07jaccardMat[, 3]) >= 0.8)
cutoff08jaccardMat <- cutoff07jaccardMat[cutoff08, ]
save(cutoff08jaccardMat, file = 'cutoff08jaccardMat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~ATP with Entrez names~~~~~~~~~~~~~~~~~~~~~
load('cutoff08jaccardMat.RData')
load('mitoPhyloData.RData')
fromList <- cutoff08jaccardMat[, 1] %in% ATPgeneList[, 2]
toList <- cutoff08jaccardMat[, 2] %in% ATPgeneList[, 2]
fullList <- fromList | toList
rm(fromList, toList)
gc()
cutoff08ATP <- cutoff08jaccardMat[fullList, ]
save(cutoff08ATP, file = 'cutoff08jaccardATP.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



###############################################################################


