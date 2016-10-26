setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

##########################hand TP/FP/TN/FN#################################
HandTF <- function(valM, bin = 10000, increasing = TRUE) {

  require('foreach')
  require('doMC')
  registerDoMC(4)

  ## INPUT: "valM" is the measure data frame, 1st col is value, 2nd col is indicator. "bin" number of thresholds. "increasing" larger --> positive?
  ## OUTPUT: TP/FP/TN/FN matrix

  allVec <- valM[, 1]
  PVec <- allVec[valM[, 2] == 'TP']
  PNum <- length(PVec)
  NVec <- allVec[valM[, 2] == 'TN']
  NNum <- length(NVec)

  maxVal <- max(allVec)
  minVal <- min(allVec)

  cutVec <- seq(minVal, maxVal, length.out = bin)

  TFMat <- foreach(i = 1:length(cutVec), .combine = rbind) %dopar% {
    eachCut <- cutVec[i]
    TP <- sum(PVec >= eachCut)
    FN <- PNum - TP
    TN <- sum(NVec <= eachCut)
    FP <- NNum - TN

    return(c(TP, FP, TN, FN))
  }

  if (!increasing) {
    TFMat <- TFMat[, 4:1]
  } else {}

  return(TFMat)
}

PR <- function(TPMat) {
   return(cbind(TPMat[, 1]/rowSums(TPMat[, 1:2]), TPMat[, 1]/rowSums(TPMat[, c(1, 4)])))
}
#####################################################################

##############################test all complex#######################
library('pROC')
library('ggplot2')
library('foreach')
library('doMC')
registerDoMC(4)

## load file
load('complexAll/simdistROCNPP70_cutInf_seed456.RData')
corMatNPP <- corMat
corRocNPP <- corRoc
load('complexAll/simdistROCSVD100_cutInf_seed456.RData')
euMatSVD100 <- euMat
euRocSVD100 <- euRoc
load('complexAll/simdistROCSVD30_cutInf_seed456.RData')
euMatSVD30 <- euMat
euRocSVD30 <- euRoc

pat <- 'ROC_cutInf_seed456'
rocDataFiles <- dir('complexAll', pattern = pat, full.names = TRUE)
for(i in rocDataFiles) {load(i)}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P <- sum(topMat[, 2] == 'TP')
stList <- list(Top = topMat[1:(P*2), ],
               Tree = LRMat[1:(P*2), ],
               Dollo = DolloMat[1:(P*2), ],
               Cor = corMat[1:(P*2), ],
               Jaccard = jacMat[1:(P*2), ],
               MI = MIMat[1:(P*2), ],
               Hamming = hamMat[1:(P*2), ],
               NPP = corMatNPP[1:(P*2), ],
               SVD100 = euMatSVD100[1:(P*2), ],
               SVD30 = euMatSVD30[1:(P*2), ])
rocList <- foreach(i = 1:length(stList)) %dopar% {
  x <- stList[[i]]
  return(roc(x[, 2], x[, 1], levels = c('TP', 'TN')))
}
rocMatList <- foreach(i = 1:length(stList)) %dopar% {
  x <- rocList[[i]]
  return(cbind(1 - x$specificities, x$sensitivities))
}

## better plot
minRow <- min(sapply(rocMatList, nrow))
rocMatList <- list(Top = rocMatList[[1]],
                   Tree = rocMatList[[2]],
                   Dollo = rocMatList[[3]][round(seq(1, nrow(rocMatList[[3]]), length.out = minRow)), ],
                   Cor = rocMatList[[4]],
                   Jaccard = rocMatList[[5]][round(seq(1, nrow(rocMatList[[5]]), length.out = minRow)), ],
                   MI = rocMatList[[6]],
                   Hamming = rocMatList[[7]][round(seq(1, nrow(rocMatList[[7]]), length.out = minRow)), ],
                   NPP = rocMatList[[8]],
                   SVD100 = rocMatList[[9]],
                   SVD30 = rocMatList[[10]][round(seq(1, nrow(rocMatList[[10]]), length.out = minRow)), ]
                   )

## plot ROC
mergedRocMat <- do.call(rbind, rocMatList)
mergedRocMat <- data.frame(FPR = mergedRocMat[, 1],
                           TPR = mergedRocMat[, 2],
                           Methods = rep(names(stList), sapply(rocMatList, nrow)))

aucAnno <- paste0(names(stList), '=', round(sapply(rocList, function(x){return(x$auc)}), 3))
colAnno <- c('#F8766D', rep('#D39200', 2), rep('#00BA38', 2), rep('#00B9E3', 2),  '#FF61C3', rep('#DB72FB', 2))
names(colAnno) <- names(stList)
lineAnno <- c('solid', 'solid', 'dashed', 'solid', 'dashed', 'solid', 'dashed', 'solid', 'solid', 'dashed')
names(lineAnno) <- names(stList)

pdf('complexAll/our_complexAll_cutInf_seed456_ROC.pdf', height = 7, width = 9)
ggplot(data = mergedRocMat, mapping = aes(x = FPR, y = TPR, colour = Methods, linetype = Methods)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "longdash") +
  scale_linetype_manual(name = 'AUC',
                        labels = aucAnno,
                        breaks = names(stList),
                        values = lineAnno) +
  scale_colour_manual(name = 'AUC',
                      labels = aucAnno,
                      breaks = names(stList),
                      values = colAnno)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~precision and recall~~~~~~~~~~~~~~~~
prMatList <- list(Top = PR(HandTF(topMat, bin = 200000, increasing = FALSE)),
               Tree = PR(HandTF(LRMat)),
               Dollo =  PR(HandTF(DolloMat, increasing = FALSE)),
               Cor = PR(HandTF(corMat)),
               Jaccard = PR(HandTF(jacMat)),
               MI = PR(HandTF(MIMat)),
               Hamming = PR(HandTF(hamMat, increasing = FALSE)),
               NPP = PR(HandTF(corMatNPP)),
               SVD100 = PR(HandTF(euMatSVD100)),
               SVD30 = PR(HandTF(euMatSVD30)))

mergedPrMat <- do.call(rbind, prMatList)
mergedPrMat <- data.frame(Precision = mergedPrMat[, 1],
                          Recall = mergedPrMat[, 2],
                          Methods = rep(names(prMatList), sapply(prMatList, nrow)))

pdf('complexAll/our_complexAll_cutInf_seed123_PR.pdf', height = 7, width = 9)
ggplot(data = mergedPrMat, mapping = aes(x = Recall, y = Precision, colour = Methods, linetype = Methods)) +
  geom_line() +
  xlab('Recall') +
  ylab('Precision') +
   scale_linetype_manual(name = 'Methods',
                        breaks = names(prMatList),
                        values = lineAnno) +
  scale_colour_manual(name = 'Methods',
                      breaks = names(prMatList),
                      values = colAnno)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################################################################


###############################plot NPP and SVD######################
library('pROC')
library('ggplot2')
library('foreach')
require('doMC')
registerDoMC(4)

## load file
load('complexAll/simdistROCSVD30_cutInf_seed123.RData')

## set TP number
P <- sum(corMat[, 2] == 'TP')

stList <- list(Cor = corMat[1:(P*2), ],
               MI = MIMat[1:(P*2), ],
               Euclidean = euMat[1:(P*2), ])
rocList <- foreach(i = 1:length(stList)) %dopar% {
  x <- stList[[i]]
  return(roc(x[, 2], x[, 1], levels = c('TP', 'TN')))
}
rocMatList <- foreach(i = 1:length(stList)) %dopar% {
  x <- rocList[[i]]
  return(cbind(1 - x$specificities, x$sensitivities))
}

## plot ROC
mergedRocMat <- do.call(rbind, rocMatList)
mergedRocMat <- data.frame(FPR = mergedRocMat[, 1],
                           TPR = mergedRocMat[, 2],
                           Methods = rep(names(stList), sapply(rocMatList, nrow)))
aucAnno <- paste0(names(stList), ' AUC=', round(sapply(rocList, function(x){return(x$auc)}), 3))

pdf('complexAll/our_complexAll_cutInf_seed123_SVD100ROC.pdf', height = 7, width = 9)
ggplot(data = mergedRocMat, mapping = aes(x = FPR, y = TPR, colour = Methods)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'Methods',
    breaks = names(stList),
    labels = aucAnno)
dev.off()

#####################################################################
