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
load('complexAll/simdistROCNPP70_cut40_seed123.RData')
corMatNPP <- corMat
corRocNPP <- corRoc
load('complexAll/simdistROCSVD100_cut40_seed123.RData')
euMatSVD100 <- euMat
euRocSVD100 <- euRoc
load('complexAll/simdistROCSVD30_cut40_seed123.RData')
euMatSVD30 <- euMat
euRocSVD30 <- euRoc

pat <- 'ROC_cut40_seed123'
rocDataFiles <- dir('complexAll', pattern = pat, full.names = TRUE)
for(i in rocDataFiles) {load(i)}

## set P and N number
P <- sum(topMat[, 2] == 'TP')
N <- sum(topMat[, 2] == 'TN')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

## plot ROC
mergedRocMat <- do.call(rbind, rocMatList)
mergedRocMat <- data.frame(FPR = mergedRocMat[, 1],
                           TPR = mergedRocMat[, 2],
                           Methods = rep(names(stList), sapply(rocMatList, nrow)))
aucAnno <- paste0(names(stList), '=', round(sapply(rocList, function(x){return(x$auc)}), 3))

pdf('complexAll/our_complexAll_cut40_seed123_ROC.pdf', height = 7, width = 9)
ggplot(data = mergedRocMat, mapping = aes(x = FPR, y = TPR, colour = Methods)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'AUC',
    breaks = names(stList),
    labels = aucAnno)
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

pdf('complexAll/our_complexAll_cut40_seed123_PR.pdf', height = 7, width = 9)
ggplot(data = mergedPrMat, mapping = aes(x = Recall, y = Precision, colour = Methods)) +
  geom_line() +
  xlab('Recall') +
  ylab('Precision')
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
load('complexAll/simdistROCSVD30_cut40_seed123.RData')

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

pdf('complexAll/our_complexAll_cut40_seed123_SVD100ROC.pdf', height = 7, width = 9)
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
############################### plot ROC with other method #############

## AUC 0.5580
load('complex/Bioinfor_mutualROCMat.RData')
## AUC 0.6064
load('complex/Bioinfor_jacsimROCMat.RData')
## AUC 0.5416
load('complex/Bioinfor_corROCMat.RData')
## AUC 0.5349
load('complex/Bioinfor_bayesROCMat.RData')
## AUC 0.6389
load('complex/Bioinfor_topROCMat.RData')
## ## load('complex/Bioinfor_ROCMat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~plot PPV and RPP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ROC2PPV <- function(rocMat, Pnum, Nnum) {
  ## USE: conver ROC matrix to PPV (positive predictive value) matrix.
  ## INPUT: 'rocMat' is a matrix, of which the first column is 'TPR', and the second column is 'FPR'. The rownames of 'rocMat' is the threshold. 'Pnum' and 'Nnum' is the total positive and negative number.
  ## OUTPUT: 'ppvMat', and the PPV = TP/(TP+FP)

  TPFPMat <- cbind(rocMat[, 1] * Pnum, rocMat[, 2] * Nnum)
  ppv <- TPFPMat[, 1] / rowSums(TPFPMat)
  ppvMat <- data.frame(Cutoff = as.numeric(rownames(rocMat)),
                       PPV = ppv)
  return(ppvMat)
}


ROC2RPPMat <- function(rocMat, Pnum, Nnum) {
  ## USE: conver ROC matrix to RPP (rate of positive predictions) and PPV (positive predictive value) matrix.
  ## INPUT: 'rocMat' is a matrix, of which the first column is 'TPR', and the second column is 'FPR'. The rownames of 'rocMat' is the threshold. 'Pnum' and 'Nnum' is the total positive and negative number.
  ## OUTPUT: 'rppMat', and the PPV = TP/(TP+FP), RPP = (TP+FP)/(TP+FP+TN+FN)

  ppv <- ROC2PPV(rocMat, Pnum, Nnum)[, 2]

  TPFPMat <- cbind(rocMat[, 1] * Pnum, rocMat[, 2] * Nnum)
  rpp <- rowSums(TPFPMat) / (Pnum + Nnum)

  ppvrppMat <- data.frame(PPV = ppv, RPP = rpp)

  return(ppvrppMat)
}

jacsimRppMat <- ROC2RPPMat(jacsimROCMat, 57114, 57114)
corRppMat <- ROC2RPPMat(corROCMat, 57114, 57114)
mutualRppMat <- ROC2RPPMat(mutualROCMat, 57114, 57114)
topRppMat <- ROC2RPPMat(topROCMat, 57114, 57114)
bayesRppMat <- ROC2RPPMat(bayesROCMat, 57114, 57114)

mergedRppMat <- rbind(jacsimRppMat, corRppMat, mutualRppMat, bayesRppMat, topRppMat)
mergedRppMat <- data.frame(PPV = mergedRppMat[, 1],
                           RPP = mergedRppMat[, 2],
                           Method = rep(c('Jaccard', 'Cor', 'Mutual', 'Tree', 'Top'), c(nrow(jacsimRppMat), nrow(corRppMat), nrow(mutualRppMat), nrow(bayesRppMat), nrow(topRppMat))))


pdf('complexAll/our_complexAll_rpp.pdf', height = 7, width = 8)
ggplot(data = mergedRppMat, mapping = aes(x = RPP, y = PPV, colour = Method)) +
  geom_line() +
  xlab('Rate of positive predictions') +
  ylab('Positive predictive value') +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') +
  scale_y_continuous(limits=c(0, 1))
dev.off()
 

jacsimPpvMat <- ROC2PPV(jacsimROCMat, 57114, 57114)
pdf('complexAll/our_complexAll_ppv_jaccard.pdf')
ggplot(data = jacsimPpvMat, mapping = aes(x = Cutoff, y = PPV)) +
  geom_line(color = '#F8766D') +
  xlab('Jaccard cut-off') +
  ylab('Positive predictive value') +
  scale_y_continuous(limits=c(0, 1)) +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') 
dev.off()

corPpvMat <- ROC2PPV(corROCMat, 57114, 57114)
pdf('complexAll/our_complexAll_ppv_cor.pdf')
ggplot(data = corPpvMat, mapping = aes(x = Cutoff, y = PPV)) +
  geom_line(color = '#D39200') +
  xlab('Cor cut-off') +
  ylab('Positive predictive value') +
  scale_y_continuous(limits=c(0, 1)) +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') 
dev.off()

mutualPpvMat <- ROC2PPV(mutualROCMat, 57114, 57114)
pdf('complexAll/our_complexAll_ppv_mutual.pdf')
ggplot(data = mutualPpvMat, mapping = aes(x = Cutoff, y = PPV)) +
  geom_line(color = '#00C19F') +
  xlab('MI cut-off') +
  ylab('Positive predictive value') +
  scale_y_continuous(limits=c(0, 1)) +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') 
dev.off()

topPpvMat <- ROC2PPV(topROCMat, 57114, 57114)
topPpvMat[, 1] <- 20127 - topPpvMat[, 1]
pdf('complexAll/our_complexAll_ppv_top.pdf')
ggplot(data = topPpvMat, mapping = aes(x = Cutoff, y = PPV)) +
  geom_line(color = '#619CFF') +
  xlab('Top cut-off') +
  ylab('Positive predictive value') +
  scale_y_continuous(limits=c(0, 1)) +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') +
  scale_x_reverse() 
dev.off()

bayesPpvMat <- ROC2PPV(bayesROCMat, 57114, 57114)
pdf('complexAll/our_complexAll_ppv_bayes.pdf')
ggplot(data = bayesPpvMat, mapping = aes(x = Cutoff, y = PPV)) +
  geom_line(color = '#DB72FB') +
  xlab('LR cut-off') +
  ylab('Positive predictive value') +
  scale_y_continuous(limits=c(0, 1)) +
  geom_hline(yintercept = 0.5, colour='grey', linetype = 'dashed') 
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################################
