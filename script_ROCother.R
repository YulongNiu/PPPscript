setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')

##############################test all complex#######################
library('pROC')
library('ggplot2')

load('complexAll/allRS_cut40_seed123.RData')
load('complexAll/simdistROC_cut40_seed123.RData')
load('complexAll/DolloROC_cut40_seed123.RData')
load('complexAll/LRROC_cut40_seed123.RData')
load('complexAll/topROC_cut40_seed123.RData')


N <- sum(allRS[, 3] == 'TP')
## AUC = 0.721
topRoc <- roc(status ~ distTop, topMat[1:N*2, ], levels = c('TP', 'TN'))
topRocMat <- cbind(1 - topRoc$specificities, topRoc$sensitivities)

## AUC = 0.609
LRRoc <- roc(status ~ simLR, LRMat[1:(N*2), ], levels = c('TP', 'TN'))
LRRocMat <- cbind(1 - LRRoc$specificities, LRRoc$sensitivities)

## AUC = 0.623
DolloRoc <- roc(status ~ distDollo, DolloMat[1:(N*2), ], levels = c('TP', 'TN'))
DolloRocMat <- cbind(1 - DolloRoc$specificities, DolloRoc$sensitivities)

## AUC = 0.702
corRoc <- roc(status ~ simcor, corMat[1:(N*2), ], levels = c('TP', 'TN'))
corRocMat <- cbind(1 - corRoc$specificities, corRoc$sensitivities)

## AUC = 0.716
jacRoc <- roc(status ~ simjac, jacMat[1:(N*2), ], levels = c('TP', 'TN'))
jacRocMat <- cbind(1 - jacRoc$specificities, jacRoc$sensitivities)

## AUC = 0.710
MIRoc <- roc(status ~ simMI, MIMat[1:(N*2), ], levels = c('TP', 'TN'))
MIRocMat <- cbind(1 - MIRoc$specificities, MIRoc$sensitivities)

## AUC = 0.666
hamRoc <- roc(status ~ distham, hamMat[1:(N*2), ], levels = c('TP', 'TN'))
hamRocMat <- cbind(1 - hamRoc$specificities, hamRoc$sensitivities)

## plot ROC
mergedRocMat <- rbind(topRocMat,
                      LRRocMat,
                      DolloRocMat,
                      corRocMat,
                      jacRocMat,
                      MIRocMat,
                      hamRocMat)
mergedRocMat <- data.frame(FPR = mergedRocMat[, 1],
                           TPR = mergedRocMat[, 2],
                           Method = rep(c('Top', 'Tree', 'Dollo', 'Cor', 'Jaccard', 'MI', 'Hamming'), c(nrow(topRocMat), nrow(LRRocMat), nrow(DolloRocMat), nrow(corRocMat), nrow(jacRocMat), nrow(MIRocMat), nrow(hamRocMat))))

pdf('complexAll/our_complexAll_ROC.pdf', height = 7, width = 9)
ggplot(data = mergedRocMat, mapping = aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'Methods',
    breaks = c('Top', 'Tree', 'Dollo', 'Cor', 'Jaccard', 'MI', 'Hamming'),
    labels = c('Top AUC=0.721',
               'Tree AUC=0.609',
               'Dollo AUC=0.623',
               'Cor AUC=0.702',
               'Jaccard AUC=0.716',
               'MI AUC=0.710',
               'Hamming AUC=0.666'))
dev.off()

#####################################################################

############################### mutual test ###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mutual information~~~~~~~~~~~~~
require('ggplot2')
load('complexAll/our_complex.RData')
## load('complexAll/our_complex_567.RData')
## load('complex/Bioinfor_complex.RData')
load('wholePhyloData.RData')

PRSJacCor <- apply(PRSJacCor, 1:2, as.character)
NRSJacCor <- apply(NRSJacCor, 1:2, as.character)
PRSmutual <- GetSimiMat(PRSJacCor[, 1:2], wholePhyloDataNet)
NRSmutual <- GetSimiMat(NRSJacCor[, 1:2], wholePhyloDataNet)

mutualVec <- c(PRSmutual[, 'simi'], NRSmutual[, 'simi'])
mutualsim <- data.frame(Mutual = as.numeric(mutualVec), Status = rep(c('TP', 'TN'), c(nrow(PRSmutual), nrow(NRSmutual))))

p <- ggplot(data = mutualsim, aes(x = Mutual))
p +
  geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = factor(Status))) +
  stat_density(geom = 'line', position = 'identity', aes(colour = factor(Status)))

save(PRSmutual, NRSmutual, file = 'complexAll/our_complex_mutual.RData')
## save(PRSmutual, NRSmutual, file = 'complexAll/our_complex_mutual_567.RData')
## save(PRSmutual, NRSmutual, file = 'complex/Bioinfor_complex_mutua.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test mutual ROC~~~~~~~~~~~~~~~~~~~~
load('complexAll/our_complex_mutual.RData')
## load('complexAll/our_complex_mutual_567.RData')
## load('complex/Bioinfor_complex_mutual.RData')
mutualVec <- as.numeric(c(PRSmutual[, 'simi'], NRSmutual[, 'simi']))
mutualMat <- data.frame(Mutual = as.numeric(mutualVec), Status = rep(c('pos', 'neg'), c(nrow(PRSmutual), nrow(NRSmutual))))

cutVec <- seq(min(mutualVec), max(mutualVec), 0.01)
mutualROCMat <- SingROCMat(cutVec, mutualMat)
save(mutualROCMat, file = 'complexAll/mutualROCMat.RData')
## save(mutualROCMat, file = 'complexAll/mutualROCMat_567.RData')
## save(mutualROCMat, file = 'complex/Bioinfor_mutualROCMat.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################################


############################### plot ROC with other method ################
library('ggplot2')
## AUC 0.6759
load('complexAll/mutualROCMat.RData')
## AUC 0.7085
load('complexAll/jacsimROCMat.RData')
## AUC 0.6635
load('complexAll/corROCMat.RData')
## AUC 0.5947
load('complexAll/bayesROCMat.RData')
## AUC 0.7332
load('complexAll/topROCMat.RData')
## load('complexAll/ROCMat.RData')


## AUC 0.6747
load('complexAll/mutualROCMat_567.RData')
## AUC 0.7066
load('complexAll/jacsimROCMat_567.RData')
## AUC 0.6623
load('complexAll/corROCMat_567.RData')
## AUC 0.5906
load('complexAll/bayesROCMat_567.RData')
## AUC 0.7322
load('complexAll/topROCMat_567.RData')
## ## load('complexAll/ROCMat_567.RData')

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~ROC plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## topROCMat <- ROCMat / 57114
## topROCMat <- ROCMat / 26525
## colnames(topROCMat) <- c('TPR', 'FPR')
mergedROCMat <- rbind(jacsimROCMat, corROCMat, mutualROCMat, bayesROCMat, topROCMat)
mergedROCMat <- cbind(mergedROCMat, Method = rep(c('Jaccard', 'Cor', 'Mutual', 'Tree', 'Top'), c(nrow(jacsimROCMat), nrow(corROCMat), nrow(mutualROCMat), nrow(bayesROCMat), nrow(topROCMat))))


pdf('complexAll/our_complexAll_ROC.pdf', height = 7, width = 9)
ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'Method',
    breaks = c('Cor', 'Jaccard', 'Mutual', 'Top', 'Tree'),
    labels = c('Cor AUC=0.6635',
                        'Jaccard AUC=0.7085',
                        'Mutual AUC=0.6759',
                        'Top AUC=0.7332',
                        'Tree AUC=0.5947'))
dev.off()

pdf('complexAll/our_complexAll_ROC_567.pdf', height = 7, width = 9)
ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'Method',
    breaks = c('Cor', 'Jaccard', 'Mutual', 'Top', 'Tree'),
    labels = c('Cor AUC=0.6623',
                        'Jaccard AUC=0.7066',
                        'Mutual AUC=0.6747',
                        'Top AUC=0.7322',
                        'Tree AUC=0.5906'))
dev.off()

pdf('complex/Bioinfor_complex_ROC.pdf', height = 7, width = 9)
ggplot(data = mergedROCMat, mapping = aes(x = FPR, y = TPR, colour = Method)) +
  geom_line() +
  xlab('False positive rate') +
  ylab('True positive rate') +
  geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "dashed") +
  scale_color_discrete(
    name = 'Method',
    breaks = c('Cor', 'Jaccard', 'Mutual', 'Top', 'Tree'),
    labels = c('Cor AUC=0.5416',
                        'Jaccard AUC=0.6064',
                        'Mutual AUC=0.5580',
                        'Top AUC=0.6389',
                        'Tree AUC=0.5349'))
dev.off()
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
