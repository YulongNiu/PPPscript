#######################test WGCNA####################
#~~~~~~~~~~~~~~~~~~~~~ choose beta~~~~~~~~~~~~~~~~~~~~~~
library('WGCNA')
load('wholePhyloData.RData')
allowWGCNAThreads()

wholePhyloDataNet <- t(wholePhyloDataNet)
gc()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(wholePhyloDataNet, powerVector = powers, verbose = 5)
save(sft, file = 'sft.RData')

##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
## 1      1 6.72e-01  1.0600          0.824    9090   10400.0  11700
## 2      2 2.05e-01  0.3320          0.837    5330    6390.0   8070
## 3      3 4.87e-02  0.1180          0.805    3550    4220.0   6070
## 4      4 6.52e-05  0.0036          0.800    2560    2910.0   4840
## 5      5 9.84e-02 -0.1130          0.831    1960    2080.0   4030
## 6      6 3.44e-01 -0.2010          0.874    1560    1530.0   3460
## 7      7 6.85e-01 -0.2960          0.906    1280    1150.0   3120
## 8      8 8.27e-01 -0.3780          0.918    1080     883.0   2860
## 9      9 8.76e-01 -0.4500          0.906     924     686.0   2640
## 10    10 9.10e-01 -0.5060          0.930     803     537.0   2460
## 11    12 9.28e-01 -0.5950          0.943     628     342.0   2170
## 12    14 9.37e-01 -0.6550          0.946     509     223.0   1940
## 13    16 9.61e-01 -0.6910          0.964     423     152.0   1760
## 14    18 9.73e-01 -0.7180          0.976     358     105.0   1610
## 15    20 9.82e-01 -0.7410          0.985     308      74.8   1490
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('WGCNA')
load('wholePhyloData.RData')
allowWGCNAThreads()

wholePhyloDataNet <- t(wholePhyloDataNet)
gc()


softPower = 10;
adjacency = adjacency(wholePhyloDataNet, power = softPower);
save(adjacency, file = 'adjacency.RData')

load('adjacency.RData')
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################
