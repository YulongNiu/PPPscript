############################ deal with MS data############################
require('stringr')
load('wholeAnno.RData')
load('mitoPhyloData.RData')

MSRaw <- read.csv('mito/MS.csv', header = FALSE)
MSRaw <- as.character(unlist(MSRaw))
## remove space
MSRaw <- sapply(MSRaw, str_trim)
## remove raw duplicated from 230 to 198
MSRaw <- MSRaw[!duplicated(MSRaw)]

MSUpper <- sapply(MSRaw, toupper)
##
MSUpperHas <- as.character(wholeAnno[wholeAnno$hsaSymbol %in% MSUpper, 6])
MSUpperHas <- unname(MSUpperHas)
MSKEGGID <- c(MSUpperHas, 'hsa:10367', 'hsa:51204', 'hsa:1632', 'hsa:51115', 'hsa:92014', 'hsa:54996', 'hsa:7416')

save(MSKEGGID, file = 'mito/MSKEGGID.RData')
##########################################################################


############################ deal with major complex########################
require('foreach')
require('doMC')
require('stringr')
registerDoMC(4)


mitoComplexFile <- dir('mito', pattern = '*.csv')
mitoComplexName <- sapply(strsplit(mitoComplexFile, split = '.', fixed = TRUE), '[[', 1)

mitoComplex <- foreach(i = 1:length(mitoComplexFile)) %dopar% {
    
    eachComplex <- read.csv(paste('mito/', mitoComplexFile[i], sep = ''))
    eachComplex <- apply(eachComplex, 1:2, function(x){
        eachEle <- as.character(x)
        eachEle <- str_trim(eachEle)
    })

    eachComplex[, 1] <- paste(rep('hsa:', nrow(eachComplex)), eachComplex[, 1], sep = '')

    return(eachComplex)
}

names(mitoComplex) <- mitoComplexName

save(mitoComplex, file = 'mitoComplex.RData')

###########################################################################


################################use mitocarta for illustration#################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heat plot data~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NEED PlotPhyloDentro() from script_phylo.R
## require('colorspace')
require('RColorBrewer')
require('foreach')
require('doMC')
registerDoMC(4)
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/mitoPhyloData.RData')
load('wholePhyloData.RData')
load('mito/mitoComplex.RData')
KEGGSpeMat <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)


# separate species by kindom
KEGGphylo <- as.character(KEGGSpeMat$Phylo)
kingdom <- foreach(i = 1:length(KEGGphylo), .combine = c) %dopar% {
  firstPhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[1]
  if (firstPhylo == 'Eukaryotes') {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
  } else {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
  }
  return(sePhylo)
}

kingdom <- c('#F8766D', '#D39200', '#93AA00', '#00C19F', '#00B9E3', '#DB72FB')[factor(kingdom)]
KEGGID <- as.character(KEGGSpeMat$KEGGID)
names(kingdom) <- KEGGID


# extract mitoproteins
mitoGenes <- as.character(mitoPhyloData[, 1])
mitophylo <- wholePhyloDataNet[rownames(wholePhyloDataNet) %in% mitoGenes, ]

## # whether in MitoCarta
## lapply(mitoComplex, function(x){
##     eachComplex <- x[!(x[, 1] %in% mitoGenes), ,drop = FALSE]
##     return(eachComplex)
## })

## # collapse complex
## mitoComplexComp <- list(complex1 = do.call(rbind, mitoComplex[1:4]),
##                         complex2 = mitoComplex$complex2,
##                         complex3 = mitoComplex$complex3,
##                         complex4 = mitoComplex$complex4,
##                         complex5 = do.call(rbind, mitoComplex[8:9]),
##                         ribo = do.call(rbind, mitoComplex[10:11]))


# color mito complex
mitoGenesVec <- rep('white', length(mitoGenes))
names(mitoGenesVec) <- mitoGenes


# set 11 complex color
## RColorBrewer::display.brewer.all()
## complexCol <- rainbow_hcl(6, c = 100, l = 65)
complex1Col <- c(brewer.pal(4,'Dark2')[3:4], brewer.pal(7, 'Set1')[7], brewer.pal(3, 'Set1')[2])
complex224Col <- c(brewer.pal(5, 'Set3')[5], brewer.pal(10, 'Paired')[10], brewer.pal(5,'Dark2')[5])
complex5Col <- c(brewer.pal(3, 'Set1')[1],  brewer.pal(4,'Dark2')[1])
complexRiboCol <- c(brewer.pal(5, 'Set1')[5], brewer.pal(6, 'Set1')[6])
complexCol <- c(complex1Col, complex224Col, complex5Col, complexRiboCol)


# complex 1 block
mitoGenesCol <- mitoGenesVec
for(i in 1:4) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

complex1Block <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
complex1Block <- complex1Block$geneBlockObj



mitoGenesCol <- mitoGenesVec
for(i in 5:5) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

complex2Block <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
complex2Block <- complex2Block$geneBlockObj


mitoGenesCol <- mitoGenesVec
for(i in 6:6) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

complex3Block <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
complex3Block <- complex3Block$geneBlockObj


mitoGenesCol <- mitoGenesVec
for(i in 7:7) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

complex4Block <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
complex4Block <- complex4Block$geneBlockObj

mitoGenesCol <- mitoGenesVec
for(i in 8:9) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

complex5Block <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
complex5Block <- complex5Block$geneBlockObj

mitoGenesCol <- mitoGenesVec
for(i in 10:11) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}

mitoRiboBlock <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)
mitoRiboBlock <- mitoRiboBlock$geneBlockObj

mitoGenesCol <- mitoGenesVec
for(i in 1:length(complexCol)) {
    eachInd <- match(mitoComplex[[i]][, 1], mitoGenes)
    # removeNA
    eachInd <- eachInd[!is.na(eachInd)]
    mitoGenesCol[eachInd] <- rep(complexCol[i], length(eachInd))
}
mitoComplexPlot <- PlotPhyloDendro(mitophylo, speCol = kingdom, geneCol = mitoGenesCol)

mitophyloPlot <- list(phyloObj = mitoComplexPlot$phyloObj,
                      speBlockObj = mitoComplexPlot$speBlockObj,
                      complex1Block = complex1Block,
                      complex2Block = complex2Block,
                      complex3Block = complex3Block,
                      complex4Block = complex4Block,
                      complex5Block = complex5Block,
                      mitoRiboBlock = mitoRiboBlock)
save(mitophyloPlot, file = 'mitophyloPlot.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mito heat plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('ggplot2')
require('grid')
load('mito/mitophyloPlot.RData')

pdf('phyloheat_mito.pdf', width = 6.5, height = 8)
grid.newpage()
print(mitophyloPlot$phyloObj, vp=viewport(0.8, 0.91, x=0.59, y=0.50))
print(mitophyloPlot$speBlockObj, vp=viewport(0.8, 0.02, x=0.59, y=0.97))
print(mitophyloPlot$complex1Block, vp=viewport(0.02, 0.91, x=0.17, y=0.50))
print(mitophyloPlot$complex2Block, vp=viewport(0.02, 0.91, x=0.14, y=0.50))
print(mitophyloPlot$complex3Block, vp=viewport(0.02, 0.91, x=0.11, y=0.50))
print(mitophyloPlot$complex4Block, vp=viewport(0.02, 0.91, x=0.08, y=0.50))
print(mitophyloPlot$complex5Block, vp=viewport(0.02, 0.91, x=0.05, y=0.50))
print(mitophyloPlot$mitoRiboBlock, vp=viewport(0.02, 0.91, x=0.02, y=0.50))
dev.off()



###############################################################################




###################mito interaction##############################
load('wholeAnno.RData')
load('mitoPhyloData.RData')
load('top400.RData')
load('mito/MSKEGGID.RData')
load('mito/mitoComplex.RData')
load('wholePhyloData.RData')


## get only mitogenes
mitoGenes <- as.character(mitoPhyloData[, 1])
top400mito <- top400[(top400[, 1] %in% mitoGenes) & (top400[, 2] %in% mitoGenes), ]

## annotation network
annoNodesMat <- apply(top400mito[, 1:2], 1:2, function(x){
  eachSymbol <- annoFirst[names(annoFirst) %in% x]
  return(eachSymbol)
})
top400mitoMat <- cbind(annoNodesMat, top400mito[, 3])
colnames(top400mitoMat)[3] <- 'jaccard'
write.table(top400mitoMat, 'mito/cytonetwork/top400mito.txt', quote = FALSE, row.names = FALSE)


## network node anno
mitoGenesComplex <- rep('other', length(mitoGenes))
mitoComplexMat <- do.call(rbind, mitoComplex)
for(i in 1:length(mitoGenesComplex)) {
    hasInx <- match(mitoGenes[i], mitoComplexMat[, 1])
    if(!is.na(hasInx)) {
        mitoGenesComplex[i] <- mitoComplexMat[hasInx, 3]
    } else {}
}

mitoGenesMS <- rep('N', length(mitoGenes))
for(i in 1:length(mitoGenesMS)) {
    hasInx <- match(mitoGenes[i], MSKEGGID)
    if(!is.na(hasInx)) {
        mitoGenesMS[i] <- 'Y'
    } else {}
}

mitoGenesAnno <- sapply(mitoGenes, function(x) {annoFirst[names(annoFirst) %in% x]})
top400mitoNode <- cbind(mitoGenesAnno, mitoGenesComplex, mitoGenesMS)
write.table(top400mitoNode, 'mito/cytonetwork/top400mitoNode.txt', quote = FALSE, row.names = FALSE, sep = '\t')


#################################################################
