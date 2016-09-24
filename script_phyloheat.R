################################Useful functions##############################
PlotPhyloDendro <- function(phyloData, geneNameSize = 3, geneNameCol = 'grey55', geneBlockCol = NA, presentCol = 'steelblue', absentCol = 'grey91', speCol, geneCol) {
  ## USE: plot the dendrogram of phylogenetic data
  ## INPUT: 'phyloData' is the phylogenetic data, of which the row is gene and column is species.
  ## 'geneNameSize' is the size of gene names label.
  ## 'geneNameCol' is the colour of gene names.
  ## 'geneBlockCol' is the space color between gene blocks, and the default value is NA meaning no space color. If the gene number is samll, less than 20, set it to 'white' is fine.
  ## 'presentCol' is the color of present 1.
  ## 'absentCol' is the color of absent 0.
  ## 'speCol' is a vector of colors with names of species, which are the same as colnames of 'phyloData' (may not in the same order).
  ## 'geneCol' is a vector of colors with names of species, which are the same as rownames of 'phyloData' (may not in the same order).
  ## OUTPUT: A List of ggplot2 object.
  ## 'phyloObj' is phylogenetic profile data.
  ## 'geneDendroObj' is the dendrogram of genes.
  ## 'speBlockObj' is a color block of species.
  ## 'geneBlockObj' is the color block of species.
  ## 'geneNamesObj' is the gene names.

  require('ggplot2')
  require('reshape')
  require('ggdendro')
  require('grid')

  # cluster genes and species
  hcGene <- hclust(dist(phyloData), method = 'average')
  rowInd <- hcGene$order
  hcSpe <- hclust(dist(t(phyloData)), method = 'average')
  colInd <- hcSpe$order
  
  # order 'phyloData'
  orderedPhyloData <- phyloData[rowInd,colInd]
  orderedColNames <- colnames(orderedPhyloData)
  orderedRowNames <- rownames(orderedPhyloData)
  breaksRow <- 1:length(orderedRowNames)

  # order 'geneCol'
  orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]
  
  # order 'speCol'
  orderedSpeCol <- speCol[match(colnames(orderedPhyloData), names(speCol))]
    
  # melt data for ggplot2
  colnames(orderedPhyloData) <- 1:ncol(orderedPhyloData)
  rownames(orderedPhyloData) <- 1:nrow(orderedPhyloData)
  orderedPhyloData <- melt(orderedPhyloData)
  colnames(orderedPhyloData) <- c('geneNames', 'speNames', 'apData')


  # plot gene names
  orderedRowNamesMat <- data.frame(x = rep(0, length(orderedRowNames)),
                        y = seq(0.5, length(orderedRowNames)-0.5, 1),
                        fillName = orderedRowNames)
  

  geneNamesObj <- ggplot(orderedRowNamesMat, aes(x, y, label = fillName)) +
    geom_text(size = geneNameSize, colour = geneNameCol) +
       labs(x = NULL, y = NULL) +
      scale_x_continuous(expand = c(0, 0), breaks = NULL) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
           theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  # plot phylogenetic matrix
  phyloObj <- ggplot(orderedPhyloData, aes(speNames, geneNames)) +
    geom_tile(aes(fill = factor(apData))) +
      scale_fill_manual(name = 'status', labels = c('absent', 'present'), values = c(absentCol, presentCol)) +
        labs(x = NULL, y = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            scale_x_continuous(expand = c(0, 0), breaks = NULL) +
              theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  # dendrogram plot for genes
  ddata <- dendro_data(hcGene, type = "rectangle")
  segData <- segment(ddata)
  segData[, c(1, 3)] <- segData[, c(1, 3)] - 0.5
  geneDendroObj <- ggplot(segData) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      labs(x = NULL, y = NULL) +
      scale_y_reverse(expand = c(0, 0), breaks = NULL) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, nrow(phyloData)), breaks = NULL) +
          coord_flip() +
            theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  # gene color block
  orderedGeneColMat <- data.frame(x = rep(0, length(orderedGeneCol)),
                        y = 1:length(orderedGeneCol),
                        fillCol = factor(orderedGeneCol))

  geneBlockObj <- ggplot(orderedGeneColMat, aes(x, y)) +
    geom_tile(aes(fill = fillCol), color = geneBlockCol) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedGeneColMat$fillCol)) +
              theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  # species color block
  orderedSpeColMat <- data.frame(y = rep(0, length(orderedSpeCol)),
                                 x = 1:length(orderedSpeCol),
                                 fillCol = factor(orderedSpeCol))

  speBlockObj <- ggplot(orderedSpeColMat, aes(x, y)) +
    geom_tile(aes(fill = fillCol)) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedSpeColMat$fillCol)) +
              theme(legend.position='none',
                    title = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.length = unit(0, "mm"),
                    axis.ticks.margin = unit(0, "mm"),
                    axis.line = element_blank(),
                    panel.margin = unit(0, 'mm'),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), 'line'),
                    legend.margin = unit(0, 'mm'))

  # return results
  plotRes <- list(phyloObj = phyloObj,
                  geneDendroObj = geneDendroObj,
                  speBlockObj = speBlockObj,
                  geneBlockObj = geneBlockObj,
                  geneNamesObj = geneNamesObj)

  return(plotRes)
}
##############################################################################




###############################plot phyloheat########################
#~~~~~~~~~~~~~~~~~~~~~~~~plot ATPsynthase~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('foreach')
library('doMC')
registerDoMC(4)
load('wholePhyloData.RData')
atpSub <- read.table('/home/Yulong/RESEARCH/neuro/Bioinfor/phylogenetic_profile_old/F1F0Gene.txt')
KEGGSpeMat <- read.csv('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholeListFile.csv', row.names = 1)

# separate species by domain
KEGGphylo <- as.character(KEGGSpeMat$Phylo)
domain <- foreach(i = 1:length(KEGGphylo), .combine = c) %dopar% {
  firstPhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[1]
  if (firstPhylo == 'Eukaryotes') {
    return(firstPhylo)
  } else {
    sePhylo <- unlist(strsplit(KEGGphylo[i], split = ';', fixed = TRUE))[2]
    return(sePhylo)
  }
}
domain <- c('#D39200', '#93AA00', '#00C19F')[factor(domain)]
KEGGID <- as.character(KEGGSpeMat$KEGGID)
names(domain) <- KEGGID

# 17 ATPsubunits phyloData
ATPphylo <- wholePhyloDataNet[rownames(wholePhyloDataNet) %in% atpSub[, 2], ]

# sort atp enzemble names according to ATPphylo
atpSubUni <- cbind(atpSub, rep(c('F1', 'Fo'), c(6, 11)))
atpSubUni <- atpSubUni[order(atpSubUni[, 2]), ]
atpSubUni <- atpSubUni[rank(rownames(ATPphylo)), ]
rownames(ATPphylo) <- atpSubUni[, 1]

atpSubVec <- c('#F8766D', '#00BA38')[factor(atpSubUni[, 3])]
names(atpSubVec) <- atpSubUni[, 1]


ATPphyloPlot <- PlotPhyloDendro(ATPphylo, geneBlockCol = 'white', speCol = domain, geneCol = atpSubVec)

pdf('phyloheat_17ATP.pdf', width = 7, height = 2.7)
grid.newpage()
print(ATPphyloPlot$phyloObj, vp=viewport(0.8, 0.8, x=0.58, y=0.45))
print(ATPphyloPlot$speBlockObj, vp=viewport(0.8, 0.05, x=0.58, y=0.88))
print(ATPphyloPlot$geneBlockObj, vp=viewport(0.03, 0.8, x=0.16, y=0.45))
print(ATPphyloPlot$geneNamesObj, vp=viewport(0.1, 0.8, x=0.11, y=0.45))
print(ATPphyloPlot$geneDendroObj, vp=viewport(0.05, 0.8, x=0.05, y=0.45))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################################








