########################### pathway function#######################
## path2matInter <- function(path, uniProtMat, isUniProt = FALSE) {
##   # USE: transfer pathway data to from-to interaction matrix.
##   # INPUT: 'path' is the pathway object from 'graphite' databases, which namely are 'reactome', 'kegg', 'biocarta', 'nci'. 'isUniProt' is whether the uniprot names. 'uniProtMat' is the KEGGID and uniProt matrix. 
##   # OUTPUT: list of mat interaction.

##   require('graphite')
  
##   pathList <- lapply(path, function(x) {
##     # to graphNEL object
##     g <- pathwayGraph(x)
##     return(g)
##   })
##   names(pathList) <- names(path)

##   # remove pathways only has 1 node
##   nodeNum <- sapply(pathList, function(x){return(length(nodes(x)))})
##   pathList <- pathList[which(nodeNum > 1)]
  
##   pathListMat <- lapply(pathList, function(x) {
    
##     # deal with a 2 nodes/1 edge bug
##     if (length(nodes(x)) == 2 & sum(sapply(edges(x), length)) == 1) {
##       ftMat <- matrix(c(nodes(x), '1'), nrow = 1)
##       colnames(ftMat) <- c('from', 'to', 'weight')
##     } else {
##       ftMat <- extractFromTo(as(x, 'graphBAM'))
##     }

##     ftMat <- ftMat[, 1:2, drop = FALSE]

##     if (!isUniProt) {
##       # KEGG IDs
##       ftMat <- apply(ftMat, 1:2, sub, pattern = 'EntrezGene', replacement = 'hsa')
##     } else {
##       ftMat <- apply(ftMat, 1:2, function(i) {
##         eachNode <- sub(pattern = 'UniProt', replacement = 'up', i)
##         eachNode <- uniProtMat[match(eachNode, uniProtMat[, 2]), 1]
##         return(eachNode)
##       })
##        # remove NA
##       fromIsNa <- is.na(ftMat[, 1])
##       toIsNa <- is.na(ftMat[, 2])
##       ftMat <- ftMat[!(fromIsNa | toIsNa), ,drop = FALSE]
##     }

##     # deal with the no nodes matrix
##     if (nrow(ftMat) == 0) {
##       ftMat <- NULL
##     } else {
##       # directed nework to undirected
##       fromNum <- sapply(strsplit(ftMat[, 1], split = ':', fixed = TRUE), '[[', 2)
##       toNum <- sapply(strsplit(ftMat[, 2], split = ':', fixed = TRUE), '[[', 2)
##       isSmall <- fromNum <= toNum
##       for (i in 1:length(isSmall)) {
##         if (!isSmall[i]) {
##           ftMat[i, 1:2] <- ftMat[i, 2:1]
##         } else {}
##       }
##       ftMat <- ftMat[!duplicated(ftMat, MARGIN = 1), , drop = FALSE]

##       # remove self directed network
##       ftMat <- ftMat[!(ftMat[, 1] == ftMat[, 2]), , drop = FALSE]
##     }

##     return(ftMat)
##   })


##   # remove NULL
##   pathListMat <- pathListMat[which(!sapply(pathListMat, is.null))]

##   return(pathListMat)
## }

path2matInter <- function(path, uniProtMat, isUniProt = FALSE, ncor = 4) {
  # USE: transfer pathway data to from-to interaction matrix.
  # INPUT: 'path' is the pathway object from 'graphite' databases, which namely are 'reactome', 'kegg', 'biocarta', 'nci'. 'isUniProt' is whether the uniprot names. 'uniProtMat' is the KEGGID and uniProt matrix. 'ncor' is the threads number.
  # OUTPUT: list of mat interaction.

  require('graph')
  require('graphite')
  require('foreach')
  require('doMC')
  registerDoMC(ncor)

  pathList <- foreach(i = 1:length(path)) %dopar% {
    print(paste('It is running ', i, ' in total of ', length(path), '.', sep = ''))
    g <- pathwayGraph(path[[i]])

    if (length(nodes(g)) == 1) {
      ## give NAs to pathways only has 1 node
      ftMatRes <- NA
    } else {
      if (length(nodes(g)) == 2 & sum(sapply(edges(g), length)) == 1) {
        ftMat <- matrix(c(nodes(g), '1'), nrow = 1)
        colnames(ftMat) <- c('from', 'to', 'weight')
      }
      else {
        ftMat <- extractFromTo(as(g, 'graphBAM'))
      }

      ftMat <- ftMat[, 1:2, drop = FALSE]

      if (!isUniProt) {
        ## KEGG IDs
        ftMat <- apply(ftMat, 1:2, sub, pattern = 'EntrezGene', replacement = 'hsa')
      } else {
        ## uniprot IDs
        ftMat <- apply(ftMat, 1:2, function(i) {
          eachNode <- sub(pattern = 'UniProt', replacement = 'up', i)
          eachNode <- uniProtMat[match(eachNode, uniProtMat[, 2]), 1]
          return(eachNode)
        })
        ## remove NA
        fromIsNa <- is.na(ftMat[, 1])
        toIsNa <- is.na(ftMat[, 2])
        ftMat <- ftMat[!(fromIsNa | toIsNa), ,drop = FALSE]
      }

      ## deal with the no nodes matrix
      if (nrow(ftMat) == 0) {
        ftMatRes <- NA
      } else {
        ## directed nework to undirected
        fromNum <- sapply(strsplit(ftMat[, 1], split = ':', fixed = TRUE), '[[', 2)
        toNum <- sapply(strsplit(ftMat[, 2], split = ':', fixed = TRUE), '[[', 2)
        isSmall <- fromNum <= toNum
        ftMat <- rbind(ftMat[isSmall, ,drop = FALSE], ftMat[!isSmall, 2:1, drop = FALSE])
        ftMat <- ftMat[!duplicated(ftMat, MARGIN = 1), , drop = FALSE]

        ## remove self directed network
        ftMatRes <- ftMat[!(ftMat[, 1] == ftMat[, 2]), , drop = FALSE]
      }
      
    }
    return(ftMatRes)
  }

  # give names
  names(pathList) <- names(path)
  
  # remove NA
  isMatNA <- sapply(pathList, function(x) {
    return(is.null(dim(x)))
  })
  pathListMat <- pathList[!isMatNA]

  return(pathListMat)
}

FilterPath <- function(path, geneSet, ncor = 4) {
  # USE: filter path and only select genes in 'geneSet'.
  # INPUT: 'path' is pathway object, which is a list containing matrix. 'geneSet' is the whole gene list. 'ncor' is the threads number. 
  # OUTPUT: filtered path object

  require('foreach')
  require('doMC')
  registerDoMC(ncor)

  filetedPath <- foreach (i = 1:length(path)) %dopar% {
    fromNodes <- path[[i]][, 1]
    toNodes <- path[[i]][, 2]

    fromHas <- fromNodes %in% geneSet
    toHas <- toNodes %in% geneSet
    wholeHas <- fromHas & toHas

    if (sum(wholeHas) > 0) {
      eachPath <- path[[i]][wholeHas, , drop = FALSE]
    } else {
      eachPath <- NA
    }

    return(eachPath)
  }

  names(filetedPath) <- names(path)

  # remove NA
  isMatNA <- sapply(filetedPath, function(x) {
    return(is.null(dim(x)))
  })
  filetedPath <- filetedPath[!isMatNA]

  return(filetedPath)
  
}

PredSummary <- function(bioSet, wholeInter, ncor = 4) {
  # USE: summary the prediction results.
  # INPUT: 'bioSet' is the biological geneset, which is a matrix with two columns that are two nodes. 'wholeInter' is the whole interaction, which is a vector. The two nodes are combined with '|'. 'ncor' is the threads number. 
  # OUTPUT: sorted prediction matrix.

  require('foreach')
  require('doMC')
  registerDoMC(ncor)

  # number of prediction
  predNum <- foreach (i = 1:length(bioSet), .combine = c) %dopar% {
    print(paste('It is running ', i, ' in total of ', length(bioSet), '.', sep = ''))
    eachPath <- apply(bioSet[[i]], 1, paste, collapse = '|')
    eachPathNum <- sum(eachPath %in% wholeInter)
    
    return(eachPathNum)
  }

  # number of original
  orgNum <- sapply(bioSet, nrow)

  # percentage
  percent <- predNum/orgNum

  predMat <- cbind(predNum, orgNum, percent)
  colnames(predMat) <- c('predNum', 'orgNum', 'percentage')

  # orderd by precentage
  predMat <- predMat[order(percent, decreasing = TRUE), ]

  return(predMat)
}

OrderHumMat <- function(humMat) {
  # USE: order the human matrix (pathways, genesets) by the first and second rows
  # INPUT: 'humMat' is the human matrix. The first and second rows should be KEGGIDs, like 'hsa:1000'.
  # OUTPUT: ordered matrix

  fromNum <- sapply(strsplit(humMat[, 1], split = ':', fixed = TRUE), '[[', 2)
  toNum <- sapply(strsplit(humMat[, 2], split = ':', fixed = TRUE), '[[', 2)
  isSmall <- fromNum <= toNum
  humMat <- rbind(humMat[isSmall, ,drop = FALSE], humMat[!isSmall, 2:1, drop = FALSE])

  return(humMat)
}

###################################################################


############################## preprocess pathway data to interaction##############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require('graphite')
load('uniProt.RData')
keggPath <- path2matInter(kegg, uniProt)
biocartaPath <- path2matInter(biocarta, uniProt)
nciPath <- path2matInter(nci, uniProt, TRUE)

require('KEGGBioCycAPI')
cutMat <- CutSeqEqu(length(reactome), 100)
for (i in 1:ncol(cutMat)) {
  reactomePath <- path2matInter(reactome[cutMat[1, i]:cutMat[2, i]], uniProt, TRUE)
  fileName <- paste('reactomePath', cutMat[1, i], '_', cutMat[2, i], '.RData', sep = '')
  save(reactomePath, file = fileName)
}

fileName <- dir(pattern = '^reactomePath\\d')
reactomePathAll <- foreach(i = 1:length(fileName), .combine = append) %dopar% {
  load(fileName[i])
  return(reactomePath)
}

reactomePath <- reactomePathAll
save(keggPath, biocartaPath, nciPath, reactomePath, file = 'path.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~select nodes only exist in human genome~~~~~~~~~~~~~~
require('foreach')
require('doMC')
registerDoMC(4)
load('pathway/path.RData')
load('wholePhyloData.RData')

keggPathFilter <- FilterPath(keggPath, names(annoFirst))
biocartaPathFilter <- FilterPath(biocartaPath, names(annoFirst))
nciPathFilter <- FilterPath(nciPath, names(annoFirst))
reactomePathFilter <- FilterPath(reactomePath, names(annoFirst))

save(keggPathFilter, biocartaPathFilter, nciPathFilter, reactomePathFilter, file = 'pathFilter.RData')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################


####################### map predicted interaction to pathways ##################
load('pathway/pathFilter.RData')
load('top400.RData')
wholeInter <- apply(top400[, 1:2], 1, paste, collapse = '|')

keggPred <- PredSummary(keggPathFilter, wholeInter)
write.csv(keggPred, 'keggPred.csv')
biocartaPred <- PredSummary(biocartaPathFilter, wholeInter)
write.csv(biocartaPred, 'biocartaPred.csv')
nciPred <- PredSummary(nciPathFilter, wholeInter)
write.csv(nciPred, 'nciPred.csv')
reactomePred <- PredSummary(reactomePathFilter, wholeInter)
write.csv(reactomePred, 'reactomePred.csv')
################################################################################


###########################deal with complexAll################################
load('complexAll/comInterJacCor.RData')
load('top400.RData')

humComp <- lapply(MIPSHumList, function(x) {
  eachHumMat <- x[[2]][, 1:2, drop = FALSE]
  # order each mat
  eachHumMat <- OrderHumMat(eachHumMat)
  return(eachHumMat)
})

wholeInter <- apply(top400[, 1:2], 1, paste, collapse = '|')
humCompPred <- PredSummary(humComp, wholeInter)

# asign names
compNames <- sapply(MIPSHumList, '[[', 4)
compMat <- cbind(names(compNames), compNames)
compMat <- compMat[order(compMat[, 1]), ]
compMat <- compMat[rank(rownames(humCompPred)), ]

rownames(humCompPred) <- compMat[, 2]
###############################################################################


