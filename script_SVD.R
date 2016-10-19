##########################SVD preprocess############################
SVDPhy <- function(input_matrix, trimming = 1, min_conservation = -0.1){

  s <- La.svd(input_matrix)
  svdm <- (s$u)

  ## filter genes
  filteredIndices <- apply(input_matrix, 1, function(x){
    counter <- sum(x > 0)
    return(counter > min_conservation)
  })
  svdm <- svdm[filteredIndices, ]

  ## trim species
  svdmTrimmed <- svdm[, 1:round(trimming * ncol(svdm))]
  resultM <- t(apply(svdmTrimmed, 1, function(x){
    return(x/sqrt(sum(x^2)))
  }))

  rownames(resultM) <- rownames(input_matrix)
  colnames(resultM) <- colnames(input_matrix)[1:round(trimming * ncol(svdm))]

   return(resultM)
}

setwd('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/phyloMito/wholenetwork0001/')
load('/home/Yulong/RESEARCH/neuro/Bioinfor/PhyloViz/wholePhyloDataHit.RData')


## ref: SVD-phy: improved prediction of protein functional associations through singular value decomposition of phylogenetic profiles
## step1: M < 60 to 1
norProfile <- apply(wholePhyloData, 1:2, function(x){
  x <- ifelse(x < 60, 1, x)
  return(x)
})

## step2: x[i, j]/max(x[i, ])
norProfile <- apply(norProfile, 1, function(x) {
  x <- x/max(x)
  return(x)
})
norProfile <- t(norProfile)

## step3: try different trimming score, 0%, 15%, 30%, 45%
norProfile100 <- SVDPhy(norProfile, trimming = 1)
norProfile15 <- SVDPhy(norProfile, trimming = 0.15)
norProfile30 <- SVDPhy(norProfile, trimming = 0.30)
norProfile45 <- SVDPhy(norProfile, trimming = 0.45)

save(norProfile100, norProfile45, norProfile30, norProfile15, file = 'SVD_profile.RData')
####################################################################
