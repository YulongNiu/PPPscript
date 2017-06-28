library('STRINGdb')

sp <- get_STRING_species(version="10", species_name=NULL)

## human
sdb <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="" )
pMat <- sdb$get_proteins()

sdb$get_homologs(pMat[1, 1], target_species_id=10090)
sdb$get_homologs(pMat[1, 1], target_species_id=sp[1:10, 1])

sdb$get_homologs_besthits(pMat[1, 1], symbets = FALSE, bitscore_threshold = 0)
lapply(1:10, function(x)return(sdb$get_homologs(pMat[1, 1], target_species_id=sp[x, 1])))
