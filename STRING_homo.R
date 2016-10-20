library('STRINGdb')

sp <- get_STRING_species(version="10", species_name=NULL)

## human
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="" )
pList <- string_db$get_proteins()

string_db$get_homologs(pList[1, 1], target_species_id=10090)
string_db$get_homologs_besthits(pList[1, 1], symbets = TRUE)
sapply(1:10, function(x)return(string_db$get_homologs(pList[1, 1], target_species_id=sp[x, 1])))
