
glycan_map <- function(prep_pdb, map_file, glycan_map_file)
{
  source('~/Desktop/Protein_analysis/scripts/LinkInfo.R')
  links <- link_info(in_file = prep_pdb, out_file = '', return_d = T)
  
  new_data <- links[rep(1:nrow(links), 11) ,]
  for (i in 1:nrow(links))
  {
    start <- as.numeric(links$hetresno[i])
    new_data$hetresno[new_data$resno == links$resno[i]] <- start:(start + 10)
  }
  
  map <- read.csv(map_file, sep = "", stringsAsFactors = F)
  
  new_data <- merge(new_data, map, by.x = 'resno', by.y = 'seq')
  write.table(new_data, glycan_map_file, row.names = F)
}
glycan_map('~/Desktop/Protein_analysis/Prep/Du156/Du156_m9_2_301_prep.pdb', 
           '~/Desktop/Protein_analysis/map/Du156_301/map.txt', 
           '~/Desktop/Protein_analysis/map/Du156_301/glycan_map.txt')

glycan_map_to_CA <- function(prep_pdb, map_file, glycan_map_file)
{
  source('~/Desktop/Protein_analysis/scripts/LinkInfo.R')
  links <- link_info(in_file = prep_pdb, out_file = '', return_d = T, cut_off = 6, NLN_elety = "CA")
  
  new_data <- links[rep(1:nrow(links), 11) ,]
  for (i in 1:nrow(links))
  {
    start <- as.numeric(links$hetresno[i])
    new_data$hetresno[new_data$resno == links$resno[i]] <- start:(start + 10)
  }
  
  map <- read.csv(map_file, sep = "", stringsAsFactors = F)
  
  new_data <- merge(new_data, map, by.x = 'resno', by.y = 'seq')
  write.table(new_data, glycan_map_file, row.names = F)
}
