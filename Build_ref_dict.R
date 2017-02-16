Build_Ref_dict <- function(in_file = '~/Desktop/Protein_analysis/sequence_files/CAP45/CAP45env_pm.fasta'){
  require(bio3d)
  in_data <- read.fasta(in_file)
  mat <- in_data[['ali']]
  num_pos <- 0
  gap_pos <- 0
  pro_pos <- 0
  temp <- data.frame(seq = numeric(), ref_n = numeric(), ref_s = character(), stringsAsFactors = F)
  for (i in 1:dim(mat)[2])
  {
    if (mat[2,i] != '-')
    {
      num_pos <- num_pos + 1
      gap_pos <- 0
    }
    if (mat[2,i] == '-')
    {
      gap_pos <- gap_pos + 1
    }
    if (mat[1,i] != '-')
    {
      pro_pos <- pro_pos + 1
      temp[i, 1] <- pro_pos
      temp[i, 2] <- as.numeric(paste(num_pos, '.', gap_pos, sep = ''))
      temp[i, 3] <- paste(num_pos, ifelse(gap_pos == 0, '', letters[gap_pos]), sep = '')
    } 
  }
  temp <- subset(temp, !is.na(seq))
  temp$region <- 'gp120'
  temp$region <- ifelse(temp$ref_n < 131, 'C1', temp$region)
  temp$region <- ifelse(temp$ref_n >= 131 & temp$ref_n < 158, 'V1', temp$region)
  temp$region <- ifelse(temp$ref_n >= 158 & temp$ref_n < 197, 'V2', temp$region)
  temp$region <- ifelse(temp$ref_n >= 197 & temp$ref_n < 296, 'C2', temp$region)
  temp$region <- ifelse(temp$ref_n >= 296 & temp$ref_n < 332, 'V3', temp$region)
  temp$region <- ifelse(temp$ref_n >= 332 & temp$ref_n < 385, 'C3', temp$region)
  temp$region <- ifelse(temp$ref_n >= 385 & temp$ref_n < 419, 'V4', temp$region)
  temp$region <- ifelse(temp$ref_n >= 419 & temp$ref_n < 460, 'C4', temp$region)
  temp$region <- ifelse(temp$ref_n >= 460 & temp$ref_n < 472, 'V5', temp$region)
  temp$region <- ifelse(temp$ref_n >= 472 & temp$ref_n < 512, 'V5', temp$region)
  temp$region <- ifelse(temp$ref_n >= 512 & temp$ref_n <= 665, 'gp41', temp$region)
  
  temp$ref_s <- gsub('_0', '', temp$ref_s)
  temp$seq_name <- in_data[['id']][1]
  n <- max(temp$seq)
  temp$mon <- 'A'
  build_temp <- temp
  mon <- c('A', 'B', 'C')
  for (i in 2:3){
    temp$seq <- temp$seq + n
    temp$mon <- mon[i]
    build_temp <- rbind(build_temp, temp)
  }
  return(build_temp)
}

write.table(Build_Ref_dict('~/Desktop/Protein_analysis/sequence_files/CAP45/CAP45env_pm.fasta'), 
            '~/Desktop/Protein_analysis/map/CAP45/map.txt', row.names = F)
write.table(Build_Ref_dict('~/Desktop/Protein_analysis/sequence_files/Du156/Du156env_pm.fasta'), 
             '~/Desktop/Protein_analysis/map/Du156/map.txt', row.names = F)
