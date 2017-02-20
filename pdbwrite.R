write_pdb <- function(in_file, out_file, p_end = c(), fails = c())
{
  require(bio3d)
  
  # Read in pdb and extract dataframe ------------------------------------- 
  pdb <- read.pdb(in_file)
  atom <- pdb[['atom']]
  
  # Remap resid to force field convention ---------------------------------
  atom$resid[atom$resid == "HIS"] <- "HIE"
  atom$resid[atom$resid == "CYS"] <- "CYX"
  atom$resid[atom$resid == "OFA"] <- "OfA"
  
  # Remap faild wigler results to ASN -------------------------------------
  atom$resid[atom$resid == "NLN" & atom$resno %in% fails] <- "ASN"
  
  # Insert TER after proteins 
  p_start <- min(atom$resno)
  
  if (length(p_end) == 0)
  {
    p_end <- atom$resno[atom$elety == "OXT"]
    if(length(p_end) == 0)
    {
      stop('p_end is missing, with no default')
    }
  }
  
  p_end <- p_end[order(p_end)]
  
  for (i in 1:length(p_end))
  {
    small_protein <- atom[atom$resno >= p_start & atom$resno <= p_end[i],]
    xyz <- c(t(small_protein[, c('x', 'y', 'z')]))
    
    write.pdb(file = out_file,
              xyz = xyz, type = "ATOM", 
              resno = small_protein$resno, 
              resid = small_protein$resid, 
              eleno = small_protein$eleno,
              elety = small_protein$elety, 
              append = TRUE, chainter = TRUE, end = FALSE)  
    p_start <- p_end[i] + 1
  }
  
  # Insert TER after glycan residues --------------------------------------
  het_start <- p_start
  hetatoms <- atom[atom$resno >= het_start, ]
  
  het_res <- unique(hetatoms$resno)
  
  for (i in 1:length(het_res))
  {
    small_glycan <- atom[atom$resno == het_res[i],]
    xyz <- c(t(small_glycan[, c('x', 'y', 'z')]))
    
    write.pdb(file = out_file,
              xyz = xyz, type = "ATOM", 
              resno = small_glycan$resno, 
              resid = small_glycan$resid, 
              eleno = small_glycan$eleno,
              elety = small_glycan$elety, 
              append = TRUE, chainter = TRUE, end = FALSE)  
  }
}

#write_pdb(in_file = '~/Desktop/Protein_analysis/Prep/Du156/Du156_m9_2_301.pdb', 
#          out_file = '~/Desktop/Protein_analysis/Prep/Du156/Du156_m9_2_301_prep.pdb', 
#          fails = c(1313,303,898,1493,357,952,1547))
