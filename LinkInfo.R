
link_info <- function(in_file, out_file, cut_off = 2, return_d = F, NLN_elety = "ND2"){
  require(bio3d)
  
  # Read in pdb and extract data frame -----------------------------------------------------------
  pdb <- read.pdb(in_file)
  atom <- pdb[['atom']]
  
  # Find cystine residues
  CYX <- atom[atom$resid == 'CYX' & atom$elety == 'SG', c('resno','x','y','z')]

  # Create links for cystine bonds ---------------------------------------------------------------
  ssbond <- data.frame(fresno = numeric(), 
                       felety = character(),
                       sresno = numeric(),
                       selety = character(), stringsAsFactors = F)
  
  if (nrow(CYX) != 0)
  {
    for (i in 1:(dim(CYX)[1]/2))
    {
      ref <- CYX[1,2:4]
      a <- apply(CYX[-1,2:4], 1, function(x) sqrt(sum((x-ref)^2)))
      dist_to_ref <- data.frame(site = CYX[-1,1],dist = a)
      min_dist <- which.min(dist_to_ref$dist)
      if (dist_to_ref$dist[min_dist] < 3)
      {
        ssbond[i,] <- c(CYX[1,1], 'SG', dist_to_ref$site[min_dist], 'SG')
        CYX <- CYX[-c(1, min_dist + 1),]
      }
    } 
  }
  
  # Find the first residue of each glycan -------------------------------------------------------
  all_4YB <- unique(atom$resno[atom$resid == '4YB'])
  first_het_res <- all_4YB[all_4YB %in% (all_4YB - 1)]
  hetatoms <- subset(atom, resno %in% first_het_res & elety == 'C1', select = c('resno','x','y','z'))
  n_glyc <- dim(hetatoms)[1]

  # Create links for glycan bonds ---------------------------------------------------------------
  het_to_het_links <- data.frame(fresno = numeric(), 
                                 felety = character(),
                                 sresno = numeric(),
                                 selety = character(), stringsAsFactors = F)
  for (i in 1:n_glyc)
  {
    het_to_het_links[i,] <- c(hetatoms$resno[i], 'O4', hetatoms$resno[i] + 1, 'C1')
    het_to_het_links[i + n_glyc,] <- c(hetatoms$resno[i] + 1, 'O4', hetatoms$resno[i] + 2, 'C1')
    het_to_het_links[i + 2*n_glyc,] <- c(hetatoms$resno[i] + 2, 'O6', hetatoms$resno[i] + 3, 'C1')
    het_to_het_links[i + 3*n_glyc,] <- c(hetatoms$resno[i] + 2, 'O3', hetatoms$resno[i] + 8, 'C1')
    het_to_het_links[i + 4*n_glyc,] <- c(hetatoms$resno[i] + 8, 'O2', hetatoms$resno[i] + 9, 'C1')
    het_to_het_links[i + 5*n_glyc,] <- c(hetatoms$resno[i] + 9, 'O2', hetatoms$resno[i] + 10, 'C1')
    het_to_het_links[i + 6*n_glyc,] <- c(hetatoms$resno[i] + 3, 'O6', hetatoms$resno[i] + 4, 'C1')
    het_to_het_links[i + 7*n_glyc,] <- c(hetatoms$resno[i] + 3, 'O3', hetatoms$resno[i] + 6, 'C1')
    het_to_het_links[i + 8*n_glyc,] <- c(hetatoms$resno[i] + 6, 'O2', hetatoms$resno[i] + 7, 'C1')
    het_to_het_links[i + 9*n_glyc,] <- c(hetatoms$resno[i] + 4, 'O2', hetatoms$resno[i] + 5, 'C1')
  }
  
  # Find NLN residues
  Linkto <- atom[atom$elety == NLN_elety & 
                   atom$resid == "NLN", 
                 c('resno','x','y','z')]
  
  links <- data.frame(resno = numeric(), elety = character(), 
                    hetresno = numeric(), hetelety = character(), stringsAsFactors = F)
  
  if (dim(Linkto)[1] != dim(hetatoms)[1]) 
  {
    stop('The number of glycans and the number of NLNs do not match.')
  }
  
  for (i in 1:n_glyc)
  {
    ref <- Linkto[i, 2:4]
    a <- apply(hetatoms[, 2:4], 1, function(x) sqrt(sum((x-ref)^2)))
    dist_to_ref <- data.frame(site = hetatoms[, 1], dist = a)
    min_dist <- which.min(dist_to_ref$dist)
    if (dist_to_ref$dist[min_dist] < cut_off)
    {
      links[i,] <- c(Linkto[i,1], 'ND2', dist_to_ref$site[min_dist], 'C1')
    }
    else
    {
      warning(paste('Check if residue', Linkto[i,1], 'is really glycosylated.',
                    'If it is set the cutoff value higher'))
    }
  }
  
  if (return_d == F)
  {
    # Write cystine bonds to output file ----------------------------------------------------------
    ssbond$out_col <- paste('bond mol.', ssbond[,1], '.', ssbond[,2], ' mol.', ssbond[,3], '.', ssbond[,4], sep = '')
    write.table(ssbond$out_col, out_file, quote = F, row.names = F, col.names = F)
    
    # Write glycan bonds to output file -------------------------------------------------------------
    het_to_het_links$out_col <- paste('bond mol.', het_to_het_links[,1], '.', het_to_het_links[,2], 
                                      ' mol.', het_to_het_links[,3], '.', het_to_het_links[,4], sep = '')
    write.table(het_to_het_links$out_col, out_file, quote = F, row.names = F, col.names = F, append = T)
    
    # Write protein to glycan bonds to output file
    links$out_col <- paste('bond mol.', links[,1], '.', links[,2], 
                           ' mol.', links[,3], '.', links[,4], sep = '')
    write.table(links$out_col, out_file, quote = F, row.names = F, col.names = F, append = T)
  }
  
  if (return_d == T)
  {
    return(links[, c('resno', 'hetresno')])
  }
}


