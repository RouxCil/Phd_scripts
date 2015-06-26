library(bio3d)
wd <- 'Desktop/Protein_analysis/Prep/Du156/'
in_file <- 'Du156_m8_prep.pdb'
out_file <- 'm8_link_info.txt'
glyc_type <- 10 # m8: 8+2=10 m9: 9+2=11
cut_off <- 2
PNGS_failed <- read.table('Desktop/Protein_analysis/Wiggler/Du156/nosolution_m8.txt', header = T)

pdb <- read.pdb(paste(wd, in_file, sep =''))
atom <- pdb[['atom']]

hetatoms <- atom[atom$eleno == 3 & atom$elety == "C1", c('resno','x','y','z')]
n_glyc <- dim(hetatoms)[1]
dat_dim <- (glyc_type - 2)*n_glyc

het_to_het_links <- data.frame(fresno = numeric(dat_dim), 
                               felety = character(dat_dim),
                               sresno = numeric(dat_dim),
                               selety = character(dat_dim), stringsAsFactors = F)
for (i in 1:n_glyc)
{
  if (glyc_type == 10)
  {
    het_to_het_links[i,] <- c(hetatoms$resno[i], 'O4', hetatoms$resno[i] + 1, 'C1')
    het_to_het_links[i + n_glyc,] <- c(hetatoms$resno[i] + 1, 'O4', hetatoms$resno[i] + 2, 'C1')
    het_to_het_links[i + 2*n_glyc,] <- c(hetatoms$resno[i] + 2, 'O6', hetatoms$resno[i] + 3, 'C1')
    het_to_het_links[i + 3*n_glyc,] <- c(hetatoms$resno[i] + 2, 'O3', hetatoms$resno[i] + 7, 'C1')
    het_to_het_links[i + 4*n_glyc,] <- c( hetatoms$resno[i] + 7, 'O2', hetatoms$resno[i] + 8, 'C1')
    het_to_het_links[i + 5*n_glyc,] <- c(hetatoms$resno[i] + 8, 'O2', hetatoms$resno[i] + 9, 'C1')
    het_to_het_links[i + 5*n_glyc,] <- c(hetatoms$resno[i] + 3, 'O6', hetatoms$resno[i] + 4, 'C1')
    het_to_het_links[i + 6*n_glyc,] <- c(hetatoms$resno[i] + 3, 'O3', hetatoms$resno[i] + 6, 'C1')
    het_to_het_links[i + 7*n_glyc,] <- c(hetatoms$resno[i] + 4, 'O2', hetatoms$resno[i] + 5, 'C1') 
  }
  if (glyc_type == 11)
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
}

write.table(het_to_het_links, paste(wd, out_file, sep = ''), quote = F, row.names = F, col.names = F)

Linkto <- atom[atom$elety == "ND2" & 
                 atom$resid == "NLN" & 
                 !(atom$resno %in% PNGS_failed$SITE), c('resno','x','y','z')]
links <- data.frame(resno = numeric(n_glyc), elety = character(n_glyc), 
                    hetresno = numeric(n_glyc), hetelety = character(n_glyc), stringsAsFactors = F)
for (i in 1:n_glyc)
{
  ref <- Linkto[i,2:4]
  a <- apply(hetatoms[,2:4], 1, function(x) sqrt(sum((x-ref)^2)))
  dist_to_ref <- data.frame(site = hetatoms[,1],dist = a)
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

write.table(links, paste(wd, out_file, sep =''), quote = F, row.names = F, col.names = F, append = T)