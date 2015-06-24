library(bio3d)
wd <- 'Desktop/Protein_analysis/Prep/Du156/'
in_file <- 'Du156_m9.pdb'
out_file <- 'Du156_m9_prep.pdb'
gly_len <- 129 - 3 #atom length - start m8: 118 m9: 129

pdb <- read.pdb(paste(wd, in_file, sep =''))
atom <- pdb[['atom']]
atom$resid[atom$resid=="HIS"] <- "HIE"
atom$resid[atom$resid=="CYS"] <- "CYX"
atom$resid[atom$resid=="OFA"] <- "OfA"

p_end <- max(as.numeric(row.names(atom[atom$elety == "OXT",])))
hetatoms <- atom[atom$eleno == 3 & atom$elety == "C1",c('resno','x','y','z')]
n_glyc <- dim(hetatoms)[1]

protein_xyz <- c(matrix(c(atom$x[1:p_end],atom$y[1:p_end],atom$z[1:p_end]),nrow = 3, byrow = TRUE))
write.pdb(file = paste(wd, out_file, sep =''), xyz = protein_xyz, 
          resno = atom$resno[1:p_end], 
          resid = atom$resid[1:p_end], 
          eleno = atom$eleno[1:p_end], 
          elety = atom$elety[1:p_end],
          end = FALSE)

het_e <- p_end
for (i in 1:n_glyc)
{
  het_s <- het_e + 1
  het_e <- het_s + gly_len
  hetatm <- c(matrix(c(atom$x[het_s:het_e],
                       atom$y[het_s:het_e],atom$z[het_s:het_e]), 
                     nrow = 3, byrow = TRUE))
  write.pdb(file = paste(wd, out_file, sep =''), 
            xyz = hetatm, type = "ATOM", 
            resno = atom$resno[het_s:het_e], 
            resid = atom$resid[het_s:het_e], 
            eleno = atom$eleno[het_s:het_e], 
            elety = atom$elety[het_s:het_e], 
            append = TRUE, chainter = TRUE, end = FALSE)
}

