library(bio3d)
wd <- 'Desktop/Protein_analysis/Prep/Du156/'
in_file <- 'Du156_m9.pdb'
out_file <- 'Du156_m9_prep.pdb'
p_end <- c()


pdb <- read.pdb(paste(wd, in_file, sep =''))
atom <- pdb[['atom']]
atom$resid[atom$resid=="HIS"] <- "HIE"
atom$resid[atom$resid=="CYS"] <- "CYX"
atom$resid[atom$resid=="OFA"] <- "OfA"

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

n_chain <- length(p_end)
for (i in 1:n_chain)
{
  small_protein <- atom[atom$resno >= p_start & atom$resno <= p_end[i],]
  atm <- c(matrix(c(small_protein$x,
                    small_protein$y,
                    small_protein$z),
                  nrow = 3, byrow = TRUE))
  write.pdb(file = paste(wd, out_file, sep =''),
            xyz = atm, type = "ATOM", 
            resno = small_protein$resno, 
            resid = small_protein$resid, 
            eleno = small_protein$eleno,
            elety = small_protein$elety, 
            append = TRUE, chainter = TRUE, end = FALSE)  
  p_start <- p_end[i] + 1
}

hetatoms <- atom[atom$resno >= p_start,]

het_res <- unique(hetatoms$resno)
n_het <- length(het_res)

for (i in 1:n_het)
{
  small_glycan <- atom[atom$resno == het_res[i],]
  hetatm <- c(matrix(c(small_glycan$x,
                       small_glycan$y,
                       small_glycan$z),
                     nrow = 3, byrow = TRUE))
  write.pdb(file = paste(wd, out_file, sep =''),
            xyz = hetatm, type = "ATOM", 
            resno = small_glycan$resno, 
            resid = small_glycan$resid, 
            eleno = small_glycan$eleno,
            elety = small_glycan$elety, 
            append = TRUE, chainter = TRUE, end = FALSE)  
}