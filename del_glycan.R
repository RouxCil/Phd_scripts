del_glycan <- function(in_file, out_file, del = del)
{
  require(bio3d)
  require(plyr)
  
  pdb <- read.pdb(in_file)
  
  atom <- pdb[['atom']]
  del_coord <- subset(atom, elety == "ND2" & resno %in% del, select = c(resno, x, y, z))
  glycan <- subset(atom, resid == "4YB" & elety == "C1", select = c(resno, x, y, z))
  
  dist_mat <- dist.xyz(del_coord[,c('x','y','z')], glycan[,c('x','y','z')])
  sel_min <- function(x,  res_del, atom = atom) {
    n <- which.min(x)
    if (x[n] < 1.5) {
      start <- res_del[n, 'resno']
      end <- start + 10#res_del[n + 2, 'resno']
      del <- start:end
      return(del)
    }
  }
  
  del_res <- c(apply(dist_mat, 1, sel_min, res_del = glycan, atom = atom))
  
  write_out <- subset(atom, !(resno %in% del_res))
  
   unq_resno <- unique(write_out$resno)
   renumber_by <- which(c(unq_resno[-1], 0) - unq_resno != 1)[1]
   het_res <- write_out$resno[write_out$resno > renumber_by] - unq_resno[renumber_by + 1] + renumber_by + 1
   write_out$resno[write_out$resno > renumber_by] <- het_res
  
  xyz_mat <- c(matrix(c(write_out$x,
                        write_out$y,
                        write_out$z),
                      nrow = 3, byrow = TRUE))
  
  write.pdb(file = out_file,
            xyz = xyz_mat, type = "ATOM", 
            resno = write_out$resno, 
            resid = write_out$resid, 
            eleno = write_out$eleno,
            elety = write_out$elety) 
}

del_glycan(in_file = '~/Desktop/Protein_analysis/naccess/replicates/Du156_334/r10.pdb', 
           out_file = '~/Desktop/Protein_analysis/naccess/replicates/Du156_334/r10_no_301.pdb', del = c(266, 834, 1402))

