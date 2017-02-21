#####################################################################################
# Function that loads data and choosese the bigest probe per time point and residue #
# defaulting to 0 if there is none.                                                 #
#####################################################################################

write_max_probe <- function(path)
{
  require(data.table)
  
  getData <- function(i, glycan)
  {
    dat <- read.csv(files[i], sep = "", stringsAsFactors = FALSE)
    probe <- as.numeric(gsub('.*md.(.*)_(.*).rsa', '\\2', files[i]))
    if (probe == 1.4)
    {
      dat <- subset(dat, !(RES %in% c('4YB', '0MA', '2MA', 'VMA', 'VMB')), select = c('RES', 'RESNO', 'all_a'))
      dat$probe <- 0
      dat$probe[dat$all_a > 0] <- probe
    }
    else
    {
      dat <- subset(dat, !(RES %in% c('4YB', '0MA', '2MA', 'VMA', 'VMB')) & all_a > 0, select = c('RES', 'RESNO', 'all_a'))
      dat$probe <- probe 
    }
    dat$time <- as.numeric(gsub('.*md.(.*)_(.*).rsa', '\\1', files[i]))
    dat$glycan <- glycan
    return(dat)
  }
  
  glycan <- paste('glycan', gsub('.*SASA_(.*)/', '\\1', path), sep = '_')
  files <- list.files(path = path, 
                      pattern = '*.rsa', full.names = T)
  
  require(plyr)
  tmp <- data.frame( rbindlist(lapply(1:length(files), function(x, glycan) getData(x, glycan), glycan = glycan )) )
  tmp <- ddply(tmp, .(RESNO, time, glycan), summarise, max_probe = max(probe))
  write.table(tmp, paste(path, 'max_probe.txt', sep = ''))
}

################################################
# Pick the xx residues affected by the glycan #
################################################

which_residue <- function(SASA_all_p, SASA_mglycan_p, map_p)
{
  SASA_all <- read.table(SASA_all_p, sep = '', stringsAsFactors = F)
  SASA_mglycan <- read.table(SASA_mglycan_p, sep = '', stringsAsFactors = F)
  map <- read.table(map_p, sep = '', stringsAsFactors = F, header = T)
  require(plyr)
  tmp <- merge(SASA_all, SASA_mglycan, by = c('RESNO', 'time'))
  tmp <- merge(tmp, map, by.x = 'RESNO', by.y = 'seq')
  tmp$shield <- ifelse(tmp$max_probe.x != tmp$max_probe.y, 'Yes', 'No')
  all_resno <- ddply(subset(tmp, max_probe.x != max_probe.y), 
                     .(RESNO, mon), summarise, len = length(RESNO))
  all_resno <- ddply(all_resno, .(mon), mutate, sep_top = sort(len, decreasing = T)[20])
  out <- subset(all_resno, len >= sep_top)$RESNO
  out <- all_resno$RESNO
  return(out)
}

