##############
# Hbond Data #
##############

####################################
# Create hbond data frame function #
####################################
Build_hbond_data <- function(hbond_file, map_file, out_file)
{
  require(data.table)
  require(reshape2)
  require(plyr)
  
  # Read in hbond data
  all_hbonds <- data.frame(fread(hbond_file))
  colnames(all_hbonds)[1] <- 'Frame'
  
  # Reshape and data
  m_all_hbond <- melt(all_hbonds, id.vars = 'Frame')
  
  m_all_hbond$A <- as.numeric(gsub('.*([0-9][0-9][0-9][0-9]).*([0-9][0-9][0-9][0-9]).*$', '\\1', m_all_hbond$variable))
  m_all_hbond$D <- as.numeric(gsub('.*([0-9][0-9][0-9][0-9]).*([0-9][0-9][0-9][0-9]).*$', '\\2', m_all_hbond$variable))
  
  # Remove residue-residue interactions
  sub_all_hbond <- subset(m_all_hbond, A != D)
  
  # Read in map
  glycan_map <- read.csv(map_file, sep = "", stringsAsFactors = FALSE)
  
  # Match hbond data and map data
  sub_all_hbond$A_ref_n <- glycan_map$ref_n[match(sub_all_hbond$A, glycan_map$hetresno)]
  sub_all_hbond$D_ref_n <- glycan_map$ref_n[match(sub_all_hbond$D, glycan_map$hetresno)]
  
  WD <- subset(sub_all_hbond, A_ref_n != D_ref_n)
  
  WD$large <- with(WD, ifelse(A_ref_n > D_ref_n, A, D))
  WD$small <- with(WD, ifelse(A_ref_n < D_ref_n, A, D))
  
  WD <- WD[, !(colnames(WD) %in% c("A_ref_n", "D_ref_n", "A", "D"))]
  
  WD$small_ref_n <- glycan_map$ref_n[match(WD$small, glycan_map$hetresno)]
  WD$small_ref_s <- glycan_map$ref_s[match(WD$small, glycan_map$hetresno)]
  WD$small_mon <- glycan_map$mon[match(WD$small, glycan_map$hetresno)]
  
  WD$large_ref_n <- glycan_map$ref_n[match(WD$large, glycan_map$hetresno)]
  WD$large_ref_s <- glycan_map$ref_s[match(WD$large, glycan_map$hetresno)]
  WD$large_mon <- glycan_map$mon[match(WD$large, glycan_map$hetresno)]
  
  hbond <- ddply(WD, .(Frame, small_ref_n, small_ref_s, small_mon, large_ref_n, large_ref_s, large_mon), summarise, value = ifelse(any(value == 1) == T, 1, 0))
  write.table(hbond, out_file, row.names = F)
}

###############################
# Plot Hbond network function #
###############################

plot_hbond_networks <- function(in_data, l, col, path = '~/Desktop/')
{
  require(igraph)
  require(plyr)
  
  net <- graph.data.frame(in_data[, c('small_ref_s', 'large_ref_s', 'mod', 'Frac')], directed = F)
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  
  color_data <- data.frame(glycan = c(in_data$small_ref_s, in_data$large_ref_s), mod = c(in_data$mod, in_data$mod))
  color_data$n_mod <- as.numeric(factor(color_data$mod))
  n_map <- unique(color_data[, c('mod', 'n_mod')])
  
  which_color <<- function(input) paste(unique(input[order(input)]), collapse = '')

  color_data <- ddply(color_data, .(glycan), summarise, string = which_color(mod), n_string = which_color(n_mod))
  
  if (length(unique(color_data$n_string)) != length(col))
  {
    temp <- unique(color_data$string)
    stop(paste("You need to specify", length(temp), 
               "colours for the combinations:",paste(temp, collapse = ', ')))
  }
  
  user_color <- data.frame(comb = unique(color_data$n_string), col = col, stringsAsFactors = F)
  color_data$color <- user_color[match(color_data$n_string, user_color$comb), 'col']
  
  V(net)$color <- color_data[match(V(net)$name, color_data$glycan), 'color']
  V(net)$color <- adjustcolor(V(net)$color, 0.8)
  #V(net)$label.cex <- 1
  V(net)$shape <- "rectangle"
  V(net)$frame.color <- NA
  
  for (i in 1:nrow(n_map))
  {
    which_black <- grep(n_map[i, 'n_mod'], color_data[match(V(net)$name, color_data$glycan), 'n_string'])
    V(net)$label.color <- ifelse(1:length(V(net)$name) %in% which_black, "black", adjustcolor('black', 0.5))
    
    E(net)$curved <- 0
    E(net)$color <- ifelse(E(net)$mod == n_map[i, 'mod'], 'black', NA)
    
    lty <- as.numeric(as.character(cut(E(net)$Frac,
                                       breaks = c(0,0.25,0.5,0.75,1),
                                       labels = c(3,4,2,1))))
    E(net)$lty <- round(lty, 0)
    
    pdf(paste(path, n_map[i, 'mod'], '.pdf', sep = ''), pointsize = 14, width = 8, height = 8)
    plot(net, layout = l, asp = 0,
         vertex.size = ifelse(V(net)$name %in% c('190_b', '464_b'), 18, 10),
         vertex.size2 = 8,
         edge.arrow.size = 0)
    dev.off()
  }
}

#############################
# Read CAP45 and Du156 data #
#############################

CAP45 <- read.csv("~/Desktop/Protein_analysis/hbonds/CAP45/hbond.dat", 
                  sep = "", stringsAsFactors = F)
CAP45$mod <- 'CAP45'

CAP45_301 <- read.csv("~/Desktop/Protein_analysis/hbonds/CAP45_301/hbond.dat", 
                      sep = "", stringsAsFactors = F)
CAP45_301$mod <- 'CAP45_301'

Du156 <- read.csv("~/Desktop/Protein_analysis/hbonds/Du156/hbond.dat", 
                  sep = "", stringsAsFactors = F)
Du156$mod <- 'Du156'

Du156_301 <- read.csv("~/Desktop/Protein_analysis/hbonds/Du156_301/hbond.dat", 
                      sep = "", stringsAsFactors = F)
Du156_301$mod <- 'Du156_301'

all_data <- rbind(CAP45, Du156, CAP45_301, Du156_301)

###################
# Generate layout #
###################
require(plyr)
net_data <- ddply(all_data, .(small_ref_s, small_mon, 
                              large_ref_s, large_mon, mod), summarise, Frac = sum(value)/length(value))

require(igraph)
net_graph <- graph.data.frame(net_data[, c('small_ref_s', 'large_ref_s', 'Frac')], directed = F)
net_graph <- simplify(net_graph, remove.multiple = T, remove.loops = T)

l <- layout_with_lgl(net_graph, maxit = 500000)
tkid <- tkplot(net_graph, canvas.width = 985, canvas.height = 540,
               vertex.color = "gray", layout = l, vertex.label.cex = 4, vertex.size = 16)

l <- tkplot.getcoords(tkid)
tk_close(tkid, window.close = T)

#write.table(l, '~/Desktop/Protein_analysis/hbonds/layout.txt', row.names = F)

l <- as.matrix(read.table('~/Desktop/Protein_analysis/hbonds/layout.txt', stringsAsFactors = F, header = T))
plot_data <- ddply(net_data, .(small_ref_s, large_ref_s, mod), summarise, Frac = sum(Frac)/3)

plot_hbond_networks(plot_data, l, col = c('#56B4E9','#E69F00','gray','#009E9B'), path = '~/Desktop/Protein_analysis/hbonds/')






