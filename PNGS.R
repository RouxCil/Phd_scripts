library(ggplot2)
library(bio3d)
PNGS_data <- data.frame(pos_f = character(29*2+31), 
                        mod = character(29*2+31), 
                        col = character(29*2+31), stringsAsFactors = F)
rec_data <- data.frame(min = numeric(10), 
                       max = numeric(10),
                       lab = character(10))

CAP45_in <- read.fasta('~/Desktop/Protein_analysis/sequence_files/CAP45/CAP45envaa.fasta')[[2]]
CAP45_pm <- paste(CAP45_in[1 ,], collapse = '')
CAP45_pngs_pos <- gregexpr('N(-*)(?=[^P](-*)[TS])', CAP45_pm, perl = T)[[1]]

pos <-  0
gap <- 0
CAP45_in <- rbind(CAP45_in, rep('', ncol(CAP45_in)), rep('', ncol(CAP45_in)))
for (i in 1:ncol(CAP45_in))
{
  if (CAP45_in[2 , i] != '-') {
    pos <- pos + 1
    gap <- 0
  }
  else gap <- gap + 1
  CAP45_in[3, i] <- paste(pos, ifelse(gap == 0, '', paste('.', gap, sep = '')) , sep = '')
  CAP45_in[4, i] <- paste(pos, ifelse(gap == 0, '', letters[gap]) , sep = '')
}

Du156_in <- read.fasta('~/Desktop/Protein_analysis/sequence_files/Du156/Du156envaa.fasta')[[2]]
Du156_pm <- paste(Du156_in[1 ,], collapse = '')
Du156_pngs_pos <- gregexpr('N(-*)(?=[^P](-*)[TS])', Du156_pm, perl = T)[[1]]

pos <-  0
gap <- 0
Du156_in <- rbind(Du156_in, rep('', ncol(Du156_in)), rep('', ncol(Du156_in)))
for (i in 1:ncol(Du156_in))
{
  if (Du156_in[2 , i] != '-') {
    pos <- pos + 1
    gap <- 0
  }
  else gap <- gap + 1
  Du156_in[3, i] <- paste(pos, ifelse(gap == 0, '', paste('.', gap, sep = '')) , sep = '')
  Du156_in[4, i] <- paste(pos, ifelse(gap == 0, '', letters[gap]) , sep = '')
}

HXB2_in <- read.fasta('~/Desktop/Protein_analysis/sequence_files/Du156/Du156envaa.fasta')[[2]]
HXB2_pm <- paste(HXB2_in[2 ,], collapse = '')
HXB2_pngs_pos <- gregexpr('N(-*)(?=[^P](-*)[TS])', HXB2_pm, perl = T)[[1]]

pos <-  0
gap <- 0
HXB2_in <- rbind(HXB2_in, rep('', ncol(HXB2_in)), rep('', ncol(HXB2_in)))
for (i in 1:ncol(HXB2_in))
{
  if (HXB2_in[2 , i] != '-') {
    pos <- pos + 1
    gap <- 0
  }
  else gap <- gap + 1
  HXB2_in[3, i] <- paste(pos, ifelse(gap == 0, '', paste('.', gap, sep = '')) , sep = '')
  HXB2_in[4, i] <- paste(pos, ifelse(gap == 0, '', letters[gap]) , sep = '')
}

levels_n <- unique(c(CAP45_in[3, CAP45_pngs_pos], Du156_in[3, Du156_pngs_pos], HXB2_in[3, HXB2_pngs_pos]))
levels_s <- unique(c(CAP45_in[4, CAP45_pngs_pos], Du156_in[4, Du156_pngs_pos], HXB2_in[4, HXB2_pngs_pos]))
ord <- order(as.numeric(levels_n))
levels_s <- levels_s[ord]
levels_n <- levels_n[ord]

PNGS_data$pos_f <- factor(c(CAP45_in[3, CAP45_pngs_pos], 
                            Du156_in[3, Du156_pngs_pos], 
                            HXB2_in[3, HXB2_pngs_pos]),
                          levels = levels_n)
PNGS_data$mod <- c(rep("CAP45", 29), rep("Du156", 29), rep("HXB2", 31))
PNGS_data$col <- c(rep("CAP45.G3", 29), rep("Du156.12", 29), rep("HXB2", 31))

rec_data$min <- as.numeric(factor(c('88','133','160','197','301','332','386','442','459b','611'), levels = levels_s, ordered = T)) - 0.5
rec_data$max <- as.numeric(factor(c('88','156','190b','295','301','356','406','448','463','816'), levels = levels_s, ordered = T)) + 0.5
rec_data$lab <- c('C1','V1','V2','C2','V3','C3','V4','C4','V5','GP41')

my_theme <- theme(plot.margin = unit(c(0,1,-5,1), "mm"), 
                  text = element_text(size = 14))

p <- ggplot(data = PNGS_data) + 
  geom_rect(aes(xmin = as.numeric(factor(pos_f)) - 0.4, 
                xmax = as.numeric(factor(pos_f)) + 0.4, 
                ymin = as.numeric(factor(mod)), 
                ymax = as.numeric(factor(mod)) + 1,
                fill = col)) + 
  geom_text(aes(x = (min+max)/2, y = 4.3, label = lab), rec_data) +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#CC79A7")) +
  scale_x_continuous('HXB2 position',
                     breaks = 1:48,
                     labels = levels_s,
                     expand = c(0,0.4)) + 
  scale_y_continuous(NULL, breaks = NULL, expand = c(0,0.4), labels = NULL) +
  geom_rect(aes(xmin = min, xmax = max, 
                ymin = -Inf, ymax = Inf), alpha = 0.2, subset(rec_data, lab %in% c('V1','V2','V3','V4','V5'))) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_blank()) +
  my_theme

ggsave(file = '~/Desktop/Protein_analysis/Graphs/PNGS_ggplot.pdf', p, width = 16, height = 5, units = "cm", dpi = 600)

