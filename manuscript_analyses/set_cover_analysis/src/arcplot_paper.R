
library(readr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggraph)
library(igraph)
library(viridis)

print(sessionInfo())

args = commandArgs(trailingOnly=TRUE)

# Rscript arcplot_generic.R infile.tsv out_prefix filter_level gene_start gene_end window plot_title

out_prefix = toString(args[2])
filt = as.numeric(args[3])
gene_start = as.numeric(args[4])
gene_end = as.numeric(args[5])
window = as.numeric(args[6])
title = toString(args[7])

# generate arcplot

counts_filt_overall = subset(read.table(file = toString(args[1]), sep='\t',header=TRUE), percent_pop_covered >= filt)

vert = unique(c(counts_filt_overall$var1,counts_filt_overall$var2)) # setting up vertices for each arc
gr = graph_from_data_frame(counts_filt_overall %>% select(var1, var2, n_inds, percent_pop_covered) %>% arrange(percent_pop_covered),
                           vertices = data.frame(name=vert,POS=vert))
lout = create_layout(gr,layout = 'linear',use.numeric=TRUE,sort.by='POS')


p = ggraph(lout) +
  geom_edge_arc(aes(colour = percent_pop_covered), edge_width=0.7) +
  scale_alpha_continuous(range=c(0.7,1)) +
  labs(colour='% population with pair',colors = inferno(256), family='Times',size=12, x='Chromosomal Position', title=title) +
  scale_x_continuous(breaks=c(gene_start - window, gene_start, gene_end, gene_end + window), limits=c(gene_start - window, gene_end + window)) +
  scale_edge_color_gradientn(colors = inferno(256,end=0.77), name='% population with\nallele-specific\nsite pair') +
  theme(panel.background = element_rect(fill = NA, colour = 'white'),
        plot.title = element_text(size=18, face="italic", family='sans'),
        legend.text=element_text(size=12, family='sans', hjust=0),
        legend.title=element_text(size=12, family='sans', hjust=0),
        axis.ticks = NULL,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family = 'sans', size=12),
        axis.text.x = element_blank(),#element_text(angle=0, margin=margin(t = 20), family='Times',size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank()) 

ggsave(paste(out_prefix,'arcplot.pdf', sep=''),plot=p, width=2.8, height=1.8, dpi=300)