

library(tidyverse)
library(ComplexUpset)
DisGeNET_data <- read.csv("DisGeNET.csv")
DisGeNET_gene <- DisGeNET_data$gene_symbol
GeneCards_data <- read.csv("GeneCards.csv") %>% filter(Relevance.score>5)
GeneCards_gene <- GeneCards_data$Gene.Symbol
OMIM_data <- read.csv("OMIM.csv") %>% filter(Approved.Symbol != '')
OMIM_gene <- OMIM_data$Approved.Symbol
disease_gene <- sort(unique(c(DisGeNET_gene,GeneCards_gene,OMIM_gene)))
df <- data.frame(
  gene   = disease_gene,
  DisGeNET = disease_gene %in% DisGeNET_gene,
  GeneCards = disease_gene %in% GeneCards_gene,
  OMIM = disease_gene %in% OMIM_gene
)

p <- upset(
  data = df,
  intersect = c("DisGeNET", "GeneCards", "OMIM"),
  base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background='brown', on_bar='yellow'))+ annotate(
    geom='text', x=Inf, y=Inf,
    label=paste('Total:', nrow(df)),
    vjust=1, hjust=1)+
      ylab('Intersection size')+
      xlab('')),
  name = '')
p
h <- 3;w <- 5
ggsave('01.upset.png', height = h, width = w,dpi = 300)
ggsave('01.upset.pdf', height = h, width = w,dpi = 300)