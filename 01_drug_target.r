
library(tidyverse)
library(ComplexUpset)
STP_data <- read.csv("SwissTargetPrediction.csv")
STP_gene <- STP_data$Common.name
STITCH_data <- read.csv("STITCH.csv")
STITCH_gene <- STITCH_data$X.node1
SEA_data <- read.csv("SEA.csv") %>% filter(Name != '')
SEA_gene <- SEA_data$Name %>% toupper()
drug_gene <- sort(unique(c(STP_gene,STITCH_gene,SEA_gene)))
df <- data.frame(
  gene   = drug_gene,
  STP = drug_gene %in% STP_gene,
  STITCH = drug_gene %in% STITCH_gene,
  SEA = drug_gene %in% SEA_gene
)

p <- upset(df,intersect = c("STP", "STITCH", "SEA"),name = "",
           base_annotations=list(annotate(geom='text', x=Inf, y=Inf,label=paste('Total:', nrow(df)),vjust=2, hjust=1))
           )

p <- upset(
  data = df,
  intersect = c("STP", "STITCH", "SEA"),
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