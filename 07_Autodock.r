


library(tidyverse)
library(pheatmap)
df <- read.csv("01.results_data.csv") %>% distinct(Gene,.keep_all = T)
mat <- matrix(df$Affinity.kcal.mol.,nrow = 2,byrow = T)

pheatplot <- pheatmap::pheatmap(mat=mat,
                   color = colorRampPalette(c("#4865A9","#FDD589","#EF8A43"))(50),
                   display_numbers = T,
                   fontsize_number = 15,
                   fontsize_row = 13,
                   labels_row = c('Affinity','Non-affinity'),
                   show_colnames = F,
                   cluster_cols = F,
                   cluster_rows = F,
                   annotation_names_row = F,
                   treeheight_row = 15,
                   treeheight_col = 15,
                   border_color = 'black',
                   legend_breaks = c(-6,-5.5,-5,0),
                   number_format = "%.3f"
) 
h <- 2.5;w <- 5
ggsave(file="03.energy_heatmap.png", height = h, width = w, pheatplot)
ggsave(file="03.energy_heatmap.pdf", height = h, width = w, pheatplot)

