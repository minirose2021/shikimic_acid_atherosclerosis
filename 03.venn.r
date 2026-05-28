
library(magrittr)
library(dplyr)
STP_data <- read.csv("../01_drug_target/SwissTargetPrediction.csv")
STP_gene <- STP_data$Common.name
STITCH_data <- read.csv("../01_drug_target/STITCH.csv")
STITCH_gene <- STITCH_data$X.node1
SEA_data <- read.csv("../01_drug_target/SEA.csv") %>% filter(Name != '')
SEA_gene <- SEA_data$Name %>% toupper()
drug_gene <- c(STP_gene,STITCH_gene,SEA_gene) %>% unique() %>% as.data.frame()#108
colnames(drug_gene) <- 'gene'
drug_gene$target <- rep('Shikimic acid',nrow(drug_gene))
write.csv(drug_gene,'01.drug_gene.csv',row.names = F,quote = F)

DisGeNET_data <- read.csv("../02_AS_target/DisGeNET.csv")
DisGeNET_gene <- DisGeNET_data$gene_symbol
GeneCards_data <- read.csv("../02_AS_target/GeneCards.csv") %>% filter(Relevance.score>5)
GeneCards_gene <- GeneCards_data$Gene.Symbol
OMIM_data <- read.csv("../02_AS_target/OMIM.csv") %>% filter(Approved.Symbol != '')
OMIM_gene <- OMIM_data$Approved.Symbol
disease_gene <- c(DisGeNET_gene,GeneCards_gene,OMIM_gene) %>% unique() %>% as.data.frame()#101
colnames(disease_gene) <- 'gene'
disease_gene$target <- rep('AS',nrow(disease_gene))
write.csv(disease_gene,'02.disease_gene.csv',row.names = F,quote = F)

rm(list = ls()); gc()
data1 <- read.csv("01.drug_gene.csv")
data2 <- read.csv("02.disease_gene.csv")
target1 <- unique(data1$gene)# 108
target2 <- unique(data2$gene)# 101
length(Reduce(intersect,list(target1,target2)))#6
library(ggvenn)
data <- list(target1,target2)
names(data) <- c('Shikimic acid','AS')
venn_plot <- ggvenn::ggvenn(data,
                            set_name_color = c("#EF8A43","#4865A9"),
                            set_name_size = 5,
                            text_size = 7,
                            fill_color = c("#EF8A43","#4865A9"),
                            stroke_color = "black",
                            stroke_size = 0.3)
venn_plot
h <- 4;w <- 5
ggsave('03.venn.png', height = h, width = w)
ggsave('03.venn.pdf', height = h, width = w)

inner_gene <- Reduce(intersect,list(target1,target2)) %>% as.data.frame()
colnames(inner_gene) <- 'symbol'
write.csv(inner_gene,"04.Final_Common_Gene.csv",row.names = F)

