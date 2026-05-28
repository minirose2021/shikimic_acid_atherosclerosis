

library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

gene_df <- read.csv("../03.venn/04.Final_Common_Gene.csv")
gene <- gene_df$symbol
gene_transform <- bitr(gene,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
dim(gene_transform)
#GO---------------------------------------------------------------
go <- enrichGO(gene = gene_transform$ENTREZID, 
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID", 
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               readable = TRUE)
saveRDS(go,'go.rds') 

method <- "p.adjust"
go_result <- go@result
go_result <- go_result[go_result[,method]<0.05,]
table(go_result$ONTOLOGY)

write.csv(go_result,"01.GO_diff_results.csv",row.names = F)
go_result_df <- go_result
go_result_df$GeneRatio <- go_result_df$Count/as.numeric(strsplit(go_result_df$GeneRatio[1],split = "/")[[1]][2])
library(ggplot2)
tmp_df_all <- data.frame(row.names = colnames(go_result_df)) %>% t %>% as.data.frame()

for (num in 1:2) {
  class <- unique(go_result_df$ONTOLOGY)[num]
  print(class)
  tmp_df <- go_result_df[go_result_df$ONTOLOGY==class,]
  tmp_df <- tmp_df[order(tmp_df[,method]),] %>% head(10)
  tmp_df[,method] <- format(tmp_df[,method],scitific=T,digits=3) %>% as.numeric()
  tmp_df$Description <- stringr::str_wrap(tmp_df$Description, width = 50)
  tmp_df$Description <- factor(tmp_df$Description,levels = unique(tmp_df$Description))
  tmp_df_all <- rbind(tmp_df_all,tmp_df)
  num <- num+1
}
tmp_df_all <- tmp_df_all%>% arrange(GeneRatio)
tmp_df_all$Description <- factor(tmp_df_all$Description,levels = unique(tmp_df_all$Description))
ggplot(data = tmp_df_all,
       mapping = aes(x=Description,y=GeneRatio,color=tmp_df_all[,method],size=Count))+
  geom_point()+
  coord_flip()+
  facet_grid(ONTOLOGY~.,scales = 'free')+
  theme_bw()+
  theme(text = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black',size = 15,
                                   lineheight = 0.7,
                                   margin = margin(t = 5)),
        axis.text.x = element_text(color = 'black',size = 15),
        axis.title.x = element_text(color = 'black',size = 18),
        strip.text = element_text(color = 'black',size = 15),
        legend.title = element_text(color = 'black',size = 18),
        legend.text = element_text(color = 'black',size = 15))+
  scale_size(limits = range(tmp_df_all$Count),breaks = unique(tmp_df_all$Count))+
  scale_color_gradient(low = "#EF8A43", high = "#4865A9")+
  labs(colour=method,x='',y='Gene Ratio')
h <- 8;w <- 9
ggsave('02.GO_dot.png',height = h,width = w)
ggsave('02.GO_dot.pdf',height = h,width = w)

# KEGG --------------------------------------------------------------------
kegg <- enrichKEGG(gene = gene_transform$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1)
saveRDS(kegg,'kegg.rds')
kegg_result <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result <- kegg_result@result
kegg_result <- kegg_result[kegg_result[,method]<0.05,]
kegg_result <- kegg_result[order(kegg_result[,method]),]
length(kegg_result$Description)
write.csv(kegg_result,"03.KEGG_diff_results.csv",row.names = F)
kegg_result_df <- kegg_result
kegg_result_df$GeneRatio <- kegg_result_df$Count/as.numeric(strsplit(kegg_result_df$GeneRatio[1],split = "/")[[1]][2])
kegg_result_df <- kegg_result_df[order(kegg_result_df[,method]),] %>% head(15)
kegg_result_df <- kegg_result_df %>% arrange(GeneRatio)
kegg_result_df$Description <- stringr::str_wrap(kegg_result_df$Description, width = 50)
kegg_result_df$Description <- factor(kegg_result_df$Description,levels = kegg_result_df$Description)

ggplot(data = kegg_result_df,
       mapping = aes(x=Description,y=GeneRatio,color=kegg_result_df[,method],size=Count))+
  geom_point()+
  coord_flip()+
  theme_bw()+
  theme(text = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black',size = 15,
                                   lineheight = 0.7,
                                   margin = margin(t = 5)),
        axis.text.x = element_text(color = 'black',size = 15),
        axis.title.x = element_text(color = 'black',size = 18),
        strip.text = element_text(color = 'black',size = 15),
        legend.title = element_text(color = 'black',size = 18),
        legend.text = element_text(color = 'black',size = 15))+
  scale_color_gradient(low = "#EF8A43", high = "#4865A9")+
  labs(colour=method,x='',y='Gene Ratio')

h <- 4;w <- 6.5
ggsave('04.KEGG_dot.png',height = h,width = w)
ggsave('04.KEGG_dot.pdf',height = h,width = w)

