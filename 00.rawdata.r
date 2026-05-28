
library(GEOquery)
library(magrittr)
set_df <- c('GSE20129-GPL6104','GSE20129-GPL10558')
for(i in 1:length(set_df)){
  geo_num <- set_df[i]
  print(paste0(i,":",geo_num))
  matrix_file <- paste0(geo_num,'_series_matrix.txt.gz')
  anno_file <- paste0(geo_num,'_anno_raw.txt')
  
  print('matrix')
  if (!file.exists(matrix_file)){
    query <- GEOquery::getGEO(geo_num, getGPL = F, destdir = "./")
    query <- query[[1]]
  }else{
    query <- getGEO(filename = matrix_file,getGPL = F)
  }
  
  print('phe')
  phenotype_raw <- pData(query)
  saveRDS(phenotype_raw,paste0(geo_num,'_phenotype_raw.rds'))
  
  print('expr')
  expr_raw <- exprs(query) %>% as.data.frame()
  saveRDS(expr_raw,paste0(geo_num,'_expr_raw.rds'))
  
  print('anno')
  if (!file.exists(anno_file)){
    gpl_num <- query@annotation
    gpl <- getGEO(gpl_num)
    anno_raw <- Table(gpl)
  }else{
    anno_raw <- read.delim(anno_file,comment.char = '#')
  }
  saveRDS(anno_raw,paste0(geo_num,'_anno_raw.rds'))
  # break
}


rm(list = ls()); gc()
library(tibble)
library(dplyr)
set_df <- c('GSE20129-GPL6104','GSE20129-GPL10558')
for(i in 1:length(set_df)){
  geo_num <- set_df[i]
  print(paste0(i,":",geo_num))
  anno_raw <- readRDS(paste0(geo_num,'_anno_raw.rds'))
  expr_raw <- readRDS(paste0(geo_num,'_expr_raw.rds'))
  phe_raw <- readRDS(paste0(geo_num,'_phenotype_raw.rds'))
  range(expr_raw)
  
  expr_raw <- na.omit(expr_raw)
  ex <- expr_raw
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  expr_raw <- log2(ex)}
  range(expr_raw)
  
  anno <- anno_raw
  anno[anno == ""] <- NA
  anno <- anno %>% filter(Species == 'Homo sapiens') %>% 
    tidyr::separate(col="Symbol", into=c("GeneSymbol","n3"), extra = "merge", fill = "right", sep=" /// ") %>%
    dplyr::select(ID,GeneSymbol) %>% na.omit()
  expr <- expr_raw
  expr <- expr %>% tibble::rownames_to_column(var = "ID")
  expr_use <- merge(anno, expr, by = "ID") %>%
    dplyr::select(-ID)
  expr_use <- limma::avereps(expr_use,ID = expr_use$GeneSymbol) %>% 
    as.data.frame 
  rownames(expr_use) <- expr_use$GeneSymbol
  expr_use <- expr_use %>%
    dplyr::select(-GeneSymbol) %>%
    dplyr::mutate_all(as.numeric)

  phe <- phe_raw
  group <- data.frame(sample=phe$geo_accession,
                      group=ifelse(!grepl('no',phe$characteristics_ch1),'AS','Control')) %>% 
    arrange(group)
  table(group$group)
  group <- group[group$sample %in% intersect(group$sample,colnames(expr_use)),]
  expr_use <- expr_use[,group$sample]
  write.csv(expr_use,paste0(geo_num,'_expr_use.csv'),row.names = T)
  saveRDS(expr_use,paste0(geo_num,'_expr_use.rds'))
  saveRDS(group,paste0(geo_num,'_group.rds'))
}

rm(list = ls()); gc()
library(sva)
expr_GPL10558 <- readRDS("GSE20129-GPL10558_expr_use.rds")
expr_GPL6104 <- readRDS("GSE20129-GPL6104_expr_use.rds")
common_genes <- intersect(rownames(expr_GPL10558), rownames(expr_GPL6104))
expr_merged <- cbind(expr_GPL10558[common_genes,], expr_GPL6104[common_genes,])
boxplot(expr_merged)
batch <- c(rep("GPL10558", ncol(expr_GPL10558)),
           rep("GPL6104", ncol(expr_GPL6104)))

group_GPL10558 <- readRDS("GSE20129-GPL10558_group.rds")
group_GPL6104 <- readRDS("GSE20129-GPL6104_group.rds")
group <- rbind(group_GPL10558,group_GPL6104)

expr_merged.corrected <- limma::removeBatchEffect(as.matrix(expr_merged), batch=factor(batch), group=group$group)
saveRDS(expr_merged.corrected,'GSE20129_expr_use.rds')
saveRDS(group,'GSE20129_group.rds')
