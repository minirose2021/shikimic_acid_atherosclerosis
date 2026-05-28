

library(tidyverse)
train_expr <- readRDS("../00.rawdata/GSE20129_expr_use.rds") %>% t %>% as.data.frame %>% tibble::rownames_to_column("sample")
train_group <- readRDS("../00.rawdata/GSE20129_group.rds")
train_group$group <- factor(train_group$group)
gene <- read.csv("../07_Autodock/01.results_data.csv") %>% filter(Affinity.kcal.mol. < -5)
ml_data <- merge(train_expr, train_group, by = "sample") %>% tibble::column_to_rownames(var = "sample")
colnames(ml_data)[which(colnames(ml_data)=='NOS2A')] <- 'NOS2'
ml_data <- ml_data[,c(unique(gene$Gene),"group")]
saveRDS(ml_data,'df_merge.rds')

library(glmnet)
lasso_data <- ml_data
i <- 99
set.seed(i)
fit <- glmnet(as.matrix(lasso_data[-ncol(lasso_data)]), lasso_data$group, family="binomial") 
saveRDS(fit,'fit.rds')
set.seed(i)
cvfit <- cv.glmnet(as.matrix(lasso_data[-ncol(lasso_data)]), lasso_data$group, 
                   family = "binomial", type.measure = "deviance",nfolds = 10)

{
  x <- coef(fit) 
  tmp <- as.data.frame(as.matrix(x)) 
  tmp <- tmp[-1,]
  tmp$coef <- row.names(tmp) 
  tmp <- reshape::melt(tmp, id = "coef") 
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
  tmp$coef <- gsub('_','-',tmp$coef) 
  tmp$lambda <- fit$lambda[tmp$variable+1] 
  tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
  fit_plot <- ggplot() + 
    geom_vline(xintercept = log(cvfit$lambda.min),
               size=0.5,color='grey60',
               alpha=0.8,linetype=2)+
    geom_line(data=tmp,mapping = aes(log(lambda),value,color = coef),size=0.5) + 
    geom_text(aes(x=-5,y=8,label=paste0('Minimum lambda: ',round(cvfit$lambda.min,6))))+
    labs(x='Log(lambda)',y='Coefficients',color='')+
    guides(color=guide_legend(nrow = 1 ))+
    theme(axis.title.y = element_text(size = 18,color = 'black'),
          axis.title.x =element_text(size = 18,color = 'black'),
          axis.text.x = element_text(size = 15,color = 'black',hjust = 1),
          axis.text.y = element_text(size = 15,color = 'black'),
          legend.position = 'top',
          legend.key.spacing.y = unit(0.1,'mm'),
          legend.key.spacing.x = unit(0.1,'mm'),
          legend.key.width = unit(3,'mm'),
          legend.key.height = unit(1,'mm'),
          panel.background = element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill=NA))
  
  xx <- data.frame(lambda=cvfit[["lambda"]],
                   cvm=cvfit[["cvm"]],
                   cvsd=cvfit[["cvsd"]], 
                   cvup=cvfit[["cvup"]],
                   cvlo=cvfit[["cvlo"]],
                   nozezo=cvfit[["nzero"]]) 
  xx$ll<- log(xx$lambda) 
  xx$NZERO<- paste0(xx$nozezo,' vars')
  
  cv_plot <- ggplot(xx,aes(ll,cvm))+ 
    geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),
                  width=0.05,size=0.3)+ 
    geom_vline(xintercept = log(cvfit$lambda.min),
               size=0.5,color='grey60',alpha=0.8,
               linetype=2)+ 
    annotate('text',x=-5,y=1.42,label=paste0('Minimum lambda: ',round(cvfit$lambda.min,6)))+
    geom_point(size=0.5,color='red')+
    labs(x='Log(lambda)',y='Binomial Deviance',color='')+
    guides(color=guide_legend(nrow = 2))+
    theme(axis.title.y = element_text(size = 18,color = 'black'),
          axis.title.x =element_text(size = 18,color = 'black'),
          axis.text.x = element_text(size = 15,color = 'black',hjust = 1),
          axis.text.y = element_text(size = 15,color = 'black'),
          panel.background = element_rect(fill=NA,color=NA),
          legend.position = 'top',
          panel.border = element_rect(fill=NA))
  library(cowplot)
  plot_grid(fit_plot,cv_plot,rel_widths = c(1,1))
  h <- 4;w <- 10
  ggsave("01.lasso_plot.png",height=h,width=w)
  ggsave("01.lasso_plot.pdf",height=h,width=w)
}

coef.min <- coef(cvfit, s = cvfit$lambda.min)
df.coef <- cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef <- subset(df.coef, coefficient != 0) %>% as.data.frame
df.coef <- df.coef[-1,]
dim(df.coef)#3
write.csv(df.coef, "02.lasso_gene.csv")
