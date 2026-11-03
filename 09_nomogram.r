library(magrittr)
hub_gene <- read.csv("../07_Autodock/01.results_data.csv") %>% filter(Affinity.kcal.mol. < -5)
train_expr <- readRDS("../00.rawdata/GSE20129_expr_use.rds")
train_group <- readRDS("../00.rawdata/GSE20129_group.rds")
rownames(train_expr)[which(rownames(train_expr)=='NOS2A')] <- 'NOS2'
hub_expr <- train_expr[unique(hub_gene$Gene),train_group$sample] %>% t %>% as.data.frame()
library(rms)
library(regplot)
d <- hub_expr
d$group <- train_group$group[match(rownames(d),train_group$sample)]
d$y <- factor(d$group,levels = c('Control','AS'))
ddist <- datadist(d)
options(datadist='ddist')
fit <- rms::lrm(y ~ .,
                data = d[,!grepl('group',colnames(d))], x = TRUE, y = TRUE,maxit=1000)
save(d,fit,file = 'fit.rda')
regplot::regplot(fit,
                 plots=c("violin", "boxes"), 
                 observation=d[3,!grepl('group',colnames(d))], 
                 title="Risk of AS", 
                 clickable=F,
                 points=T,
                 droplines=T,
                 center = T,
                 dencol="#95A2BF",
                 boxcol="#ED7474",
                 showP = F,
                 odd=F)
nomo <- recordPlot()
h <- 5;w <- 8
pdf(file = "nomogram.pdf",width = w,height = h)
nomo
dev.off()
png(file = "nomogram.png",width = w,height = h, units="in", res=300) 
nomo
dev.off()
graphics.off()

beta <- coef(fit)

vcov_mat <- vcov(fit)

SE <- sqrt(diag(vcov_mat))

beta_lower <- beta - 1.96 * SE
beta_upper <- beta + 1.96 * SE

OR <- exp(beta)
OR_lower <- exp(beta_lower)
OR_upper <- exp(beta_upper)
Z <- beta / SE
P_value <- 2 * pnorm(abs(Z), lower.tail = FALSE)
results <- data.frame(
  Variable = names(beta),
  Beta = round(beta, 3),
  # SE = round(SE, 3),
  Beta_95CI = sprintf(
    "%.3f - %.3f",
    beta_lower, beta_upper
  ),
  # OR = round(OR, 3),
  # OR_95CI = sprintf(
  #   "%.3f - %.3f",
  #   OR_lower, OR_upper
  # ),
  # P_value = signif(P_value, 3),
  row.names = NULL
)

results <- results[-1,]
print(results)
write.csv(results,"Three_gene_Logistic_coefficients_95CI.csv",row.names = FALSE)
