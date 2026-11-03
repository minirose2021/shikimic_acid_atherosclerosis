library(rms)
library(pROC)
library(ggplot2)

load("../09_nomogram/fit.rda")
genes <- c("ALOX5", "IL6", "NOS2")
stopifnot(all(genes %in% colnames(d)))
d$y <- factor(d$y, levels = c("Control", "AS"))
table(d$y)
dat <- d[, c(genes, "y")]
dat <- na.omit(dat)

dd <- datadist(dat)
options(datadist = "dd")
fit3 <- lrm(
  y ~ ALOX5 + IL6 + NOS2,
  data = dat,
  x = TRUE,
  y = TRUE
)

pred <- plogis(fit3$linear.predictors)
roc_combined <- roc(
  response = dat$y,
  predictor = pred,
  levels = c("Control", "AS"),
  direction = "<",
  quiet = TRUE,
  ci=T
)
cat("Combined AUC:", as.numeric(auc(roc_combined)), "\n")
print(ci.auc(roc_combined))
      
set.seed(123)
val <- validate(fit3, method = "boot", B = 1000)
print(val)
auc_corrected <- (val["Dxy", "index.corrected"] + 1) / 2
cat("Bootstrap-corrected AUC:",round(auc_corrected, 3), "\n")

p <- ggroc(roc_combined, legacy.axes = T,size = 0.8,color = "firebrick",linewidth = 1.8,alpha=0.9) + 
  geom_abline(slope = 1,linetype=2)+
  annotate("text",x=0.7,y=0.2,
           label=paste0('AUC: ',round(as.numeric(auc(roc_combined)), 3)," (",round(roc_combined$ci[1],3),'~',round(roc_combined$ci[2],3),")"),size=5)+
  annotate("text",x = 0.65,y = 0.1,label = paste0("Bootstrap-corrected AUC: ",round(auc_corrected, 3)),size = 5) +
  labs(x = "Specificity",y = "Sensitivity",title = "Model ROC (GSE20129)")+
  theme(axis.title.y = element_text(size = 18,color = 'black'),
        axis.title.x =element_text(size = 18,color = 'black'),
        axis.text.x = element_text(size = 15,color = 'black'),
        axis.text.y = element_text(size = 15,color = 'black'),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA),
        panel.background = element_rect(fill=NA,color = NA),
        plot.title = element_text(size = 20,hjust = 0.5))
print(p)
ggsave(filename = "03.ROC.pdf", width = 5, height = 5)
ggsave(filename = "03.ROC.png", width = 5, height = 5)

library(rms)
load("../09_nomogram/fit.rda")
hl2 <- stats::resid(fit,"gof")
hl2
pval <- signif(hl2[5], 3)
pval
set.seed(3)
cal1 <- rms::calibrate(fit,method='boot', B=1000) 
cali_plot <- function(){
  par(mar = c(6,5,2,2))
  plot(cal1, lwd=2, lty=1, 
       cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, 
       xlim=c(0, 1), ylim= c(0, 1), 
       xlab="Nomogram-Predicted Probability of AS Risk", 
       ylab="Actual Probability of AS Risk", 
       legend=FALSE)
  lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
  abline(0, 1, lty=3, lwd=2) 
  legend(x=0.6, y=0.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
         lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
  text(x = 0.15, y = 0.8, paste0("Hosmer-Lemeshow "))
  text(x = 0.4, y = 0.8, as.expression(bquote(italic('P')==.(pval))))
}
cali_plot()
h <- 5;w <- 7
pdf("01.calibrate.pdf",width=w,height=h)
cali_plot()
dev.off()
png("01.calibrate.png", width=w,height=h, units = "in", res = 300)
cali_plot()
dev.off()

rm(list = ls()); gc()
library(magrittr)
library(rms)
load("../09_nomogram/fit.rda")
library(rmda)
d$y <- as.numeric(d$y)-1
set.seed(123)
IL6 <- rmda::decision_curve(y~IL6,data = d, family = binomial(link ='logit'),
                              thresholds= seq(0,1, by = 0.01),
                              confidence.intervals =0.95,study.design = 'case-control',
                              population.prevalence = 0.3)
set.seed(123)
NOS2 <- rmda::decision_curve(y~NOS2,data = d, family = binomial(link ='logit'), 
                            thresholds= seq(0,1, by = 0.01),
                            confidence.intervals =0.95,study.design = 'case-control',
                            population.prevalence = 0.3)
set.seed(123)
ALOX5 <- rmda::decision_curve(y~ALOX5,data = d, family = binomial(link ='logit'), 
                              thresholds= seq(0,1, by = 0.01),
                              confidence.intervals =0.95,study.design = 'case-control',
                              population.prevalence = 0.3)
set.seed(123)
nomogram <- rmda::decision_curve(y~IL6+NOS2+ALOX5,data = d,
                                 family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                                 confidence.intervals= 0.95,study.design = 'case-control',
                                 population.prevalence= 0.3)
List <- list(IL6,NOS2,ALOX5,nomogram)

dca_plot <- function(){
  rmda::plot_decision_curve(List,curve.names= c('IL6','NOS2','ALOX5','nomogram'),
                            cost.benefit.axis =FALSE,
                            legend.position = 'topright',
                            col = c("orange","#4865A9",'pink','red'),
                            confidence.intervals =FALSE,
                            standardize = T,
                            xlim = c(0,0.5),
                            lwd = 2.5,
                            xlab='',
                            ylab='',
                            cex.axis = 1.5
  )
  mtext("Risk threshold", side = 1, line = 2.5, cex = 1.8)
  mtext("Net benefit", side = 2, line = 2.5, cex = 1.8)
}
h <- 6;w <- 6
png("02.DCA.png",width = w,height = h,res = 300,units = "in")
dca_plot()
dev.off()
pdf("02.DCA.pdf",width = w,height = h)
dca_plot()
dev.off()
