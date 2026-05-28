

library(Seurat)
library(ggplot2)
samples <- list.files('../00_rawdata/GSE159677_RAW/')
sceList <- lapply(samples,function(smp){
  print(smp)
  sce <- CreateSeuratObject(counts =  Read10X(file.path('../00_rawdata/GSE159677_RAW',smp,'outs/filtered_feature_bc_matrix')),
                            project =  smp,
                            min.cells = 5,
                            min.features = 300)
  return(sce)
})
sce.all <- merge(x=sceList[[1]],y=sceList[-1])

sce <- sce.all
group <- readRDS("../00_rawdata/GSE159677_group.rds")
sce$group <- group$group[match(sce$orig.ident,group$sample)]
table(sce$group)
saveRDS(sce,"scRNA_raw.rds")

rm(list=ls()); gc()
scRNA <- readRDS("scRNA_raw.rds")
scRNA[['percent.mt']] <- PercentageFeatureSet(scRNA,pattern = "^MT-")
scRNA <- subset(scRNA,subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA<30000)
length(colnames(scRNA))
length(rownames(scRNA))
qc_plot <- VlnPlot(scRNA, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                   group.by = "group",pt.size=0,
                   cols = c("#4865A9","#EF8A43"))
w <- 8;h <- 4
ggsave('01.QC.png',width = w,height = h,plot = qc_plot)
ggsave('01.QC.pdf',width = w,height = h,plot = qc_plot)

scRNA.norm <- Seurat::NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)

scRNA.norm <- Seurat::FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
scRNA.nor.sca <- Seurat::ScaleData(scRNA.norm)
plot1 <- Seurat::VariableFeaturePlot(scRNA.nor.sca,cols = c("black", "#EF8A43"))
plot1$data <- plot1$data %>% na.omit()
library(dplyr)
top10 <- plot1$data %>% arrange(variance.standardized) %>% tail(10) %>% rownames
var_plot <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)+
  theme(axis.title =element_text(size = 18,color = 'black'),
        axis.text = element_text(size = 15,color = 'black'),
        legend.text = element_text(size = 13,color = 'black'),
        panel.background = element_rect(fill=NA,color=NA),
        panel.border = element_rect(fill=NA),
        legend.position = 'top')
w <- 8;h <- 7
ggsave('02.var_gene10.png',width = w,height = h,plot = var_plot)
ggsave('02.var_gene10.pdf',width = w,height = h,plot = var_plot)

scRNA.norm.pca <- Seurat::RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca))
library(harmony)
scRNA.norm.pca <- harmony::RunHarmony(scRNA.norm.pca,group.by.vars='orig.ident')
saveRDS(scRNA.norm.pca,'scRNA.norm.pca.rds')

pdf('test.pdf',width = 15,height = 15)
DimHeatmap(scRNA.norm.pca, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

Elbow_plot <- Seurat::ElbowPlot(scRNA.norm.pca, ndims = 50)
w <- 6;h <- 3
ggsave('03.Elbow_plot.png',width = w,height = h,plot = Elbow_plot)
ggsave('03.Elbow_plot.pdf',width = w,height = h,plot = Elbow_plot)

scRNA.norm.pca.c <- Seurat::FindNeighbors(scRNA.norm.pca,reduction = 'harmony',dims = 1:30)
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.8)#
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.7)#
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.6)#
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.5)#
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.4)#
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.3)#17
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.2)#14
scRNA.norm.pca.c <- Seurat::FindClusters(scRNA.norm.pca.c,resolution = 0.1)#9
saveRDS(scRNA.norm.pca.c,"scRNA_noAnno.rds")

UMAP <- RunUMAP(scRNA.norm.pca.c, check_duplicates = FALSE,dims = 1:30)
saveRDS(UMAP,"UMAP.rds")
sel.clust <- "RNA_snn_res.0.2"
UMAP <- SetIdent(UMAP, value = sel.clust)
table(UMAP@active.ident,UMAP$group)

gene <- c(
  'PRDM1', 'XBP1', 'SDC1',
  'CD2','CD7','CD3D','CD3E','CD3G',
  'CD19','CD79A', 'CD79B', 'MS4A1','PAX5','BLK','VPREB3','BANK1',
  'CD5','NKG7','GNLY','KLRF1',
  'TAGLN','MYH11','ACTA2','MYL9',
  'DCN','GSN','COL1A2','COL3A1',
  'CD14','PTPRC','HAMP','CD163',
  'CLDN5','PECAM1','VWF','AQP1','CDH5'
)

Seurat::DotPlot(object = UMAP,features = unique(gene))+
  theme(axis.text.x = element_text(angle = 90))

UMAP$seurat_clusters <- UMAP$RNA_snn_res.0.2
celltype <- data.frame(ClusterID=0:(length(unique(UMAP$seurat_clusters))-1),celltype= 0:(length(unique(UMAP$seurat_clusters))-1)) 

{
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(11),"Plasma cells",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(0,9),"T cells",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(6),"B cells",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(7,8),"NK",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(2,10),"SMCs",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(5),"Fibroblast",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(3,4,12),"Macrophages",celltype$celltype)
  celltype$celltype <- ifelse(celltype$ClusterID %in% c(1,13),"ECs",celltype$celltype)
  
}

UMAP@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  UMAP@meta.data[which(UMAP@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(UMAP@meta.data$celltype)

cluster_plot <- DimPlot(UMAP, reduction = "umap",
                        group.by = "RNA_snn_res.0.2",
                        label = T,label.box = T,raster = T,
                        label.color = "black",label.size = 2.5,
                        pt.size = 1.5,alpha = 1)+
  guides(color=guide_legend(nrow = 2))+
  labs(x='UMAP-1',y='UMAP-2',title = 'Clusters (resolution 0.2)')+
  theme(panel.border = element_rect(fill = NA,size = 0.5),
        legend.position = 'top',
        axis.text = element_text(size = 15,color = 'black'))
cluster_plot
anno_plot <- DimPlot(UMAP, reduction = "umap", 
                     group.by = "celltype", 
                     label = T,pt.size = 1.5,raster = T,
                     label.box = F,repel = T)+
  labs(x='UMAP-1',y='UMAP-2',title = "Cell type")+
  guides(color=guide_legend(nrow = 2,ncol = 4,byrow = T))+
  theme(panel.border = element_rect(fill = NA,size = 1),
        legend.position = 'top',
        legend.key.spacing.x = unit(0.1,'mm'),
        legend.key.spacing.y = unit(0.1,'mm'),
        axis.text = element_text(size = 15,color = 'black'))+
  scale_color_manual(values = c(`T cells`="#00BFFF",
                                `SMCs`="#8FBC8F",
                                Macrophages="#EF8A43",
                                ECs='#FFD700',
                                `Plasma cells`="#FF5F90",
                                `NK`="#9370DB",
                                `B cells`="#6495ED",
                                Fibroblast='#FFBFFF'
  ))
library(patchwork)
cell_plot <- cluster_plot+anno_plot+plot_layout(design = "AB")
w <- 10.5;h <- 6
ggsave('04.UMAP_anno.png',width = w,height = h,plot = cell_plot)
ggsave('04.UMAP_anno.pdf',width = w,height = h,plot = cell_plot)

saveRDS(UMAP,"UMAP_anno.rds")

UMAP$celltype <- factor(UMAP$celltype,levels = c('B cells',
                                                 'T cells',
                                                 'SMCs',
                                                 'Macrophages',
                                                 'ECs',
                                                 'NK',
                                                 'Fibroblast',
                                                 'Plasma cells'))

marker_gene <- c(
  'CD79A', 'CD79B', 'MS4A1','BANK1',
  'CD2','CD3D','CD3E','CD3G',
  'TAGLN','MYH11','ACTA2','MYL9',
  'CD14','CD163',
  'PECAM1','VWF','AQP1',
  'NKG7','GNLY','KLRF1',
  'DCN','COL1A2','COL3A1',
  'XBP1'
)
Seurat::DotPlot(object = UMAP,features = unique(marker_gene),group.by = "celltype",cols = c("lightgrey","#EF8A43"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))+
  labs(x='',y='')
h <- 4;w <- 8
ggsave("05.marker_plot.pdf",width = w,height = h)
ggsave("05.marker_plot.png",width = w,height = h)


UMAP <- readRDS("UMAP_anno.rds")
hub_gene <- read.csv("../07_Autodock/01.results_data.csv") %>% 
  distinct(Gene,.keep_all = T) %>% 
  filter(Affinity.kcal.mol. < -5)
Seurat::DotPlot(object = UMAP,features = hub_gene$Gene,group.by = "celltype",cols = c("lightgrey","#EF8A43"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))+
  labs(x='',y='')

h <- 4;w <- 6
ggsave("06.hub_plot.pdf",width = w,height = h)
ggsave("06.hub_plot.png",width = w,height = h)

gene_plot <- FeaturePlot(UMAP,features = hub_gene$Gene,
                         cols = c("gray","#EF8A43"),ncol = 3)+
  labs(x='UMAP-1',y='UMAP-2')+
  theme(panel.border = element_rect(fill = NA,size = 1),
        axis.text = element_text(size = 15,color = 'black'))
gene_plot[[1]]$labels$x <- 'UMAP-1'
gene_plot[[1]]$labels$y <- 'UMAP-2'
gene_plot[[2]]$labels$x <- 'UMAP-1'
gene_plot[[2]]$labels$y <- 'UMAP-2'
gene_plot
h <- 3;w <- 9
ggsave("07.hub_dim.pdf",width = w,height = h)
ggsave("07.hub_dim.png",width = w,height = h)
