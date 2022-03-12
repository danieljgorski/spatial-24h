## Libraries, functions and colors----
library(Seurat) #>=4.0.1
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggrepel)
library(stringr)
library(RColorBrewer)
library(scales)
default_colors <- (hue_pal()(7))
replicants <- c("#5E4FA2",
                "#388FBA", 
                "#A7DBA4",
                "#E0F298",
                "#F5FBB0",
                "#FDD27F",
                "#F78850",
                "#D8434D",
                "#AB1044")

## Load in clustered object----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

## Genes of interest----
GOI_AL <- c("Cd40", "Cd40lg", "Ccr7", "Ccl5", "Il4ra", "Il12b", "Cd3e", "Cd4", "Cd8a", "Cd44")
activ_DCs <- c("Fscn1", "Ccl5", "Ccr7", "Fabp5", "Tmem123", "Il4i1", "Ccl17", "Samsn1", "Ccl22",
               "Tbc1d4", "Traf1", "Relb", "Socs2")
Inflammatory <- c("Il1b", "Ifna2", "Ifng", "Tnf", "Ccl2", "Il6", "Il18", "Il23a", "Il33", "Il17a")
CD40_complex <- c("Cd40", "Traf1", "Traf2", "Traf3", "Traf5", "Traf6")
all <- c(GOI_AL, activ_DCs, Inflammatory, CD40_complex)
spatial_assay_only_GOIs <- c("Cd40lg", "Il12b", "Ifna2", "Ifng", "Il17a")
all <- setdiff(all, spatial_assay_only_GOIs)

## SpatialFeature Plots----
for (i in all) {
  Sham_max <- max(GetAssayData(object = hearts_24Sham, slot = "data")[i,])
  IR_max <- max(GetAssayData(object = hearts_24IR, slot = "data")[i,])
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max))))
  print(p)
  dev.off()
}
for (i in spatial_assay_only_GOIs) {
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) 
  print(p)
  dev.off()
}

## Signatures----
hearts <- AddModuleScore(hearts, features = list(activ_DCs), assay = "Spatial", name= "activ_DCs")
hearts <- AddModuleScore(hearts, features = list(Inflammatory), assay = "Spatial", name= "Inflammatory")
hearts <- AddModuleScore(hearts, features = list(CD40_complex), assay = "Spatial", name= "CD40_complex")
length(colnames(hearts@meta.data))
colnames(hearts@meta.data)[13] <- "activ_DCs"
colnames(hearts@meta.data)[14] <- "Inflammatory"
colnames(hearts@meta.data)[15] <- "CD40_complex"
signatures <- c("activ_DCs", "Inflammatory", "CD40_complex")
Idents(hearts) <- "surgery"
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"
for (i in signatures) {
  Sham_max <-max(hearts_24Sham[[as.character(i)]])
  IR_max <-max(hearts_24IR[[as.character(i)]])
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max)))) 
  print(p)
  dev.off()
}

## CD40 high spot markers 24 h IR----
DefaultAssay(hearts_24IR) <- "SCT"
Idents(hearts_24IR) <- "surgery"
VlnPlot(hearts_24IR, features = "Cd40")
cd40_high_whichcells <- WhichCells(hearts_24IR, expression = Cd40 > 0.5)
hearts_24IR@meta.data$Cd40_expression <- ifelse(rownames(hearts_24IR@meta.data) %in% 
                                                  cd40_high_whichcells, "High", "Low")
head(hearts_24IR@meta.data)
Cd40_high_spot_markers_24_IR <- FindMarkers(hearts_24IR, 
                                            group.by = "Cd40_expression",
                                            ident.1 = "High",
                                            ident.2 = "Low",
                                            assay = "Spatial", 
                                            logfc.threshold = 0.1,
                                            only.pos= T)
Cd40_high_spot_markers_24_IR$gene <- rownames(Cd40_high_spot_markers_24_IR)
write.xlsx(Cd40_high_spot_markers_24_IR, file="shared_alexander_lang/Cd40_high_spot_markers_24_IR.xlsx", 
           row.names = F)
ten <- Cd40_high_spot_markers_24_IR$gene[1:10]
for (i in ten) {
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, 
       width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) 
  print(p)
  dev.off()
}

## CD40 high spot markers 24 h Sham----
DefaultAssay(hearts_24Sham) <- "SCT"
Idents(hearts_24Sham) <- "surgery"
VlnPlot(hearts_24Sham, features = "Cd40")
cd40_high_whichcells <- WhichCells(hearts_24Sham, expression = Cd40 > 0.5)
hearts_24Sham@meta.data$Cd40_expression <- ifelse(rownames(hearts_24Sham@meta.data) %in% 
                                                    cd40_high_whichcells, "High", "Low")
head(hearts_24Sham@meta.data)
Cd40_high_spot_markers_24_Sham <- FindMarkers(hearts_24Sham, 
                                              group.by = "Cd40_expression",
                                              ident.1 = "High",
                                              ident.2 = "Low",
                                              assay = "Spatial", 
                                              logfc.threshold = 0.1,
                                              only.pos= T)
Cd40_high_spot_markers_24_Sham$gene <- rownames(Cd40_high_spot_markers_24_Sham)
write.xlsx(Cd40_high_spot_markers_24_Sham, file="shared_alexander_lang/Cd40_high_spot_markers_24_Sham.xlsx", 
           row.names = F)
ten <- Cd40_high_spot_markers_24_Sham$gene[1:10]
for (i in ten) {
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, 
       width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) 
  print(p)
  dev.off()
}

## Highly correlated genes  24 h IR----
matrix<-hearts_24IR@assays$Spatial@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Cd40",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations <- as.data.frame(correlations)
Cd40_correlated_24_IR <- correlations %>% dplyr::arrange(dplyr::desc(correlations))
head(Cd40_correlated_24_IR, 10)
Cd40_correlated_24_IR$gene <- rownames(Cd40_correlated_24_IR)
write.xlsx(Cd40_correlated_24_IR, file = "shared_alexander_lang/Cd40_correlated_24_IR.xlsx", row.names = F)
ten <- Cd40_correlated_24_IR$gene[1:10]
for (i in ten) {
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, 
       width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) 
  print(p)
  dev.off()
}

## Highly correlated genes  24 h Sham----
matrix<-hearts_24Sham@assays$Spatial@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Cd40",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations <- as.data.frame(correlations)
Cd40_correlated_24_Sham <- correlations %>% dplyr::arrange(dplyr::desc(correlations))
head(Cd40_correlated_24_Sham, 10)
Cd40_correlated_24_Sham$gene <- rownames(Cd40_correlated_24_Sham)
write.xlsx(Cd40_correlated_24_Sham, file = "shared_alexander_lang/Cd40_correlated_24_Sham.xlsx", row.names = F)
ten <- Cd40_correlated_24_Sham$gene[1:10]
for (i in ten) {
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, 
       width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) 
  print(p)
  dev.off()
}

## scRNAseq cluster signatures----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
clus_sig <- read_excel("shared_alexander_lang/scRNAseq_cluster_signatures.xlsx")
for (i in colnames(clus_sig)) {
  genes <- clus_sig[i] 
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts, features = list(genes), name = i, assay = "Spatial")
}
signatures <- tail(colnames(hearts@meta.data), length(colnames(clus_sig)))
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

# SpatialFeature Plot of scRNAseq cluster signatures
for (i in signatures) {
  Sham_max <-max(hearts_24Sham[[as.character(i)]])
  IR_max <-max(hearts_24IR[[as.character(i)]])
  jpeg(file = paste0("shared_alexander_lang/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max)))) 
  print(p)
  dev.off()
}
## scRNAseq cluster markers all----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
marker_all4 <- read_csv("shared_alexander_lang/cluster_markers/marker_all4.csv")
marker_all4$cluster <- str_replace_all(marker_all4$cluster , "/", "-")
marker_all4$cluster <- str_replace_all(marker_all4$cluster , "_", "-")

for (i in unique(marker_all4$cluster)) {
  df <- marker_all4[marker_all4$cluster==i,]
  genes <- df$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts, features = list(genes), name = i, assay = "Spatial")
}

signatures <- tail(colnames(hearts@meta.data), length(unique(marker_all4$cluster)))
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

# SpatialFeature Plot of cluster markers
for (i in signatures) {
  Sham_max <-max(hearts_24Sham[[as.character(i)]])
  IR_max <-max(hearts_24IR[[as.character(i)]])
  jpeg(file = paste0("shared_alexander_lang/cluster_markers/all/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max)))) 
  print(p)
  dev.off()
}

## scRNAseq cluster markers top 10----
rm(hearts, hearts_24IR, hearts_24Sham) # to keep meta data clean
load("hearts.Rdata")
Idents(hearts) <- "surgery"
top10_all4 <- read_csv("shared_alexander_lang/cluster_markers/top10_all4.csv")
top10_all4$cluster <- str_replace_all(top10_all4$cluster , "/", "-")
top10_all4$cluster <- str_replace_all(top10_all4$cluster , "_", "-")

for (i in unique(top10_all4$cluster)) {
  df <- top10_all4[top10_all4$cluster==i,]
  genes <- df$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts, features = list(genes), name = i, assay = "Spatial")
}

signatures <- tail(colnames(hearts@meta.data), length(unique(top10_all4$cluster)))
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

# SpatialFeature Plot of cluster markers
for (i in signatures) {
  Sham_max <-max(hearts_24Sham[[as.character(i)]])
  IR_max <-max(hearts_24IR[[as.character(i)]])
  jpeg(file = paste0("shared_alexander_lang/cluster_markers/top10/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max)))) 
  print(p)
  dev.off()
}
