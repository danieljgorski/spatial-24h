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
library(xlsx)
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

## Analysis for A01 Kooperationen----
GOI <- c("Cacna1d", "Cacna1c")
for (i in GOI) {
  Sham_max <- max(GetAssayData(object = hearts_24Sham, slot = "data")[i,])
  IR_max <- max(GetAssayData(object = hearts_24IR, slot = "data")[i,])
  jpeg(file = paste0("shared_A01_kooperationen/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1)) & 
    scale_fill_gradientn(colours=replicants,  limits = c(0, max(c(Sham_max, IR_max))))
  print(p)
  dev.off()
}