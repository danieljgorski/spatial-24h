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

## Genes of interest----
GOI <- c("Plin5", "Plin2", "Pdk4", "Pdp2", "Cpt1a", "Acacb",
         "Angptl4", "Slc2a4", "Cd36", "Lpl", "Acadm",
         "Ppargc1a", "Has1", "Has2", "Has3", "Vcan", "Acan",
         "Ncan", "Bcan", "Dcn", "Bgn", "Col1a1", "Fn1", "Il2", 
         "Cemip", "Hyal1", "Hyal2", "Hmmr", "Adamts5", "Dmkn", 
         "Tbx18", "Wt1", "Aqp1", "Lcn2", "Serpinb2", "Krt19",
         "Krt8", "Upk3b", "Clu", "Gpm6a", "Bnc1", "Saa3",
         "Tcf21", "Pdgfra", "Postn", "Cthrc1", "Pecam1",
         "Adgre1", "Cx3cr1", "S100a8", "S100a9", "Ccr1", "Csf3r",
         "Cd3e", "Cd79a", "Ppard", "Pparg", "Lipa", "Plin3", 
         "Slc2a1", "Myh7", "Map1lc3b")

## SpatialFeature Plots----
for (i in GOI) {
  tryCatch(
    expr = {
      pdf(file = paste0("genes_of_interest/SpatialFeaturePlot_", i,".pdf"), height = 6, width = 12.75, useDingbats = F)
      Sham_max <- max(GetAssayData(object = hearts_24Sham, slot = "data")[i,])
      IR_max <- max(GetAssayData(object = hearts_24IR, slot = "data")[i,])
      suppressMessages(
        a <- SpatialFeaturePlot(hearts, features = i, images = "Sham", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
          ggtitle("Sham - 24 h") +
          theme(legend.position = "right", 
                legend.title = element_text(face = "italic", size = 16),
                title = element_text(size = 16)) & 
          scale_fill_gradientn(colours=replicants, name=i, limits = c(0, max(c(Sham_max, IR_max))))
      )
      suppressMessages(
        b <- SpatialFeaturePlot(hearts, features = i, images = "IR", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
          ggtitle("I/R - 24 h") +
          theme(legend.position = "right", 
                legend.title = element_text(face = "italic", size = 16),
                title = element_text(size = 16)) & 
          scale_fill_gradientn(colours=replicants,name=i, limits = c(0, max(c(Sham_max, IR_max))))
      )
      print(a + b + plot_layout(guides = "collect"))
      dev.off()
      print(paste0(i, " - done"))
    },
    error = function(e) {
      pdf(file = paste0("genes_of_interest/SpatialFeaturePlot_", i,".pdf"), height = 6, width = 12.75, useDingbats = F)
      Sham_max <- max(GetAssayData(object = hearts_24Sham, slot = "data", assay = "Spatial")[i,])
      IR_max <- max(GetAssayData(object = hearts_24IR, slot = "data", assay = "Spatial")[i,])
      DefaultAssay(hearts) <- "Spatial"
      suppressMessages(
        a <- SpatialFeaturePlot(hearts, features = i, images = "Sham", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
          ggtitle("Sham - 24 h") +
          theme(legend.position = "right", 
                legend.title = element_text(face = "italic", size = 16),
                title = element_text(size = 16)) & 
          scale_fill_gradientn(colours=replicants, name=i, limits = c(0, max(c(Sham_max, IR_max))))
      )
      suppressMessages(
        b <- SpatialFeaturePlot(hearts, features = i, images = "IR", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
          ggtitle("I/R - 24 h") +
          theme(legend.position = "right", 
                legend.title = element_text(face = "italic", size = 16),
                title = element_text(size = 16)) & 
          scale_fill_gradientn(colours=replicants,name=i, limits = c(0, max(c(Sham_max, IR_max))))
      )
      DefaultAssay(hearts) <- "SCT"
      print(a + b + plot_layout(guides = "collect"))
      dev.off()
      print(paste0(i, " - done using 'Spatial' assay"))
      dev.off()
  }
  )
}

## VlnPlots----
for (i in GOI) {
  pdf(file = paste0("genes_of_interest/VlnPlot_", i,".pdf"), height = 6, width = 8, useDingbats = F)
  p <- VlnPlot(hearts, 
               features = i, 
               split.by = "surgery",
               assay = "Spatial",
               split.plot=TRUE) +
    theme(plot.title = element_text(face = "italic", size = 20),
          axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "right")
  print(p)
  dev.off()
  print(paste0(i, " - done"))
}

