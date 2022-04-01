## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

## Genes of interest----

# genes
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

# SpatialFeaturePlots
for (i in GOI) {
  pdf(file = paste0("results/genes-of-interest/SpatialFeaturePlot_", i,".pdf"),
      height = 6,
      width = 12.75,
      useDingbats = F)
  SpatialFeaturePlotScaled(object = hearts,
                           group = "surgery",
                           group.1 = "Sham",
                           group.2 = "IR",
                           feature_of_interest = i,
                           from.meta.data = FALSE,
                           group.1.title = "Sham - 24 h",
                           group.2.title = "IR - 24 h")
  dev.off()
}

# VlnPlots
for (i in GOI) {
  pdf(file = paste0("results/genes-of-interest/VlnPlot_", i,".pdf"),
      height = 6,
      width = 8,
      useDingbats = F)
  p <- VlnPlot(hearts, 
               features = i, 
               split.by = "surgery",
               assay = "Spatial",
               split.plot = TRUE) +
    theme(plot.title = element_text(face = "italic", size = 20),
          axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "right")
  print(p)
  dev.off()
  print(paste0(i, " - done."))
}
