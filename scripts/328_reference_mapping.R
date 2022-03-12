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
## Map to day 5 Schrader data set----
# load data, re-run SCT to update model list
day5 <- readRDS("reference_maps/all4_20211215.rds")
head(day5@meta.data)
day5 <- SCTransform(day5, verbose = T) 

#Replace "/" and "_" and with "-"
day5$idents <- str_replace_all(day5$idents, "/", "-")
day5$idents <- str_replace_all(day5$idents, "_", "-")

# transfer anchors
DefaultAssay(hearts) <- "SCT"
anchors <- FindTransferAnchors(reference = day5, 
                               query = hearts, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = day5$idents, 
                                  prediction.assay = TRUE,
                                  weight.reduction = hearts[["pca"]], 
                                  dims = 1:30)
hearts[["day5_predictions"]] <- predictions.assay
DefaultAssay(hearts) <- "day5_predictions"

# plot
for (i in unique(day5$idents)) {
  jpeg(file = paste0("shared_alexander_lang/reference_mapping/day5/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1))
  print(p)
  dev.off()
}



## Map to Farbehi et al data set----
# load data, re-run SCT to update model list
load("reference_maps/tip.Rdata")
head(tip@meta.data)
tip <- SCTransform(tip, verbose = T) 

# transfer anchors
DefaultAssay(hearts) <- "SCT"
anchors <- FindTransferAnchors(reference = tip, 
                               query = hearts, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = tip$published_annotations, 
                                  prediction.assay = TRUE,
                                  weight.reduction = hearts[["pca"]], 
                                  dims = 1:30)
hearts[["farbehi_predictions"]] <- predictions.assay
DefaultAssay(hearts) <- "farbehi_predictions"

# plot
for (i in unique(tip$published_annotations)) {
  jpeg(file = paste0("shared_alexander_lang/reference_mapping/farbehi/SpatialFeaturePlot_", i,".jpg"), height = 2000, width = 2200, res = 300, units = "px")
  p <- SpatialFeaturePlot(hearts, features = i, pt.size.factor = 2, alpha = c(0.01,1))
  print(p)
  dev.off()
}


