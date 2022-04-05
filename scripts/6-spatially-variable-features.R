## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")

## Spatially variable features of 24 h IR heart----

# subset 24 h IR heart, re-run SCT to normalize, scale, center data in 24h IR
# alone
Idents(hearts) <- "surgery"
hearts_ir <- subset(hearts, idents = "IR")
hearts_ir <- SCTransform(hearts_ir, assay = "Spatial")

# find svf in 24 IR heart
hearts_ir <- FindSpatiallyVariableFeatures(hearts_ir,
  assay = "SCT",
  selection.method = "markvariogram",
  features = VariableFeatures(hearts_ir),
  image = "IR",
  verbose = T,
  nfeatures = 20
)
top_20_svf_ir <- head(
  SpatiallyVariableFeatures(
    hearts_ir,
    selection.method = "markvariogram"
  ),
  20
)
write.csv(top_20_svf_ir,
  file = "results/spatially-variable-features/top_20_svf_ir.csv",
  row.names = F
)

# SpatialFeaturePlots, plotted with original integrated object
for (i in top_20_svf_ir) {
  pdf(
    file = paste0(
      "results/spatially-variable-features/SpatialFeaturePlot_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 12.75,
    useDingbats = F
  )
  SpatialFeaturePlotScaled(hearts,
    group = "surgery",
    group.1 = "Sham",
    group.2 = "IR",
    feature_of_interest = i,
    from.meta.data = FALSE,
    group.1.title = "24 h - Sham",
    group.2.title = "24 h - IR",
    legend.title = i
  )
  dev.off()
}
