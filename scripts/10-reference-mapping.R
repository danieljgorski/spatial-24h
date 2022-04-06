## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")

## Map to Farbehi et al. 2019 total interstitial population (TIP) data set----

# load data, re-run SCT to update model list
load("data/reference-maps/tip.Rdata")
tip <- SCTransform(tip, verbose = T)

# transfer anchors
DefaultAssay(hearts) <- "SCT"
anchors <- FindTransferAnchors(
  reference = tip,
  query = hearts,
  normalization.method = "SCT"
)
predictions_assay <- TransferData(
  anchorset = anchors,
  refdata = tip$published_annotations,
  prediction.assay = TRUE,
  weight.reduction = hearts[["pca"]],
  dims = 1:30
)
hearts[["farbehi_predictions"]] <- predictions_assay
DefaultAssay(hearts) <- "farbehi_predictions"

# plot
for (i in unique(tip$published_annotations)) {
  pdf(
    file = paste0(
      "results/reference-mapping/SpatialFeaturePlot_farbehi_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 12.75,
    useDingbats = F
  )
  SpatialFeaturePlotScaled(
    object = hearts,
    group = "surgery",
    group.1 = "Sham",
    group.2 = "IR",
    feature_of_interest = i,
    from.meta.data = FALSE,
    group.1.title = "Sham - 24 h",
    group.2.title = "IR - 24 h"
  )
  dev.off()
}
