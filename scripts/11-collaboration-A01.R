## Load libraries, functions and object----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")

## Analysis for A01 Kooperationen----
genes <- c("Cacna1d", "Cacna1c")
for (i in genes) {
  jpeg(
    file = paste0(
      "results/collaborations/A01/SpatialFeaturePlot_",
      i,
      ".jpg"
    ),
    height = 1500,
    width = 3200,
    res = 300,
    units = "px"
  )
  SpatialFeaturePlotScaled(
    object = hearts,
    group = "surgery",
    group.1 = "Sham",
    group.2 = "IR",
    feature_of_interest = i,
    from.meta.data = FALSE,
    group.1.title = "Sham - 24 h",
    group.2.title = "IR - 24 h",
    legend.title = i
  )
  dev.off()
}
