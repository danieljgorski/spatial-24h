## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggrepel)
source("scripts/SpatialFeaturePlotScaled.R")
source("scripts/GOBP_fold_enrich.R")
source("scripts/VolcanoPlot.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

## Differential gene expression, each cluster, Sham v IR----

# No significance threshold
de_genes <- list()
for (i in levels(Idents(hearts))) {
  results <- FindMarkers(hearts,
    group.by = "surgery",
    subset.ident = i,
    ident.1 = "IR",
    assay = "Spatial",
    logfc.threshold = 0
  )
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  de_genes[[i]] <- results
}
deg_clusters <- do.call(rbind, de_genes)
rownames(deg_clusters) <- NULL
write.csv(deg_clusters,
  file = "results/differential-gene-expression/deg_clusters.csv",
  row.names = F
)
remove(de_genes)

# Volcano plot
for (i in unique(deg_clusters$cluster)) {
  pdf(
    file = paste0(
      "results/differential-gene-expression/VolcanoPlot_cluster_",
      i,
      "_Sham_vs_IR.pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = deg_clusters,
    identity = i,
    title = paste0("DEG cluster ", i, " - IR/Sham")
  )
  dev.off()
}

## Differential gene expression, remote clusters, Sham v IR----

# DEG testing no threshold
deg_remote <- FindMarkers(hearts,
  subset.ident = c("2", "4"),
  group.by = "surgery",
  ident.1 = "IR",
  assay = "Spatial",
  logfc.threshold = 0
)
deg_remote$cluster <- "remote_zone"
deg_remote$gene <- rownames(deg_remote)
deg_remote$regulation <- ifelse(deg_remote$avg_log2FC > 0, "Up", "Down")
rownames(deg_remote) <- NULL
write.csv(deg_remote,
  file = "results/differential-gene-expression/remote-clusters/deg_remote.csv",
  row.names = F
)

# VolcanoPlot of deg
pdf(
  file = paste0("results/differential-gene-expression/remote-clusters/VolcanoPlot_remote_Sham_vs_IR.pdf"),
  height = 6,
  width = 8,
  useDingbats = F
)
VolcanoPlot(
  df = deg_remote,
  identity = "remote_zone",
  title = "DEG remote zone - IR/Sham"
)
dev.off()
