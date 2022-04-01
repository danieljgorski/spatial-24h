## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

## Export basic figures----

# VlnPlot of UMI counts
pdf(
  file = "results/basic-figures/VlnPlot_UMI_count.pdf",
  height = 6,
  width = 6,
  useDingbats = F
)
VlnPlot(hearts,
  features = "nCount_Spatial",
  pt.size = 0.25,
  group.by = "surgery"
) +
  NoLegend() +
  ggtitle("") +
  ylab("UMI count") +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(size = 16)
  )
dev.off()

# SpatialFeaturePlot of UMI counts
pdf(
  file = "results/basic-figures/SpatialFeaturePlot_UMI_count.pdf",
  height = 6,
  width = 12.75,
  useDingbats = F
)
SpatialFeaturePlotScaled(
  object = hearts,
  group = "surgery",
  group.1 = "Sham",
  group.2 = "IR",
  feature_of_interest = "nCount_Spatial",
  from.meta.data = TRUE,
  group.1.title = "Sham - 24 h",
  group.2.title = "I/R - 24 h",
  legend.title = "UMI count"
)
dev.off()

# SpatialDimPlot
a <- SpatialDimPlot(hearts,
  images = "Sham",
  pt.size.factor = 1.8
) +
  labs(fill = "Cluster") +
  ggtitle("Sham - 24 h") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5))) + NoLegend()
b <- SpatialDimPlot(hearts,
  images = "IR",
  pt.size.factor = 1.8
) +
  labs(fill = "Cluster") +
  ggtitle("I/R - 24 h") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
pdf(
  file = "results/basic-figures/SpatialDimPlot.pdf",
  height = 6,
  width = 12.75,
  useDingbats = F
)
print(a + b + plot_layout(guides = "collect"))
dev.off()

# Normal DimPlot
p <- DimPlot(hearts, pt.size = 1.5) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Cluster") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 16)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p2 <- LabelClusters(
  plot = p,
  id = "ident",
  repel = F,
  box = T,
  fill = alpha("white", 0.25),
  size = 8,
  label.r = unit(0.75, "lines"),
  label.size = 0
)
pdf(
  file = "results/basic-figures/DimPlot.pdf",
  height = 6,
  width = 8,
  useDingbats = F
)
print(p2)
dev.off()

# DimPlot grouped by surgery
pdf(
  file = "results/basic-figures/DimPlot_surgery.pdf",
  height = 6,
  width = 8,
  useDingbats = F
)
DimPlot(hearts, pt.size = 1.5, group.by = "surgery") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Surgery") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 16),
    plot.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

## Single cluster highlight figures----
for (i in levels(Idents(hearts))) {
  coi <- default_colors[(as.numeric(i) + 1)]
  a <- SpatialDimPlot(hearts,
    images = "Sham",
    cells.highlight = WhichCells(hearts, idents = i),
    cols.highlight = c(coi, "#ff00ff00"),
    stroke = NA
  ) +
    scale_fill_manual(values = c(coi, "#ff00ff00"), labels = c(i, "rest")) +
    NoLegend() +
    labs(title = paste0("Cluster ", i), subtitle = "Sham - 24 h") +
    theme(
      plot.title = element_text(color = coi, size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  b <- SpatialDimPlot(hearts,
    images = "IR",
    cells.highlight = WhichCells(hearts, idents = i),
    cols.highlight = c(coi, "#ff00ff00"),
    stroke = NA
  ) +
    scale_fill_manual(values = c(coi, "#ff00ff00"), labels = c(i, "rest")) +
    NoLegend() +
    labs(title = "", subtitle = "I/R - 24 h") +
    theme(
      plot.title = element_text(color = coi, size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  pdf(
    file = paste0("results/basic-figures/SpatialDimPlot_Cluster_", i, ".pdf"),
    height = 6,
    width = 11.25,
    useDingbats = F
  )
  print(a + b)
  dev.off()
}
