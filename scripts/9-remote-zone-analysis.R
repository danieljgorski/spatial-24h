## Load libraries, functions and objects----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
load("results/objects/hearts.Rdata")

# Remote zone analyses----

# combine clusters 2+4 as remote zone
hearts <- RenameIdents(hearts,
  "0" = "0",
  "1" = "1",
  "2" = "Remote",
  "3" = "3",
  "4" = "Remote",
  "5" = "5",
  "6" = "6"
)
hearts@meta.data$remote_zone_annotation <- Idents(hearts)

# genes of interest
genes <- c("Fabp3", "Dgat2", "Cpt2")

# VlnPlots
for (i in genes) {
  pdf(
    file = paste0("results/remote-zone-analysis/VlnPlot_", i, ".pdf"),
    height = 6,
    width = 8,
    useDingbats = F
  )
  p <- VlnPlot(hearts,
    features = i,
    idents = "Remote",
    split.by = "surgery",
    assay = "Spatial",
    split.plot = F
  ) +
    theme(
      plot.title = element_text(face = "italic", size = 20),
      axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "right"
    )
  print(p)
  dev.off()
}

# DimPlot of highlighted remote zone
a <- SpatialDimPlot(hearts,
  images = "Sham",
  cells.highlight = WhichCells(hearts, idents = "Remote"),
  cols.highlight = c("#F4BA27", "#ff00ff00"),
  stroke = NA
) +
  scale_fill_manual(
    values = c("#F4BA27", "#ff00ff00"),
    labels = c("Remote", "rest")
  ) +
  NoLegend() +
  labs(title = "Remote zone", subtitle = "Sham - 24 h") +
  theme(
    plot.title = element_text(
      color = "#F4BA27",
      size = 16,
      face = "bold"
    ),
    plot.subtitle = element_text(size = 14)
  )
b <- SpatialDimPlot(hearts,
  images = "IR",
  cells.highlight = WhichCells(hearts, idents = "Remote"),
  cols.highlight = c("#F4BA27", "#ff00ff00"),
  stroke = NA
) +
  scale_fill_manual(
    values = c("#F4BA27", "#ff00ff00"),
    labels = c("Remote", "rest")
  ) +
  NoLegend() +
  labs(title = "", subtitle = "I/R - 24 h") +
  theme(
    plot.title = element_text(
      color = "#F4BA27",
      size = 16,
      face = "bold"
    ),
    plot.subtitle = element_text(size = 14)
  )
pdf(
  file = "results/remote-zone-analysis/SpatialDimPlot_remote_zone.pdf",
  height = 6,
  width = 11.25,
  useDingbats = F
)
print(a + b)
dev.off()

# reading in gene sets, score cells
gene_ontology <- read.csv("data/gene-signatures/gene_signatures.csv")
for (i in colnames(gene_ontology)) {
  genes <- gene_ontology[i]
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts,
    features = list(genes),
    name = i,
    assay = "Spatial"
  )
}

# removing extra "1" added by Seurat
meta_data_columns <- colnames(hearts@meta.data)
all <- as.numeric(length(meta_data_columns))
new <- as.numeric(length(colnames(gene_ontology)))
original <- all - new
left <- meta_data_columns[1:original]
right <- gsub(".{1}$", "", meta_data_columns)[(original + 1):all]
new_col_names <- c(left, right)
colnames(hearts@meta.data) <- new_col_names
signatures <- tail(colnames(hearts@meta.data), length(colnames(gene_ontology)))

# wilcoxon rank sum test of lipid catabolic process score in remote zone
data <- hearts@meta.data
remote_zone_data <- data[data$remote_zone_annotation == "Remote", ]
res <- wilcox.test(lipid_catabolic_process_GO.0016042 ~ surgery,
  data = remote_zone_data
)

# plot lipid catabolic score
pdf(
  file = paste0("results/remote-zone-analysis/BoxPlot_lipid_catabolic_process.pdf"),
  height = 6,
  width = 8,
  useDingbats = F
)
ggplot(remote_zone_data, aes(
  x = surgery,
  y = lipid_catabolic_process_GO.0016042,
  fill = surgery
)) +
  geom_violin() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.3, width = 0.2) +
  ggtitle("Lipid catabolic process GO:0016042",
    subtitle = paste0("p-value: ", res$p.value)
  ) +
  ylab("Gene signature score") +
  theme_bw()
dev.off()
