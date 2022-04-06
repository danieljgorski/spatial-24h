## Load libraries, functions and object----
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(readr)
library(stringr)
source("scripts/SpatialFeaturePlotScaled.R")
source("scripts/SpatialFeaturePlotScaledSig.R")
load("results/objects/hearts.Rdata")

## Genes and signatures of interest for AG Gerdes/Alexander Lang----

# genes of interest
genes_of_interest <- c(
  "Cd40", "Cd40lg", "Ccr7", "Ccl5", "Il4ra",
  "Il12b", "Cd3e", "Cd4", "Cd8a", "Cd44"
)
activ_DCs <- c(
  "Fscn1", "Ccl5", "Ccr7", "Fabp5",
  "Tmem123", "Il4i1", "Ccl17", "Samsn1", "Ccl22",
  "Tbc1d4", "Traf1", "Relb", "Socs2"
)
Inflammatory <- c(
  "Il1b", "Ifna2", "Ifng", "Tnf",
  "Ccl2", "Il6", "Il18", "Il23a", "Il33", "Il17a"
)
CD40_complex <- c(
  "Cd40", "Traf1", "Traf2", "Traf3",
  "Traf5", "Traf6"
)
all <- c(genes_of_interest, activ_DCs, Inflammatory, CD40_complex)

# SpatialFeaturePlots of genes of interest
for (i in all) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
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

# signatures
hearts <- AddModuleScore(hearts,
  features = list(activ_DCs),
  assay = "Spatial",
  name = "activ_DCs"
)
hearts <- AddModuleScore(hearts,
  features = list(Inflammatory),
  assay = "Spatial",
  name = "Inflammatory"
)
hearts <- AddModuleScore(hearts,
  features = list(CD40_complex),
  assay = "Spatial",
  name = "CD40_complex"
)
sigs <- length(colnames(hearts@meta.data))
colnames(hearts@meta.data)[sigs - 2] <- "activ_DCs"
colnames(hearts@meta.data)[sigs - 1] <- "Inflammatory"
colnames(hearts@meta.data)[sigs] <- "CD40_complex"
signatures <- c("activ_DCs", "Inflammatory", "CD40_complex")

# SpatialFeaturePlots for signatures
for (i in signatures) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
      i,
      ".jpg"
    ),
    height = 1500,
    width = 2700,
    res = 300,
    units = "px"
  )
  SpatialFeaturePlotScaledSig(
    object = hearts,
    group = "surgery",
    group.1 = "Sham",
    group.2 = "IR",
    feature_of_interest = i,
    group.1.title = "Sham - 24 h",
    group.2.title = "IR - 24 h",
    legend.title = i
  )
  dev.off()
}

## Finding markers of Cd40+ spots----

# split hearts object by surgery
Idents(hearts) <- "surgery"
hearts_sham <- subset(hearts, idents = "Sham")
hearts_ir <- subset(hearts, idents = "IR")
Idents(hearts) <- "seurat_clusters"

# finding CD40 high spots 24 h IR
VlnPlot(hearts_ir, features = "Cd40")
cd40_high_cells <- WhichCells(hearts_ir, expression = Cd40 > 0.5)
hearts_ir@meta.data$Cd40_expression <- ifelse(rownames(hearts_ir@meta.data) %in%
  cd40_high_cells, "High", "Low")
head(hearts_ir@meta.data)

# finding markers of CD40 high v low spots 24 h IR
Cd40_high_spot_markers_ir <- FindMarkers(hearts_ir,
  group.by = "Cd40_expression",
  ident.1 = "High",
  ident.2 = "Low",
  assay = "Spatial",
  logfc.threshold = 0.1,
  only.pos = T
)
Cd40_high_spot_markers_ir$gene <- rownames(Cd40_high_spot_markers_ir)
write.csv(Cd40_high_spot_markers_ir,
  file = "results/collaborations/Alexander-Lang/Cd40_high_spot_markers_ir.csv",
  row.names = F
)

# SpatialFeaturePlots of top ten markers of CD40 high spots 24 h IR
genes <- Cd40_high_spot_markers_ir$gene[1:10]
for (i in genes) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
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

# finding CD40 high spots 24 h Sham
VlnPlot(hearts_sham, features = "Cd40")
cd40_high_cells <- WhichCells(hearts_sham, expression = Cd40 > 0.5)
hearts_sham@meta.data$Cd40_expression <- ifelse(rownames(hearts_sham@meta.data) %in%
  cd40_high_cells, "High", "Low")
head(hearts_sham@meta.data)

# finding markers of CD40 high v low spots 24 h Sham
Cd40_high_spot_markers_sham <- FindMarkers(hearts_sham,
  group.by = "Cd40_expression",
  ident.1 = "High",
  ident.2 = "Low",
  assay = "Spatial",
  logfc.threshold = 0.1,
  only.pos = T
)
Cd40_high_spot_markers_sham$gene <- rownames(Cd40_high_spot_markers_sham)
write.csv(Cd40_high_spot_markers_sham,
  file = "results/collaborations/Alexander-Lang/Cd40_high_spot_markers_sham.csv",
  row.names = F
)

# SpatialFeaturePlots of top ten markers of CD40 high spots 24 h Sham
genes <- Cd40_high_spot_markers_sham$gene[1:10]
for (i in genes) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
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

## Expression correlated with Cd40----

# finding Cd40 correlated genes 24 h IR
matrix <- hearts_ir@assays$Spatial@data
matrix_mod <- as.matrix(matrix)
gene <- as.numeric(matrix_mod["Cd40", ])
correlations <- apply(matrix_mod, 1, function(x) {
  cor(gene, x)
})
correlations <- as.data.frame(correlations)
Cd40_correlated_ir <- correlations %>% dplyr::arrange(dplyr::desc(correlations))
head(Cd40_correlated_ir, 10)
Cd40_correlated_ir$gene <- rownames(Cd40_correlated_ir)
write.csv(Cd40_correlated_ir,
  file = "results/collaborations/Alexander-Lang/Cd40_correlated_ir.csv",
  row.names = F
)

# SpatialFeaturePlots of top 10 Cd40 correlated genes 24 h IR
genes <- Cd40_correlated_ir$gene[1:10]
for (i in genes) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
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

# finding Cd40 correlated genes 24 h Sham
matrix <- hearts_sham@assays$Spatial@data
matrix_mod <- as.matrix(matrix)
gene <- as.numeric(matrix_mod["Cd40", ])
correlations <- apply(matrix_mod, 1, function(x) {
  cor(gene, x)
})
correlations <- as.data.frame(correlations)
Cd40_correlated_sham <- correlations %>% dplyr::arrange(dplyr::desc(correlations))
head(Cd40_correlated_sham, 10)
Cd40_correlated_sham$gene <- rownames(Cd40_correlated_sham)
write.csv(Cd40_correlated_sham,
  file = "results/collaborations/Alexander-Lang/Cd40_correlated_sham.csv",
  row.names = F
)

# SpatialFeaturePlots of top 10 Cd40 correlated genes 24 h Sham
genes <- Cd40_correlated_sham$gene[1:10]
for (i in genes) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/SpatialFeaturePlot_",
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



## Making signatures from scRNAseq cluster markers from AG Schrader data----

# read in cluster markers and score spots
marker_all4 <- read_csv("data/cluster-markers/marker_all4.csv")
marker_all4$cluster <- str_replace_all(marker_all4$cluster, "/", "-")
marker_all4$cluster <- str_replace_all(marker_all4$cluster, "_", "-")
for (i in unique(marker_all4$cluster)) {
  df <- marker_all4[marker_all4$cluster == i, ]
  genes <- df$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts,
    features = list(genes),
    name = i,
    assay = "Spatial"
  )
}

# renaming signatures correctly
meta_data_columns <- colnames(hearts@meta.data)
all <- as.numeric(length(meta_data_columns))
new <- as.numeric(length(unique(marker_all4$cluster)))
original <- all - new
left <- meta_data_columns[1:original]
right <- unique(marker_all4$cluster)
new_col_names <- c(left, right)
colnames(hearts@meta.data) <- new_col_names

# SpatialFeature Plot of scRNAseq clusters marker signatures
for (i in unique(marker_all4$cluster)) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/cluster-marker-signatures/SpatialFeaturePlot_",
      i,
      ".jpg"
    ),
    height = 1500,
    width = 2700,
    res = 300,
    units = "px"
  )
  SpatialFeaturePlotScaledSig(
    object = hearts,
    group = "surgery",
    group.1 = "Sham",
    group.2 = "IR",
    feature_of_interest = i,
    group.1.title = "Sham - 24 h",
    group.2.title = "IR - 24 h",
    legend.title = i
  )
  dev.off()
}

# Reference mapping to day 5 Schrader data set----

# load seurat object, re-run SCT to update model list
day5 <- readRDS("data/reference-maps/all4_20211215.rds")
day5 <- SCTransform(day5, verbose = T)

# replace "/" and "_" and with "-"
day5$idents <- str_replace_all(day5$idents, "/", "-")
day5$idents <- str_replace_all(day5$idents, "_", "-")

# transfer anchors
anchors <- FindTransferAnchors(
  reference = day5,
  query = hearts,
  normalization.method = "SCT"
)
predictions_assay <- TransferData(
  anchorset = anchors,
  refdata = day5$idents,
  prediction.assay = TRUE,
  weight.reduction = hearts[["pca"]],
  dims = 1:30
)
hearts[["day5_predictions"]] <- predictions_assay
DefaultAssay(hearts) <- "day5_predictions"

# plot
for (i in unique(day5$idents)) {
  jpeg(
    file = paste0(
      "results/collaborations/Alexander-Lang/reference-mapping/day5/SpatialFeaturePlot_",
      i,
      ".jpg"
    ),
    height = 1500,
    width = 3500,
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
