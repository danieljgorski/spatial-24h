## Load library----
library(Seurat) # v4.0.1

## Load in data and add meta data----
c1_sham <- Load10X_Spatial(
  data.dir = "data/spaceranger/328-4_C1_spaceranger_count/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Sham",
  filter.matrix = T,
  to.upper = F
)
c1_sham@meta.data$sample <- "328-4_C1"
c1_sham@meta.data$timepoint <- "24 h"
c1_sham@meta.data$surgery <- "Sham"
d1_ir <- Load10X_Spatial(
  data.dir = "data/spaceranger/328-5_D1_spaceranger_count/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "IR",
  filter.matrix = T,
  to.upper = F
)
d1_ir@meta.data$sample <- "328-5_D1"
d1_ir@meta.data$timepoint <- "24 h"
d1_ir@meta.data$surgery <- "IR"

## SCTransform, integration and clustering----

# SCTransform
c1_sham <- SCTransform(c1_sham, assay = "Spatial")
d1_ir <- SCTransform(d1_ir, assay = "Spatial")

# Using direct list instead of merge,
# otherwise not all genes are present in the SCT model:
# https://github.com/satijalab/seurat/issues/3198
hearts_list <- list(c1_sham, d1_ir)

# integrate
features <- SelectIntegrationFeatures(
  object.list = hearts_list,
  nfeatures = 3000
)
hearts_list <- PrepSCTIntegration(
  object.list = hearts_list,
  anchor.features = features
)
hearts_anchors <- FindIntegrationAnchors(
  object.list = hearts_list,
  normalization.method = "SCT",
  anchor.features = features
)
hearts <- IntegrateData(
  anchorset = hearts_anchors,
  normalization.method = "SCT"
)
DefaultAssay(hearts) <- "integrated"

# cluster
hearts <- RunPCA(hearts, verbose = T)
ElbowPlot(hearts, ndims = 50)
hearts <- FindNeighbors(hearts, dims = 1:20)
hearts <- FindClusters(hearts, verbose = T, res = 0.4)
hearts <- RunUMAP(hearts, reduction = "pca", dims = 1:20)
DefaultAssay(hearts) <- "Spatial"

# Normalizing and scaling Spatial assay for DEG testing----
hearts <- NormalizeData(hearts)
all_genes <- rownames(hearts)
hearts <- ScaleData(hearts, features = all_genes)

# Set variable features in SCT assay, because merge wasn't used
DefaultAssay(hearts) <- "SCT"
VariableFeatures(hearts) <- c(VariableFeatures(c1_sham), VariableFeatures(d1_ir))

## Factor surgery variable, save clustered object----
hearts@meta.data$surgery <- factor(hearts@meta.data$surgery,
  levels = c("Sham", "IR")
)
save(hearts, file = "results/objects/hearts.Rdata")
