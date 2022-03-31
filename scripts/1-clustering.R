## Load library----
library(Seurat) # v4.0.1

## Load in data and add meta data----
C1 <- Load10X_Spatial(
  data.dir = "data/spaceranger/328-4_C1_spaceranger_count/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Sham",
  filter.matrix = T,
  to.upper = F)
C1@meta.data$sample <- "328-4_C1"
C1@meta.data$timepoint <- "24 h"
C1@meta.data$surgery <- "Sham"
D1 <- Load10X_Spatial(
  data.dir = "data/spaceranger/328-5_D1_spaceranger_count/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "IR",
  filter.matrix = T,
  to.upper = F)
D1@meta.data$sample <- "328-5_D1"
D1@meta.data$timepoint <- "24 h"
D1@meta.data$surgery <- "IR"

## SCTransform, integration and clustering----

# SCTransform
C1 <- SCTransform(C1, assay = "Spatial")
D1 <- SCTransform(D1, assay = "Spatial")

# Using direct list instead of merge, 
# otherwise not all genes are present in the SCT model:
# https://github.com/satijalab/seurat/issues/3198 
hearts_list <- list(C1,D1) 

# integrate
features <- SelectIntegrationFeatures(object.list = hearts_list, 
                                      nfeatures = 3000)
hearts_list <- PrepSCTIntegration(object.list = hearts_list, 
                                  anchor.features = features)
hearts_anchors <- FindIntegrationAnchors(object.list = hearts_list, 
                                         normalization.method = "SCT", 
                                         anchor.features = features)
hearts <- IntegrateData(anchorset = hearts_anchors, 
                        normalization.method = "SCT")
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
VariableFeatures(hearts) <- c(VariableFeatures(C1), VariableFeatures(D1)) 

## Factor surgery variable, save clustered object----
hearts@meta.data$surgery <- factor(hearts@meta.data$surgery, 
                                   levels = c("Sham", "IR"))
save(hearts, file = "results/objects/hearts.Rdata")
