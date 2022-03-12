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
library(xlsx)
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
GOBP_fold_enrich <- function(background_genes, genes_of_interest, plot_title){
  GO_BP_over_rep <- enrichGO(gene = genes_of_interest,
                             universe = background_genes,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             keyType = "SYMBOL",
                             pAdjustMethod = "bonferroni", 
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)
  df <- GO_BP_over_rep@result
  df$mapped <- sub(".*/", "", df$GeneRatio)
  df$mapped <- as.numeric(df$mapped)
  df$GO_in_background <- sub("/.*", "", df$BgRatio)
  df$GO_in_background <- as.numeric(df$GO_in_background)
  df$background_mapped <- sub(".*/", "", df$BgRatio)
  df$background_mapped <- as.numeric(df$background_mapped)
  df$Fold_enrichment <- df$Count/(df$mapped*(df$GO_in_background/df$background_mapped))
  GO_BP_over_rep@result <- df
  p <- barplot(GO_BP_over_rep, 
               showCategory=10, 
               x = "Fold_enrichment", 
               color = "p.adjust", 
               order=T,
               title = plot_title,
               font.size = 10) +
    theme(plot.title = element_text(size = 10)) +
    xlab("Fold enrichment")
  print(p)
}
GOBP_count <- function(background_genes, genes_of_interest, plot_title){
  GO_BP_over_rep <- enrichGO(gene = genes_of_interest,
                             universe = background_genes,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP",
                             keyType = "SYMBOL",
                             pAdjustMethod = "bonferroni", 
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)
  df <- GO_BP_over_rep@result
  df$mapped <- sub(".*/", "", df$GeneRatio)
  df$mapped <- as.numeric(df$mapped)
  df$GO_in_background <- sub("/.*", "", df$BgRatio)
  df$GO_in_background <- as.numeric(df$GO_in_background)
  df$background_mapped <- sub(".*/", "", df$BgRatio)
  df$background_mapped <- as.numeric(df$background_mapped)
  df$Fold_enrichment <- df$Count/(df$mapped*(df$GO_in_background/df$background_mapped))
  GO_BP_over_rep@result <- df
  p <- barplot(GO_BP_over_rep, 
               showCategory=10, 
               x = "Count", #could use x = "Fold_enrichment" here if desired
               color = "p.adjust", 
               order=T,
               title = plot_title,
               font.size = 10) +
    theme(plot.title = element_text(size = 10)) +
    ylab("Number of genes") #change to "Fold enrichment" if needed
  print(p)
}

## Load in data and add meta data----
C1 <- Load10X_Spatial(data.dir = "spaceranger/328-4_C1_spaceranger_count/outs",
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial",
                      slice = "Sham",
                      filter.matrix = T,
                      to.upper = F)
C1@meta.data$sample <- "328-4_C1"
C1@meta.data$timepoint <- "24 h"
C1@meta.data$surgery <- "Sham"
D1 <- Load10X_Spatial(data.dir = "spaceranger/328-5_D1_spaceranger_count/outs",
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial",
                      slice = "IR",
                      filter.matrix = T,
                      to.upper = F)
D1@meta.data$sample <- "328-5_D1"
D1@meta.data$timepoint <- "24 h"
D1@meta.data$surgery <- "IR"

## SCTransform, integration and clustering----
C1 <- SCTransform(C1, assay = "Spatial")
D1 <- SCTransform(D1, assay = "Spatial")
hearts.list <- list(C1,D1) #Using direct list instead of merge, otherwise not all genes are present in the SCT model: https://github.com/satijalab/seurat/issues/3198 
features <- SelectIntegrationFeatures(object.list = hearts.list, nfeatures = 3000)
hearts.list <- PrepSCTIntegration(object.list = hearts.list, anchor.features = features)
hearts.anchors <- FindIntegrationAnchors(object.list = hearts.list, normalization.method = "SCT", 
                                         anchor.features = features)
hearts <- IntegrateData(anchorset = hearts.anchors, normalization.method = "SCT")
DefaultAssay(hearts) <- "integrated"
hearts <- RunPCA(hearts, verbose = T)
ElbowPlot(hearts, ndims = 50)
hearts <- FindNeighbors(hearts, dims = 1:20) #Used more dimensions than usual because SCTransform allows for more
hearts <- FindClusters(hearts, verbose = T, res = 0.4) #Lowered resolution from 0.8 to 0.4
hearts <- RunUMAP(hearts, reduction = "pca", dims = 1:20)
DefaultAssay(hearts) <- "Spatial"
hearts <- NormalizeData(hearts) #Normalizing Spatial assay for DEG testing
all_genes <- rownames(hearts)
hearts <- ScaleData(hearts, features = all_genes) #Scaling Spatial assay for visualization
DefaultAssay(hearts) <- "SCT"
VariableFeatures(hearts) <- c(VariableFeatures(C1), VariableFeatures(D1)) #Have to set this after integration because I didn't use merge above
SpatialDimPlot(hearts) + DimPlot(hearts, label = T) + DimPlot(hearts, group.by = "surgery")

## Save objects----
hearts@meta.data$surgery <- factor(hearts@meta.data$surgery, levels = c("Sham", "IR"))
save(hearts, file = "hearts.Rdata")

## Load in clustered object----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

## UMI count figures----
#VlnPlot
hearts@meta.data$surgery <- factor(hearts@meta.data$surgery, levels = c("Sham", "IR"))
pdf(file = "figures/VlnPlot_UMI_count.pdf", height =6, width = 6, useDingbats = F)
VlnPlot(hearts, features = "nCount_Spatial", pt.size = 0.25, group.by = "surgery") + 
  NoLegend() + 
  ggtitle("") + 
  ylab("UMI count") + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size = 16))
dev.off()

#SpatialFeaturePlot
Sham_max <-max(hearts_24Sham[["nCount_Spatial"]])
IR_max <-max(hearts_24IR[["nCount_Spatial"]])
pdf(file = "figures/SpatialFeaturePlot_UMI_count.pdf", height = 6, width = 12.75, useDingbats = F)
a <- SpatialFeaturePlot(hearts, features = "nCount_Spatial", images = "Sham", pt.size.factor = 1.8) + 
  ggtitle("Sham - 24 h") +
  theme(legend.position = "right", 
        title = element_text(size = 16)) & 
  scale_fill_gradientn(colours=replicants, name="UMI count", limits = c(0, max(c(Sham_max, IR_max))))
b <- SpatialFeaturePlot(hearts, features = "nCount_Spatial", images = "IR", pt.size.factor = 1.8) + 
  ggtitle("I/R - 24 h") +
  theme(legend.position = "right", 
        title = element_text(size = 16)) & 
  scale_fill_gradientn(colours=replicants, name="UMI count", limits = c(0, max(c(Sham_max, IR_max))))
a + b + plot_layout(guides = "collect")
dev.off()

## DimPlot figures----
#SpatialDimPlot
pdf(file = "figures/SpatialDimPlot.pdf", height = 6, width = 12.75, useDingbats = F)
a <- SpatialDimPlot(hearts, images = "Sham", pt.size.factor = 1.8) +
  labs(fill="Cluster") +
  ggtitle("Sham - 24 h") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3,"point"),
        title = element_text(size = 16)) +
  guides(fill = guide_legend(override.aes = list(size=5))) + NoLegend()

b <- SpatialDimPlot(hearts, images = "IR", pt.size.factor = 1.8) +
  labs(fill="Cluster") +
  ggtitle("I/R - 24 h") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3,"point"),
        title = element_text(size = 16)) +
  guides(fill = guide_legend(override.aes = list(size=5)))
a + b + plot_layout(guides = "collect")
dev.off()

#Normal DimPlot
pdf(file = "figures/DimPlot.pdf", height = 6, width = 8, useDingbats = F)
p <- DimPlot(hearts, pt.size = 1.5) + 
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color="Cluster") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3,"point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size=5)))
LabelClusters(plot = p, id = "ident", 
              repel = F, 
              box = T, 
              fill = alpha("white", 0.25),
              size = 8,
              label.r = unit(0.75, "lines"),
              label.size = 0)
dev.off()

#DimPlot grouped by surgery
pdf(file = "figures/DimPlot_surgery.pdf", height = 6, width = 8, useDingbats = F)
DimPlot(hearts, pt.size = 1.5, group.by = "surgery") + 
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color="Surgery") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.justification = "top",
        legend.key.size = unit(3,"point"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        plot.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5)))
dev.off()

## Single cluster highlight figures----
for (i in levels(Idents(hearts))) {
  coi <- default_colors[(as.numeric(i)+1)]
  a <- SpatialDimPlot(hearts, 
                 images = "Sham",
                 cells.highlight = WhichCells(hearts, idents = i),
                 cols.highlight = c(coi, "#ff00ff00"),
                 stroke = NA) + 
    scale_fill_manual(values= c(coi, "#ff00ff00"), labels=c(i, "rest")) +
    NoLegend() +
    labs(title = paste0("Cluster ", i), subtitle = "Sham - 24 h") +
    theme(plot.title = element_text(color = coi, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14))
  b <- SpatialDimPlot(hearts, 
                      images = "IR",
                      cells.highlight = WhichCells(hearts, idents = i),
                      cols.highlight = c(coi, "#ff00ff00"),
                      stroke = NA) + 
    scale_fill_manual(values= c(coi, "#ff00ff00"), labels=c(i, "rest")) +
    NoLegend() +
    labs(title = "", subtitle = "I/R - 24 h") +
    theme(plot.title = element_text(color = coi, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14))
  pdf(file = paste0("figures/SpatialDimPlot_Cluster_", i, ".pdf"), 
      height = 6, 
      width = 11.25, 
      useDingbats = F)
  print(a + b)
  dev.off()
}
## Cluster markers----
#Finding markers
cluster_markers <- FindAllMarkers(hearts, assay = "SCT", only.pos = T)
top_10_cluster_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_100_cluster_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(cluster_markers, file = "cluster_markers/cluster_markers.csv", row.names = F)
write.csv(top_10_cluster_markers, file = "cluster_markers/top_10_cluster_markers.csv", row.names = F)
write.csv(top_100_cluster_markers, file = "cluster_markers/top_100_cluster_markers.csv", row.names = F)

#SpatialFeaturePlots of top 10 cluster markers
GOI <- top_10_cluster_markers$gene
for (i in GOI) {
  pdf(file = paste0("cluster_markers/SpatialFeaturePlot_", i,".pdf"), height = 6, width = 12.75, useDingbats = F)
  Sham_max <- max(GetAssayData(object = hearts_24Sham, slot = "data")[i,])
  IR_max <- max(GetAssayData(object = hearts_24IR, slot = "data")[i,])
  a <- SpatialFeaturePlot(hearts, features = i, images = "Sham", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
    ggtitle("Sham - 24 h") +
    theme(legend.position = "right", 
          legend.title = element_text(face = "italic", size = 16),
          title = element_text(size = 16)) & 
    scale_fill_gradientn(colours=replicants, name=i, limits = c(0, max(c(Sham_max, IR_max))))
  b <- SpatialFeaturePlot(hearts, features = i, images = "IR", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
    ggtitle("I/R - 24 h") +
    theme(legend.position = "right", 
          legend.title = element_text(face = "italic", size = 16),
          title = element_text(size = 16)) & 
    scale_fill_gradientn(colours=replicants,name=i, limits = c(0, max(c(Sham_max, IR_max))))
  print(a + b + plot_layout(guides = "collect"))
  dev.off()
}

#DotPlot
top_5_cluster_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf(file = "cluster_markers/DotPlot_top_5_cluster_markers.pdf", height = 6, width = 12, useDingbats = F)
DotPlot(hearts, features = unique(top_5_cluster_markers$gene)) + 
  RotatedAxis() +
  ylab("Cluster") +
  xlab("Marker genes") +
  theme(axis.text.x = element_text(face = "italic")) 
dev.off()

#Heatmap
pdf(file = "cluster_markers/Heatmap_top_30_cluster_markers.pdf", height = 8, width = 6, useDingbats = F)
top_30_cluster_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
heatmap_genes <- unique(top_30_cluster_markers$gene)
DoHeatmap(hearts, 
          features = heatmap_genes, 
          assay = "Spatial",
          angle = 0,
          size = 3.5) +
  ylab("Marker genes") +
  scale_color_manual(values = default_colors) +
  theme(axis.text.y = element_blank()) +
  labs(color="Cluster") 
dev.off()

#Gene ontology biological processes overrepresentation test
for (i in unique(cluster_markers$cluster)) {
  GOI <- cluster_markers[cluster_markers$cluster == i,]$gene
  p <- GOBP_fold_enrich(rownames(hearts), GOI, paste0("GO Biological Processes Overrepresentation: Cluster ", i, " Markers"))
  pdf(file = paste0("cluster_markers/GOBP_Overrep_Cluster_", i, ".pdf"), height = 6, width = 10, useDingbats = F)
  print(p)
  dev.off()
}
