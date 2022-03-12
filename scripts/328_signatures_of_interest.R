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

## Load in clustered object, score with gene signatures, then split----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
gene_ontology <- read_excel("gene_signatures/gene_signatures.xlsx", sheet = "Gene_ontology")
for (i in colnames(gene_ontology)) {
  genes <- gene_ontology[i] 
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts, features = list(genes), name = i, assay = "Spatial")
}
signatures <- tail(colnames(hearts@meta.data), length(colnames(gene_ontology)))
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"



length(gene_ontology$`lipid_catabolic_process_GO:0016042`)
lip_cata <- unique(gene_ontology$`lipid_catabolic_process_GO:0016042`)
hearts <- AddModuleScore(hearts, features = list(lip_cata), name = "lip_cata", assay = "Spatial")
VlnPlot(hearts, features = "lip_cata1", group.by = "seurat_clusters", split.by = "surgery")




cluster2_degs <- cluster_deg_IR_v_Sham[cluster_deg_IR_v_Sham$cluster=="2",]
intersect(lip_cata, cluster2_degs)

?AddModuleScore
## SpatialFeature Plot of signatures----
for (i in signatures) {
  Sham_max <-max(hearts_24Sham[[i]])
  IR_max <-max(hearts_24IR[[i]])
  pdf(file = paste0("gene_signatures/SpatialFeaturePlot_", i, ".pdf"), height = 7, width = 12.5, useDingbats = F)
  
  suppressMessages(
  a <- SpatialFeaturePlot(hearts, features = i, images = "Sham", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
    ggtitle("Sham - 24 h") +
    theme(legend.position = "bottom", 
          title = element_text(size = 16)) & 
    scale_fill_gradientn(colours=replicants, name = i, limits = c(0, max(c(Sham_max, IR_max))))
  )
  
  suppressMessages(
  b <- SpatialFeaturePlot(hearts, features = i, images = "IR", pt.size.factor = 1.8, alpha = c(0.01,1)) + 
   ggtitle("I/R - 24 h") +
   theme(legend.position = "bottom", 
          title = element_text(size = 16)) & 
    scale_fill_gradientn(colours=replicants, name = i, limits = c(0, max(c(Sham_max, IR_max))))
  )
  
  print(a + b + plot_layout(guides = "collect") & theme(legend.position = 'bottom') )
  dev.off()
  print(paste0(i, " - done"))
}