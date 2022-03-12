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

## Load in clustered object, combining clusters 2+4----
load("hearts.Rdata")
Idents(hearts) <- "seurat_clusters"
hearts <- RenameIdents(hearts,
                       "0" = "0",
                       "1" = "1",
                       "2" = "Remote",
                       "3" = "3",
                       "4" = "Remote",
                       "5" = "5",
                       "6" = "6")
hearts@meta.data$remote_combination_2_4_annotation <- Idents(hearts)

## Genes of interest----
GOI <- c("Fabp3", "Dgat2", "Cpt2")

## VlnPlots----
for (i in GOI) {
  pdf(file = paste0("remote_combination_2+4/VlnPlot_", i,".pdf"), height = 6, width = 8, useDingbats = F)
  p <- VlnPlot(hearts, 
               features = i, 
               idents = "Remote",
               split.by = "surgery",
               assay = "Spatial",
               split.plot=F) +
    theme(plot.title = element_text(face = "italic", size = 20),
          axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "right")
  print(p)
  dev.off()
  print(paste0(i, " - done"))
}

## DEG testing 0.01 significance threshold remote zone 2+4 Sham v IR----
Idents(hearts) <- "remote_combination_2_4_annotation"
results <- FindMarkers(hearts, 
            assay = "Spatial",
            subset.ident = "Remote", 
            group.by = "surgery",
            ident.1 = "IR")
results$gene <- row.names(results)
row.names(results) <- NULL
results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
results <- results[results$p_val_adj < 0.01, ]
remote_2_4_deg_IR_v_Sham <- results
write.xlsx(remote_2_4_deg_IR_v_Sham, file = "remote_combination_2+4/remote_2_4_deg_IR_v_Sham.xlsx", row.names = F)

## DEG testing no significance threshold remote zone 2+4 Sham v IR----
Idents(hearts) <- "remote_combination_2_4_annotation"
results <- FindMarkers(hearts, 
                       assay = "Spatial",
                       subset.ident = "Remote", 
                       group.by = "surgery",
                       ident.1 = "IR",
                       logfc.threshold = 0)
results$gene <- row.names(results)
row.names(results) <- NULL
results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
remote_2_4_deg_IR_v_Sham_no_threshold <- results
write.xlsx(remote_2_4_deg_IR_v_Sham_no_threshold, file = "remote_combination_2+4/remote_2_4_deg_IR_v_Sham_no_threshold.xlsx", row.names = F)

## Remote zone highlight----
coi <- "#F4BA27"
i <- "Remote"
a <- SpatialDimPlot(hearts, 
                      images = "Sham",
                      cells.highlight = WhichCells(hearts, idents = i),
                      cols.highlight = c(coi, "#ff00ff00"),
                      stroke = NA) + 
    scale_fill_manual(values= c(coi, "#ff00ff00"), labels=c(i, "rest")) +
    NoLegend() +
    labs(title = "Remote zone clusters 2+4", subtitle = "Sham - 24 h") +
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
pdf(file = "remote_combination_2+4/SpatialDimPlot_remote_zone_clusters_2+4.pdf", 
      height = 6, 
      width = 11.25, 
      useDingbats = F)
print(a + b)
dev.off()



## Lipid catabolic process----
#Scoring
gene_ontology <- read_excel("gene_signatures/gene_signatures.xlsx", sheet = "Gene_ontology")
gene_ontology$`lipid_catabolic_process_GO:0016042`
hearts <- AddModuleScore(hearts, features = list(unique(gene_ontology$`lipid_catabolic_process_GO:0016042`)), 
                         name = "lipid_catabolic_process_GO", 
                         assay = "Spatial")
#Wilcoxon rank sum test
data <- hearts@meta.data
ggplot(data, aes(x=surgery, y=lipid_catabolic_process_GO1)) + geom_boxplot() 
res <- wilcox.test(lipid_catabolic_process_GO1 ~ surgery, data = data)
res 
#VlnPlot
i <- "lipid_catabolic_process_GO1"
pdf(file = paste0("remote_combination_2+4/VlnPlot_", i,".pdf"), height = 6, width = 8, useDingbats = F)
p <- VlnPlot(hearts, 
               features = i, 
               idents = "Remote",
               split.by = "surgery",
               assay = "Spatial",
               split.plot=F) +
  ggtitle(label="Lipid catabolic process GO:0016042",
          subtitle = "Wilcoxon rank sum test p-value: < 2.2e-16") +
    theme(plot.title = element_text(size = 16),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "right") 
print(p)
dev.off()





