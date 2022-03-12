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

## Load in clustered object----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

## Detected background genes for overrepresentation testing reference----
detected_genes <- rownames(hearts)
write.table(detected_genes, file = "detected_genes/detected_genes.txt", row.names = F, col.names = F, quote = F)

## Differential gene expression, each cluster, Sham v IR----

#DEG testing 0.01 significance threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
DEgenes <- list()
for (i in levels(Idents(hearts))) {
  results <- FindMarkers(hearts, 
                         subset.ident = i, 
                         ident.1 = "IR", 
                         group.by = "surgery", 
                         assay = "Spatial")
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  results <- results[results$p_val_adj < 0.01, ]
  DEgenes[[i]] <- results
}
cluster_deg_IR_v_Sham <- do.call(rbind, DEgenes)
rownames(cluster_deg_IR_v_Sham) <- NULL
write.csv(cluster_deg_IR_v_Sham, file = "differential_gene_expression/cluster_deg_IR_v_Sham.csv", row.names = F)
remove(DEgenes)
DefaultAssay(hearts) <- "SCT"

#No significance threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
DEgenes <- list()
for (i in levels(Idents(hearts))) {
  results <- FindMarkers(hearts, 
                         subset.ident = i, 
                         ident.1 = "IR", 
                         group.by = "surgery", 
                         assay = "Spatial", 
                         logfc.threshold = 0)
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  DEgenes[[i]] <- results
}
cluster_deg_IR_v_Sham_no_threshold <- do.call(rbind, DEgenes)
rownames(cluster_deg_IR_v_Sham_no_threshold) <- NULL
write.csv(cluster_deg_IR_v_Sham_no_threshold, file = "differential_gene_expression/cluster_deg_IR_v_Sham_no_threshold.csv", row.names = F)
remove(DEgenes)
DefaultAssay(hearts) <- "SCT"

#Volcano plot
df <- cluster_deg_IR_v_Sham_no_threshold
for (i in unique(df$cluster)) {
  deg <- df[df$cluster == i,]
  deg$neglog10p <- -(log10(deg$p_val_adj))
  deg_up <- deg[deg$regulation == "Up",]
  deg_up_sig <- deg_up[deg_up$p_val_adj < 0.001,]
  deg_up_sig <- deg_up_sig[deg_up_sig$avg_log2FC > 0.25,]
  deg_up_top20 <- deg_up_sig[1:20,]
  deg_down <- deg[deg$regulation == "Down",]
  deg_down_sig <- deg_down[deg_down$p_val_adj < 0.001,]
  deg_down_sig <- deg_down_sig[deg_down_sig$avg_log2FC < -0.25,]
  deg_down_top20 <- deg_down_sig[1:20,]
  p <- ggplot(deg,
              aes(x = avg_log2FC, 
                  y = neglog10p)) 
  e <- p + geom_point(colour = "grey", 
                      size = 1) +
    xlab("avg_log2FC") +
    ylab("-log10(adjusted p-value)") +
    xlim(-6, 6) +
    geom_hline(yintercept = -log10(0.001), 
               linetype="dashed", 
               color = "grey", 
               size = 0.1) +
    geom_vline(xintercept = -0.25, 
               linetype="dashed", 
               color = "grey", 
               size = 0.1) +
    geom_vline(xintercept = 0.25, 
               linetype="dashed", 
               color = "grey", 
               size = 0.1) +
    geom_text_repel(data = deg_up_top20, 
                    aes(label=gene), 
                    colour = "#eb2a0e", 
                    fontface = "italic", 
                    size = 3,
                    box.padding = 0.5, 
                    nudge_y = 8, 
                    nudge_x = 1, 
                    segment.alpha = .3, 
                    force = 5,
                    max.overlaps = 20,
                    segment.size = 0.2) +
    geom_point(data = deg_up_sig, 
               colour="#eb2a0e",
               size = 1) + 
    geom_text_repel(data = deg_down_top20, 
                    aes(label=gene), 
                    colour = "#2664ad", 
                    fontface = "italic", 
                    size = 3,
                    box.padding = 0.5, 
                    nudge_y = 2, 
                    nudge_x = -1, 
                    segment.alpha = .3, 
                    force = 5,
                    max.overlaps = 20,
                    segment.size = 0.2) +
    geom_point(data = deg_down_sig, 
               colour="#2664ad",
               size = 1) +
    theme_bw() + 
    theme(panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14)) +
    ggtitle(paste0("Differentially expressed genes - Cluster ", i, ": Regulation in IR"))
  pdf(file = paste0("differential_gene_expression/VolcanoPlot_cluster_", i, "_deg_Sham_vs_IR.pdf"), 
      height = 6, 
      width = 8, 
      useDingbats = F)
  print(e)
  dev.off()
}


## Differential gene expression, remote clusters, Sham v IR----

#DEG testing 0.01 significance threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
remote_deg <- FindMarkers(hearts,
                          subset.ident = c("2", "3", "4"),
                          ident.1 = "IR",
                          group.by = "surgery",
                          assay = "Spatial")
remote_deg$gene <- rownames(remote_deg)
remote_deg$regulation <- ifelse(remote_deg$avg_log2FC > 0, "Up", "Down")
remote_deg <- remote_deg[remote_deg$p_val_adj < 0.01, ]
rownames(remote_deg) <- NULL
write.csv(remote_deg, file = "differential_gene_expression/remote_clusters/remote_deg.csv", row.names = F)
DefaultAssay(hearts) <- "SCT"

#DEG testing no threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
remote_deg_no_threshold <- FindMarkers(hearts,
                          subset.ident = c("2", "3", "4"),
                          ident.1 = "IR",
                          group.by = "surgery",
                          assay = "Spatial",
                          logfc.threshold = 0)
remote_deg_no_threshold$gene <- rownames(remote_deg_no_threshold)
remote_deg_no_threshold$regulation <- ifelse(remote_deg_no_threshold$avg_log2FC > 0, "Up", "Down")
rownames(remote_deg_no_threshold) <- NULL
write.csv(remote_deg_no_threshold, file = "differential_gene_expression/remote_clusters/remote_deg_no_threshold.csv", row.names = F)
DefaultAssay(hearts) <- "SCT"

#Volcano plot
deg <- remote_deg_no_threshold
deg$neglog10p <- -(log10(deg$p_val_adj))
deg_up <- deg[deg$regulation == "Up",]
deg_up_sig <- deg_up[deg_up$p_val_adj < 0.001,]
deg_up_sig <- deg_up_sig[deg_up_sig$avg_log2FC > 0.25,]
deg_up_top20 <- deg_up_sig[1:20,]
deg_down <- deg[deg$regulation == "Down",]
deg_down_sig <- deg_down[deg_down$p_val_adj < 0.001,]
deg_down_sig <- deg_down_sig[deg_down_sig$avg_log2FC < -0.25,]
deg_down_top20 <- deg_down_sig[1:20,]
p <- ggplot(deg,
            aes(x = avg_log2FC, 
                 y = neglog10p)) 
e <- p + geom_point(colour = "grey", 
                    size = 1) +
  xlab("avg_log2FC") +
  ylab("-log10(adjusted p-value)") +
  xlim(-6, 6) +
  geom_hline(yintercept = -log10(0.001), 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_vline(xintercept = -0.25, 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_vline(xintercept = 0.25, 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_text_repel(data = deg_up_top20, 
                  aes(label=gene), 
                  colour = "#eb2a0e", 
                  fontface = "italic", 
                  size = 3,
                  box.padding = 0.5, 
                  nudge_y = 8, 
                  nudge_x = 1, 
                  segment.alpha = .3, 
                  force = 5,
                  max.overlaps = 20,
                  segment.size = 0.2) +
  geom_point(data = deg_up_sig, 
             colour="#eb2a0e",
             size = 1) + 
  geom_text_repel(data = deg_down_top20, 
                  aes(label=gene), 
                  colour = "#2664ad", 
                  fontface = "italic", 
                  size = 3,
                  box.padding = 0.5, 
                  nudge_y = 2, 
                  nudge_x = -1, 
                  segment.alpha = .3, 
                  force = 5,
                  max.overlaps = 20,
                  segment.size = 0.2) +
  geom_point(data = deg_down_sig, 
             colour="#2664ad",
             size = 1) +
  theme_bw() + 
  theme(panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14)) +
  ggtitle(paste0("Differentially expressed genes - Remote zone clusters : Regulation in IR"))
pdf(file = paste0("differential_gene_expression/remote_clusters/VolcanoPlot_remote_clusters_deg_Sham_vs_IR.pdf"),
    height = 6, 
    width = 8, 
    useDingbats = F)
print(e)
dev.off()
## Differential gene expression, remote clusters (2,4), Sham v IR----

#DEG testing 0.01 significance threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
remote_deg <- FindMarkers(hearts,
                          subset.ident = c("2", "4"),
                          ident.1 = "IR",
                          group.by = "surgery",
                          assay = "Spatial")
remote_deg$gene <- rownames(remote_deg)
remote_deg$regulation <- ifelse(remote_deg$avg_log2FC > 0, "Up", "Down")
remote_deg <- remote_deg[remote_deg$p_val_adj < 0.01, ]
rownames(remote_deg) <- NULL
write.csv(remote_deg, file = "differential_gene_expression/remote_clusters_2_4/remote_deg.csv", row.names = F)
DefaultAssay(hearts) <- "SCT"

#DEG testing no threshold
Idents(hearts) <- "seurat_clusters"
DefaultAssay(hearts) <- "Spatial"
remote_deg_no_threshold <- FindMarkers(hearts,
                                       subset.ident = c("2", "4"),
                                       ident.1 = "IR",
                                       group.by = "surgery",
                                       assay = "Spatial",
                                       logfc.threshold = 0)
remote_deg_no_threshold$gene <- rownames(remote_deg_no_threshold)
remote_deg_no_threshold$regulation <- ifelse(remote_deg_no_threshold$avg_log2FC > 0, "Up", "Down")
rownames(remote_deg_no_threshold) <- NULL
write.csv(remote_deg_no_threshold, file = "differential_gene_expression/remote_clusters_2_4/remote_deg_no_threshold.csv", row.names = F)
DefaultAssay(hearts) <- "SCT"

#Volcano plot
deg <- remote_deg_no_threshold
deg$neglog10p <- -(log10(deg$p_val_adj))
deg_up <- deg[deg$regulation == "Up",]
deg_up_sig <- deg_up[deg_up$p_val_adj < 0.001,]
deg_up_sig <- deg_up_sig[deg_up_sig$avg_log2FC > 0.25,]
deg_up_top20 <- deg_up_sig[1:20,]
deg_down <- deg[deg$regulation == "Down",]
deg_down_sig <- deg_down[deg_down$p_val_adj < 0.001,]
deg_down_sig <- deg_down_sig[deg_down_sig$avg_log2FC < -0.25,]
deg_down_top20 <- deg_down_sig[1:20,]
p <- ggplot(deg,
            aes(x = avg_log2FC, 
                y = neglog10p)) 
e <- p + geom_point(colour = "grey", 
                    size = 1) +
  xlab("avg_log2FC") +
  ylab("-log10(adjusted p-value)") +
  xlim(-6, 6) +
  geom_hline(yintercept = -log10(0.001), 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_vline(xintercept = -0.25, 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_vline(xintercept = 0.25, 
             linetype="dashed", 
             color = "grey", 
             size = 0.1) +
  geom_text_repel(data = deg_up_top20, 
                  aes(label=gene), 
                  colour = "#eb2a0e", 
                  fontface = "italic", 
                  size = 3,
                  box.padding = 0.5, 
                  nudge_y = 8, 
                  nudge_x = 1, 
                  segment.alpha = .3, 
                  force = 5,
                  max.overlaps = 20,
                  segment.size = 0.2) +
  geom_point(data = deg_up_sig, 
             colour="#eb2a0e",
             size = 1) + 
  geom_text_repel(data = deg_down_top20, 
                  aes(label=gene), 
                  colour = "#2664ad", 
                  fontface = "italic", 
                  size = 3,
                  box.padding = 0.5, 
                  nudge_y = 2, 
                  nudge_x = -1, 
                  segment.alpha = .3, 
                  force = 5,
                  max.overlaps = 20,
                  segment.size = 0.2) +
  geom_point(data = deg_down_sig, 
             colour="#2664ad",
             size = 1) +
  theme_bw() + 
  theme(panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14)) +
  ggtitle(paste0("Differentially expressed genes - Remote zone clusters : Regulation in IR"))
pdf(file = paste0("differential_gene_expression/remote_clusters_2_4/VolcanoPlot_remote_clusters_2_4_deg_Sham_vs_IR.pdf"),
    height = 6, 
    width = 8, 
    useDingbats = F)
print(e)
dev.off()