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

## Load in clustered object----
load("hearts.Rdata")
Idents(hearts) <- "surgery"
hearts_24Sham <- subset(hearts, idents="Sham")
hearts_24IR <- subset(hearts, idents="IR")
Idents(hearts) <- "seurat_clusters"

## Differential proportion----
differential_proportion <- (prop.table(table(Idents(hearts), hearts$surgery), margin = 2))
differential_proportion <- as.data.frame(differential_proportion)
colnames(differential_proportion) <- c("Cluster", "Sample", "Fraction")
write.csv(differential_proportion, file = "differential_proportion/differential_proportion.csv", row.names = F)
