## Libraries, functions and colors----
library(Seurat) # v4.0.1
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
replicants <- c(
  "#5E4FA2",
  "#388FBA",
  "#A7DBA4",
  "#E0F298",
  "#F5FBB0",
  "#FDD27F",
  "#F78850",
  "#D8434D",
  "#AB1044"
)
source("scripts/GOBP_fold_enrich.R")
