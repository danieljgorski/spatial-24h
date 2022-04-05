## Load libraries, functions and objects----
library(Seurat) # v4.0.1
load("results/objects/hearts.Rdata")

## Differential abundance----
differential_abundance <- (prop.table(table(
  Idents(hearts),
  hearts$surgery
),
margin = 2
))
differential_abundance <- as.data.frame(differential_abundance)
colnames(differential_abundance) <- c("Cluster", "Sample", "Fraction")
write.csv(differential_abundance,
  file = "results/differential-abundance/differential_abundance.csv",
  row.names = F
)
