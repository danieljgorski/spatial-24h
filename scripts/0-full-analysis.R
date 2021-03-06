# make directories
output_dirs <- c(
  "results",
  "results/basic-figures",
  "results/cluster-markers",
  "results/differential-abundance",
  "results/differential-gene-expression",
  "results/differential-gene-expression/remote-clusters",
  "results/genes-of-interest",
  "results/objects",
  "results/reference-mapping",
  "results/remote-zone-analysis",
  "results/signatures-of-interest",
  "results/spatially-variable-features"
)

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# run analysis in order
source("scripts/1-clustering.R")
source("scripts/2-basic-figures.R")
source("scripts/3-cluster-markers.R")
source("scripts/4-genes-of-interest.R")
source("scripts/5-signatures-of-interest.R")
source("scripts/6-spatially-variable-features.R")
source("scripts/7-differential-abundance.R")
source("scripts/8-differential-gene-expression.R")
source("scripts/9-remote-zone-analysis.R")
source("scripts/10-reference-mapping.R")
print("Full analysis complete")
