```r
# ============================================================================
# SPARK — Required Packages
# ============================================================================

required_packages <- c(

  "tidyverse",
  "moments",
  "mclust",
  "BiocParallel",
  "diptest",
  "umap",
  "WGCNA",
  "igraph",
  "purrr",
  "ComplexHeatmap",
  "circlize"

)

installed <- required_packages %in% installed.packages()

if (any(!installed)) {

  install.packages(
    required_packages[!installed]
  )

}

invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)

message("All SPARK dependencies loaded successfully.")
```
