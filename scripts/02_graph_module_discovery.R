# ============================================================================
# SPARK — Graph-based Module Discovery
# ============================================================================
#
# Description:
# Construction of pathway co-activity networks and identification
# of stable transcriptomic pathway modules using Leiden community
# detection with resolution optimisation.
#
# Inputs:
#   example_output/gmm_filtering/filtered_gsva_scores.csv
#
# Outputs:
#   example_output/module_discovery/
#
# ============================================================================

# ============================================================================
# Load Libraries
# ============================================================================

suppressPackageStartupMessages({

  library(tidyverse)
  library(WGCNA)
  library(igraph)
  library(mclust)
  library(purrr)
  library(ComplexHeatmap)
  library(circlize)

})

# ============================================================================
# Setup
# ============================================================================

output_dir <- "example_output/module_discovery"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================================
# Load Filtered GSVA Matrix
# ============================================================================

gsva_scores <- read.csv(
  "example_output/gmm_filtering/filtered_gsva_scores.csv",
  row.names = 1,
  check.names = FALSE
)

datExpr <- as.data.frame(
  t(gsva_scores)
)

message(
  sprintf(
    "Loaded filtered GSVA matrix: %d pathways x %d samples",
    ncol(datExpr),
    nrow(datExpr)
  )
)

# ============================================================================
# Soft-threshold Selection
# ============================================================================

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 0
)

fit_df <- sft$fitIndices

optimal_power <- fit_df$Power[
  which(fit_df$SFT.R.sq > 0.8)[1]
]

if (is.na(optimal_power)) {

  optimal_power <- fit_df$Power[
    which.max(fit_df$SFT.R.sq)
  ]

}

message(
  paste(
    "Selected soft-threshold power:",
    optimal_power
  )
)

# ============================================================================
# Adjacency Construction
# ============================================================================

adjacency_matrix <- adjacency(
  datExpr,
  power = optimal_power,
  type = "signed"
)

graph <- graph_from_adjacency_matrix(
  adjacency_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

message(
  sprintf(
    "Constructed graph with %d nodes and %d edges",
    vcount(graph),
    ecount(graph)
  )
)

# ============================================================================
# Resolution-wise Leiden Clustering
# ============================================================================

resolution_grid <- seq(
  0.5,
  2,
  by = 0.1
)

n_boot <- 100

results_df <- map_dfr(

  resolution_grid,

  function(resolution) {

    message(
      paste(
        "Testing resolution:",
        resolution
      )
    )

    leiden_res <- cluster_leiden(
      graph,
      resolution = resolution,
      objective_function = "modularity",
      n_iterations = -1
    )

    membership <- setNames(
      leiden_res$membership,
      leiden_res$names
    )

    module_df <- tibble(
      Pathway = names(membership),
      Module = paste0(
        "Module_",
        membership
      )
    ) %>%
      group_by(Module) %>%
      filter(n() >= 2) %>%
      ungroup()

    num_modules <- length(
      unique(module_df$Module)
    )

    # ------------------------------------------------------------------------
    # MAPC Calculation
    # ------------------------------------------------------------------------

    mapc_df <- module_df %>%

      group_by(Module) %>%

      summarise(
        pathways = list(Pathway),
        n_pathways = n(),
        .groups = "drop"
      ) %>%

      mutate(

        MAPC = map_dbl(

          pathways,

          function(pathways) {

            if (length(pathways) < 2) {
              return(NA)
            }

            sub_cor <- cor(
              datExpr[, pathways]
            )

            mean(
              abs(
                sub_cor[
                  lower.tri(sub_cor)
                ]
              )
            )

          }

        )

      )

    weighted_mapc <- weighted.mean(
      mapc_df$MAPC,
      mapc_df$n_pathways,
      na.rm = TRUE
    )

    # ------------------------------------------------------------------------
    # Bootstrap Stability
    # ------------------------------------------------------------------------

    ari_scores <- numeric(n_boot)

    for (b in seq_len(n_boot)) {

      boot_pathways <- sample(
        colnames(datExpr),
        size = floor(
          0.8 * ncol(datExpr)
        )
      )

      graph_boot <- induced_subgraph(
        graph,
        vids = boot_pathways
      )

      leiden_boot <- cluster_leiden(
        graph_boot,
        resolution = resolution,
        objective_function = "modularity",
        n_iterations = -1
      )

      overlap <- intersect(
        boot_pathways,
        module_df$Pathway
      )

      if (length(overlap) > 2) {

        ari_scores[b] <- adjustedRandIndex(

          module_df$Module[
            match(overlap, module_df$Pathway)
          ],

          paste0(
            "Module_",
            leiden_boot$membership[
              match(overlap, leiden_boot$names)
            ]
          )

        )

      }

    }

    mean_ari <- mean(
      ari_scores,
      na.rm = TRUE
    )

    tibble(
      Resolution = resolution,
      Num_Modules = num_modules,
      Weighted_MAPC = weighted_mapc,
      Stability = mean_ari
    )

  }

)

# ============================================================================
# Composite Score Optimisation
# ============================================================================

scale_metric <- function(x) {

  rng <- range(
    x,
    na.rm = TRUE
  )

  if (diff(rng) == 0) {

    return(
      rep(
        0.5,
        length(x)
      )
    )

  }

  (x - rng[1]) / diff(rng)

}

results_df <- results_df %>%

  mutate(
    MAPC_Norm = scale_metric(Weighted_MAPC),
    Stability_Norm = scale_metric(Stability),
    Composite_Score = (
      MAPC_Norm +
      Stability_Norm
    ) / 2
  )

best_resolution <- results_df$Resolution[
  which.max(results_df$Composite_Score)
]

message(
  paste(
    "Optimal resolution selected:",
    best_resolution
  )
)

write.csv(
  results_df,
  file.path(
    output_dir,
    "resolution_metrics.csv"
  ),
  row.names = FALSE
)

# ============================================================================
# Final Leiden Clustering
# ============================================================================

final_leiden <- cluster_leiden(
  graph,
  resolution = best_resolution,
  objective_function = "modularity",
  n_iterations = -1
)

final_membership <- setNames(
  final_leiden$membership,
  final_leiden$names
)

module_assignments <- tibble(
  Pathway = names(final_membership),
  Module = paste0(
    "Module_",
    final_membership
  )
) %>%
  group_by(Module) %>%
  filter(n() >= 2) %>%
  ungroup()

write.csv(
  module_assignments,
  file.path(
    output_dir,
    "module_assignments.csv"
  ),
  row.names = FALSE
)

# ============================================================================
# Module Score Generation
# ============================================================================

module_scores <- list()

pc1_loadings <- list()

for (module_name in unique(module_assignments$Module)) {

  module_pathways <- module_assignments$Pathway[
    module_assignments$Module == module_name
  ]

  module_matrix <- datExpr[
    ,
    module_pathways,
    drop = FALSE
  ]

  if (ncol(module_matrix) == 1) {

    score <- module_matrix[, 1]

    loading_df <- tibble(
      Pathway = module_pathways,
      Module = module_name,
      PC1_Loading = 1
    )

  } else {

    pca <- prcomp(
      module_matrix,
      scale. = TRUE
    )

    score <- pca$x[, 1]

    loading_df <- tibble(
      Pathway = rownames(pca$rotation),
      Module = module_name,
      PC1_Loading = pca$rotation[, 1]
    )

  }

  module_scores[[module_name]] <- score

  pc1_loadings[[module_name]] <- loading_df

}

module_scores_df <- as.data.frame(
  module_scores
)

pc1_loadings_df <- bind_rows(
  pc1_loadings
)

write.csv(
  module_scores_df,
  file.path(
    output_dir,
    "module_scores.csv"
  )
)

write.csv(
  pc1_loadings_df,
  file.path(
    output_dir,
    "module_pc1_loadings.csv"
  ),
  row.names = FALSE
)

# ============================================================================
# Resolution Optimisation Plot
# ============================================================================

p_resolution <- ggplot(
  results_df,
  aes(
    x = Resolution,
    y = Composite_Score
  )
) +

  geom_line(
    linewidth = 1,
    color = "#1F77B4"
  ) +

  geom_point(
    size = 2,
    color = "#D62728"
  ) +

  geom_vline(
    xintercept = best_resolution,
    linetype = "dashed"
  ) +

  theme_bw() +

  labs(
    title = "Composite Score Across Leiden Resolutions",
    x = "Resolution",
    y = "Composite Score"
  )

ggsave(
  file.path(
    output_dir,
    "resolution_optimization.svg"
  ),
  p_resolution,
  width = 7,
  height = 5,
  bg = "white"
)

# ============================================================================
# Module Correlation Heatmap
# ============================================================================

if (ncol(module_scores_df) > 1) {

  module_cor <- cor(
    module_scores_df,
    use = "pairwise.complete.obs"
  )

  color_fun <- colorRamp2(
    c(-1, 0, 1),
    c("blue", "white", "red")
  )

  heatmap_plot <- Heatmap(
    module_cor,
    name = "Correlation",
    col = color_fun
  )

  svg(
    file.path(
      output_dir,
      "module_correlation_heatmap.svg"
    ),
    width = 7,
    height = 7
  )

  draw(heatmap_plot)

  dev.off()

}

message("SPARK graph-based module discovery completed.")