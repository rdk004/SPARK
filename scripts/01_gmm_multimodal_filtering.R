# ============================================================================
# SPARK — Multimodal Pathway Filtering
# ============================================================================
#
# Description:
# Identification of biologically informative pathways exhibiting
# multimodal GSVA score distributions using Gaussian Mixture Models.
#
# Workflow:
#   1. GSVA matrix import
#   2. Variance quality control
#   3. Gaussian mixture modelling
#   4. ΔBIC-based multimodality filtering
#   5. Distribution diagnostics
#   6. Export of retained pathway matrix
#
# Input:
#   example_data/GSVA_scores.csv
#
# Output:
#   example_output/gmm_filtering/
#
#
# ============================================================================

# ============================================================================
# Load Libraries
# ============================================================================

suppressPackageStartupMessages({

  library(mclust)
  library(tidyverse)
  library(BiocParallel)
  library(ggplot2)
  library(moments)
  library(diptest)
  library(ggridges)
  library(umap)

})

# ============================================================================
# Setup
# ============================================================================

output_dir <- "example_output/gmm_filtering"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_dir, "pipeline_log.txt")

log_message <- function(text) {

  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

  message_line <- paste(timestamp, text)

  cat(message_line, "\n")

  write(
    message_line,
    file = log_file,
    append = TRUE
  )

}

log_message("SPARK multimodal pathway filtering started.")

# ============================================================================
# Load GSVA Matrix
# ============================================================================

input_file <- "example_data/GSVA_scores.csv"

gsva_scores <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

log_message(
  sprintf(
    "Loaded GSVA matrix: %d pathways x %d samples",
    nrow(gsva_scores),
    ncol(gsva_scores)
  )
)

# ============================================================================
# Variance Quality Control
# ============================================================================

log_message("Running variance quality control.")

variance_df <- tibble(
  Pathway = rownames(gsva_scores),
  Variance = apply(gsva_scores, 1, var, na.rm = TRUE)
) %>%
  arrange(desc(Variance))

write.csv(
  variance_df,
  file.path(output_dir, "variance_qc.csv"),
  row.names = FALSE
)

p_variance <- ggplot(
  variance_df,
  aes(Variance)
) +

  geom_histogram(
    bins = 50,
    fill = "#4C78A8",
    color = "white"
  ) +

  theme_bw() +

  labs(
    title = "Variance Distribution Across Pathways",
    x = "Variance",
    y = "Pathway Count"
  )

ggsave(
  file.path(output_dir, "variance_distribution.svg"),
  p_variance,
  width = 7,
  height = 5,
  bg = "white"
)

# ============================================================================
# ΔBIC Multimodality Filtering
# ============================================================================

log_message("Running Gaussian mixture modelling.")

num_cores <- max(
  1,
  parallel::detectCores() - 1
)

param <- MulticoreParam(
  workers = num_cores
)

gmm_results <- bplapply(

  seq_len(nrow(gsva_scores)),

  function(i) {

    pathway_name <- rownames(gsva_scores)[i]

    y <- as.numeric(gsva_scores[i, ])

    y <- y[is.finite(y)]

    if (length(unique(y)) < 4) {

      return(
        tibble(
          Pathway = pathway_name,
          Delta_BIC = NA,
          Components = NA
        )
      )

    }

    bic_model <- tryCatch(

      mclustBIC(
        y,
        verbose = FALSE
      ),

      error = function(e) NULL

    )

    if (is.null(bic_model)) {

      return(
        tibble(
          Pathway = pathway_name,
          Delta_BIC = NA,
          Components = NA
        )
      )

    }

    model_summary <- tryCatch(

      summary(
        bic_model,
        data = y
      ),

      error = function(e) NULL

    )

    if (is.null(model_summary)) {

      return(
        tibble(
          Pathway = pathway_name,
          Delta_BIC = NA,
          Components = NA
        )
      )

    }

    bic_best <- suppressWarnings(
      as.numeric(model_summary$bic[1])
    )

    bic_single <- if (
      "1" %in% rownames(bic_model)
    ) {
      suppressWarnings(
        max(bic_model["1", ], na.rm = TRUE)
      )
    } else {
      NA
    }

    tibble(
      Pathway = pathway_name,
      Delta_BIC = bic_best - bic_single,
      Components = model_summary$G[1]
    )

  },

  BPPARAM = param

)

gmm_df <- bind_rows(gmm_results) %>%

  filter(is.finite(Delta_BIC)) %>%

  mutate(
    Evidence = case_when(
      Delta_BIC < 2 ~ "Weak",
      Delta_BIC < 6 ~ "Positive",
      Delta_BIC < 10 ~ "Strong",
      Delta_BIC >= 10 ~ "Very Strong"
    )
  ) %>%

  arrange(desc(Delta_BIC))

retained_pathways <- gmm_df %>%
  filter(Delta_BIC > 10) %>%
  pull(Pathway)

log_message(
  sprintf(
    "Retained %d multimodal pathways.",
    length(retained_pathways)
  )
)

write.csv(
  gmm_df,
  file.path(output_dir, "bic_results.csv"),
  row.names = FALSE
)

write.csv(
  tibble(Pathway = retained_pathways),
  file.path(output_dir, "retained_pathways.csv"),
  row.names = FALSE
)

# ============================================================================
# ΔBIC Distribution Plot
# ============================================================================

plot_df <- gmm_df %>%

  mutate(
    Group = ifelse(
      Delta_BIC > 10,
      "Retained",
      "Filtered"
    )
  )

retained_n <- sum(plot_df$Delta_BIC > 10)

total_n <- nrow(plot_df)

p_bic <- ggplot(
  plot_df,
  aes(
    x = Delta_BIC,
    fill = Group
  )
) +

  geom_histogram(
    bins = 40,
    alpha = 0.95,
    color = "white"
  ) +

  geom_vline(
    xintercept = 10,
    linetype = "dashed",
    linewidth = 1,
    color = "#D62728"
  ) +

  scale_fill_manual(
    values = c(
      "Filtered" = "#9ECAE1",
      "Retained" = "#1F4E79"
    )
  ) +

  theme_classic(base_size = 13) +

  labs(
    title = "ΔBIC Distribution Across Pathways",
    x = expression(Delta*"BIC"),
    y = "Pathway Count",
    fill = NULL
  )

ggsave(
  file.path(output_dir, "bic_distribution.svg"),
  p_bic,
  width = 8,
  height = 6,
  bg = "white"
)

# ============================================================================
# Dip Test Quality Control
# ============================================================================

log_message("Running dip test quality control.")

dip_results <- bplapply(

  seq_len(nrow(gsva_scores)),

  function(i) {

    pathway_name <- rownames(gsva_scores)[i]

    y <- as.numeric(gsva_scores[i, ])

    y <- y[is.finite(y)]

    if (length(unique(y)) < 4) {

      return(
        tibble(
          Pathway = pathway_name,
          Dip_p = NA
        )
      )

    }

    dip_res <- tryCatch(

      dip.test(y),

      error = function(e) NULL

    )

    if (is.null(dip_res)) {

      return(
        tibble(
          Pathway = pathway_name,
          Dip_p = NA
        )
      )

    }

    tibble(
      Pathway = pathway_name,
      Dip_p = dip_res$p.value
    )

  },

  BPPARAM = param

)

dip_df <- bind_rows(dip_results) %>%

  mutate(
    FDR = p.adjust(
      Dip_p,
      method = "BH"
    ),
    Significant = FDR < 0.05
  )

write.csv(
  dip_df,
  file.path(output_dir, "dip_test_results.csv"),
  row.names = FALSE
)

bic_dip_df <- left_join(
  gmm_df,
  dip_df,
  by = "Pathway"
)

p_dip <- ggplot(
  bic_dip_df,
  aes(
    x = Delta_BIC,
    y = -log10(FDR),
    color = Significant
  )
) +

  geom_point(alpha = 0.7) +

  scale_color_manual(
    values = c(
      "TRUE" = "#D62728",
      "FALSE" = "grey70"
    )
  ) +

  theme_bw() +

  labs(
    title = "ΔBIC vs Dip Test Significance",
    x = expression(Delta*"BIC"),
    y = expression(-log[10]*"(FDR)")
  )

ggsave(
  file.path(output_dir, "dip_vs_bic.svg"),
  p_dip,
  width = 7,
  height = 5,
  bg = "white"
)

# ============================================================================
# Distribution Shape Diagnostics
# ============================================================================

log_message("Running distribution diagnostics.")

distribution_df <- t(

  apply(

    gsva_scores,

    1,

    function(y) {

      y <- y[is.finite(y)]

      c(
        Skewness = skewness(y),
        Kurtosis = kurtosis(y)
      )

    }

  )

) %>%

  as.data.frame()

distribution_df$Pathway <- rownames(distribution_df)

distribution_df$Group <- ifelse(
  distribution_df$Pathway %in% retained_pathways,
  "Retained",
  "Filtered"
)

write.csv(
  distribution_df,
  file.path(output_dir, "distribution_summary.csv"),
  row.names = FALSE
)

# ============================================================================
# Example Density Plots
# ============================================================================

set.seed(123)

high_bic_sample <- sample(
  retained_pathways,
  min(6, length(retained_pathways))
)

low_bic_sample <- sample(
  setdiff(
    rownames(gsva_scores),
    retained_pathways
  ),
  min(
    6,
    length(
      setdiff(
        rownames(gsva_scores),
        retained_pathways
      )
    )
  )
)

plot_pathways <- c(
  high_bic_sample,
  low_bic_sample
)

density_df <- gsva_scores[
  plot_pathways,
  ,
  drop = FALSE
] %>%

  as.data.frame() %>%

  rownames_to_column("Pathway") %>%

  pivot_longer(
    -Pathway,
    names_to = "Sample",
    values_to = "Score"
  ) %>%

  mutate(
    Group = ifelse(
      Pathway %in% high_bic_sample,
      "Retained",
      "Filtered"
    )
  )

p_density <- ggplot(
  density_df,
  aes(
    x = Score,
    fill = Group
  )
) +

  geom_density(alpha = 0.4) +

  facet_wrap(
    ~ Pathway,
    scales = "free"
  ) +

  theme_bw() +

  labs(
    title = "Representative GSVA Score Distributions",
    x = "GSVA Score",
    y = "Density"
  )

ggsave(
  file.path(output_dir, "example_density_plots.svg"),
  p_density,
  width = 14,
  height = 8,
  bg = "white"
)

# ============================================================================
# UMAP Projection
# ============================================================================

log_message("Running UMAP projection.")

if (length(retained_pathways) > 5) {

  umap_res <- umap(
    t(
      gsva_scores[
        retained_pathways,
        ,
        drop = FALSE
      ]
    )
  )

  umap_df <- as.data.frame(
    umap_res$layout
  )

  colnames(umap_df) <- c(
    "UMAP1",
    "UMAP2"
  )

  p_umap <- ggplot(
    umap_df,
    aes(
      x = UMAP1,
      y = UMAP2
    )
  ) +

    geom_point(
      alpha = 0.7,
      color = "#1F77B4"
    ) +

    theme_bw() +

    labs(
      title = "UMAP Projection Using Retained Pathways"
    )

  ggsave(
    file.path(output_dir, "patient_umap.svg"),
    p_umap,
    width = 7,
    height = 5,
    bg = "white"
  )

}

# ============================================================================
# Export Filtered GSVA Matrix
# ============================================================================

log_message("Exporting filtered GSVA matrix.")

gsva_filtered <- gsva_scores[
  retained_pathways,
  ,
  drop = FALSE
]

write.csv(
  gsva_filtered,
  file.path(output_dir, "filtered_gsva_scores.csv")
)

# ============================================================================
# Summary Statistics
# ============================================================================

summary_df <- tibble(

  Metric = c(
    "Total pathways",
    "Retained pathways",
    "Retention proportion",
    "Median retained ΔBIC",
    "Median retained skewness",
    "Median retained kurtosis"
  ),

  Value = c(
    nrow(gsva_scores),
    length(retained_pathways),
    round(
      length(retained_pathways) / nrow(gsva_scores),
      3
    ),
    median(
      gmm_df$Delta_BIC[
        gmm_df$Pathway %in% retained_pathways
      ],
      na.rm = TRUE
    ),
    median(
      distribution_df$Skewness[
        distribution_df$Pathway %in% retained_pathways
      ],
      na.rm = TRUE
    ),
    median(
      distribution_df$Kurtosis[
        distribution_df$Pathway %in% retained_pathways
      ],
      na.rm = TRUE
    )

  )

)

write.csv(
  summary_df,
  file.path(output_dir, "qc_summary.csv"),
  row.names = FALSE
)

log_message("SPARK multimodal pathway filtering completed.")
# ============================================================================
