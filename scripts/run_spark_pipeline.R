# ============================================================================
# SPARK — Full Pipeline Runner
# ============================================================================
#
# Description:
# Executes the complete SPARK workflow:
#
#   1. Multimodal pathway filtering
#   2. Graph-based module discovery
#
# Required Input:
#   example_data/GSVA_scores.csv
#
# Outputs:
#   example_output/
#
# ============================================================================

# ============================================================================
# Setup
# ============================================================================

message("========================================")
message("Starting SPARK pipeline")
message("========================================")

pipeline_start_time <- Sys.time()

# ============================================================================
# Check Input Files
# ============================================================================

required_input <- "example_data/GSVA_scores.csv"

if (!file.exists(required_input)) {

  stop(
    paste(
      "Required input file not found:",
      required_input
    )
  )

}

# ============================================================================
# Create Output Directories
# ============================================================================

dir.create(
  "example_output",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================================
# Run GMM Multimodal Filtering
# ============================================================================

message("----------------------------------------")
message("Step 1: Multimodal pathway filtering")
message("----------------------------------------")

source(
  "scripts/01_gmm_multimodal_filtering.R"
)

message("Multimodal pathway filtering completed.")

# ============================================================================
# Verify Filtered GSVA Output
# ============================================================================

filtered_matrix <- file.path(
  "example_output",
  "gmm_filtering",
  "filtered_gsva_scores.csv"
)

if (!file.exists(filtered_matrix)) {

  stop(
    "Filtered GSVA matrix was not generated."
  )

}

# ============================================================================
# Run Graph-based Module Discovery
# ============================================================================

message("----------------------------------------")
message("Step 2: Graph-based module discovery")
message("----------------------------------------")

source(
  "scripts/02_graph_module_discovery.R"
)

message("Graph-based module discovery completed.")

# ============================================================================
# Final Summary
# ============================================================================

pipeline_end_time <- Sys.time()

runtime_minutes <- round(
  as.numeric(
    difftime(
      pipeline_end_time,
      pipeline_start_time,
      units = "mins"
    )
  ),
  2
)

message("========================================")
message("SPARK pipeline completed successfully")
message(
  paste(
    "Total runtime:",
    runtime_minutes,
    "minutes"
  )
)
message("Outputs saved in: example_output/")
message("========================================")
