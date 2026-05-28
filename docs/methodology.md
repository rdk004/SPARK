# SPARK Methodology

## Overview

SPARK (Stability-optimised Program Architecture Reconstruction Framework)
is a graph-based systems biology framework designed to identify stable,
coordinated transcriptomic pathway modules from pathway activity matrices.

## Workflow

### 1. GSVA Pathway Activity Generation
Pathway activity scores are generated using GSVA or equivalent methods.

### 2. Multimodal Pathway Filtering
Gaussian Mixture Models (GMMs) are fitted to each pathway distribution.
Pathways with ΔBIC > 10 are retained for downstream analysis.

### 3. Correlation Network Construction
A signed pathway co-activity network is constructed using WGCNA-style adjacency transformation.

### 4. Leiden Community Detection
Pathway modules are identified using Leiden graph partitioning across multiple resolutions.

### 5. Stability Optimisation
Bootstrap-adjusted Rand Index (ARI) and pathway coherence metrics are combined to identify the optimal graph resolution.

### 6. Module Score Generation
Principal component analysis (PCA) is used to derive module-level transcriptomic scores.

## Outputs

SPARK generates:
- retained multimodal pathways
- graph-derived pathway modules
- module score matrices
- pathway loading tables
- diagnostic optimisation plots
