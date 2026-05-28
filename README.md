###### # SPARK

###### 

###### \## Stability-optimised Program Architecture Reconstruction Framework

###### 

###### SPARK is a graph-based transcriptomic framework for identifying stable, coordinated pathway modules from pathway activity landscapes using multimodality filtering, network-based organisation, Leiden community detection, and stability optimisation.

###### 

###### The framework reconstructs higher-order transcriptomic programs from GSVA-derived pathway activity matrices and is designed for robust pathway-level systems biology analyses in cancer and complex biological systems.

###### 

###### \---

###### 

###### \# Workflow Overview

###### 

###### !\[SPARK Workflow](docs/workflow\_schematic.png)

###### 

###### \---

###### 

###### \# Overview

###### 

###### SPARK combines:

###### 

###### \- Gaussian mixture model (GMM)-based multimodal pathway filtering

###### \- Bayesian Information Criterion (ΔBIC)-driven pathway selection

###### \- Signed pathway co-activity network construction

###### \- WGCNA-inspired soft-thresholding

###### \- Leiden graph community detection

###### \- Bootstrap stability optimisation

###### \- Principal component-derived module eigengenes

###### 

###### to identify stable transcriptomic pathway modules representing coordinated biological programs.

###### 

###### \---

###### 

###### \# Features

###### 

###### \- Multimodal pathway filtering using Gaussian mixture models

###### \- ΔBIC-based identification of biologically informative pathways

###### \- Graph-based pathway organisation

###### \- Resolution-wise Leiden community detection

###### \- Stability-aware module optimisation

###### \- Bootstrap-adjusted Rand Index (ARI) evaluation

###### \- Principal component-based module summarisation

###### \- Publication-ready module-level outputs and diagnostics

###### 

###### \---

###### 

###### \# Repository Structure

###### 

###### ```text

###### SPARK\_repository/

###### │

###### ├── README.md

###### ├── LICENSE

###### ├── .gitignore

###### │

###### ├── scripts/

###### │   ├── 01\_gmm\_multimodal\_filtering.R

###### │   ├── 02\_graph\_module\_discovery.R

###### │   └── run\_spark\_pipeline.R

###### │

###### ├── example\_data/

###### │   └── GSVA\_scores.csv

###### │

###### ├── example\_output/

###### │   ├── gmm\_filtering/

###### │   └── module\_discovery/

###### │

###### ├── docs/

###### │   ├── workflow\_schematic.png

###### │   ├── mathematical\_overview.pdf

###### │   ├── methodology.md

###### │   └── parameter\_descriptions.md

###### │

###### └── environment/

###### &#x20;   ├── required\_packages.R

###### &#x20;   └── sessionInfo.txt

###### ```

###### 

###### \---

###### 

###### \# Installation

###### 

###### \## 1. Clone Repository

###### 

###### ```bash

###### git clone https://github.com/YOUR\_USERNAME/SPARK.git

###### cd SPARK

###### ```

###### 

###### \---

###### 

###### \## 2. Install Dependencies

###### 

###### Open R and run:

###### 

###### ```r

###### source("environment/required\_packages.R")

###### ```

###### 

###### \---

###### 

###### \# Input Data

###### 

###### SPARK requires a GSVA pathway activity matrix formatted as:

###### 

###### | Pathway | Sample1 | Sample2 | Sample3 |

###### |---|---|---|---|

###### | HALLMARK\_MYC\_TARGETS | 0.42 | -0.15 | 0.91 |

###### | HALLMARK\_E2F\_TARGETS | -0.33 | 0.71 | -0.08 |

###### 

###### \- Rows represent pathways

###### \- Columns represent samples

###### \- Values represent GSVA enrichment scores

###### 

###### Place the matrix inside:

###### 

###### ```text

###### example\_data/GSVA\_scores.csv

###### ```

###### 

###### \---

###### 

###### \# Running SPARK

###### 

###### Run the complete pipeline using:

###### 

###### ```r

###### source("scripts/run\_spark\_pipeline.R")

###### ```

###### 

###### The workflow sequentially performs:

###### 

###### 1\. Multimodal pathway filtering

###### 2\. Correlation network construction

###### 3\. Leiden-based module discovery

###### 4\. Stability optimisation

###### 5\. Module eigengene generation

###### 

###### \---

###### 

###### \# Workflow Components

###### 

###### \## 1. Multimodal Pathway Filtering

###### 

###### Pathways are evaluated using Gaussian mixture models (GMMs).

###### 

###### Model selection is performed using the Bayesian Information Criterion (BIC):

###### 

###### ```math

###### \\Delta BIC = BIC\_{best} - BIC\_{G=1}

###### ```

###### 

###### Pathways with:

###### 

###### ```math

###### \\Delta BIC > 10

###### ```

###### 

###### are retained for downstream graph analysis.

###### 

###### \---

###### 

###### \## 2. Graph-based Module Discovery

###### 

###### A signed pathway co-activity network is constructed using pathway-pathway correlations across samples.

###### 

###### Soft-thresholding is performed using WGCNA-inspired adjacency transformation:

###### 

###### ```math

###### A\_{ij} = \\left(\\frac{1 + cor(i,j)}{2}\\right)^\\beta

###### ```

###### 

###### Leiden community detection is then performed across multiple resolutions.

###### 

###### \---

###### 

###### \## 3. Stability Optimisation

###### 

###### For each resolution:

###### 

###### \- bootstrap pathway resampling is performed

###### \- bootstrap sample resampling is performed

###### \- Adjusted Rand Index (ARI) stability is computed

###### 

###### Resolution selection is based on a composite score integrating:

###### 

###### \- module coherence (MAPC)

###### \- graph stability

###### 

###### \---

###### 

###### \## 4. Module Eigengene Generation

###### 

###### For each identified module:

###### 

###### \- principal component analysis (PCA) is performed

###### \- PC1 is used as the module eigengene

###### \- pathway PC1 loadings are exported

###### 

###### These eigengenes represent coordinated transcriptomic pathway programs.

###### 

###### \---

###### 

###### \# Outputs

###### 

###### \## GMM Filtering Outputs

###### 

###### Located in:

###### 

###### ```text

###### example\_output/gmm\_filtering/

###### ```

###### 

###### Includes:

###### 

###### \- ΔBIC statistics

###### \- retained pathway lists

###### \- filtered GSVA matrices

###### \- multimodality diagnostics

###### \- UMAP projections

###### 

###### \---

###### 

###### \## Module Discovery Outputs

###### 

###### Located in:

###### 

###### ```text

###### example\_output/module\_discovery/

###### ```

###### 

###### Includes:

###### 

###### \- module assignments

###### \- module eigengene matrices

###### \- PC1 loading matrices

###### \- resolution optimisation metrics

###### \- module correlation heatmaps

###### 

###### \---

###### 

###### \# Example Output Files

###### 

###### \## GMM Filtering

###### 

###### ```text

###### bic\_results.csv

###### retained\_pathways.csv

###### filtered\_gsva\_scores.csv

###### bic\_distribution.svg

###### patient\_umap.svg

###### ```

###### 

###### \---

###### 

###### \## Module Discovery

###### 

###### ```text

###### module\_assignments.csv

###### module\_scores.csv

###### module\_pc1\_loadings.csv

###### resolution\_metrics.csv

###### resolution\_optimization.svg

###### module\_correlation\_heatmap.svg

###### ```

###### 

###### \---

###### 

###### \# Documentation

###### 

###### Additional methodological details are available in:

###### 

###### ```text

###### docs/

###### ```

###### 

###### Including:

###### 

###### \- workflow schematic

###### \- mathematical overview

###### \- parameter descriptions

###### \- methodological summaries

###### 

###### \---

###### 

###### \# Applications

###### 

###### SPARK is suitable for:

###### 

###### \- transcriptomic module discovery

###### \- pathway-level systems biology

###### \- cancer transcriptomics

###### \- pathway architecture reconstruction

###### \- tumour subtype characterisation

###### \- multimodal pathway organisation analyses

###### 

###### \---

###### 

###### \# Reproducibility

###### 

###### The repository includes:

###### 

###### \- dependency installation scripts

###### \- session information

###### \- example input data

###### \- example output files

###### 

###### for reproducible execution.

###### 

###### Environment information is provided in:

###### 

###### ```text

###### environment/sessionInfo.txt

###### ```

###### 

###### \---

###### 

###### \# Citation

###### 

###### If you use SPARK in your work, please cite the associated manuscript.

###### 

###### ```text

###### Rishabh Kulkarni et al.

###### SPARK: Stability-optimised Program Architecture Reconstruction Framework

###### (Manuscript in preparation)

###### ```

###### 

###### \---

###### 

###### \# Contact

###### 

###### Rishabh Kulkarni  

###### Indian Institute of Science Education and Research (IISER) Pune

###### 

###### \---

###### 

###### \# License

###### 

###### This project is released under the MIT License.

