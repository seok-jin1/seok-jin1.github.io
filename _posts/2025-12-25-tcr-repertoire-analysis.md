---
layout: post
title: "TCR/BCR Repertoire Analysis: scRepertoire and Dandelion for Immune Profiling"
date: 2025-12-25
permalink: /blog/tcr-repertoire-analysis/
published: true
categories: [tutorial]
tags:
  - immunology
  - single-cell
  - bioinformatics
  - python
  - tutorial
---

## Introduction: Why Repertoire Analysis Matters

Every adaptive immune response begins with a molecular lottery. During V(D)J recombination, developing T cells and B cells randomly rearrange variable (V), diversity (D), and joining (J) gene segments to produce unique antigen receptors --- T cell receptors (TCRs) and B cell receptors (BCRs). The combinatorial diversity of segment selection, combined with junctional diversity from random nucleotide additions and deletions at the recombination junctions, generates an estimated 10^15 to 10^18 possible TCR sequences in humans. This extraordinary diversity is what allows the immune system to recognize virtually any pathogen it encounters.

The **immune repertoire** --- the complete set of TCR and BCR sequences in an individual --- encodes the history and current state of adaptive immunity. A naive T cell carries a unique TCR. When that T cell encounters its cognate antigen, it proliferates, producing a **clonal expansion** --- a population of cells sharing the identical TCR sequence. By sequencing these receptors at single-cell resolution, we can:

- **Track clonal expansions** that mark active immune responses against tumors or infections
- **Measure repertoire diversity** to assess immune competence or reconstitution after transplant
- **Identify shared clonotypes** across tissues, revealing immune cell trafficking
- **Link receptor sequences to transcriptomic phenotypes**, connecting what a cell recognizes to what it does

With the advent of 10x Genomics 5' V(D)J sequencing, we can now pair TCR/BCR sequences with gene expression data from the same single cell. This tutorial covers two complementary tools for analyzing this data: **scRepertoire** (R/Bioconductor) and **Dandelion** (Python/scverse), each offering distinct strengths for immune repertoire analysis.

---

## Data: 10x Genomics V(D)J Output

Both tools start from the output of the 10x Genomics Cell Ranger `vdj` pipeline. The key input file is `filtered_contig_annotations.csv`, which contains one row per contig (reconstructed receptor chain) with the following critical columns:

| Column             | Description                                |
| ------------------ | ------------------------------------------ |
| `barcode`          | Cell barcode (matches GEX data)            |
| `is_cell`          | Whether the barcode is called as a cell    |
| `contig_id`        | Unique identifier for the assembled contig |
| `chain`            | Chain type (TRA, TRB, IGH, IGK, IGL)       |
| `v_gene`           | Assigned V gene (e.g., TRAV12-1)           |
| `d_gene`           | Assigned D gene (TRB and IGH only)         |
| `j_gene`           | Assigned J gene                            |
| `c_gene`           | Assigned constant gene                     |
| `cdr3`             | CDR3 amino acid sequence                   |
| `cdr3_nt`          | CDR3 nucleotide sequence                   |
| `raw_clonotype_id` | Cell Ranger clonotype assignment           |

For a typical paired GEX + VDJ experiment, you will have:

```
sample_1/
├── filtered_feature_bc_matrix/    # Gene expression
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── vdj_t/                         # TCR V(D)J
    ├── filtered_contig_annotations.csv
    ├── clonotypes.csv
    └── consensus_annotations.csv
```

For this tutorial, we will use a multi-sample tumor immunology dataset. Assume we have matched GEX and VDJ data from tumor, adjacent normal tissue, and peripheral blood of the same patient.

---

## Part 1: scRepertoire (R)

scRepertoire is the most widely adopted R package for immune repertoire analysis from single-cell data. It provides a comprehensive suite of functions for clonotype quantification, diversity analysis, and integration with Seurat or SingleCellExperiment objects.

### Installation and Setup

```r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scRepertoire")

# Load libraries
library(scRepertoire)
library(Seurat)
library(ggplot2)
library(dplyr)
```

### Loading V(D)J Data

The first step is reading the Cell Ranger output and combining contigs into paired-chain clonotypes using `combineTCR()` or `combineBCR()`.

```r
# Read filtered_contig_annotations.csv for each sample
tumor_tcr <- read.csv("data/tumor/vdj_t/filtered_contig_annotations.csv")
normal_tcr <- read.csv("data/normal/vdj_t/filtered_contig_annotations.csv")
blood_tcr <- read.csv("data/blood/vdj_t/filtered_contig_annotations.csv")

# Combine into a list
contig_list <- list(tumor_tcr, normal_tcr, blood_tcr)

# Create combined TCR object with paired alpha-beta chains
combined_tcr <- combineTCR(
  contig_list,
  samples = c("Tumor", "Normal", "Blood"),
  removeNA = TRUE,       # Remove cells without paired chains
  removeMulti = TRUE     # Remove cells with >2 chains (doublets)
)

# Inspect the structure
str(combined_tcr[[1]])
# 'data.frame': ~2500 obs. of 10 variables
#  $ barcode   : chr "Tumor_AAACCTGCAGTATGCT-1" ...
#  $ TCR1      : chr "TRAV12-1.TRAJ33.CVVNMGDSSYKLIF" ...
#  $ TCR2      : chr "TRBV6-1.None.TRBJ2-1.CASSEGQGANEQFF" ...
#  $ CTgene    : chr "TRAV12-1.TRAJ33_TRBV6-1.None.TRBJ2-1" ...
#  $ CTnt      : chr "TGTGTGGTGAATATG..._TGCGCCAGCAGTGA..." ...
#  $ CTaa      : chr "CVVNMGDSSYKLIF_CASSEGQGANEQFF" ...
#  $ CTstrict  : chr "TRAV12-1.TRAJ33.CVVNMGDSSYKLIF_..." ...
#  $ sample    : chr "Tumor" ...
#  $ cloneSize : num ...
```

The `combineTCR()` function pairs alpha (TRA) and beta (TRB) chains for each cell and creates multiple clonotype definitions at different stringency levels:

- **CTgene** --- clonotype defined by V and J gene usage only
- **CTnt** --- clonotype defined by CDR3 nucleotide sequence
- **CTaa** --- clonotype defined by CDR3 amino acid sequence (most commonly used)
- **CTstrict** --- clonotype defined by V gene + J gene + CDR3 nucleotide sequence

For BCR analysis, the equivalent function handles heavy and light chains:

```r
# For B cell receptor data
combined_bcr <- combineBCR(
  bcr_contig_list,
  samples = c("Tumor", "Normal", "Blood"),
  threshold = 0.85  # Hamming distance threshold for clonal grouping
)
```

The `threshold` parameter in `combineBCR()` accounts for somatic hypermutation --- unlike TCRs, BCR sequences undergo affinity maturation, so clonally related B cells may have slightly different CDR3 sequences. A threshold of 0.85 means sequences with >= 85% identity are grouped together.

### Clonotype Quantification

How expanded is the repertoire? `clonalQuant()` visualizes the number of clonotypes per sample, while controlling for differences in cell recovery.

```r
# Quantify clonotypes per sample
clonalQuant(
  combined_tcr,
  cloneCall = "CTaa",       # Use amino acid-level clonotype definition
  chain = "both",           # Require both alpha and beta chains
  scale = TRUE,             # Scale by total number of clonotypes
  exportTable = FALSE
)
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure1-clonal-quantification.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 1. Clonotype quantification across tissue compartments. Tumor-infiltrating T cells show fewer unique clonotypes but greater clonal expansion compared to blood, reflecting antigen-driven selection in the tumor microenvironment.
</div>

To visualize the distribution of clonal expansion sizes:

```r
# Clonal homeostasis --- proportion of cells in different expansion categories
clonalHomeostasis(
  combined_tcr,
  cloneCall = "CTaa",
  cloneSize = c(Rare = 1e-4, Small = 0.001, Medium = 0.01,
                Large = 0.1, Hyperexpanded = 1)
)
```

This categorizes clonotypes by their relative frequency: rare singletons, small clones, medium expansions, and hyperexpanded clones. In tumors with active anti-tumor immunity, you expect a shift toward large and hyperexpanded clones.

```r
# Clonal proportion --- top clonotypes by abundance
clonalProportion(
  combined_tcr,
  cloneCall = "CTaa",
  split = c(10, 50, 100, 500, 1000)
)
```

### Diversity Metrics

Repertoire diversity captures how evenly TCR sequences are distributed. A highly diverse repertoire (many unique clonotypes, each at low frequency) suggests a naive or polyclonal state. A skewed repertoire (dominated by a few expanded clones) indicates antigen-driven selection.

```r
# Calculate multiple diversity indices
clonalDiversity(
  combined_tcr,
  cloneCall = "CTaa",
  metrics = c("shannon", "simpson", "chao1", "ACE"),
  group.by = "sample",
  n.boots = 100          # Bootstrap for confidence intervals
)
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure2-diversity-metrics.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 2. Diversity metrics across compartments. Shannon entropy and inverse Simpson index are lower in tumor compared to blood, indicating oligoclonal expansion of tumor-reactive T cells. Chao1 estimates total species richness including unobserved clonotypes.
</div>

Key diversity indices and their interpretation:

| Metric      | Formula concept                          | Interpretation                                                                 |
| ----------- | ---------------------------------------- | ------------------------------------------------------------------------------ |
| **Shannon** | $$H = -\sum p_i \ln(p_i)$$               | Sensitive to rare clonotypes; higher = more diverse                            |
| **Simpson** | $$D = 1 - \sum p_i^2$$                   | Probability two random cells have different TCRs; dominated by abundant clones |
| **Chao1**   | Estimates total richness from singletons | Predicts unobserved clonotypes; useful for under-sampled repertoires           |
| **ACE**     | Abundance-based coverage estimator       | Similar to Chao1 but uses more rare-species information                        |

Rarefaction analysis is critical when comparing samples with different cell numbers:

```r
# Rarefaction/extrapolation curves
clonalRarefaction(
  combined_tcr,
  cloneCall = "CTaa",
  plot.type = "curve",    # "curve" for interpolation/extrapolation
  n.boots = 50
)
```

### Clonal Overlap Between Samples

Shared clonotypes between tissues reveal immune cell migration and systemic immune responses.

```r
# Pairwise clonal overlap
clonalOverlap(
  combined_tcr,
  cloneCall = "CTaa",
  method = "morisita"     # Options: "overlap", "morisita", "jaccard", "cosine"
)
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure3-clonal-overlap.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 3. Morisita overlap index between tissue compartments. Significant overlap between tumor and blood TCR repertoires suggests active trafficking of tumor-reactive T cells, while lower overlap with adjacent normal tissue indicates tumor-specific recruitment.
</div>

The Morisita index is preferred over simple Jaccard overlap because it accounts for clone abundance, not just presence/absence. Two samples sharing a single hyperexpanded clone are more biologically similar than two samples sharing many rare singletons.

```r
# Visualize specific shared clonotypes
clonalNetwork(
  combined_tcr,
  cloneCall = "CTaa",
  filter.clones = 3,     # Show clonotypes present in 2+ samples
  filter.identity = 0,
  exportClones = "Tumor"  # Focus on tumor-origin clones
)
```

### Integration with Seurat

The real power of paired VDJ + GEX data emerges when you overlay clonotype information onto transcriptomic clusters. `combineExpression()` adds clonotype metadata to a Seurat object.

```r
# Load and process the matched GEX data (standard Seurat workflow)
gex <- Read10X("data/tumor/filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = gex, min.cells = 3, min.features = 200)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Add clonotype information to Seurat object
seurat_obj <- combineExpression(
  combined_tcr,
  seurat_obj,
  cloneCall = "CTaa",
  group.by = "sample",
  proportion = TRUE,       # Use proportional abundance
  cloneSize = c(
    Single = c(0, 1e-4),
    Small = c(1e-4, 0.001),
    Medium = c(0.001, 0.01),
    Large = c(0.01, 0.1),
    Hyperexpanded = c(0.1, 1)
  )
)

# Check added metadata columns
head(seurat_obj@meta.data[, c("CTaa", "cloneSize", "Frequency")])
#                              CTaa       cloneSize  Frequency
# AAACCTGCAGTATGCT-1  CVVN..._CASS...  Large        0.032
# AAACCTGGTCTTGTCC-1  CAVN..._CASS...  Single       0.0004
```

Now visualize clonal expansion on the UMAP embedding:

```r
# Color UMAP by clonal expansion size
DimPlot(seurat_obj, group.by = "cloneSize") +
  scale_color_manual(
    values = c("grey90", "#2166AC", "#67A9CF", "#FDDBC7", "#EF8A62", "#B2182B"),
    na.value = "grey90"
  ) +
  theme_minimal() +
  ggtitle("Clonal Expansion on UMAP")
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure4-umap-clonal-expansion.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 4. UMAP embedding colored by clonal expansion category. Hyperexpanded clones (red) concentrate in cytotoxic CD8+ T cell clusters, while naive and regulatory T cell clusters are predominantly composed of singleton clonotypes (grey), consistent with antigen-experienced effector expansion.
</div>

You can also highlight specific clonotypes of interest:

```r
# Highlight the top 5 most expanded clonotypes
clonalOverlay(
  seurat_obj,
  reduction = "umap",
  freq.cutpoint = 0.01,    # Minimum frequency to highlight
  bins = 25,
  facet.by = "sample"
)

# Alluvial plot: track clonotypes across clusters
alluvialClonotypes(
  seurat_obj,
  cloneCall = "CTaa",
  y.axes = c("sample", "seurat_clusters"),
  color = "CTaa",
  facet.by = NULL
)
```

### Comparing Clonal Expansion Across Clusters

Which cell clusters harbor the most expanded clones? This directly addresses the question of which T cell phenotypes are engaged in the anti-tumor response.

```r
# Clonal diversity per cluster
clonalDiversity(
  seurat_obj,
  cloneCall = "CTaa",
  group.by = "seurat_clusters",
  metrics = c("shannon", "simpson")
)

# Scatter plot: clonal expansion vs. gene expression
clonalScatter(
  seurat_obj,
  cloneCall = "CTaa",
  x.axis = "Tumor",
  y.axis = "Blood",
  graph = "proportion"
)
```

### Gene Usage Analysis

V and J gene usage patterns reveal biases in the TCR repertoire that may reflect antigen-driven selection.

```r
# Visualize V gene usage
vizGenes(
  combined_tcr,
  x = "Vgene",
  y = NULL,
  plot = "barplot",
  chain = "TRB",
  scale = TRUE
)

# Paired V-J gene usage as a chord diagram
vizGenes(
  combined_tcr,
  x = "Vgene",
  y = "Jgene",
  plot = "heatmap",
  chain = "TRB"
)
```

---

## Part 2: Dandelion (Python)

Dandelion is a Python package within the scverse ecosystem designed for single-cell BCR/TCR analysis. It integrates natively with Scanpy and AnnData objects and provides unique capabilities including network-based clonotype analysis and mutation profiling for BCR data.

### Installation and Setup

```python
# Install dandelion and dependencies
# pip install sc-dandelion

import dandelion as ddl
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)
```

### Loading V(D)J Data

Dandelion reads Cell Ranger output into its own `Dandelion` object, which stores per-contig and per-cell information.

```python
# Load filtered_contig_annotations.csv
vdj_tumor = ddl.read_10x_vdj("data/tumor/vdj_t/filtered_contig_annotations.csv")
vdj_normal = ddl.read_10x_vdj("data/normal/vdj_t/filtered_contig_annotations.csv")
vdj_blood = ddl.read_10x_vdj("data/blood/vdj_t/filtered_contig_annotations.csv")

print(vdj_tumor)
# Dandelion object with n_contigs = 5200
#   data: 5200 x 33
#   metadata: 2800 x 15
```

The `Dandelion` object has two core DataFrames:

- **`vdj.data`** --- per-contig information (one row per chain)
- **`vdj.metadata`** --- per-cell information (one row per cell, with paired chain info)

```python
# Examine the contig-level data
print(vdj_tumor.data.columns.tolist())
# ['cell_id', 'contig_id', 'locus', 'v_call', 'd_call', 'j_call',
#  'c_call', 'junction', 'junction_aa', 'productive', ...]

# Examine cell-level metadata
print(vdj_tumor.metadata.head())
#                          clone_id  locus_VDJ  locus_VJ  ...
# cell_id
# AAACCTGCAGTATGCT-1      clone_42   TRB        TRA      ...
# AAACCTGGTCTTGTCC-1      clone_187  TRB        TRA      ...
```

### Preprocessing and Filtering

```python
# Filter to productive contigs only
ddl.tl.filter_contigs(vdj_tumor)

# Quantify clonal groups
# Dandelion uses the AIRR-standard junction (CDR3 + flanking conserved residues)
ddl.tl.find_clones(vdj_tumor, key="junction_aa")

print(f"Unique clonotypes: {vdj_tumor.metadata['clone_id'].nunique()}")
print(f"Total cells with TCR: {len(vdj_tumor.metadata)}")
```

### Gene Usage Visualization

Dandelion provides built-in plotting functions for V/J gene usage analysis.

```python
# V gene usage barplot
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
ddl.pl.barplot(
    vdj_tumor,
    color="v_call_VDJ",    # V gene for TRB (VDJ locus)
    figsize=(12, 5)
)
plt.title("TRB V Gene Usage - Tumor")
plt.tight_layout()
plt.savefig("assets/img/blog/tcr-repertoire/figure5-gene-usage.png", dpi=150)
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure5-gene-usage.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 5. TRB V gene usage in tumor-infiltrating T cells. Skewed V gene usage (e.g., over-representation of TRBV6-1 or TRBV20-1) can indicate antigen-driven selection, where specific V gene-encoded CDR1/CDR2 loops contribute to pMHC recognition.
</div>

### Integration with Scanpy AnnData

Dandelion integrates directly with Scanpy's AnnData object, adding VDJ metadata to `adata.obs`.

```python
# Load matched GEX data
adata = sc.read_10x_mtx("data/tumor/filtered_feature_bc_matrix/")

# Standard scanpy preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs["pct_counts_mt"] < 15, :].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6)

# Transfer VDJ data to AnnData
ddl.tl.transfer(adata, vdj_tumor)

# Check what was added
vdj_columns = [c for c in adata.obs.columns if c.startswith(("clone", "v_call", "j_call",
                                                               "junction", "locus"))]
print(f"VDJ columns added: {vdj_columns}")
# ['clone_id', 'clone_id_size', 'v_call_VDJ', 'v_call_VJ',
#  'j_call_VDJ', 'j_call_VJ', 'junction_aa_VDJ', 'junction_aa_VJ',
#  'locus_VDJ', 'locus_VJ']
```

Now you can use standard Scanpy plotting with VDJ metadata:

```python
# Visualize clone size on UMAP
sc.pl.umap(adata, color="clone_id_size", cmap="Reds",
           title="Clonal Expansion (clone size)")

# Visualize V gene usage on UMAP
sc.pl.umap(adata, color="v_call_VDJ",
           title="TRB V Gene on UMAP")
```

### Network Analysis of Clonotypes

One of Dandelion's unique features is network-based analysis of clonotype relationships. This constructs a graph where cells are nodes and edges connect cells with similar receptor sequences, enabling community detection of clonal families.

```python
# Build a sequence similarity network
ddl.tl.generate_network(
    vdj_tumor,
    key="junction_aa",     # Use CDR3 amino acid sequences
    min_size=2             # Minimum clone size to include
)

# Visualize the clonotype network
ddl.pl.clone_network(
    vdj_tumor,
    color="v_call_VDJ",
    layout="fr",           # Fruchterman-Reingold layout
    figsize=(10, 10)
)
```

{% include figure.liquid loading="eager" path="assets/img/blog/tcr-repertoire/figure6-clonotype-network.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 6. Clonotype similarity network from Dandelion. Each node represents a cell, edges connect cells with similar CDR3 sequences. Dense clusters indicate clonal expansions. Coloring by V gene reveals that expanded clones tend to share V gene usage, consistent with convergent selection toward common antigens.
</div>

For BCR data, Dandelion provides additional mutation analysis capabilities:

```python
# BCR-specific: analyze somatic hypermutation
# Requires IMGT-gapped germline alignment (from IgBLAST or IMGT/HighV-QUEST)
ddl.tl.calculate_threshold(vdj_bcr)        # Determine clonal distance threshold
ddl.tl.find_clones(vdj_bcr, by_alleles=True)

# Mutation analysis
ddl.tl.clone_centrality(vdj_bcr)           # Identify founder sequences
ddl.tl.clone_diversity(vdj_bcr)            # Intra-clonal diversity
```

### Combining Multiple Samples

```python
# Concatenate Dandelion objects from multiple samples
vdj_tumor.metadata["sample"] = "Tumor"
vdj_normal.metadata["sample"] = "Normal"
vdj_blood.metadata["sample"] = "Blood"

# Merge for cross-sample analysis
vdj_combined = ddl.concat([vdj_tumor, vdj_normal, vdj_blood])

# Find shared clonotypes across samples
shared = vdj_combined.metadata.groupby("clone_id")["sample"].nunique()
shared_clones = shared[shared > 1].index.tolist()
print(f"Clonotypes shared across samples: {len(shared_clones)}")

# Annotate shared status
vdj_combined.metadata["is_shared"] = vdj_combined.metadata["clone_id"].isin(shared_clones)
```

---

## Part 3: Immunological Insights from Repertoire Analysis

### Clonal Expansion as a Marker of Immune Response

The degree of clonal expansion directly reflects the magnitude of an adaptive immune response. In the context of tumor immunology, clonal expansion patterns tell a rich biological story:

**Tumor-infiltrating lymphocytes (TILs)** in immunologically "hot" tumors --- those responsive to checkpoint immunotherapy --- typically show oligoclonal expansion. A handful of T cell clones dominate the infiltrate, each presumably recognizing a tumor-associated antigen. In contrast, immunologically "cold" tumors may show either minimal T cell infiltration or a polyclonal repertoire of bystander T cells not engaged in anti-tumor immunity.

```r
# Quantify the fraction of expanded clones per sample
# A higher proportion of expanded clones suggests active immune engagement
expansion_summary <- combined_tcr %>%
  lapply(function(x) {
    total_cells <- nrow(x)
    clone_counts <- table(x$CTaa)
    expanded <- sum(clone_counts > 1)
    singleton <- sum(clone_counts == 1)
    data.frame(
      expanded_clones = expanded,
      singleton_clones = singleton,
      expansion_ratio = expanded / (expanded + singleton),
      top_clone_freq = max(clone_counts) / total_cells
    )
  }) %>%
  bind_rows(.id = "sample")

print(expansion_summary)
#   sample expanded_clones singleton_clones expansion_ratio top_clone_freq
# 1 Tumor            145              320           0.312          0.045
# 2 Normal            42              580           0.068          0.008
# 3 Blood             89             1200           0.069          0.003
```

### Shared Clonotypes Across Tissues

Clonotypes found in multiple tissue compartments provide evidence of active immune surveillance and T cell trafficking. A tumor-reactive T cell clone that is also detected in blood may represent a systemically circulating anti-tumor response --- these are precisely the clones that could be expanded ex vivo for adoptive cell therapy.

```r
# Identify shared clonotypes and their tissue distribution
clonalOverlap(
  combined_tcr,
  cloneCall = "CTaa",
  method = "raw"          # Return raw overlap counts
)

# Extract specific shared clones for downstream analysis
shared_clones <- intersect(
  combined_tcr[["Tumor"]]$CTaa,
  combined_tcr[["Blood"]]$CTaa
)
cat(sprintf("Tumor-Blood shared clonotypes: %d\n", length(shared_clones)))

# These shared clonotypes can be prioritized for:
# 1. Antigen specificity testing (peptide-MHC tetramer staining)
# 2. TCR cloning for adoptive T cell therapy
# 3. Tracking treatment response longitudinally
```

### TCR-Antigen Specificity Prediction

Knowing that a T cell clone is expanded is only half the story --- the critical question is _what antigen does it recognize_? Several computational tools attempt to predict TCR-antigen specificity from sequence alone:

**GLIPH2** (Grouping of Lymphocyte Interactions by Paratope Hotspots) clusters TCR sequences by shared CDR3 motifs that may recognize the same antigen:

```bash
# GLIPH2 command-line usage
# Input: tab-separated file with CDR3b, TRBV, TRBJ, patient, condition
gliph2 \
  --tcr=tcr_input.tsv \
  --refdb=ref_CD4_v2.txt \
  --local_min_pvalue=0.001 \
  --global_min_pvalue=0.001 \
  --output_prefix=gliph2_results
```

**TCRdist3** computes pairwise distances between TCR sequences based on biochemical similarity of CDR loops, enabling hierarchical clustering of TCRs by predicted specificity:

```python
# TCRdist3 example
from tcrdist.repertoire import TCRrep

tr = TCRrep(
    cell_df=tcr_df,           # DataFrame with cdr3_b_aa, v_b_gene, etc.
    organism="human",
    chains=["beta"],
    compute_distances=True
)

# Access the pairwise distance matrix
print(tr.pw_beta.shape)       # (n_clones, n_clones)

# Cluster TCRs by sequence similarity
from scipy.cluster.hierarchy import linkage, fcluster
Z = linkage(tr.pw_beta, method="average")
clusters = fcluster(Z, t=50, criterion="distance")
```

These tools are complementary: GLIPH2 identifies groups of TCRs that likely share antigen specificity based on convergent sequence motifs, while TCRdist3 provides a quantitative distance metric for more nuanced clustering.

### Relevance to CAR-T and Adoptive Cell Therapy

Repertoire analysis from single-cell data has direct translational applications in cellular immunotherapy:

1. **Identifying tumor-reactive clonotypes**: Expanded TIL clones, especially those expressing exhaustion markers (PD-1, TIM-3, LAG-3) alongside cytotoxic effector genes (GZMB, PRF1, IFNG), are strong candidates for tumor reactivity. Their TCR sequences can be cloned into viral vectors for TCR-engineered T cell therapy.

2. **Monitoring CAR-T persistence**: After CAR-T infusion, single-cell VDJ sequencing of the patient's blood can track the endogenous T cell repertoire for signs of epitope spreading --- new clonal expansions against non-CAR-targeted tumor antigens, which indicate a broader anti-tumor immune response.

3. **Predicting treatment response**: Pre-treatment TIL repertoire features (higher clonality, presence of specific V gene biases, greater overlap with blood) have been associated with better responses to anti-PD-1 checkpoint immunotherapy.

```r
# Example: identifying candidate tumor-reactive clones
# Criteria: expanded in tumor, express exhaustion + effector signature
tumor_cells <- subset(seurat_obj, sample == "Tumor")
tumor_cells$is_expanded <- tumor_cells$cloneSize %in% c("Large", "Hyperexpanded")

# Score cells for exhaustion and cytotoxicity signatures
exhaustion_genes <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX")
cytotoxic_genes <- c("GZMB", "GZMA", "PRF1", "IFNG", "NKG7", "GNLY")

tumor_cells <- AddModuleScore(tumor_cells,
                               features = list(exhaustion_genes),
                               name = "exhaustion_score")
tumor_cells <- AddModuleScore(tumor_cells,
                               features = list(cytotoxic_genes),
                               name = "cytotoxic_score")

# Filter for candidate tumor-reactive clones
candidate_cells <- subset(tumor_cells,
                           is_expanded == TRUE &
                           exhaustion_score1 > 0.5 &
                           cytotoxic_score1 > 0.3)

# Extract their TCR sequences for cloning
candidate_tcrs <- candidate_cells@meta.data %>%
  select(CTaa, Frequency, exhaustion_score1, cytotoxic_score1) %>%
  distinct(CTaa, .keep_all = TRUE) %>%
  arrange(desc(Frequency))

print(head(candidate_tcrs, 10))
# Top tumor-reactive TCR candidates ranked by expansion
```

---

## Part 4: Practical Considerations and Best Practices

### Quality Control for V(D)J Data

Before diving into analysis, apply these QC filters:

```r
# Check for doublets: cells with >2 productive chains of the same type
# combineTCR handles this with removeMulti = TRUE, but verify:
qc_summary <- lapply(contig_list, function(x) {
  productive <- x %>% filter(productive == "True", is_cell == "True")
  cells_per_chain <- productive %>%
    group_by(barcode, chain) %>%
    summarise(n_chains = n(), .groups = "drop")
  multi_chain <- cells_per_chain %>% filter(n_chains > 2)
  cat(sprintf("Multi-chain cells (potential doublets): %d / %d\n",
              n_distinct(multi_chain$barcode),
              n_distinct(productive$barcode)))
})

# Check chain pairing rates
pairing_rate <- sapply(combined_tcr, function(x) {
  has_alpha <- !is.na(x$TCR1) & x$TCR1 != "NA"
  has_beta <- !is.na(x$TCR2) & x$TCR2 != "NA"
  mean(has_alpha & has_beta)
})
cat(sprintf("Alpha-beta pairing rates: %s\n",
            paste(sprintf("%.1f%%", pairing_rate * 100), collapse = ", ")))
# Typical pairing rates: 50-70% for 10x 5' chemistry
```

### Choosing Clonotype Definitions

The choice of clonotype definition affects your results significantly:

| Definition | Stringency | Best for                                              |
| ---------- | ---------- | ----------------------------------------------------- |
| CTgene     | Lowest     | V/J gene usage analysis, broad patterns               |
| CTaa       | Medium     | Most analyses; groups convergent TCRs                 |
| CTnt       | High       | Distinguishing truly independent recombination events |
| CTstrict   | Highest    | When V gene identity matters for specificity          |

For most immunological questions, **CTaa** (CDR3 amino acid) is the recommended default. It captures convergent recombination --- different nucleotide sequences encoding the same amino acid CDR3 --- which is a hallmark of antigen-driven selection.

### Downsampling and Statistical Considerations

Comparing repertoire metrics between samples with vastly different cell numbers requires careful normalization:

```python
# In Python/Dandelion: bootstrap diversity estimation
import numpy as np

def bootstrap_diversity(clone_counts, n_iter=1000, subsample_size=None):
    """Bootstrap Shannon diversity with subsampling."""
    if subsample_size is None:
        subsample_size = min(len(clone_counts), 500)

    diversities = []
    cells = np.repeat(np.arange(len(clone_counts)), clone_counts)

    for _ in range(n_iter):
        sample = np.random.choice(cells, size=subsample_size, replace=False)
        _, counts = np.unique(sample, return_counts=True)
        freqs = counts / counts.sum()
        shannon = -np.sum(freqs * np.log(freqs))
        diversities.append(shannon)

    return np.mean(diversities), np.std(diversities)

# Apply to each sample
for sample_name, vdj in [("Tumor", vdj_tumor), ("Normal", vdj_normal), ("Blood", vdj_blood)]:
    clone_sizes = vdj.metadata["clone_id"].value_counts().values
    mean_div, std_div = bootstrap_diversity(clone_sizes)
    print(f"{sample_name}: Shannon = {mean_div:.3f} +/- {std_div:.3f}")
```

---

## Key Takeaways

1. **V(D)J recombination creates the foundation of adaptive immunity.** Each unique TCR/BCR sequence is a molecular barcode for a clonal lineage, and sequencing these receptors at single-cell resolution connects immune recognition to cell state.

2. **scRepertoire excels at clonotype quantification and Seurat integration.** Its comprehensive set of diversity metrics (Shannon, Simpson, Chao1), overlap measures (Morisita, Jaccard), and visualization functions make it the go-to R package for repertoire analysis. The `combineExpression()` function seamlessly bridges VDJ and GEX data.

3. **Dandelion brings Python/scverse integration and network analysis.** Native AnnData compatibility means VDJ metadata flows directly into Scanpy workflows. The network-based clonotype analysis is particularly powerful for BCR data, where somatic hypermutation creates clonal families with graded sequence similarity.

4. **Clonal expansion patterns are biologically informative.** Oligoclonal expansion in tumors marks antigen-driven immune responses. Shared clonotypes across tissues reveal immune cell trafficking. The degree of expansion correlates with effector differentiation.

5. **Repertoire analysis has direct translational applications.** Identifying expanded, exhausted TIL clonotypes provides TCR sequences for engineered T cell therapy. Repertoire features predict checkpoint immunotherapy response. Monitoring clonal dynamics tracks treatment efficacy.

6. **Choose your clonotype definition thoughtfully.** CDR3 amino acid sequence (CTaa) captures convergent selection and is appropriate for most analyses. Use stricter definitions (CTstrict) when V gene identity matters for specificity prediction. For BCR data, account for somatic hypermutation with distance-based clustering.

The combination of these tools with the growing databases of known TCR-antigen pairs (VDJdb, McPAS-TCR, IEDB) is steadily closing the gap between observing an immune response and understanding what drives it --- bringing us closer to truly personalized immunotherapy.

---

## References and Resources

- **scRepertoire**: [Borcherding et al., F1000Research 2020](https://doi.org/10.12688/f1000research.22139.2) --- [Documentation](https://www.borch.dev/uploads/screpertoire/)
- **Dandelion**: [Suo et al., Nature Biotechnology 2024](https://doi.org/10.1038/s41587-023-01734-7) --- [Documentation](https://sc-dandelion.readthedocs.io/)
- **GLIPH2**: [Huang et al., Nature Biotechnology 2020](https://doi.org/10.1038/s41587-020-0505-4)
- **TCRdist3**: [Mayer-Blackwell et al., eLife 2021](https://doi.org/10.7554/eLife.68605)
- **10x Genomics V(D)J**: [Support documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest)
- **AIRR Community Standards**: [Vander Heiden et al., Frontiers in Immunology 2018](https://doi.org/10.3389/fimmu.2018.02206)
