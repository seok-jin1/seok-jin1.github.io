---
layout: post
title: "CellChat and CellPhoneDB: Decoding Cell-Cell Communication in the Immune Microenvironment"
date: 2025-12-28
permalink: /blog/cellchat-cell-communication/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - single-cell
  - immunology
  - python
  - tutorial
description: "A practical guide to cell-cell communication inference with CellChat and CellPhoneDB — building ligand-receptor networks, visualising pathway activity, and comparing signalling between tumour and normal tissue in immune microenvironments."
---

## Introduction

Single-cell RNA sequencing tells us what each cell is expressing, but cells do not exist in isolation. In the tumor microenvironment (TME), immune cells, stromal cells, and cancer cells constantly exchange signals through ligand-receptor interactions --- checkpoint molecules like PD-1/PD-L1 that suppress T cell killing, chemokines like CXCL9/CXCR3 that recruit effector cells into the tumor bed, and costimulatory signals like CD28/B7 that determine whether a T cell mounts a productive response or becomes anergic.

Computational tools that infer cell-cell communication from single-cell transcriptomic data allow us to reconstruct these signaling networks and identify the molecular axes that drive immune evasion. In this tutorial, we walk through **CellChat** (R) and **CellPhoneDB** (Python), covering the full workflow from input preparation to publication-quality visualizations.

---

## Environment Setup

### R Environment (CellChat)

```r
# Install CellChat and dependencies
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("jinworks/CellChat")

# Additional packages
install.packages(c("NMF", "circlize", "ComplexHeatmap"))
BiocManager::install("BiocNeighbors")

# Load libraries
library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
library(circlize)
```

### Python Environment (CellPhoneDB)

```python
# Install CellPhoneDB
# pip install cellphonedb scanpy pandas matplotlib

import cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
```

---

## Part 1: CellChat (R)

CellChat uses a curated database of ligand-receptor interactions (CellChatDB) that includes secreted signaling, ECM-receptor interactions, and cell-cell contact. It models communication probability using the law of mass action, accounting for the expression of ligands, receptors, and their cofactors.

### 1.1 Creating a CellChat Object from Seurat

Assume you have a Seurat object from a tumor scRNA-seq experiment with cell type annotations stored in `seurat_obj$cell_type`.

```r
# Load your Seurat object (example: tumor microenvironment dataset)
seurat_obj <- readRDS("tme_seurat.rds")

# Verify cell type annotations
table(seurat_obj$cell_type)
# CD8_T   CD4_T   Treg   Macrophage   DC   NK   B_cell   Tumor
# 1823    1456    387    2104          512  634  891      3421

# Extract normalized expression matrix and metadata
data.input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
meta <- data.frame(labels = seurat_obj$cell_type, row.names = colnames(seurat_obj))

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human  # Use CellChatDB.mouse for mouse data
showDatabaseCategory(CellChatDB)
# [1] "Secreted Signaling"   "ECM-Receptor"   "Cell-Cell Contact"

# Use all interaction categories (or subset to specific ones)
cellchat@DB <- CellChatDB
```

### 1.2 Identifying Overexpressed Ligands and Receptors

CellChat identifies genes overexpressed in each cell group, focusing the analysis on biologically meaningful interactions.

```r
# Identify overexpressed signaling genes
cellchat <- subsetData(cellchat)

# Identify overexpressed ligand-receptor interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
```

### 1.3 Inference of Communication Probability

This is the core computational step. CellChat computes the communication probability for each ligand-receptor pair between each pair of cell groups using a triMean function (a robust measure that reduces the effect of outlier expression values).

```r
# Compute communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Filter communications with fewer than 10 cells in a group
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred communication network as a data frame
df.net <- subsetCommunication(cellchat)
head(df.net[, c("source", "target", "ligand", "receptor", "prob", "pval")])
#       source     target  ligand   receptor       prob     pval
# 1     Tumor  Macrophage  CSF1     CSF1R        0.0234   0.001
# 2     Tumor      CD8_T   PVR     TIGIT        0.0189   0.003
# 3  Macrophage  CD8_T     CXCL9    CXCR3        0.0156   0.002
# 4     Tumor      CD8_T   PDL1     PDCD1        0.0312   0.001
# 5        DC      CD4_T   CD80     CD28         0.0278   0.001
```

### 1.4 Signaling Pathway Analysis

CellChat aggregates individual ligand-receptor pairs into biologically meaningful signaling pathways, making it easier to interpret the results.

```r
# Compute communication probability at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the cell-cell communication network
cellchat <- aggregateNet(cellchat)

# View the number of interactions and interaction strength
groupSize <- as.numeric(table(cellchat@idents))
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure1-aggregated-network.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 1. Aggregated cell-cell communication network showing the number of interactions (left) and total interaction strength (right) between cell types in the tumor microenvironment.
</div>

### 1.5 Visualization

CellChat provides a rich set of visualization functions. Here we demonstrate the most useful ones for immunology studies.

#### Hierarchy Plot

The hierarchy plot shows the inferred communication network for a specific signaling pathway, with the cell types arranged hierarchically.

```r
# Hierarchy plot for PD-L1 signaling
pathways.show <- "PD-L1"
vertex.receiver <- seq(1, 4)  # Index of target cell groups
netVisual_aggregate(cellchat,
                    signaling = pathways.show,
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy")
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure2-hierarchy-plot.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 2. Hierarchy plot of PD-L1 signaling. Tumor cells and macrophages express PD-L1 that engages PD-1 on CD8+ T cells, representing a key immune checkpoint axis in the TME.
</div>

#### Circle Plot

```r
# Circle plot for CXCL signaling
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "circle")

# Chord diagram for all significant interactions
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "chord")
```

#### Bubble Plot

The bubble plot is one of the most informative visualizations, showing communication probability (color) and statistical significance (size) for each ligand-receptor pair.

```r
# Bubble plot: interactions from Macrophage and DC to T cells
netVisual_bubble(cellchat,
                 sources.use = c("Macrophage", "DC"),
                 targets.use = c("CD8_T", "CD4_T", "Treg"),
                 remove.isolate = FALSE)
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure3-bubble-plot.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 3. Bubble plot showing ligand-receptor interactions from myeloid cells (macrophages and DCs) to T cell subsets. Color encodes the communication probability (redder = stronger predicted signalling); dot size encodes the statistical significance (p-value) of the interaction.
</div>

### 1.6 Comparing Conditions: Tumor vs Normal

One of CellChat's most powerful features is the ability to compare communication networks across conditions (e.g., tumor vs adjacent normal tissue).

```r
# Create CellChat objects for each condition
cellchat_tumor <- createCellChat(
  object = GetAssayData(seurat_tumor, assay = "RNA", layer = "data"),
  meta = data.frame(labels = seurat_tumor$cell_type,
                    row.names = colnames(seurat_tumor)),
  group.by = "labels"
)

cellchat_normal <- createCellChat(
  object = GetAssayData(seurat_normal, assay = "RNA", layer = "data"),
  meta = data.frame(labels = seurat_normal$cell_type,
                    row.names = colnames(seurat_normal)),
  group.by = "labels"
)

# Process each object through the full pipeline
process_cellchat <- function(cc) {
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  return(cc)
}

cellchat_tumor <- process_cellchat(cellchat_tumor)
cellchat_normal <- process_cellchat(cellchat_normal)

# Merge the two CellChat objects
object.list <- list(Normal = cellchat_normal, Tumor = cellchat_tumor)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
```

```r
# Compare interaction counts, differential network, and pathway ranking
compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure4-condition-comparison.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 4. Comparison of cell-cell communication between normal and tumor tissue. Red edges indicate increased interactions in tumor; blue edges indicate decreased interactions. Note the strong upregulation of immune checkpoint signaling in the tumor condition.
</div>

---

## Part 2: CellPhoneDB (Python)

CellPhoneDB uses a curated database of heteromeric receptor complexes and employs permutation testing to identify significant ligand-receptor interactions.

### 2.1 Input Preparation from AnnData

```python
import scanpy as sc
import pandas as pd

# Load your AnnData object
adata = sc.read_h5ad("tme_annotated.h5ad")

# Verify cell type annotations
print(adata.obs["cell_type"].value_counts())
# Tumor         3421
# Macrophage    2104
# CD8_T         1823
# CD4_T         1456
# B_cell         891
# NK             634
# DC             512
# Treg           387

# Ensure the data is normalized (log1p-transformed counts)
# If raw counts, normalize first:
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# 1. Count matrix (genes x cells, normalized)
count_df = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
).T  # Transpose: genes as rows, cells as columns

count_df.to_csv("cellphonedb_counts.txt", sep="\t")

# 2. Metadata file (cell barcode -> cell type)
meta_df = pd.DataFrame({
    "Cell": adata.obs_names,
    "cell_type": adata.obs["cell_type"].values
})
meta_df.to_csv("cellphonedb_meta.txt", sep="\t", index=False)
```

### 2.2 Running the Statistical Analysis

```python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Run CellPhoneDB statistical analysis
deconvoluted, means, pvalues, significant_means = \
    cpdb_statistical_analysis_method.call(
        cpdb_file_path="cellphonedb.zip",       # CellPhoneDB database
        meta_file_path="cellphonedb_meta.txt",
        counts_file_path="cellphonedb_counts.txt",
        counts_data="hgnc_symbol",               # Gene identifier type
        threshold=0.1,                            # Min % of cells expressing gene
        iterations=1000,                          # Number of permutations
        threads=4,                                # Parallel threads
        output_path="cpdb_output/"
    )
```

### 2.3 Interpreting Results

CellPhoneDB produces four key output files:

```python
import pandas as pd

# Load results
means = pd.read_csv("cpdb_output/means.txt", sep="\t")
pvalues = pd.read_csv("cpdb_output/pvalues.txt", sep="\t")
significant_means = pd.read_csv("cpdb_output/significant_means.txt", sep="\t")
deconvoluted = pd.read_csv("cpdb_output/deconvoluted.txt", sep="\t")

# The columns represent cell type pairs (e.g., "CD8_T|Macrophage")
# Each row is a ligand-receptor interaction

# Filter for significant interactions (p < 0.05) involving T cells
interaction_cols = [c for c in pvalues.columns if "CD8_T" in c]
sig_interactions = pvalues[
    pvalues[interaction_cols].min(axis=1) < 0.05
][["interacting_pair"] + interaction_cols]

print(sig_interactions.head(10))
# interacting_pair           CD8_T|Macrophage  Macrophage|CD8_T  CD8_T|Tumor
# PDCD1_CD274                0.001             NaN               0.002
# TIGIT_PVR                  NaN               NaN               0.003
# CXCR3_CXCL9               NaN               0.001             NaN
# CD28_CD86                  NaN               0.012             NaN
# CTLA4_CD80                 NaN               0.008             NaN
# IFNG_IFNGR1_IFNGR2         0.004             NaN               0.011
# TNF_TNFRSF1A               0.002             NaN               NaN

# Get the mean expression values for significant interactions
sig_means_filtered = significant_means[
    significant_means["interacting_pair"].isin(sig_interactions["interacting_pair"])
]
```

The four output files contain: **means.txt** (average ligand-receptor expression per cell pair), **pvalues.txt** (permutation test p-values), **significant_means.txt** (means only for significant interactions), and **deconvoluted.txt** (per-gene expression in multi-subunit complexes).

### 2.4 Dot Plot Visualization

```python
import matplotlib.pyplot as plt
import numpy as np

# Select interactions and cell pairs of interest
interactions_of_interest = [
    "PDCD1_CD274", "CTLA4_CD80", "CTLA4_CD86", "CD28_CD86",
    "TIGIT_PVR", "CXCR3_CXCL9", "CXCR3_CXCL10", "CCR5_CCL5",
]
cell_pairs_of_interest = [
    "CD8_T|Tumor", "CD8_T|Macrophage", "Macrophage|CD8_T",
    "DC|CD4_T", "DC|CD8_T", "Treg|CD8_T", "Tumor|Macrophage"
]

# Build dot plot data
mask = means["interacting_pair"].isin(interactions_of_interest)
plot_means = means.loc[mask, ["interacting_pair"] + cell_pairs_of_interest].set_index("interacting_pair")
plot_pvals = pvalues.loc[mask, ["interacting_pair"] + cell_pairs_of_interest].set_index("interacting_pair")

# Create dot plot
fig, ax = plt.subplots(figsize=(12, 8))
for i, interaction in enumerate(plot_means.index):
    for j, pair in enumerate(plot_means.columns):
        mean_val = plot_means.loc[interaction, pair]
        p_val = plot_pvals.loc[interaction, pair]
        if pd.notna(mean_val) and pd.notna(p_val):
            size = max(20, -np.log10(p_val + 1e-4) * 30)
            color = plt.cm.YlOrRd(min(mean_val / plot_means.max().max(), 1.0))
            ax.scatter(j, i, s=size, c=[color], edgecolors="black", linewidth=0.5)

ax.set_xticks(range(len(cell_pairs_of_interest)))
ax.set_xticklabels(cell_pairs_of_interest, rotation=45, ha="right")
ax.set_yticks(range(len(plot_means.index)))
ax.set_yticklabels(plot_means.index)
ax.set_title("CellPhoneDB: Significant Interactions in TME")
plt.tight_layout()
plt.savefig("cpdb_dotplot.png", dpi=150, bbox_inches="tight")
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure5-cellphonedb-dotplot.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 5. CellPhoneDB dot plot showing significant ligand-receptor interactions between immune cell types in the TME. Dot size represents statistical significance (-log10 p-value); color intensity indicates mean expression level.
</div>

---

## Part 3: CellChat vs CellPhoneDB vs LIANA

Choosing between tools depends on your analysis goals. Here is a practical comparison:

| Feature                     | CellChat                                              | CellPhoneDB                              | LIANA                                |
| --------------------------- | ----------------------------------------------------- | ---------------------------------------- | ------------------------------------ |
| **Language**                | R                                                     | Python                                   | R / Python                           |
| **Database**                | CellChatDB (secreted, ECM, cell-cell contact)         | CellPhoneDB DB (multi-subunit complexes) | Consensus of multiple DBs            |
| **Statistical method**      | Permutation test + triMean                            | Permutation test on means                | Multiple methods (aggregate ranking) |
| **Multi-subunit receptors** | Cofactor modeling                                     | Explicit complex definition              | Depends on method                    |
| **Condition comparison**    | Built-in (mergeCellChat)                              | Manual                                   | Built-in                             |
| **Visualization**           | Extensive (hierarchy, circle, chord, bubble, heatmap) | Basic (dot plot, heatmap)                | Standardized output                  |
| **Spatial support**         | CellChat v2 (spatial transcriptomics)                 | CellPhoneDB v5 (spatial)                 | Via spatialDM                        |
| **Speed**                   | Moderate                                              | Fast                                     | Varies by method                     |

**CellChat** is best for comprehensive pathway-level analysis with rich visualizations and built-in condition comparison. **CellPhoneDB** excels at fast permutation-based statistics with careful handling of multi-subunit receptors, integrating well with Scanpy workflows. **LIANA** provides consensus rankings across multiple methods and databases, reducing method-specific biases.

---

## Part 4: Immunology Case Study --- T Cell, Macrophage, and DC Interactions in the TME

Now let us apply these tools to dissect three critical immunological axes in the tumor microenvironment.

### 4.1 The PD-1/PD-L1 Checkpoint Axis

The PD-1 (PDCD1) / PD-L1 (CD274) axis is the most therapeutically important immune checkpoint in oncology. Tumor cells and tumor-associated macrophages (TAMs) express PD-L1 to engage PD-1 on exhausted CD8+ T cells, delivering an inhibitory signal that suppresses cytotoxic activity.

```r
# Visualize PD-L1 signaling network
netVisual_aggregate(cellchat, signaling = "PD-L1", layout = "circle",
                    edge.width.max = 10)

# Which cell pairs contribute most to PD-L1 signaling?
netAnalysis_contribution(cellchat, signaling = "PD-L1")

# Gene expression of PD-L1 axis components
plotGeneExpression(cellchat, signaling = "PD-L1",
                  enriched.only = FALSE, type = "violin")

# In the tumor condition, we expect:
# - High CD274 (PD-L1) expression in Tumor and Macrophage populations
# - High PDCD1 (PD-1) expression in CD8_T (especially exhausted subset)
# - Significant Tumor -> CD8_T and Macrophage -> CD8_T interactions
```

In many solid tumors, TAMs are the dominant source of PD-L1, not tumor cells. CellChat reveals this by showing relative sender contributions --- a finding with direct implications for anti-PD-1 therapy response prediction.

### 4.2 CXCL/CXCR Chemokine Signaling

The CXCL9/CXCL10/CXCL11-CXCR3 axis is critical for recruiting effector T cells and NK cells into the tumor. Macrophages and DCs in the TME produce these chemokines in response to IFN-gamma signaling, creating a positive feedback loop with incoming T cells.

```r
# Visualize CXCL chemokine signaling
netVisual_aggregate(cellchat, signaling = "CXCL", layout = "chord")

# Bubble plot focusing on chemokine interactions
netVisual_bubble(cellchat,
                 sources.use = c("Macrophage", "DC"),
                 targets.use = c("CD8_T", "NK", "CD4_T"),
                 signaling = "CXCL",
                 remove.isolate = FALSE)

# Compare CXCL signaling between tumor and normal
netVisual_bubble(cellchat_merged,
                 sources.use = c("Macrophage", "DC"),
                 targets.use = c("CD8_T", "NK"),
                 signaling = "CXCL",
                 comparison = c(1, 2),
                 angle.x = 45)
```

In "hot" (T cell-inflamed) tumors, you expect strong CXCL9/10 expression from myeloid cells and corresponding CXCR3 on CD8+ T cells. In "cold" tumors, this axis is often silenced by immunosuppressive cytokines like TGF-beta.

### 4.3 CD28/B7 Costimulatory Signaling

T cell activation requires two signals: TCR engagement (signal 1) and costimulation through CD28 binding to CD80/CD86 on antigen-presenting cells (signal 2). In the TME, this costimulatory axis competes with the inhibitory CTLA-4 pathway, as CTLA-4 binds the same B7 ligands with higher affinity.

```r
# Examine costimulatory vs inhibitory balance
# CD28 (costimulatory) and CTLA-4 (inhibitory) compete for CD80/CD86
netVisual_bubble(cellchat,
                 sources.use = c("DC", "Macrophage"),
                 targets.use = c("CD8_T", "CD4_T", "Treg"),
                 pairLR.use = data.frame(
                   ligand = c("CD80", "CD86", "CD80", "CD86"),
                   receptor = c("CD28", "CD28", "CTLA4", "CTLA4")
                 ),
                 remove.isolate = FALSE)
```

{% include figure.liquid loading="eager" path="assets/img/blog/cellchat/figure6-costimulation-balance.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 6. Costimulatory (CD28) versus inhibitory (CTLA-4) signaling balance. DCs provide stronger CD28 costimulation to conventional T cells, while Tregs preferentially engage the CTLA-4 inhibitory pathway, reflecting their role in maintaining immune tolerance in the TME.
</div>

### 4.4 Integrating CellPhoneDB Results

We can validate key findings using CellPhoneDB for orthogonal confirmation:

```python
# Filter CellPhoneDB results for the three axes
checkpoint_interactions = [
    "PDCD1_CD274", "CTLA4_CD80", "CTLA4_CD86",
    "TIGIT_PVR", "LAG3_LGALS3"
]

chemokine_interactions = [
    "CXCR3_CXCL9", "CXCR3_CXCL10", "CXCR3_CXCL11",
    "CCR5_CCL5", "CCR5_CCL4"
]

costim_interactions = [
    "CD28_CD86", "CD28_CD80",
    "ICOS_ICOSLG"
]

all_axes = checkpoint_interactions + chemokine_interactions + costim_interactions

# Extract significant interactions
sig_mask = significant_means["interacting_pair"].isin(all_axes)
tme_interactions = significant_means[sig_mask].copy()

# Summarize by axis
for axis_name, interactions in [
    ("Checkpoint", checkpoint_interactions),
    ("Chemokine", chemokine_interactions),
    ("Costimulatory", costim_interactions)
]:
    subset = tme_interactions[tme_interactions["interacting_pair"].isin(interactions)]
    n_sig = subset.notna().sum(axis=1).sum()
    print(f"{axis_name}: {len(subset)} interactions, {n_sig} significant cell-pair entries")

# Checkpoint: 4 interactions, 12 significant cell-pair entries
# Chemokine: 4 interactions, 8 significant cell-pair entries
# Costimulatory: 3 interactions, 6 significant cell-pair entries
```

---

## Key Takeaways

**1. CellChat and CellPhoneDB answer different but complementary questions.** CellChat excels at pathway-level interpretation and condition comparison with rich visualizations. CellPhoneDB provides rigorous permutation-based statistics with careful handling of multi-subunit receptor complexes. Using both gives you higher confidence in your findings.

**2. The database matters as much as the method.** Each tool uses a different curated database of ligand-receptor interactions. Interactions missing from the database will never be detected. Consider supplementing with LIANA, which aggregates results across multiple databases and methods.

**3. Communication analysis is hypothesis-generating, not definitive.** These tools infer _potential_ communication based on co-expression of ligands and receptors. They cannot confirm that signaling actually occurs. Validation through protein-level assays (flow cytometry, imaging), functional experiments (co-culture, blocking antibodies), or spatial transcriptomics is essential.

**4. For immunology, focus on biologically interpretable axes.** Rather than reporting hundreds of significant interactions, organize your results around known immunological circuits: checkpoint inhibition, costimulation, chemokine recruitment, and cytokine polarization. This makes the analysis actionable for downstream experiments.

**5. Condition comparison reveals therapeutic opportunities.** Comparing tumor vs normal (or responder vs non-responder) communication networks can identify the specific signaling axes that are rewired in disease. These represent candidate targets for immunotherapy --- whether through checkpoint blockade, chemokine modulation, or costimulatory agonists.

**6. Consider spatial context.** Both CellChat v2 and CellPhoneDB v5 now support spatial transcriptomics data, allowing you to infer communication only between physically proximal cells --- especially important in the TME where immune cell localization determines functional outcomes.
