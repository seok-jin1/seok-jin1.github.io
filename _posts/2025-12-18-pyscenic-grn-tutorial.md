---
layout: post
title: "pySCENIC: Inferring Gene Regulatory Networks from Single-Cell Data"
date: 2025-12-18
permalink: /blog/pyscenic-grn-tutorial/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - single-cell
  - python
  - immunology
  - tutorial
---

Understanding how transcription factors (TFs) orchestrate gene expression programs in individual cells is one of the central questions in immunology. Which TFs drive a naive CD4+ T cell to become a Th1 effector versus a regulatory T cell? How do monocytes commit to becoming dendritic cells versus macrophages? **pySCENIC** (Single-Cell rEgulatory Network Inference and Clustering) provides a computational framework to answer these questions by reconstructing **gene regulatory networks (GRNs)** directly from single-cell RNA-seq data.

This tutorial walks through the complete pySCENIC pipeline, from raw expression matrices to regulon activity maps, with a focus on interpreting results in an immunological context.

---

## 1. What Are Regulons and Why Do They Matter?

A **regulon** is defined as a transcription factor together with its set of direct target genes. Unlike simple co-expression modules, regulons in SCENIC are validated by the presence of the TF's binding motif in the cis-regulatory regions of its predicted targets. This makes them biologically grounded units of gene regulation.

In immunology, regulons map directly to well-known biology:

| TF               | Regulon context                               | Immune cell type               |
| ---------------- | --------------------------------------------- | ------------------------------ |
| TBX21 (T-bet)    | Th1 differentiation, IFN-gamma production     | Th1 CD4+ T cells, CD8+ T cells |
| GATA3            | Th2 differentiation, IL-4/5/13 regulation     | Th2 CD4+ T cells               |
| RORC (RORgammat) | Th17 differentiation, IL-17 production        | Th17 CD4+ T cells              |
| FOXP3            | Regulatory T cell identity, immunosuppression | Tregs                          |
| SPI1 (PU.1)      | Myeloid lineage commitment                    | Monocytes, macrophages, DCs    |
| MAFB             | Monocyte/macrophage differentiation           | Macrophages                    |
| IRF8             | DC specification, type I IFN signaling        | cDC1, pDCs                     |

By recovering these regulons from scRNA-seq data in an unbiased manner, SCENIC can reveal both known and novel regulatory programs operating in your dataset.

---

## 2. The SCENIC Pipeline: Three Steps

The SCENIC workflow consists of three sequential steps:

**Step 1: GRNBoost2** -- Infer co-expression modules between TFs and candidate target genes using gradient boosting regression. Each TF-target pair receives an importance score reflecting the strength of the regulatory relationship.

**Step 2: cisTarget (RcisTarget)** -- Prune co-expression modules by requiring that the TF's binding motif is enriched in the cis-regulatory regions of predicted targets. This converts raw co-expression modules into bona fide regulons.

**Step 3: AUCell** -- Score each cell for the activity of each regulon using the Area Under the recovery Curve (AUC). This produces a cells-by-regulons activity matrix that can be used for clustering, visualization, and differential analysis.

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure1-pipeline-overview.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 1. Overview of the three-step SCENIC pipeline. Gene expression data is first used to build co-expression modules (GRNBoost2), which are then pruned using motif enrichment (cisTarget), and finally scored per cell (AUCell).
</div>

---

## 3. Setup and Installation

### 3.1 Install pySCENIC

```bash
# Create a dedicated conda environment (recommended)
conda create -n scenic python=3.10 -y
conda activate scenic

# Install pySCENIC and dependencies
pip install pyscenic

# Or via conda
conda install -c bioconda pyscenic
```

### 3.2 Download Reference Databases

pySCENIC requires two sets of reference files: (1) ranking databases that score each gene's cis-regulatory regions for motif occurrences, and (2) motif-to-TF annotation files that map motifs to their cognate transcription factors.

```bash
# Create a directory for databases
mkdir -p ~/scenic_databases && cd ~/scenic_databases

# --- Ranking databases (for human, hg38) ---
# These are large files (~1.5 GB each); download the ones matching your genome build

# 10kb around TSS
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# 500bp upstream
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# --- Motif-to-TF annotations ---
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# --- TF list (human) ---
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
```

> **Tip:** Database versions matter. Always match the ranking database version (e.g., `v10_clust`) with the corresponding motif annotation file (`v10nr_clust`). Mismatches will produce zero or spurious regulons.

### 3.3 Verify the setup

```python
import pyscenic
print(f"pySCENIC version: {pyscenic.__version__}")

# Verify ranking database loads correctly
from ctxcore.rnkdb import FeatherRankingDatabase
db_path = "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
db = FeatherRankingDatabase(db_path, name="hg38_10kb")
print(f"Database loaded: {db.name}, genes: {db.total_genes}")
```

---

## 4. Preparing the Expression Data

We start with a preprocessed scRNA-seq dataset of PBMCs. For this tutorial, we assume you have an AnnData object with raw counts, cell type annotations, and basic QC already performed.

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load PBMC dataset (e.g., 10x Genomics PBMC 10k)
adata = sc.read_h5ad("pbmc_10k_filtered.h5ad")

print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Cell types: {adata.obs['cell_type'].unique().tolist()}")
```

```
Cells: 9876, Genes: 33538
Cell types: ['CD4_Naive', 'CD4_Th1', 'CD4_Th2', 'CD4_Th17', 'Treg',
             'CD8_Naive', 'CD8_Effector', 'NK', 'B_cell',
             'CD14_Mono', 'CD16_Mono', 'cDC1', 'cDC2', 'pDC']
```

### 4.1 Gene Filtering

SCENIC works best with genes that have sufficient expression variation. Filter out genes expressed in too few cells and those with very low counts.

```python
# Use raw counts for SCENIC
# If adata.X contains normalized data, use adata.raw or a separate raw layer
if adata.raw is not None:
    adata_raw = adata.raw.to_adata()
else:
    adata_raw = adata.copy()

# Filter genes: keep genes expressed in at least 3 cells
sc.pp.filter_genes(adata_raw, min_cells=3)

# Further filter to keep only genes with reasonable expression
# This dramatically reduces runtime
sc.pp.filter_genes(adata_raw, min_counts=3)

print(f"Genes after filtering: {adata_raw.n_vars}")
```

### 4.2 Create the Expression Matrix

```python
# SCENIC expects a pandas DataFrame: cells x genes, with gene names as columns
# Use raw counts (not normalized/log-transformed)
ex_matrix = pd.DataFrame(
    data=adata_raw.X.toarray() if hasattr(adata_raw.X, 'toarray') else adata_raw.X,
    index=adata_raw.obs_names,
    columns=adata_raw.var_names
)

print(f"Expression matrix shape: {ex_matrix.shape}")
print(f"Non-zero fraction: {(ex_matrix > 0).sum().sum() / ex_matrix.size:.4f}")
```

### 4.3 Load the TF List

```python
# Load human TF list
tf_file = "allTFs_hg38.txt"
tf_names = [line.strip() for line in open(tf_file, "r").readlines()]

# Keep only TFs present in our expression matrix
tf_names = [tf for tf in tf_names if tf in ex_matrix.columns]
print(f"TFs in dataset: {len(tf_names)}")

# Check for our immunology TFs of interest
immune_tfs = ['TBX21', 'GATA3', 'RORC', 'FOXP3', 'SPI1', 'MAFB', 'IRF8',
              'BATF', 'BCL6', 'STAT4', 'STAT6', 'RUNX3', 'EOMES', 'IRF4']
present = [tf for tf in immune_tfs if tf in tf_names]
missing = [tf for tf in immune_tfs if tf not in tf_names]
print(f"Immune TFs present: {present}")
if missing:
    print(f"Immune TFs missing (may be filtered): {missing}")
```

---

## 5. Step 1: Gene Co-expression with GRNBoost2

GRNBoost2 uses gradient boosting regression trees to infer regulatory links between TFs and potential target genes. For each gene, it fits a model predicting that gene's expression from the expression of all TFs, then extracts feature importances as regulatory scores.

### 5.1 Python API

```python
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster

# Create a Dask distributed client for parallel computation
local_cluster = LocalCluster(n_workers=8, threads_per_worker=1)
client = Client(local_cluster)
print(f"Dask dashboard: {client.dashboard_link}")

# Run GRNBoost2
# This is the most computationally expensive step
adjacencies = grnboost2(
    expression_data=ex_matrix,
    tf_names=tf_names,
    verbose=True,
    client_or_address=client,
    seed=42
)

# Shut down the Dask client
client.close()
local_cluster.close()

# Inspect the output
print(f"Number of TF-target adjacencies: {len(adjacencies)}")
print(adjacencies.head(10))
```

```
   TF       target    importance
0  SPI1     CD14      234.56
1  SPI1     TYROBP    198.34
2  TBX21    IFNG      187.21
3  IRF8     BATF3     165.89
4  GATA3    IL4       152.77
5  FOXP3    IL2RA     148.92
6  SPI1     CSF1R     142.55
7  RORC     IL17A     138.44
8  TBX21    CXCR3     131.20
9  MAFB     APOE      125.67
```

### 5.2 CLI Alternative

For large datasets, the CLI is often more convenient and can be submitted as a batch job:

```bash
pyscenic grn \
    --num_workers 8 \
    --method grnboost2 \
    --output adjacencies.tsv \
    --seed 42 \
    pbmc_expression.loom \
    allTFs_hg38.txt
```

> **Note on input formats:** The CLI accepts `.loom` files. Convert your expression matrix:
>
> ```python
> import loompy
> # Create loom file from expression matrix
> row_attrs = {"Gene": np.array(ex_matrix.columns)}
> col_attrs = {"CellID": np.array(ex_matrix.index),
>              "cell_type": np.array(adata.obs["cell_type"])}
> loompy.create("pbmc_expression.loom",
>               ex_matrix.values.T,  # genes x cells
>               row_attrs, col_attrs)
> ```

### 5.3 Saving and Loading Adjacencies

```python
# Save adjacencies for later use
adjacencies.to_csv("adjacencies.tsv", sep="\t", index=False)

# Load later
adjacencies = pd.read_csv("adjacencies.tsv", sep="\t")
```

> **Tip:** GRNBoost2 runtime scales roughly linearly with the number of genes and cells. For 10k cells x 20k genes with 8 cores, expect 15-30 minutes. For 100k+ cells, consider subsampling to ~10k-20k representative cells per cell type.

---

## 6. Step 2: Regulon Prediction with cisTarget

This step prunes the co-expression modules by verifying that the TF binding motif is enriched among the cis-regulatory regions of predicted target genes. Only TF-target links supported by motif evidence are retained, converting raw modules into regulons.

### 6.1 Python API

```python
from pyscenic.prune import prune2df, df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase

# Load ranking databases
db_10kb = FeatherRankingDatabase(
    "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    name="hg38_10kb"
)
db_500bp = FeatherRankingDatabase(
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    name="hg38_500bp"
)
dbs = [db_10kb, db_500bp]

# Load motif annotations
motif_annotations_fname = "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# Derive modules from adjacencies (group targets by TF)
from pyscenic.utils import modules_from_adjacencies

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
print(f"Number of initial modules: {len(modules)}")
```

```python
# Prune modules using cisTarget (motif enrichment analysis)
df = prune2df(
    dbs,
    modules,
    motif_annotations_fname,
    num_workers=8
)

print(f"Pruned dataframe shape: {df.shape}")

# Convert to regulon objects
regulons = df2regulons(df)
print(f"Number of regulons: {len(regulons)}")

# Inspect regulons
for reg in sorted(regulons, key=lambda r: -len(r))[:10]:
    print(f"  {reg.name}: {len(reg)} target genes")
```

```
Number of regulons: 247

  SPI1(+): 423 target genes
  IRF8(+): 312 target genes
  TBX21(+): 287 target genes
  FOXP3(+): 245 target genes
  GATA3(+): 198 target genes
  MAFB(+): 176 target genes
  RORC(+): 142 target genes
  RUNX3(+): 138 target genes
  BATF(+): 127 target genes
  BCL6(+): 119 target genes
```

The `(+)` suffix indicates activating regulons (positive correlation between TF and targets). SCENIC can also detect `(-)` repressive regulons.

### 6.2 CLI Alternative

```bash
pyscenic ctx \
    adjacencies.tsv \
    hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname pbmc_expression.loom \
    --output regulons.csv \
    --num_workers 8 \
    --mask_dropouts
```

### 6.3 Examining Regulon Content

```python
# Find specific immune regulons
def find_regulon(regulons, tf_name):
    """Find a regulon by TF name."""
    for reg in regulons:
        if reg.transcription_factor == tf_name:
            return reg
    return None

# Examine the TBX21 (T-bet) regulon
tbet_reg = find_regulon(regulons, "TBX21")
if tbet_reg:
    targets = list(tbet_reg.gene2weight.keys())
    print(f"TBX21 regulon: {len(targets)} targets")
    print(f"Top targets: {targets[:20]}")
    # Expected: IFNG, CXCR3, CCL3, CCL4, IL12RB2, STAT4, etc.

# Examine the FOXP3 regulon
foxp3_reg = find_regulon(regulons, "FOXP3")
if foxp3_reg:
    targets = list(foxp3_reg.gene2weight.keys())
    print(f"FOXP3 regulon: {len(targets)} targets")
    print(f"Top targets: {targets[:20]}")
    # Expected: IL2RA (CD25), CTLA4, TNFRSF18 (GITR), IKZF2 (Helios), etc.
```

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure2-regulon-targets.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 2. Target gene networks for key immune TFs. Each node is a target gene; edge width reflects the regulatory weight. Only the top 20 targets per regulon are shown.
</div>

---

## 7. Step 3: Regulon Activity Scoring with AUCell

AUCell scores each cell for the activity of each regulon. The method ranks all genes in a cell by expression level, then computes the AUC for the recovery of regulon target genes among highly expressed genes. The result is a regulon activity score between 0 and 1 for each cell-regulon pair.

### 7.1 Python API

```python
from pyscenic.aucell import aucell

# Compute AUCell scores
# ex_matrix: cells x genes DataFrame with raw or normalized counts
auc_mtx = aucell(ex_matrix, regulons, num_workers=8)

print(f"AUCell matrix shape: {auc_mtx.shape}")
# (n_cells, n_regulons)
print(f"Columns (regulons): {auc_mtx.columns.tolist()[:10]}")
```

```
AUCell matrix shape: (9876, 247)
Columns (regulons): ['ARID5A(+)', 'BATF(+)', 'BCL6(+)', 'CEBPB(+)',
                      'ETS1(+)', 'FOXP3(+)', 'GATA3(+)', 'IRF4(+)',
                      'IRF8(+)', 'MAFB(+)']
```

### 7.2 CLI Alternative

```bash
pyscenic aucell \
    pbmc_expression.loom \
    regulons.csv \
    --output pbmc_scenic_output.loom \
    --num_workers 8
```

### 7.3 Store Results in AnnData

```python
# Add AUCell scores to AnnData for integrated analysis
adata_scenic = adata.copy()

# Ensure cells match
common_cells = adata_scenic.obs_names.intersection(auc_mtx.index)
auc_mtx_aligned = auc_mtx.loc[common_cells]
adata_scenic = adata_scenic[common_cells].copy()

# Store as obsm
adata_scenic.obsm["X_aucell"] = auc_mtx_aligned.values
adata_scenic.uns["regulon_names"] = auc_mtx_aligned.columns.tolist()

# Also add individual regulon scores to obs for easy plotting
for col in auc_mtx_aligned.columns:
    adata_scenic.obs[col] = auc_mtx_aligned[col].values

# Save
adata_scenic.write_h5ad("pbmc_scenic_results.h5ad")
print("Saved SCENIC results to AnnData.")
```

---

## 8. Visualization

With the AUCell matrix in hand, we can produce several informative visualizations.

### 8.1 Regulon Activity on UMAP

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

# Immune TF regulons to visualize
immune_regulons = ['TBX21(+)', 'GATA3(+)', 'RORC(+)', 'FOXP3(+)',
                   'SPI1(+)', 'MAFB(+)', 'IRF8(+)', 'EOMES(+)']

fig, axes = plt.subplots(2, 4, figsize=(24, 12))
axes = axes.flatten()

for idx, reg_name in enumerate(immune_regulons):
    ax = axes[idx]
    if reg_name in adata_scenic.obs.columns:
        sc.pl.umap(
            adata_scenic,
            color=reg_name,
            ax=ax,
            show=False,
            title=reg_name,
            color_map='viridis',
            vmin=0,
            frameon=False
        )
    else:
        ax.set_title(f"{reg_name}\n(not found)")
        ax.axis('off')

plt.tight_layout()
plt.savefig("regulon_umap.png", dpi=200, bbox_inches="tight")
plt.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure3-regulon-umap.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 3. Regulon activity scores projected onto UMAP embeddings. Each panel shows the AUCell score for a key immune TF regulon. TBX21 activity localizes to Th1 and CD8 effector clusters; FOXP3 activity is restricted to the Treg cluster; SPI1 and MAFB mark myeloid populations.
</div>

### 8.2 Regulon Activity Heatmap Across Cell Types

```python
import seaborn as sns

# Select regulons of interest
selected_regulons = [
    'TBX21(+)', 'EOMES(+)', 'RUNX3(+)',     # Th1/CTL
    'GATA3(+)', 'STAT6(+)',                   # Th2
    'RORC(+)', 'BATF(+)',                     # Th17
    'FOXP3(+)',                                # Treg
    'BCL6(+)',                                 # Tfh / B cells
    'SPI1(+)', 'MAFB(+)', 'CEBPB(+)',        # Myeloid
    'IRF8(+)', 'IRF4(+)',                     # DC
    'TCF7(+)', 'LEF1(+)',                     # Naive T
    'PAX5(+)', 'EBF1(+)',                     # B cells
]

# Filter to regulons actually present in our results
selected_regulons = [r for r in selected_regulons if r in auc_mtx.columns]

# Compute mean regulon activity per cell type
cell_types = adata_scenic.obs['cell_type']
mean_activity = pd.DataFrame(index=cell_types.unique())

for reg in selected_regulons:
    mean_activity[reg] = adata_scenic.obs.groupby('cell_type')[reg].mean()

# Z-score normalize across cell types for visualization
mean_activity_z = (mean_activity - mean_activity.mean()) / mean_activity.std()

# Plot
plt.figure(figsize=(14, 8))
sns.heatmap(
    mean_activity_z.T,
    cmap='RdBu_r',
    center=0,
    xticklabels=True,
    yticklabels=True,
    linewidths=0.5,
    cbar_kws={'label': 'Z-scored AUCell activity'}
)
plt.title("Regulon Activity Across Immune Cell Types", fontsize=14)
plt.xlabel("Cell Type", fontsize=12)
plt.ylabel("Regulon", fontsize=12)
plt.tight_layout()
plt.savefig("regulon_heatmap.png", dpi=200, bbox_inches="tight")
plt.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure4-regulon-heatmap.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 4. Heatmap of mean regulon activity (Z-scored) across immune cell types. Columns are cell types; rows are regulons. The pattern recapitulates known biology: TBX21 is active in Th1 cells, FOXP3 in Tregs, SPI1 and MAFB in monocytes/macrophages, and IRF8 in cDC1 cells.
</div>

### 8.3 Regulon Specificity Score (RSS)

The Regulon Specificity Score quantifies how specifically a regulon is active in a given cell type compared to all other cell types. This is particularly useful for identifying the "master regulators" of each population.

```python
from pyscenic.rss import regulon_specificity_scores

# Compute RSS
rss = regulon_specificity_scores(auc_mtx, cell_types)
print(f"RSS matrix shape: {rss.shape}")

# Plot top regulons per cell type
def plot_rss(rss, cell_type, top_n=10, ax=None):
    """Plot top regulons for a given cell type ranked by RSS."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 5))

    data = rss[cell_type].sort_values(ascending=False).head(top_n)
    ax.barh(range(len(data)), data.values, color='steelblue')
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data.index)
    ax.invert_yaxis()
    ax.set_xlabel("Regulon Specificity Score")
    ax.set_title(f"Top regulons: {cell_type}")
    return ax

# Plot RSS for selected immune cell types
cell_types_to_plot = ['CD4_Th1', 'Treg', 'CD14_Mono', 'cDC1']
fig, axes = plt.subplots(1, 4, figsize=(20, 5))

for ax, ct in zip(axes, cell_types_to_plot):
    plot_rss(rss, ct, top_n=8, ax=ax)

plt.tight_layout()
plt.savefig("regulon_rss.png", dpi=200, bbox_inches="tight")
plt.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure5-rss-barplot.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 5. Regulon Specificity Scores for selected immune cell types. TBX21 and EOMES rank highest in Th1 cells, FOXP3 dominates in Tregs, SPI1 and CEBPB lead in monocytes, and IRF8 is the top regulator in cDC1 cells.
</div>

### 8.4 Binary Regulon Activity and Clustering

SCENIC can also binarize regulon activity (on/off) using an AUC threshold, enabling regulon-based cell clustering that is independent of gene expression.

```python
from pyscenic.binarize import binarize

# Binarize AUCell matrix
binary_mtx, thresholds = binarize(auc_mtx)

print(f"Binary matrix shape: {binary_mtx.shape}")
print(f"Example thresholds:")
for reg in ['TBX21(+)', 'FOXP3(+)', 'SPI1(+)']:
    if reg in thresholds.index:
        print(f"  {reg}: {thresholds.loc[reg, 'threshold']:.4f}")

# Use binary regulon activity for UMAP
from sklearn.manifold import TSNE
import umap

# Compute UMAP on binary regulon activity
reducer = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=42)
embedding = reducer.fit_transform(binary_mtx.values)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Color by cell type
for ct in cell_types.unique():
    mask = cell_types == ct
    axes[0].scatter(embedding[mask, 0], embedding[mask, 1],
                    label=ct, s=3, alpha=0.5)
axes[0].set_title("Regulon-based UMAP (colored by cell type)")
axes[0].legend(markerscale=5, bbox_to_anchor=(1.05, 1), fontsize=8)
axes[0].set_xlabel("UMAP1")
axes[0].set_ylabel("UMAP2")

# Color by a specific regulon activity
reg = 'TBX21(+)'
if reg in auc_mtx.columns:
    sc_plot = axes[1].scatter(embedding[:, 0], embedding[:, 1],
                               c=auc_mtx[reg].values, s=3, cmap='viridis', alpha=0.5)
    plt.colorbar(sc_plot, ax=axes[1], label='AUCell score')
    axes[1].set_title(f"Regulon-based UMAP ({reg} activity)")
    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig("regulon_umap_binary.png", dpi=200, bbox_inches="tight")
plt.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/pyscenic/figure6-binary-umap.png" class="img-fluid rounded z-depth-1" zoomable=true %}

<div class="caption">
    Figure 6. UMAP embedding computed on binary regulon activity rather than gene expression. Left: cells colored by annotated cell type. Right: cells colored by TBX21 regulon activity. Regulon-based dimensionality reduction often produces cleaner separation between functionally distinct populations.
</div>

---

## 9. Immunological Interpretation

One of the most valuable outputs of SCENIC is the ability to map regulon activity to immune cell biology in a data-driven manner. Here we summarize expected findings and how to interpret them.

### 9.1 T Cell Compartment

```python
# Compare regulon activity between T cell subsets
t_cell_types = ['CD4_Naive', 'CD4_Th1', 'CD4_Th2', 'CD4_Th17', 'Treg',
                'CD8_Naive', 'CD8_Effector']
t_cell_regulons = ['TBX21(+)', 'GATA3(+)', 'RORC(+)', 'FOXP3(+)',
                   'EOMES(+)', 'RUNX3(+)', 'TCF7(+)', 'LEF1(+)',
                   'STAT4(+)', 'BATF(+)', 'BCL6(+)']

t_cell_regulons = [r for r in t_cell_regulons if r in auc_mtx.columns]

# Subset to T cells
t_mask = cell_types.isin(t_cell_types)
t_auc = auc_mtx.loc[t_mask, t_cell_regulons]
t_types = cell_types[t_mask]

# Statistical test: regulon activity differences
from scipy.stats import mannwhitneyu

print("Key regulon-cell type associations (Mann-Whitney U):")
print("-" * 60)

tests = [
    ('TBX21(+)', 'CD4_Th1'),
    ('GATA3(+)', 'CD4_Th2'),
    ('RORC(+)', 'CD4_Th17'),
    ('FOXP3(+)', 'Treg'),
    ('EOMES(+)', 'CD8_Effector'),
]

for reg, ct in tests:
    if reg in t_auc.columns:
        in_group = t_auc.loc[t_types == ct, reg]
        out_group = t_auc.loc[t_types != ct, reg]
        stat, pval = mannwhitneyu(in_group, out_group, alternative='greater')
        effect_size = in_group.mean() - out_group.mean()
        print(f"  {reg:15s} in {ct:15s}: "
              f"mean_diff={effect_size:.4f}, p={pval:.2e}")
```

Expected output:

```
Key regulon-cell type associations (Mann-Whitney U):
------------------------------------------------------------
  TBX21(+)        in CD4_Th1       : mean_diff=0.0823, p=1.23e-45
  GATA3(+)        in CD4_Th2       : mean_diff=0.0654, p=3.45e-38
  RORC(+)         in CD4_Th17      : mean_diff=0.0712, p=8.91e-29
  FOXP3(+)        in Treg          : mean_diff=0.0934, p=2.67e-52
  EOMES(+)        in CD8_Effector  : mean_diff=0.0567, p=5.12e-31
```

### 9.2 Myeloid Compartment

```python
# Myeloid-specific analysis
myeloid_types = ['CD14_Mono', 'CD16_Mono', 'cDC1', 'cDC2', 'pDC']
myeloid_regulons = ['SPI1(+)', 'MAFB(+)', 'CEBPB(+)', 'CEBPA(+)',
                    'IRF8(+)', 'IRF4(+)', 'STAT1(+)', 'STAT2(+)']

myeloid_regulons = [r for r in myeloid_regulons if r in auc_mtx.columns]

m_mask = cell_types.isin(myeloid_types)
m_auc = auc_mtx.loc[m_mask, myeloid_regulons]
m_types = cell_types[m_mask]

# Heatmap for myeloid compartment
mean_m = pd.DataFrame()
for reg in myeloid_regulons:
    mean_m[reg] = m_auc.groupby(m_types)[reg].mean()

mean_m_z = (mean_m - mean_m.mean()) / mean_m.std()

plt.figure(figsize=(8, 5))
sns.heatmap(mean_m_z.T, cmap='RdBu_r', center=0,
            annot=True, fmt='.2f', linewidths=0.5)
plt.title("Myeloid Regulon Activity")
plt.tight_layout()
plt.savefig("myeloid_regulon_heatmap.png", dpi=200, bbox_inches="tight")
plt.show()
```

**Key findings to expect:**

- **SPI1 (PU.1):** Broadly active across monocytes and DCs, reflecting its role as a master myeloid TF. Highest in CD14+ monocytes.
- **MAFB:** Preferentially active in monocytes/macrophages, consistent with its role in suppressing DC differentiation and promoting macrophage identity.
- **IRF8:** Highly specific to cDC1 cells. IRF8 is essential for cDC1 development and cross-presentation.
- **IRF4:** Active in cDC2 cells. IRF4 cooperates with BATF to drive the cDC2 program.
- **CEBPB:** Active in monocytes; C/EBPbeta drives emergency myelopoiesis and inflammatory macrophage programs.

### 9.3 Discovering Novel Regulons

Beyond confirming known biology, SCENIC frequently reveals unexpected TF-cell type associations:

```python
# Find the top 5 most cell-type-specific regulons for each population
for ct in cell_types.unique():
    top = rss[ct].sort_values(ascending=False).head(5)
    print(f"\n{ct}:")
    for reg, score in top.items():
        print(f"  {reg:20s} RSS={score:.4f}")
```

Novel findings might include:

- **BHLHE40** in tissue-resident memory T cells -- known to regulate cytokine production but not always included in standard TF lists for T cell subset annotations.
- **NR4A1/NR4A2** in recently activated T cells -- the NR4A family marks TCR signaling and is involved in T cell exhaustion programs.
- **MAF** in Th2 and Tfh cells -- c-Maf cooperates with GATA3 in Th2 and with BCL6 in Tfh differentiation.

---

## 10. Tips and Best Practices

### 10.1 Database Versions

```python
# Always verify your database and annotation versions match
import os

# Check files
db_files = [
    "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
]

for f in db_files:
    if os.path.exists(f):
        size_mb = os.path.getsize(f) / (1024 ** 2)
        print(f"  {f}: {size_mb:.0f} MB")
    else:
        print(f"  {f}: MISSING!")
```

**Version compatibility table:**

| Database version | Motif annotation | Species    |
| ---------------- | ---------------- | ---------- |
| `v10_clust`      | `v10nr_clust`    | hg38, mm10 |
| `v9`             | `v9nr`           | hg19, mm9  |

### 10.2 Gene Filtering Guidelines

```python
# Recommended filtering thresholds
# Too lenient: many noisy regulons, long runtime
# Too strict: miss lowly-expressed TFs (e.g., FOXP3)

# Conservative (faster, cleaner results)
sc.pp.filter_genes(adata_raw, min_cells=10)

# Liberal (retains rare TF programs, slower)
sc.pp.filter_genes(adata_raw, min_cells=3)

# Verify key TFs survived filtering
for tf in ['TBX21', 'GATA3', 'RORC', 'FOXP3', 'SPI1', 'MAFB', 'IRF8']:
    present = tf in adata_raw.var_names
    if not present:
        print(f"WARNING: {tf} was filtered out! Consider relaxing thresholds.")
```

### 10.3 Computational Requirements

| Dataset size          | GRNBoost2 | cisTarget | AUCell  | Total (8 cores) |
| --------------------- | --------- | --------- | ------- | --------------- |
| 5k cells, 15k genes   | ~10 min   | ~15 min   | ~2 min  | ~30 min         |
| 10k cells, 20k genes  | ~25 min   | ~20 min   | ~5 min  | ~50 min         |
| 50k cells, 20k genes  | ~2 hr     | ~20 min   | ~15 min | ~2.5 hr         |
| 100k cells, 25k genes | ~6 hr     | ~25 min   | ~30 min | ~7 hr           |

> **Memory:** Ranking databases require ~3-4 GB RAM each. Budget at least 16 GB for a typical run, 32 GB+ for large datasets. GRNBoost2 memory usage scales with `n_cells x n_TFs`.

### 10.4 Subsampling Strategy for Large Datasets

```python
# For datasets with >20k cells, subsample while preserving cell type proportions
def balanced_subsample(adata, cell_type_col='cell_type', max_cells_per_type=2000,
                       seed=42):
    """Subsample cells proportionally, capping each cell type."""
    np.random.seed(seed)
    indices = []
    for ct in adata.obs[cell_type_col].unique():
        ct_idx = np.where(adata.obs[cell_type_col] == ct)[0]
        n_sample = min(len(ct_idx), max_cells_per_type)
        sampled = np.random.choice(ct_idx, n_sample, replace=False)
        indices.extend(sampled)
    return adata[sorted(indices)].copy()

adata_sub = balanced_subsample(adata_raw, max_cells_per_type=2000)
print(f"Subsampled: {adata_sub.n_obs} cells")
print(adata_sub.obs['cell_type'].value_counts())
```

### 10.5 Running the Full Pipeline as a Script

For reproducibility, here is the entire pipeline condensed into a single script:

```python
#!/usr/bin/env python
"""
pySCENIC pipeline: GRNBoost2 -> cisTarget -> AUCell
Usage: python run_scenic.py --input data.h5ad --output scenic_results.h5ad
"""
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

def main(args):
    # 1. Load data
    adata = sc.read_h5ad(args.input)
    ex_matrix = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs_names, columns=adata.var_names
    )

    # 2. Load TFs
    tf_names = [l.strip() for l in open(args.tf_list)]
    tf_names = [t for t in tf_names if t in ex_matrix.columns]

    # 3. GRNBoost2
    print("Step 1/3: Running GRNBoost2...")
    cluster = LocalCluster(n_workers=args.n_workers, threads_per_worker=1)
    client = Client(cluster)
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names,
                            client_or_address=client, seed=42)
    client.close(); cluster.close()
    adjacencies.to_csv(f"{args.output_prefix}_adjacencies.tsv",
                       sep="\t", index=False)

    # 4. cisTarget
    print("Step 2/3: Running cisTarget...")
    dbs = [FeatherRankingDatabase(db, name=f"db_{i}")
           for i, db in enumerate(args.databases)]
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    df = prune2df(dbs, modules, args.motif_annotations,
                  num_workers=args.n_workers)
    regulons = df2regulons(df)
    print(f"  Found {len(regulons)} regulons")

    # 5. AUCell
    print("Step 3/3: Running AUCell...")
    auc_mtx = aucell(ex_matrix, regulons, num_workers=args.n_workers)

    # 6. Save results
    adata.obsm["X_aucell"] = auc_mtx.loc[adata.obs_names].values
    adata.uns["regulon_names"] = auc_mtx.columns.tolist()
    adata.write_h5ad(f"{args.output_prefix}_scenic.h5ad")
    auc_mtx.to_csv(f"{args.output_prefix}_aucell.csv")
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output_prefix", default="scenic")
    parser.add_argument("--tf_list", default="allTFs_hg38.txt")
    parser.add_argument("--databases", nargs="+", required=True)
    parser.add_argument("--motif_annotations", required=True)
    parser.add_argument("--n_workers", type=int, default=8)
    args = parser.parse_args()
    main(args)
```

```bash
# Run the full pipeline
python run_scenic.py \
    --input pbmc_10k_filtered.h5ad \
    --output_prefix pbmc \
    --tf_list allTFs_hg38.txt \
    --databases \
        hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --motif_annotations motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --n_workers 8
```

---

## 11. Key Takeaways

1. **SCENIC recovers biologically meaningful regulons.** The three-step pipeline (GRNBoost2, cisTarget, AUCell) combines expression correlation with cis-regulatory evidence, producing regulons that are more reliable than co-expression alone.

2. **Regulon activity is a powerful alternative to gene expression for cell identity.** AUCell scores capture the coordinated activity of entire gene programs, making them more robust to dropout and technical noise than individual gene expression values.

3. **Immunological validation is straightforward.** Known TF-cell type relationships (TBX21/Th1, GATA3/Th2, FOXP3/Treg, SPI1/myeloid, IRF8/cDC1) serve as positive controls. If these are not recovered, check your database versions and gene filtering.

4. **Database version matching is critical.** The most common pitfall is using mismatched ranking databases and motif annotations. Always verify that the version strings match (e.g., `v10_clust` databases with `v10nr_clust` annotations).

5. **Subsampling is acceptable for large datasets.** GRNBoost2 scales linearly with cell number, and the regulon inference step (cisTarget) operates on modules, not individual cells. Subsampling to 10-20k cells typically recovers the same regulons with dramatically reduced compute time.

6. **RSS (Regulon Specificity Score) is the most interpretable output.** While UMAP and heatmaps are useful for exploration, RSS directly ranks which regulons define each cell type, making it the go-to metric for biological interpretation.

7. **Consider SCENIC+ for multimodal data.** If you have paired scRNA-seq and scATAC-seq data, the SCENIC+ framework extends this pipeline with chromatin accessibility evidence, further improving regulon prediction.

---

## References

- Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. _Nature Methods_ 14, 1083-1086 (2017).
- Van de Sande, B. et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. _Nature Protocols_ 15, 2247-2276 (2020).
- Bravo Gonzalez-Blas, C. et al. SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. _Nature Methods_ 20, 1355-1367 (2023).
- Moerman, T. et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. _Bioinformatics_ 35, 2159-2161 (2019).
