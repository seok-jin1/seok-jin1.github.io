---
layout: post
title: "Single-Cell RNA-seq Analysis with Scanpy: A Practical Guide for Immunologists"
date: 2025-12-01
permalink: /blog/scrna-seq-scanpy-tutorial/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - single-cell
  - python
  - immunology
  - tutorial
---

## Introduction

If you have spent years gating flow cytometry plots and sorting immune cell populations, single-cell RNA sequencing (scRNA-seq) can feel like a paradigm shift. Instead of defining cell types by a handful of surface markers, you suddenly have the transcriptome-wide profile of every individual cell. For immunologists, this means you can capture the full heterogeneity of an immune response --- rare subsets, transitional states, and activation signatures --- without deciding in advance which markers to stain for.

In this tutorial, we will build a complete scRNA-seq analysis pipeline from raw count matrices to annotated cell types using **Scanpy**, the most widely used Python framework for single-cell analysis. We will work with the classic **PBMC 3k dataset** from 10x Genomics (2,700 peripheral blood mononuclear cells from a healthy donor), which is an ideal starting point because PBMCs contain well-characterized immune populations that most immunologists already know.

By the end of this post, you will have a working pipeline that covers quality control, normalization, dimensionality reduction, clustering, cell type annotation, and differential expression analysis.

---

## Environment Setup

First, install the required packages. Scanpy depends on several scientific Python libraries that will be installed automatically.

```python
# Install packages (run in terminal or notebook with !)
# pip install scanpy leidenalg

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configure scanpy
sc.settings.verbosity = 3  # hints, 0=errors, 1=warnings, 2=info, 3=hints
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False)
sc.logging.print_header()
```

You should see output like:

```
scanpy==1.10.x anndata==0.10.x ...
```

The `verbosity = 3` setting gives us detailed progress information, which is helpful for understanding what each function does under the hood.

---

## Loading Data

Scanpy provides a convenience function to load the PBMC 3k dataset directly. In practice, you would use `sc.read_10x_mtx()` to load Cell Ranger output from your own experiments.

```python
# Download and load the PBMC 3k dataset
adata = sc.datasets.pbmc3k()

print(adata)
# AnnData object with n_obs x n_vars = 2700 x 32738
#     var: 'gene_ids', 'feature_types'
```

The `AnnData` object is the central data structure in Scanpy. Think of it as a richly annotated matrix:

- **`adata.X`** --- the expression matrix (cells x genes), analogous to your raw count data
- **`adata.obs`** --- per-cell metadata (like a sample annotation sheet)
- **`adata.var`** --- per-gene metadata (gene names, IDs)
- **`adata.uns`** --- unstructured annotations (color maps, analysis parameters)
- **`adata.obsm`** --- cell embeddings (PCA, UMAP coordinates)

```python
# Inspect the raw counts
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Sparsity: {1 - (adata.X.nnz / (adata.n_obs * adata.n_vars)):.3f}")
# Number of cells: 2700
# Number of genes: 32738
# Sparsity: 0.974
```

That sparsity value (around 97%) is typical for droplet-based scRNA-seq. Most genes are not detected in most cells --- this is a fundamental feature of the technology, not a bug.

---

## Quality Control

Quality control in scRNA-seq serves the same purpose as your live/dead staining in flow cytometry: removing dead cells, doublets, and other artifacts before analysis. We use three key metrics:

1. **Number of genes detected per cell (`n_genes_by_counts`)** --- too few means empty droplets or debris; too many suggests doublets
2. **Total counts per cell (`total_counts`)** --- related to sequencing depth and cell size
3. **Percentage of mitochondrial reads (`pct_counts_mt`)** --- this is the critical one for immunologists

### Why mitochondrial percentage matters

When a cell is stressed or dying (think about what happens during tissue dissociation or FACS sorting), the cytoplasmic mRNA degrades while mitochondrial transcripts, protected within the double membrane, are relatively preserved. A high fraction of mitochondrial reads therefore indicates a damaged cell. This is the computational equivalent of your live/dead stain.

```python
# Flag mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())
```

### Visualize QC metrics

```python
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc_violin.png",
)
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure1-qc-violin.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. Violin plots of QC metrics. Left: number of genes detected per cell. Center: total UMI counts per cell. Right: percentage of mitochondrial reads. Cells with extreme values in any of these metrics are likely low-quality and should be filtered out.
</div>

```python
# Scatter plots to see relationships between QC metrics
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_mt_vs_counts.png")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_genes_vs_counts.png")
```

### Apply filters

```python
# Filter cells
adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
adata = adata[adata.obs.n_genes_by_counts > 200, :].copy()
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

print(f"Cells after filtering: {adata.n_obs}")
# Cells after filtering: ~2638
```

These thresholds are dataset-specific. For tissues with high metabolic activity (e.g., cardiac tissue), you may need a more lenient mitochondrial cutoff. For PBMCs, 5% is a reasonable threshold. Always look at the distributions before choosing cutoffs.

```python
# Filter genes: keep genes expressed in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)
print(f"Genes after filtering: {adata.n_vars}")
# Genes after filtering: ~13714
```

---

## Normalization and Feature Selection

### Library size normalization

Each cell is sequenced to a different depth. Without normalization, a cell with 10,000 total counts would appear to express every gene at higher levels than a cell with 1,000 counts, regardless of actual biology. We normalize so that every cell has the same total count, then log-transform to reduce the skewness of count data.

```python
# Store raw counts for later use (e.g., differential expression)
adata.layers["counts"] = adata.X.copy()

# Normalize to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform
sc.pp.log1p(adata)
```

The `target_sum=1e4` is a convention (counts per 10,000). After log-transformation, the data are on a more interpretable scale where fold-changes are roughly additive.

### Highly variable gene selection

Of the ~13,000 genes remaining, most show little variation across cells and contribute mainly noise. We select **highly variable genes (HVGs)** --- genes whose expression varies more than expected given their mean expression level. This is conceptually similar to choosing informative markers for a flow panel, except here we let the data decide.

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(f"Number of highly variable genes: {adata.var.highly_variable.sum()}")
# Number of highly variable genes: ~1838

sc.pl.highly_variable_genes(adata, save="_hvg.png")
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure2-highly-variable-genes.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. Highly variable gene selection. Black dots represent genes selected as highly variable. These genes show higher dispersion (variance) than expected for their mean expression level, indicating they capture meaningful biological variation rather than technical noise.
</div>

```python
# Keep the full data in .raw for later use, then subset to HVGs
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()

# Regress out effects of total counts and mitochondrial percentage
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

# Scale each gene to unit variance (clip at max value 10 to reduce outlier effects)
sc.pp.scale(adata, max_value=10)
```

---

## Dimensionality Reduction

Even after selecting HVGs, we still have ~1,800 dimensions. We need to reduce this for visualization and clustering. The standard workflow uses PCA first (linear reduction), then UMAP (nonlinear reduction for visualization).

### PCA

```python
sc.tl.pca(adata, svd_solver="arpack")

# Inspect variance explained
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save="_elbow.png")
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure3-pca-elbow.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. PCA variance ratio (elbow plot). Each point shows the proportion of variance explained by the corresponding principal component. The "elbow" where the curve flattens indicates the point of diminishing returns. For this dataset, approximately 10 PCs capture the major axes of variation.
</div>

The elbow plot helps us decide how many PCs to use downstream. For PBMCs, somewhere around 10 PCs is usually sufficient. Using too many PCs introduces noise; too few discards real signal.

### Compute neighbor graph and UMAP

```python
# Build the neighbor graph using the first 10 PCs
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)

# Compute UMAP embedding
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
           save="_qc_umap.png")
```

A note on UMAP interpretation: UMAP preserves local structure (nearby cells in UMAP space are truly similar), but distances between distant clusters and cluster sizes are not meaningful. Do not over-interpret the "gaps" between clusters or how spread out a cluster appears. Think of UMAP as a way to visualize neighborhoods, not geography.

t-SNE is an alternative that is better at preserving local structure but worse at preserving global relationships. In practice, UMAP has largely replaced t-SNE in the field due to better scalability and more interpretable global structure.

```python
# Optionally compute t-SNE for comparison
# sc.tl.tsne(adata)
# sc.pl.tsne(adata, color=["n_genes_by_counts"])
```

---

## Clustering

We use the **Leiden algorithm**, a community detection method that partitions the neighbor graph into groups of cells with dense internal connections. If you are familiar with graph theory, Leiden optimizes a modularity-like objective function. If not, think of it as finding groups of cells that are more similar to each other than to cells outside the group.

```python
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)

sc.pl.umap(adata, color=["leiden"], save="_leiden_clusters.png")
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure4-umap-clusters.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. UMAP visualization colored by Leiden clusters. Each color represents a distinct cluster identified by the algorithm. At resolution 1.0, the algorithm identifies major immune cell populations that we will annotate in the next section.
</div>

### Resolution parameter tuning

The `resolution` parameter controls the granularity of clustering. Higher resolution produces more clusters, lower resolution produces fewer. There is no universally correct resolution --- it depends on your biological question.

```python
# Compare different resolutions
for res in [0.4, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res{res}",
                 flavor="igraph", n_iterations=2)

sc.pl.umap(
    adata,
    color=["leiden_res0.4", "leiden_res0.8", "leiden_res1.0", "leiden_res1.5"],
    ncols=2,
    save="_resolution_comparison.png",
)
```

For PBMC data, resolution 1.0 typically separates the major immune populations well. If you are looking for fine subtypes (e.g., Th1 vs. Th2 vs. Th17), you may need higher resolution or subclustering.

```python
# Use resolution 1.0 as our primary clustering
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)
print(adata.obs["leiden"].value_counts())
```

---

## Cell Type Annotation

This is where your immunology knowledge becomes invaluable. Computational biologists can cluster cells, but correctly naming those clusters requires understanding the biology. We will use canonical marker genes that you already know from flow cytometry and immunohistochemistry.

### Define marker genes

```python
marker_genes = {
    "T cells (general)": ["CD3D", "CD3E", "CD3G"],
    "CD4+ T cells": ["CD4", "IL7R"],
    "CD8+ T cells": ["CD8A", "CD8B"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "CD14+ Monocytes": ["CD14", "LYZ", "S100A8", "S100A9"],
    "CD16+ Monocytes": ["FCGR3A", "MS4A7"],
    "NK cells": ["GNLY", "NKG7", "KLRD1"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A"],
    "Platelets": ["PPBP", "PF4"],
}

# Flatten for plotting
marker_genes_flat = [g for genes in marker_genes.values() for g in genes]
```

### Dot plot visualization

The dot plot is the single most useful visualization for cell type annotation. The dot size encodes the fraction of cells expressing each gene (similar to percent positive in flow), while the color encodes mean expression level (similar to MFI).

```python
sc.pl.dotplot(
    adata,
    var_names=marker_genes,
    groupby="leiden",
    standard_scale="var",
    save="_marker_dotplot.png",
)
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure5-marker-dotplot.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Dot plot of canonical marker genes across Leiden clusters. Dot size represents the fraction of cells in each cluster expressing the gene (analogous to percent positive in flow cytometry). Color intensity represents the mean expression level (analogous to MFI). This plot allows systematic mapping of clusters to known immune cell types.
</div>

```python
# Stacked violin plot as an alternative view
sc.pl.stacked_violin(
    adata,
    var_names=marker_genes_flat,
    groupby="leiden",
    swap_axes=False,
    save="_marker_violin.png",
)
```

### Manual annotation

Based on the dot plot, we can now assign cell type labels. The exact cluster numbers may vary between runs, so always check the marker expression pattern rather than relying on fixed cluster IDs.

```python
# Map cluster IDs to cell type names
# IMPORTANT: Verify these mappings against your own dot plot output.
# Cluster numbers can change between runs due to stochastic graph construction.
cluster_to_celltype = {
    "0": "CD4+ T cells",
    "1": "CD14+ Monocytes",
    "2": "CD4+ T cells",
    "3": "B cells",
    "4": "CD8+ T cells",
    "5": "CD16+ Monocytes",
    "6": "NK cells",
    "7": "Dendritic cells",
    "8": "Platelets",
}

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype).astype("category")

sc.pl.umap(adata, color=["cell_type"], save="_annotated_umap.png")
```

{% include figure.liquid path="assets/img/blog/scrna-seq/figure6-annotated-umap.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 6. UMAP visualization with manual cell type annotations. Clusters have been labeled based on canonical marker gene expression. The major PBMC populations --- T cells, B cells, monocytes, NK cells, dendritic cells, and platelets --- are clearly resolved.
</div>

A few things to notice:

- **CD4+ and CD8+ T cells** form distinct but nearby clusters, reflecting their shared T cell program with lineage-specific differences.
- **CD14+ and CD16+ monocytes** are well separated, consistent with their known transcriptional and functional differences.
- **NK cells** sit between T cells and monocytes, which makes biological sense given their shared cytotoxic programs with CD8+ T cells and innate immune features.

```python
# Check the distribution of cell types
print(adata.obs["cell_type"].value_counts())
```

---

## Differential Expression Analysis

Now that we have annotated cell types, we can identify genes that distinguish each population. This is analogous to comparing sorted populations by bulk RNA-seq, except we get the data for free from our single-cell experiment.

```python
# Use the Wilcoxon rank-sum test (recommended for scRNA-seq)
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon")

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save="_de_genes.png")
```

The Wilcoxon rank-sum test is a non-parametric test that works well with the zero-inflated, non-normal distributions typical of scRNA-seq data. It is preferred over the t-test for most applications.

```python
# Extract results as a DataFrame for a specific comparison
de_results = sc.get.rank_genes_groups_df(adata, group="NK cells")
print(de_results.head(15))
```

Expected output (top genes for NK cells):

```
         names     scores       pvals   pvals_adj  logfoldchanges
0         GNLY  28.xxxxx   0.000e+00   0.000e+00        6.xx
1         NKG7  27.xxxxx   0.000e+00   0.000e+00        5.xx
2        GZMB   24.xxxxx   0.000e+00   0.000e+00        5.xx
3        PRF1   21.xxxxx   0.000e+00   0.000e+00        4.xx
4        FGFBP2 20.xxxxx   0.000e+00   0.000e+00        5.xx
...
```

These results make perfect biological sense: GNLY (granulysin), NKG7, GZMB (granzyme B), and PRF1 (perforin) are all canonical cytotoxic effector molecules that define NK cell function.

### Filtering significant DE genes

```python
# Filter for significantly upregulated genes
de_nk = sc.get.rank_genes_groups_df(adata, group="NK cells")
de_nk_sig = de_nk[(de_nk["pvals_adj"] < 0.05) & (de_nk["logfoldchanges"] > 1.0)]

print(f"Number of significant upregulated genes in NK cells: {len(de_nk_sig)}")
print(de_nk_sig.head(20)[["names", "logfoldchanges", "pvals_adj"]])
```

### Comparing specific groups

You can also run pairwise comparisons, for example CD4+ vs. CD8+ T cells:

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="cell_type",
    groups=["CD8+ T cells"],
    reference="CD4+ T cells",
    method="wilcoxon",
)

de_cd8_vs_cd4 = sc.get.rank_genes_groups_df(adata, group="CD8+ T cells")
print(de_cd8_vs_cd4.head(10)[["names", "logfoldchanges", "pvals_adj"]])
```

You should see genes like CD8A, CD8B, GZMK, and GZMH upregulated in CD8+ T cells, while CD4+ T cell markers like IL7R appear at the bottom of the list. This type of pairwise comparison is especially useful when you have a specific biological question in mind.

```python
# Restore the full comparison for downstream use
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon")
```

---

## Trajectory Analysis

Not all biological processes fit neatly into discrete clusters. Differentiation, activation, and state transitions are continuous processes. Diffusion pseudotime (DPT) can order cells along a trajectory, which is particularly useful for studying:

- Naive to effector T cell differentiation
- Monocyte to macrophage/dendritic cell differentiation
- B cell maturation stages

```python
# Compute diffusion map
sc.tl.diffmap(adata)

# Choose a root cell (e.g., a naive T cell)
# First, find the cluster that likely contains naive T cells
# We look for cells with high IL7R expression (naive/memory T cells)
root_cell = adata[adata.obs["cell_type"] == "CD4+ T cells"].obs_names[0]
adata.uns["iroot"] = np.flatnonzero(adata.obs_names == root_cell)[0]

# Compute diffusion pseudotime
sc.tl.dpt(adata)

sc.pl.umap(adata, color=["dpt_pseudotime"], save="_pseudotime.png")
```

The pseudotime values represent a relative ordering of cells along a trajectory starting from the root cell. Cells with similar pseudotime values are at similar stages of the process. Keep in mind that DPT works best when the data contains a clear continuous trajectory; for discrete populations like resting PBMCs, the results may not be as informative as they would be for, say, a thymocyte differentiation dataset.

For more sophisticated trajectory analysis, consider tools like **scVelo** (RNA velocity), **CellRank**, or **Monocle 3**, which can infer directionality and branching dynamics.

---

## Exporting Results

### Save the AnnData object

```python
# Save the complete analyzed object (can be reloaded later)
adata.write("pbmc3k_analyzed.h5ad")

# Reload later with:
# adata = sc.read_h5ad("pbmc3k_analyzed.h5ad")
```

### Export to CSV for downstream analysis

```python
# Export cell metadata (including cluster assignments and cell types)
adata.obs.to_csv("pbmc3k_cell_metadata.csv")

# Export UMAP coordinates
umap_df = pd.DataFrame(
    adata.obsm["X_umap"],
    columns=["UMAP1", "UMAP2"],
    index=adata.obs_names,
)
umap_df["cell_type"] = adata.obs["cell_type"].values
umap_df.to_csv("pbmc3k_umap_coordinates.csv")

# Export DE results for all groups
for ct in adata.obs["cell_type"].cat.categories:
    de_df = sc.get.rank_genes_groups_df(adata, group=ct)
    de_df.to_csv(f"de_genes_{ct.replace(' ', '_').replace('+', 'pos')}.csv", index=False)

print("All results exported successfully.")
```

---

## Key Takeaways

### Pipeline summary

The standard scRNA-seq analysis workflow with Scanpy follows a clear sequence:

1. **Load data** --- `sc.read_10x_mtx()` or `sc.datasets`
2. **QC and filtering** --- remove dead cells (high mito%), empty droplets, doublets
3. **Normalize** --- library size normalization + log transformation
4. **Feature selection** --- identify highly variable genes
5. **Dimensionality reduction** --- PCA, then neighbor graph, then UMAP
6. **Clustering** --- Leiden algorithm on the neighbor graph
7. **Annotation** --- marker genes + domain knowledge
8. **Differential expression** --- Wilcoxon rank-sum test
9. **Export** --- save for downstream analysis or sharing

### Common pitfalls

- **Over-filtering**: Being too aggressive with QC thresholds can remove rare but real cell populations. Always visualize distributions before setting cutoffs.
- **Resolution obsession**: There is no single correct clustering resolution. Use biological knowledge to decide what level of granularity is appropriate for your question.
- **Over-interpreting UMAP**: UMAP distances between distant clusters are not meaningful. Two clusters far apart on UMAP are not necessarily more different than two clusters that are close.
- **Ignoring batch effects**: If your data comes from multiple samples, donors, or experiments, you need batch correction (Harmony, scVI, BBKNN) before clustering. Skipping this step is the most common source of artifacts.
- **Circular reasoning**: Do not use the same marker genes for both annotation and differential expression validation. If you annotated a cluster as "NK cells" based on NKG7, finding NKG7 as a top DE gene for that cluster is not independent validation.

### Next steps

Once you are comfortable with this basic pipeline, consider exploring:

- **Batch integration**: Harmony (`sc.external.pp.harmony_integrate`) or scVI for multi-sample experiments
- **RNA velocity**: scVelo for inferring transcriptional dynamics and cell fate
- **Cell-cell communication**: CellChat or LIANA for receptor-ligand analysis
- **Gene regulatory networks**: pySCENIC for transcription factor activity inference
- **Spatial transcriptomics**: Squidpy for analyzing spatially resolved data in the Scanpy ecosystem
- **Reference-based annotation**: CellTypist or scArches for automated cell type labeling

For further reading, the [Scanpy documentation](https://scanpy.readthedocs.io/) and the [single-cell best practices book](https://www.sc-best-practices.org/) are excellent resources.

---

*The complete code from this tutorial is available as a single script. All analysis was performed using the publicly available PBMC 3k dataset from 10x Genomics.*
