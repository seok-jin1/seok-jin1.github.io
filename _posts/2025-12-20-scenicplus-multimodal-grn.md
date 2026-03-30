---
layout: post
title: "SCENIC+: Multimodal Gene Regulatory Networks from RNA + ATAC"
date: 2025-12-20
permalink: /blog/scenicplus-multimodal-grn/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - single-cell
  - multi-omics
  - epigenomics
  - tutorial
---

## Introduction

Gene regulatory networks (GRNs) describe how transcription factors (TFs) control target gene expression. Tools like **pySCENIC** infer these from scRNA-seq alone by correlating TF expression with gene modules and validating through promoter motif enrichment. This works well, but it cannot see the chromatin landscape that actually mediates regulation.

**SCENIC+** closes this gap by jointly analyzing paired scRNA-seq and scATAC-seq data --- typically from 10x Multiome --- to build **enhancer-driven gene regulatory networks**. It identifies accessible enhancer regions per cell, links them to target genes, scans for TF binding motifs, and assembles **eRegulons**: genes regulated by a TF through specific enhancer elements.

This tutorial walks through the complete SCENIC+ pipeline: from preprocessing paired multimodal data with Scanpy and pycisTopic, through building enhancer-gene links and running motif enrichment, to constructing and visualizing eRegulons. Along the way, we compare results with RNA-only pySCENIC to illustrate what the chromatin layer adds.

**Reference**: Bravo Gonzalez-Blas, C. et al. *SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks.* Nature Methods 20, 1355--1367 (2023).

---

## What SCENIC+ Adds Over pySCENIC

**pySCENIC (RNA-only)**: co-expression analysis (GRNBoost2) then promoter motif enrichment (cisTarget) then regulon scoring (AUCell).

**SCENIC+** extends this with chromatin accessibility:
1. Identify accessible regions per cell via **pycisTopic** (LDA topic modeling)
2. Link enhancers to genes using accessibility-expression correlation
3. Scan accessible regions for TF motifs via **cisTarget** databases
4. Correlate TF expression with both enhancer accessibility and target gene expression
5. Assemble **eRegulons** --- TF + enhancers + target genes with multi-layered evidence

The key insight: eRegulons tell you *which enhancers* a TF uses to regulate *which genes* in *which cell states*.

{% include figure.liquid loading="eager" path="assets/img/blog/scenicplus/figure1-scenic-vs-scenicplus.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. pySCENIC (RNA-only) versus SCENIC+ (multimodal) workflows. SCENIC+ adds enhancer-gene links and TF-to-enhancer binding to produce eRegulons.
</div>

---

## Setup and Data Requirements

SCENIC+ requires **paired** scRNA-seq and scATAC-seq from the same cells (e.g., 10x Multiome). You need:

- A gene expression count matrix (cells x genes) from the RNA modality
- An ATAC fragment file (`fragments.tsv.gz`) with cell barcodes matching the RNA data
- A reference genome annotation (e.g., hg38 for human, mm10 for mouse)
- cisTarget motif ranking databases for your genome assembly

The paired requirement is critical: each cell barcode must appear in both the RNA and ATAC data, so that expression and accessibility measurements can be directly linked. If you have separate scRNA-seq and scATAC-seq experiments, you would need to first integrate them (e.g., via bridge integration or label transfer), but this introduces noise compared to true multiome data.

### Installation

```python
# Create a dedicated conda environment
# conda create -n scenicplus python=3.11
# conda activate scenicplus

# pip install scenicplus pycisTopic pyscenic scanpy muon

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

work_dir = "/path/to/analysis"
tmp_dir = os.path.join(work_dir, "tmp")
os.makedirs(tmp_dir, exist_ok=True)

# cisTarget databases (download from https://resources.aertslab.org/cistarget/)
cistarget_dir = "/path/to/cistarget_databases"
rankings_db = os.path.join(cistarget_dir, "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather")
scores_db = os.path.join(cistarget_dir, "hg38_screen_v10_clust.regions_vs_motifs.scores.feather")
motif_annotation = os.path.join(cistarget_dir, "motifs-v10-nr.hgnc-m0.00001-o0.0.tbl")
```

---

## Preprocessing: scRNA-seq with Scanpy

The RNA component follows a standard Scanpy workflow. The critical point is that **cell barcodes must match between RNA and ATAC modalities** after quality control filtering, so keep track of which cells survive QC.

```python
adata_rna = sc.read_10x_h5(os.path.join(work_dir, "filtered_feature_bc_matrix.h5"), gex_only=True)
adata_rna.var_names_make_unique()

# QC and filtering
adata_rna.var["mt"] = adata_rna.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata_rna, min_genes=200)
sc.pp.filter_genes(adata_rna, min_cells=3)
adata_rna = adata_rna[adata_rna.obs["pct_counts_mt"] < 20, :].copy()

# Store raw counts for SCENIC+ before normalizing
adata_rna.raw = adata_rna.copy()

sc.pp.normalize_total(adata_rna, target_sum=1e4)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, n_top_genes=3000)
sc.pp.scale(adata_rna, max_value=10)
sc.tl.pca(adata_rna, n_comps=50)
sc.pp.neighbors(adata_rna, n_pcs=30)
sc.tl.umap(adata_rna)
sc.tl.leiden(adata_rna, resolution=0.8)

# Annotate cell types (example)
cell_type_map = {"0": "CD4 T", "1": "CD14 Mono", "2": "NK", "3": "B", "4": "CD8 T", "5": "CD16 Mono", "6": "DC"}
adata_rna.obs["cell_type"] = adata_rna.obs["leiden"].map(cell_type_map).astype("category")
```

{% include figure.liquid loading="eager" path="assets/img/blog/scenicplus/figure2-rna-umap.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. UMAP of scRNA-seq data colored by Leiden cluster (left) and annotated cell type (right).
</div>

---

## Preprocessing: scATAC-seq with pycisTopic

The ATAC component is processed using **pycisTopic**, which applies Latent Dirichlet Allocation (LDA) topic modeling to the binary cell-by-region matrix. Unlike simple peak calling, topic modeling discovers groups of co-accessible regions (topics) that collectively define cell states. Each cell gets a distribution over topics, and each topic gets a distribution over regions --- analogous to how LDA discovers topics in text corpora.

```python
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import impute_accessibility, normalize_scores, find_diff_features

fragments_file = os.path.join(work_dir, "atac_fragments.tsv.gz")
valid_barcodes = list(adata_rna.obs_names)

# Create cisTopic object from fragments
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=fragments_file,
    path_to_regions=os.path.join(work_dir, "consensus_peaks.bed"),
    valid_bc=valid_barcodes, n_cpu=8, project="multiome",
)

# Run LDA topic modeling
models = run_cgs_models(
    cistopic_obj, n_topics=[10, 20, 30, 40, 50],
    n_cpu=8, n_iter=300, random_state=42,
    alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False,
    save_path=os.path.join(tmp_dir, "lda_models"),
)

# Select best model (typically at the elbow of log-likelihood)
best_model = models[30]
cistopic_obj.add_LDA_model(best_model)

# Binarize topics and impute accessibility
region_bin_topics = binarize_topics(cistopic_obj, method="otsu", ntop=3000)
imputed_acc_obj = impute_accessibility(cistopic_obj, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

# Differentially accessible regions per cell type
cistopic_obj.cell_data["cell_type"] = adata_rna.obs.loc[cistopic_obj.cell_names, "cell_type"]
markers_dict = find_diff_features(
    cistopic_obj, imputed_acc_obj, variable="cell_type",
    adjpval_thr=0.05, log2fc_thr=1.0,
)
```

---

## Building Enhancer-to-Gene Links

A core innovation of SCENIC+ is establishing which enhancer regions regulate which genes. Rather than simply assigning each region to the nearest gene (which misses distal regulatory connections), SCENIC+ computes correlations between region accessibility and gene expression across cells. A region is linked to a gene if: (1) it falls within a defined genomic distance window around the gene's TSS, and (2) its accessibility is significantly correlated with the gene's expression across cells.

```python
from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships

# Define search space: potential region-gene pairs by distance
search_space = get_search_space(
    adata_rna, cistopic_obj, species="hsapiens", assembly="hg38",
    upstream=[1000, 150000], downstream=[1000, 150000],
)

# Calculate region-to-gene correlations
region_to_gene_df = calculate_regions_to_genes_relationships(
    adata_rna, cistopic_obj, search_space, imputed_acc_obj,
    temp_dir=tmp_dir, importance_scoring_method="GBM",
    correlation_scoring_method="SR", n_cpu=8,
)

# Filter for significant positive links
positive_links = region_to_gene_df[
    (region_to_gene_df["rho"] > 0.03) &
    (region_to_gene_df["importance"] > 0.005) &
    (region_to_gene_df["adj_pval"] < 0.05)
]
print(f"Significant enhancer-gene links: {len(positive_links)}")
print(f"Unique genes: {positive_links['Gene'].nunique()}, "
      f"Unique regions: {positive_links['Region'].nunique()}")
```

{% include figure.liquid loading="eager" path="assets/img/blog/scenicplus/figure3-enhancer-gene-links.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. Enhancer-to-gene link statistics. Left: distance distribution between linked regions and target gene TSSs. Right: correlation between region accessibility and gene expression for significant links.
</div>

---

## TF Motif Enrichment on Accessible Regions

With accessible regions identified per cell type and enhancer-gene links established, the next step determines which TFs are likely binding to these regions. SCENIC+ uses the cisTarget motif databases to scan regions for known TF binding motifs. A **cistrome** is the complete set of genomic regions bound by a given TF --- by intersecting cistromes with enhancer-gene links, SCENIC+ determines which TF binds which enhancer to regulate which gene.

```python
from scenicplus.TF_to_gene import calculate_TFs_to_genes_relationships
from scenicplus.motif_enrichment import run_cistarget, get_cistromes_per_region_set

# TF-to-gene expression correlations
tf_to_gene_df = calculate_TFs_to_genes_relationships(
    adata_rna, temp_dir=tmp_dir,
    importance_scoring_method="GBM", correlation_scoring_method="SR", n_cpu=8,
)

# Prepare region sets and run motif enrichment
region_sets = {**region_bin_topics, **markers_dict}
menr = run_cistarget(
    region_sets=region_sets, rankings_db_path=rankings_db,
    scores_db_path=scores_db, motif_annotation_path=motif_annotation,
    species="homo_sapiens", auc_threshold=0.005, nes_threshold=3.0,
    rank_threshold=0.05, n_cpu=8, temp_dir=tmp_dir,
)

# Extract cistromes: regions bound by each TF
cistromes = get_cistromes_per_region_set(
    menr, adata_rna, annotation_col="Direct_annot", annotation_type="TF",
)
print(f"Cistromes identified: {len(cistromes)}")
```

---

## Constructing eRegulons

The **eRegulon** is the central output of SCENIC+. Each eRegulon consists of:

- A **transcription factor** (TF)
- A set of **enhancer regions** containing the TF's binding motif, showing accessibility correlated with TF expression
- A set of **target genes** linked to those enhancers, showing expression correlated with TF expression

SCENIC+ constructs eRegulons by intersecting three data layers: TF-gene expression correlations, region-gene accessibility-expression links, and TF motif cistromes. Only regulatory connections supported by all three lines of evidence survive this integration.

```python
from scenicplus.eregulon_enrichment import build_eregulon_df, score_eregulons

# Build eRegulons integrating TF-gene, region-gene, and cistrome data
eregulon_df = build_eregulon_df(
    tf_to_gene=tf_to_gene_df, region_to_gene=region_to_gene_df,
    cistromes=cistromes, adata_rna=adata_rna,
    cistopic_obj=cistopic_obj, imputed_acc_obj=imputed_acc_obj,
    min_target_genes=5, min_regions_per_gene=1,
    adj_pval_thr=0.05, correlation_thr=0.03,
)
print(f"eRegulons: {eregulon_df['TF'].nunique()}, "
      f"TF-region-gene triplets: {len(eregulon_df)}")

# Score eRegulon activity per cell (RNA-based and ATAC-based)
eregulon_auc_rna = score_eregulons(
    adata_rna, eregulon_df, scoring_method="gene_based", auc_threshold=0.05, n_cpu=8,
)
eregulon_auc_atac = score_eregulons(
    adata_rna, eregulon_df, scoring_method="region_based",
    cistopic_obj=cistopic_obj, imputed_acc_obj=imputed_acc_obj,
    auc_threshold=0.05, n_cpu=8,
)

adata_rna.obsm["eRegulon_AUC_RNA"] = eregulon_auc_rna
adata_rna.obsm["eRegulon_AUC_ATAC"] = eregulon_auc_atac

print(f"RNA-based eRegulon scores: {eregulon_auc_rna.shape}")
print(f"ATAC-based eRegulon scores: {eregulon_auc_atac.shape}")
```

The dual scoring is a major advantage of SCENIC+. When a TF's eRegulon shows high activity in both the RNA-based and ATAC-based scores for a given cell type, you can be confident that the TF is genuinely active: its target genes are expressed *and* its enhancers are accessible. Discordance between the two scores can also be informative --- it may indicate primed but not yet active regulatory states.

---

## Visualization

SCENIC+ provides multiple ways to visualize eRegulon results. The most informative are heatmaps of eRegulon activity across cell types, UMAP overlays of individual eRegulon scores, and network graphs showing TF-enhancer-gene connections.

### eRegulon Activity Heatmap

A heatmap of eRegulon activity across cell types reveals which TFs drive each cell state from both expression and chromatin perspectives.

```python
from scenicplus.plotting import plot_eregulon_heatmap
from scenicplus.utils import get_top_eregulons_per_group

top_eregulons = get_top_eregulons_per_group(
    eregulon_auc_rna, adata_rna.obs["cell_type"], n_top=5, method="wilcoxon",
)

fig, axes = plt.subplots(1, 2, figsize=(20, 8))
plot_eregulon_heatmap(eregulon_auc_rna, adata_rna.obs["cell_type"],
    selected_eregulons=top_eregulons, ax=axes[0],
    title="eRegulon Activity (RNA)", vmin=-2, vmax=2, cmap="RdBu_r")
plot_eregulon_heatmap(eregulon_auc_atac, adata_rna.obs["cell_type"],
    selected_eregulons=top_eregulons, ax=axes[1],
    title="eRegulon Activity (ATAC)", vmin=-2, vmax=2, cmap="RdBu_r")
plt.tight_layout()
plt.savefig(os.path.join(work_dir, "eregulon_heatmap.png"), dpi=300, bbox_inches="tight")
```

{% include figure.liquid loading="eager" path="assets/img/blog/scenicplus/figure4-eregulon-heatmap.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. eRegulon activity heatmap. Left: scored from target gene expression. Right: scored from enhancer accessibility. Concordance between modalities strengthens regulatory assignments.
</div>

### eRegulon Activity on UMAP

Projecting individual eRegulon activity scores onto the UMAP embedding shows the spatial distribution of regulatory programs across cell populations. This is particularly useful for identifying TFs that are active in specific clusters or along differentiation trajectories.

```python
# UMAP overlay of selected eRegulon activities
selected_tfs = ["PAX5", "SPI1", "GATA3", "TBX21", "EOMES", "IRF4"]
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
for idx, tf in enumerate(selected_tfs):
    eregulon_name = f"{tf}_+"
    if eregulon_name in eregulon_auc_rna.columns:
        adata_rna.obs["_score"] = eregulon_auc_rna[eregulon_name].values
        sc.pl.umap(adata_rna, color="_score", title=f"{tf} eRegulon",
                   ax=axes.flat[idx], show=False, frameon=False, color_map="viridis")
plt.tight_layout()
plt.savefig(os.path.join(work_dir, "eregulon_umap.png"), dpi=300, bbox_inches="tight")
```

### Enhancer-Gene-TF Network Visualization

One of the most informative outputs of SCENIC+ is the network view showing how a specific TF connects to its enhancer regions and target genes. This visualization makes the multi-layered regulatory logic explicit.

```python
# Network visualization for a specific TF
import networkx as nx

tf_of_interest = "SPI1"
spi1_ereg = eregulon_df[eregulon_df["TF"] == tf_of_interest]
G = nx.DiGraph()
G.add_node(tf_of_interest, node_type="TF")
for _, row in spi1_ereg.iterrows():
    G.add_node(row["Gene"], node_type="gene")
    G.add_node(row["Region"], node_type="region")
    G.add_edge(tf_of_interest, row["Region"], edge_type="binding")
    G.add_edge(row["Region"], row["Gene"], edge_type="regulation")

fig, ax = plt.subplots(figsize=(14, 14))
pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
colors = {"TF": "#e74c3c", "region": "#3498db", "gene": "#2ecc71"}
sizes = {"TF": 800, "region": 100, "gene": 200}
nc = [colors[G.nodes[n].get("node_type", "gene")] for n in G.nodes()]
ns = [sizes[G.nodes[n].get("node_type", "gene")] for n in G.nodes()]
nx.draw(G, pos, ax=ax, node_color=nc, node_size=ns, edge_color="#cccccc",
        with_labels=False, arrows=True, width=0.5, alpha=0.8)
labels = {n: n for n in G.nodes() if G.nodes[n]["node_type"] != "region"}
nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
ax.set_title(f"{tf_of_interest} eRegulon Network", fontsize=16)
plt.savefig(os.path.join(work_dir, "eregulon_network.png"), dpi=300, bbox_inches="tight")
```

{% include figure.liquid loading="eager" path="assets/img/blog/scenicplus/figure5-eregulon-network.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. SPI1 (PU.1) eRegulon network. Red: TF. Blue: enhancer regions with SPI1 motifs. Green: target genes linked to those enhancers.
</div>

---

## Comparison with pySCENIC (RNA-Only)

To appreciate what the chromatin layer adds, it is useful to run pySCENIC on the same dataset and compare. The standard pySCENIC pipeline uses only gene expression data, relying on promoter motif enrichment rather than actual chromatin accessibility.

```python
from arboreto.algo import grnboost2
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

# Standard pySCENIC pipeline
adjacencies = grnboost2(adata_rna.to_df(), verbose=True, seed=42)
modules = list(modules_from_adjacencies(adjacencies, adata_rna.to_df()))
prune_df = prune2df(modules=modules, dbs=[rankings_db], motif_annotations=motif_annotation)
regulons = df2regulons(prune_df)
auc_mtx = aucell(adata_rna.to_df(), regulons, num_workers=8)

# Compare TF overlap
pyscenic_tfs = set([r.name.split("(")[0] for r in regulons])
scenicplus_tfs = set(eregulon_df["TF"].unique())
print(f"Shared TFs: {len(pyscenic_tfs & scenicplus_tfs)}")
print(f"pySCENIC only: {len(pyscenic_tfs - scenicplus_tfs)}")
print(f"SCENIC+ only: {len(scenicplus_tfs - pyscenic_tfs)}")

# Compare target genes for a shared TF
tf_compare = "SPI1"
pyscenic_targets = set(next(r.genes for r in regulons if r.name.startswith(tf_compare)))
scenicplus_targets = set(eregulon_df[eregulon_df["TF"] == tf_compare]["Gene"].unique())
print(f"\n{tf_compare} targets --- pySCENIC: {len(pyscenic_targets)}, "
      f"SCENIC+: {len(scenicplus_targets)}, Overlap: {len(pyscenic_targets & scenicplus_targets)}")
```

The comparison typically reveals several important differences:

- **SCENIC+ yields fewer but higher-confidence regulons** --- every link requires both expression and chromatin evidence, filtering out false positives from co-expression alone.
- **SCENIC+ recovers distal enhancer targets** missed by promoter-only motif enrichment. Genes regulated through enhancers 50--150 kb away from the TSS appear only in the multimodal analysis.
- **Cell-type specificity is sharper** because the same TF can use different enhancers in different cell types, and SCENIC+ captures this enhancer-level specificity.
- **eRegulons enable direct experimental follow-up** --- knowing the exact enhancer region guides CRISPR-based enhancer perturbation, reporter assays, or ChIP-seq validation.

---

## Practical Tips and Common Pitfalls

Before running SCENIC+ on your own data, keep these practical considerations in mind:

**Memory management**: The LDA topic modeling step loads the full cell-by-region binary matrix into memory. For datasets with more than 50,000 cells and 200,000 peaks, you may need 64+ GB RAM. Consider downsampling cells or filtering low-quality peaks before topic modeling.

**Choosing the number of topics**: There is no single correct answer. Run multiple models (e.g., 10 to 60 topics in steps of 10) and evaluate using log-likelihood, coherence, and biological interpretability. Too few topics merge distinct cell states; too many create redundant or noisy topics.

**cisTarget database selection**: Use region-based (not gene-based) cisTarget databases for SCENIC+. The database must match your genome assembly (hg38, mm10, etc.) and the type of regions you are analyzing (screen regions for broad peaks, encode for narrow peaks).

**Filtering thresholds**: The correlation and importance thresholds for enhancer-gene links and eRegulon construction significantly affect results. Start with the defaults from the SCENIC+ tutorials, then adjust based on the number of eRegulons recovered and their biological plausibility.

**Runtime expectations**: For a typical 10x Multiome dataset with 10,000--20,000 cells, expect the full pipeline to take 4--8 hours on a machine with 8 CPU cores and 64 GB RAM. The LDA modeling and region-to-gene correlation steps dominate runtime.

---

## Key Takeaways

1. **SCENIC+ extends pySCENIC with chromatin accessibility**, moving GRN inference from correlation-based prediction to mechanism-supported regulatory networks.

2. **eRegulons link TFs to specific enhancers and target genes**, providing three layers of evidence: TF expression, enhancer accessibility, and target gene expression.

3. **pycisTopic is essential** for the ATAC component. LDA topic modeling identifies co-accessible region sets more effectively than simple peak calling.

4. **Enhancer-gene links are computed, not assumed.** Accessibility-expression correlation captures distal regulatory connections that nearest-gene assignment would miss.

5. **Dual scoring (RNA + ATAC) increases confidence.** Concordant eRegulon activity across both modalities strengthens regulatory assignments.

6. **Paired multimodal data is required.** 10x Multiome is ideal; integrating separate scRNA-seq and scATAC-seq datasets is possible but noisier.

7. **Computational cost is substantial.** Plan for 32+ GB RAM and multiple CPU cores. LDA modeling and region-gene correlation are the bottlenecks.

8. **The API is actively evolving.** Always check the latest [SCENIC+ GitHub](https://github.com/aertslab/scenicplus) and [pycisTopic docs](https://pycistopic.readthedocs.io/) for current function signatures and module paths.
