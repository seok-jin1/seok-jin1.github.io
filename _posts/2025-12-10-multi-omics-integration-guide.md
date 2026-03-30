---
layout: post
title: "Multi-Omics Integration for Immunology: MOFA+, scVI, and Beyond"
date: 2025-12-10
permalink: /blog/multi-omics-integration-guide/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - multi-omics
  - python
  - immunology
  - machine-learning
---

No single molecular measurement can fully capture the complexity of the immune system. Transcriptomics tells us what genes are being expressed, but not whether those transcripts become functional proteins. Proteomics reveals the effector molecules on the cell surface, but misses the regulatory logic encoded in chromatin accessibility. Epigenomics maps the landscape of gene regulation, but cannot directly tell us which genes are actively transcribed. Each modality offers a partial view -- like examining a sculpture from only one angle.

Multi-omics integration aims to combine these complementary perspectives into a unified picture. For immunology, this is especially powerful: immune responses are orchestrated across multiple molecular layers. A T cell's decision to become exhausted, for example, involves coordinated changes in gene expression (upregulation of *PDCD1*, *LAG3*, *HAVCR2*), surface protein display (PD-1, LAG-3, TIM-3), and epigenetic remodeling (chromatin closing at effector gene loci). Understanding these processes -- and ultimately developing precision medicine strategies based on multi-omics profiling -- requires computational tools that can jointly analyze data from multiple modalities.

In this tutorial, we will work through practical implementations of two major integration frameworks: **MOFA+** (Multi-Omics Factor Analysis) for latent factor discovery, and **scVI/totalVI** for deep generative modeling of single-cell multi-modal data. Along the way, we will discuss when to use each approach, walk through a realistic immunology case study, and compare the broader landscape of integration methods.

---

## Overview of Integration Strategies

Before diving into code, it is worth understanding the taxonomy of multi-omics integration approaches. The choice of method depends on your data type, sample size, and biological question.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure1-integration-strategies.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 1. Taxonomy of multi-omics integration strategies. Early integration concatenates raw features; intermediate integration learns shared latent spaces; late integration combines results from independent analyses; deep learning approaches use neural networks to learn joint representations.
</div>

- **Early integration (concatenation):** Concatenate features from all modalities into a single matrix. Simple but high-dimensional modalities (RNA, 20K genes) can dominate lower-dimensional ones (200 proteins).
- **Intermediate integration (latent factor models):** Methods like **MOFA+** learn a shared low-dimensional latent space, decomposing shared and modality-specific signals. The sweet spot for exploratory analysis.
- **Late integration (meta-analysis):** Analyze each modality independently, combine results at the interpretation level. Robust to technical differences but cannot discover cross-modal interactions.
- **Deep learning approaches:** Methods like **scVI**, **totalVI**, and **MultiVI** learn nonlinear latent representations via variational autoencoders. Scale well, handle batch effects naturally, but require more data and compute.

**When to use which:**

| Scenario | Recommended Approach |
|---|---|
| Small sample size (n < 50), bulk data | MOFA+ or late integration |
| CITE-seq (RNA + protein) | totalVI |
| scRNA-seq with batch effects | scVI |
| RNA + ATAC (multiome) | MultiVI or MOFA+ |
| Multiple bulk modalities, exploratory | MOFA+ |
| Large atlas-scale single-cell | scVI family (GPU recommended) |

---

## Part 1: MOFA+ (Multi-Omics Factor Analysis)

MOFA+ is a Bayesian factor analysis method that identifies latent factors of variation across multiple data modalities and, optionally, multiple sample groups. It decomposes multi-omics data into a set of factors, each with associated feature weights per modality, allowing you to ask: *which sources of variation are shared across modalities, and which are specific to one?*

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure2-mofa-model.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 2. MOFA+ model overview. The model decomposes multi-omics data matrices into shared latent factors (Z) and modality-specific weight matrices (W), capturing both shared and view-specific variation.
</div>

### Installation

```python
# Install MOFA+ Python engine and muon for data handling
pip install mofapy2 muon anndata scanpy
```

### Loading and Preparing Multi-Omics Data

We will use a paired RNA + protein (CITE-seq) dataset as our example. The `muon` package provides the `MuData` container for multi-modal data.

```python
import muon as mu
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load a CITE-seq dataset
# muon can read 10x Genomics multiome or CITE-seq outputs
# Here we use a preprocessed MuData object
mdata = mu.read("cite_seq_pbmc.h5mu")

# Inspect the MuData structure
print(mdata)
# MuData object with n_obs x n_vars
#   2 modalities
#     rna: AnnData object with n_obs x n_vars (RNA)
#     protein: AnnData object with n_obs x n_vars (protein)

# Access individual modalities
rna = mdata.mod["rna"]
protein = mdata.mod["protein"]

print(f"RNA features: {rna.n_vars}")
print(f"Protein features: {protein.n_vars}")
```

### Preprocessing Each Modality

Each modality requires its own normalization strategy:

```python
# --- RNA preprocessing ---
sc.pp.filter_genes(rna, min_cells=10)
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=3000)

# Subset to highly variable genes for MOFA+
rna_hvg = rna[:, rna.var["highly_variable"]].copy()

# --- Protein preprocessing ---
# CLR (centered log-ratio) normalization for protein data
from scipy.sparse import issparse

protein_data = protein.X.toarray() if issparse(protein.X) else protein.X.copy()
# Add pseudocount and apply CLR
protein_data = np.log1p(protein_data)
protein_data = protein_data - protein_data.mean(axis=1, keepdims=True)
protein.X = protein_data

# Update the MuData object
mdata.mod["rna"] = rna_hvg
mdata.mod["protein"] = protein
mu.pp.intersect_obs(mdata)  # Keep only shared cells

print(f"Cells after intersection: {mdata.n_obs}")
```

### Setting Up and Training MOFA+

```python
from muon import atac as ac
import muon

# Prepare MOFA+ model using muon's interface
mu.tl.mofa(
    mdata,
    n_factors=15,
    convergence_mode="slow",    # More iterations for better convergence
    seed=42,
    outfile="mofa_model.hdf5",  # Save trained model
    use_obs_names=True,
    gpu_mode=False,             # Set True if GPU available
)

# The trained model stores factors in mdata.obsm["X_mofa"]
print(f"Latent factors shape: {mdata.obsm['X_mofa'].shape}")
```

### Interpreting MOFA+ Results

The power of MOFA+ lies in its interpretability. Each factor has a clear biological interpretation through its feature weights and variance explained.

```python
# --- Variance explained per factor per modality ---
# This tells us which factors are shared vs modality-specific

# Extract variance explained from the trained model
from mofapy2.run.entry_point import entry_point

# Alternatively, use muon's plotting utilities
mu.pl.mofa(mdata, color="cell_type", frameon=False)
plt.savefig("mofa_factors_by_celltype.png", dpi=150, bbox_inches="tight")
plt.show()
```

```python
# --- Get feature weights for each factor ---
# Weights tell us which genes/proteins drive each factor

# Access weights from the MOFA model
weights = mdata.varm["LFs"]  # Loadings stored by muon

# For RNA modality - get top genes for Factor 1
rna_weights = pd.DataFrame(
    mdata.mod["rna"].varm["LFs"],
    index=mdata.mod["rna"].var_names,
    columns=[f"Factor{i+1}" for i in range(mdata.mod["rna"].varm["LFs"].shape[1])],
)

# Top positive and negative weights for Factor 1
factor1_weights = rna_weights["Factor1"].sort_values()
print("Top negative weights (Factor 1):")
print(factor1_weights.head(10))
print("\nTop positive weights (Factor 1):")
print(factor1_weights.tail(10))
```

```python
# --- Protein weights for the same factor ---
prot_weights = pd.DataFrame(
    mdata.mod["protein"].varm["LFs"],
    index=mdata.mod["protein"].var_names,
    columns=[f"Factor{i+1}" for i in range(mdata.mod["protein"].varm["LFs"].shape[1])],
)
print("Protein weights for Factor 1:")
print(prot_weights["Factor1"].sort_values())

# --- Visualize variance decomposition ---
Z = mdata.obsm["X_mofa"]
n_factors = Z.shape[1]

fig, ax = plt.subplots(figsize=(10, 4))
r2 = {}
for mod_name in ["rna", "protein"]:
    mod_X = mdata.mod[mod_name].X
    if issparse(mod_X):
        mod_X = mod_X.toarray()
    total_var = np.var(mod_X, axis=0).sum()
    r2[mod_name] = [
        np.var(Z[:, k:k+1] @ mdata.mod[mod_name].varm["LFs"][:, k:k+1].T, axis=0).sum()
        / total_var * 100
        for k in range(n_factors)
    ]

x = np.arange(n_factors)
ax.bar(x - 0.175, r2["rna"], 0.35, label="RNA", color="#4C72B0")
ax.bar(x + 0.175, r2["protein"], 0.35, label="Protein", color="#DD8452")
ax.set_xlabel("Factor"); ax.set_ylabel("Variance Explained (%)")
ax.set_title("MOFA+ Variance Decomposition")
ax.set_xticks(x); ax.set_xticklabels([f"F{i+1}" for i in range(n_factors)])
ax.legend(); plt.tight_layout()
plt.savefig("mofa_variance_decomposition.png", dpi=150, bbox_inches="tight")
plt.show()
```

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure3-mofa-variance.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 3. MOFA+ variance decomposition across modalities. Factors that explain variance in both RNA and protein represent shared biology; factors explaining variance in only one modality capture modality-specific signals.
</div>

---

## Part 2: scVI and totalVI

While MOFA+ provides interpretable linear factors, the **scvi-tools** ecosystem offers deep generative models for large-scale single-cell data:

- **scVI**: VAE for scRNA-seq -- learns latent representations while modeling library size, batch effects, and zero-inflation.
- **totalVI**: Extends scVI to jointly model RNA + surface protein (CITE-seq), learning a shared latent space from both modalities.

### Installation

```python
# Install scvi-tools (includes scVI, totalVI, MultiVI, and more)
pip install scvi-tools

# GPU support (recommended for large datasets):
# Ensure PyTorch is installed with CUDA support
# pip install torch --index-url https://download.pytorch.org/whl/cu121
```

> **Note on compute requirements:** scVI and totalVI train neural networks. For datasets under 50,000 cells, CPU training is feasible (10-30 minutes). For larger datasets, a GPU (NVIDIA with CUDA support) reduces training time from hours to minutes. Cloud platforms like Google Colab provide free GPU access.

### Loading CITE-seq Data

```python
import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the PBMC CITE-seq dataset from scvi-tools
# This is a real public dataset with ~5,000 PBMCs measured with
# both scRNA-seq and 14 surface protein markers
adata = scvi.data.pbmc_seurat_v4_cite_seq()

print(f"Cells: {adata.n_obs}")
print(f"RNA features: {adata.shape[1]}")
print(f"Protein features stored in adata.obsm['protein_expression']")

# Inspect protein names
protein_names = adata.obsm["protein_expression"].columns.tolist()
print(f"Proteins measured: {protein_names}")
```

### Preprocessing for totalVI

```python
# Filter genes
sc.pp.filter_genes(adata, min_counts=3)

# Select highly variable genes (totalVI requires raw counts)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    flavor="seurat_v3",   # Works with raw counts
    subset=True,
)

print(f"After HVG selection: {adata.n_vars} genes")
```

### Setting Up and Training totalVI

```python
# Register the AnnData object with totalVI
# This tells totalVI where to find RNA counts and protein data
scvi.model.TOTALVI.setup_anndata(
    adata,
    protein_expression_obsm_key="protein_expression",
    layer=None,           # Use adata.X for RNA counts
    batch_key=None,       # Add batch_key="batch" if you have batches
)

# Initialize the model
model = scvi.model.TOTALVI(
    adata,
    latent_distribution="normal",
    n_latent=20,          # Dimensionality of latent space
    n_layers_encoder=2,
    n_layers_decoder=2,
)

# Train the model
model.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_patience=15,
    train_size=0.9,       # 90% train, 10% validation
    batch_size=256,
    plan_kwargs={"lr": 4e-3},
)

# Check training convergence
model.history["elbo_train"].plot(label="Train")
model.history["elbo_validation"].plot(label="Validation")
plt.xlabel("Epoch"); plt.ylabel("ELBO"); plt.legend()
plt.savefig("totalvi_training.png", dpi=150, bbox_inches="tight"); plt.show()
```

### Extracting Latent Representations and Downstream Analysis

```python
# Get the latent representation (integrates both RNA and protein)
latent = model.get_latent_representation()
adata.obsm["X_totalVI"] = latent

# Compute neighbors and UMAP on the integrated latent space
sc.pp.neighbors(adata, use_rep="X_totalVI", n_neighbors=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)

# Visualize clusters
sc.pl.umap(adata, color=["leiden"], frameon=False, save="_totalvi_clusters.png")
```

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure4-totalvi-umap.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 4. UMAP visualization of PBMC CITE-seq data in the totalVI latent space. Clusters reflect joint RNA and protein information, often yielding finer resolution than RNA-only analysis.
</div>

### Differential Expression Across Modalities

```python
# DE between clusters -- totalVI tests RNA and protein jointly
de_results = model.differential_expression(groupby="leiden", group1="0", group2="3")

# Filter significant results (Bayes factor > 3 = strong evidence)
de_sig = de_results[(de_results["bayes_factor"] > 3.0) & (de_results["non_zeros_proportion1"] > 0.1)]
print(f"Significant DE features: {len(de_sig)}")
print(de_sig.sort_values("lfc_mean", ascending=False).head(10)[["lfc_mean", "bayes_factor"]])
```

### Protein Imputation from RNA

One of totalVI's most powerful features is the ability to impute (denoise) protein expression:

```python
# Get denoised protein expression -- leverages the joint model
# Foreground probability: 1 = true signal, 0 = background noise
protein_fg_prob = model.get_protein_foreground_probability(adata, protein_list=protein_names)

# Get normalized (denoised) protein values
_, protein_denoised = model.get_normalized_expression(
    adata, n_samples=25, return_mean=True,
)

# Compare raw vs denoised for CD4
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, data, title in zip(axes,
    [adata.obsm["protein_expression"]["CD4-TotalSeqB"], protein_denoised["CD4-TotalSeqB"]],
    ["CD4 Protein (Raw)", "CD4 Protein (totalVI Denoised)"]):
    ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
               c=data, s=1, cmap="viridis", alpha=0.5)
    ax.set_title(title); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
plt.tight_layout()
plt.savefig("totalvi_denoised_protein.png", dpi=150, bbox_inches="tight")
plt.show()
```

### Using scVI for RNA-Only Batch Correction

For scRNA-seq data from multiple batches without protein:

```python
# Setup with batch key -- scVI expects raw counts
scvi.model.SCVI.setup_anndata(adata_batched, layer="counts", batch_key="batch")

vae = scvi.model.SCVI(adata_batched, n_latent=30, n_layers=2, gene_likelihood="nb")
vae.train(max_epochs=200, early_stopping=True)

# Batch-corrected latent space for downstream analysis
adata_batched.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata_batched, use_rep="X_scVI")
sc.tl.umap(adata_batched)
sc.tl.leiden(adata_batched, resolution=1.0)
```

---

## Part 3: Practical Integration Workflow

Regardless of which method you choose, the integration workflow follows a consistent pipeline. Here is a step-by-step guide.

### Step 1: QC Each Modality Separately

```python
# RNA QC
sc.pp.calculate_qc_metrics(rna, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_cells(rna, max_genes=6000)
rna = rna[rna.obs["pct_counts_mt"] < 15, :].copy()

# Protein QC -- check isotype controls (should be low)
isotype_cols = [c for c in protein.var_names if "isotype" in c.lower()]
if isotype_cols:
    iso_counts = protein[:, isotype_cols].X
    if issparse(iso_counts): iso_counts = iso_counts.toarray()
    print(f"Median isotype control counts: {np.median(iso_counts, axis=0)}")
```

### Step 2: Normalize Per Modality

```python
# RNA: log-normalization for MOFA+; raw counts for scVI/totalVI
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, n_top_genes=3000)

# Protein: CLR for MOFA+; raw counts for totalVI (handles internally)
# ATAC: TF-IDF or binary peak matrix + LSI
```

### Steps 3-5: Integrate and Analyze

```python
# After choosing your method (see guidance above) and training:
sc.pp.neighbors(adata, use_rep="X_totalVI", n_neighbors=20, metric="cosine")
sc.tl.leiden(adata, resolution=0.8)
sc.tl.umap(adata)
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=10, save="_de_genes.png")
```

### Cross-Modality Correlation Analysis

A key validation step: do RNA and protein markers agree?

```python
from scipy.stats import spearmanr

gene_protein_pairs = {
    "CD4": "CD4-TotalSeqB", "CD8A": "CD8-TotalSeqB",
    "CD19": "CD19-TotalSeqB", "NCAM1": "CD56-TotalSeqB", "CD14": "CD14-TotalSeqB",
}

for gene, prot_name in gene_protein_pairs.items():
    if gene in adata.var_names:
        rna_expr = adata[:, gene].X.toarray().flatten() if issparse(adata.X) else adata[:, gene].X.flatten()
        prot_expr = adata.obsm["protein_expression"][prot_name].values
        rho, pval = spearmanr(rna_expr, prot_expr)
        print(f"{gene} vs {prot_name}: rho={rho:.3f}, p={pval:.2e}")

# Typical findings:
# CD4 vs CD4-protein: moderate (rho ~0.3-0.5)
# CD14 vs CD14-protein: strong (rho ~0.6-0.8)
# Discrepancy is biological -- protein half-life differs from mRNA
```

---

## Immunology Case Study: Immune Cell Profiling

To illustrate the power of multi-omics integration, let us walk through how it resolves biological questions that single-modality analysis cannot answer.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure5-immunology-case-study.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 5. Multi-omics integration reveals immune cell states invisible to single-modality analysis. Top: T cell exhaustion requires concordant evidence from transcriptome, surface proteome, and epigenome. Bottom: Myeloid polarization markers across RNA and protein.
</div>

### T Cell Exhaustion: A Multi-Layer Phenomenon

T cell exhaustion is a progressive loss of effector function during chronic infection and in the tumor microenvironment -- a prime example of why multi-omics matters:

- **Transcriptomics** detects *PDCD1*, *LAG3*, *HAVCR2*, *TOX* upregulation, but transcript levels do not always match surface protein abundance.
- **Surface proteomics** (CITE-seq) directly measures PD-1, LAG-3, TIM-3 protein -- what matters for checkpoint blockade -- but panels are limited to ~200 markers.
- **Epigenomics** (scATAC-seq) reveals stable *PDCD1* locus accessibility, distinguishing truly exhausted T cells from transiently activated ones.

With multi-omics integration:

```python
# Conceptual analysis -- identifying exhausted T cells
# using totalVI latent space from CITE-seq data

# After clustering in totalVI latent space:
# 1. Identify clusters with high PD-1 protein
pd1_expr = adata.obsm["protein_expression"]["PD1-TotalSeqB"]
adata.obs["PD1_protein_high"] = pd1_expr > np.percentile(pd1_expr, 75)

# 2. Check transcriptional exhaustion signature
exhaustion_genes = ["PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT", "CTLA4"]
sc.tl.score_genes(adata, gene_list=exhaustion_genes, score_name="exhaustion_score")

# 3. Cross-reference: cells with BOTH high PD-1 protein AND
#    high exhaustion transcriptional score are truly exhausted
adata.obs["exhausted"] = (
    (adata.obs["exhaustion_score"] > 0.5) &
    (adata.obs["PD1_protein_high"])
)

# This combined criterion is more specific than either alone
print(f"High PD-1 protein only: {adata.obs['PD1_protein_high'].sum()} cells")
print(f"High exhaustion RNA only: {(adata.obs['exhaustion_score'] > 0.5).sum()} cells")
print(f"Both (truly exhausted): {adata.obs['exhausted'].sum()} cells")
```

### Myeloid Polarization and Treg Identification

Macrophage M1/M2 polarization is another case where RNA markers (*TNF*, *IL1B* vs *MRC1*, *CD163*) combined with surface proteins (CD80, CD86 vs CD206) give cleaner separation than either modality alone.

Similarly, regulatory T cells (Tregs) are classically defined by a multi-modal triad:

**FOXP3** (transcription factor, RNA) + **CD25-high** (IL-2Ra, protein) + **CD127-low** (IL-7Ra, protein). This definition inherently requires multi-modal data:

```python
# Treg identification in CITE-seq data
cd25_expr = adata.obsm["protein_expression"]["CD25-TotalSeqB"]
cd127_expr = adata.obsm["protein_expression"]["CD127-TotalSeqB"]

# FOXP3 from RNA
foxp3_expr = np.zeros(adata.n_obs)
if "FOXP3" in adata.var_names:
    foxp3_raw = adata[:, "FOXP3"].X
    foxp3_expr = foxp3_raw.toarray().flatten() if issparse(foxp3_raw) else foxp3_raw.flatten()

# Classical Treg definition
adata.obs["Treg"] = (
    (foxp3_expr > 0) &
    (cd25_expr > np.percentile(cd25_expr, 70)) &
    (cd127_expr < np.percentile(cd127_expr, 30))
)

print(f"Tregs identified: {adata.obs['Treg'].sum()} cells")
```

---

## Comparison of Methods

The multi-omics integration field is growing rapidly. Here is a practical comparison of the major methods:

| Method | Input Data | Approach | Batch Correction | GPU Required | Key Strength |
|--------|-----------|----------|-------------------|--------------|--------------|
| **MOFA+** | Any multi-omics | Bayesian factor model | Via groups | No | Interpretable factors |
| **scVI** | scRNA-seq | VAE | Yes | Recommended | Scalable, probabilistic |
| **totalVI** | RNA + protein | VAE | Yes | Recommended | Joint RNA-protein model |
| **MultiVI** | RNA + ATAC | VAE | Yes | Recommended | Multiome integration |
| **Seurat WNN** | Any paired | Weighted nearest neighbors | Separate step | No | Easy to use (R) |
| **GLUE** | Any unpaired | Graph-linked VAE | Yes | Yes | Unpaired data integration |
| **scGLUE** | RNA + ATAC (unpaired) | Knowledge-guided VAE | Yes | Yes | Uses regulatory knowledge |

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/multi-omics/figure6-method-comparison.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 6. Decision flowchart for selecting a multi-omics integration method based on data type, sample size, and whether modalities are paired or unpaired.
</div>

**Quick guide:** CITE-seq --> totalVI. Bulk multi-omics, small n --> MOFA+. Large single-cell atlas --> scVI/MultiVI with GPU. Unpaired modalities --> GLUE.

---

## Limitations and Future Directions

### Current Challenges

- **Data availability.** True multi-modal single-cell measurements remain expensive. CITE-seq panels are limited to ~200 proteins. Multiome (RNA + ATAC) has higher dropout rates than either modality alone.
- **Batch effects across modalities.** Technical variation can be confounded with biological variation, especially when modalities are measured on different platforms. No method handles this perfectly.
- **Scalability.** As atlases grow to millions of cells, GPU memory becomes a bottleneck. Approximate or streaming methods are an active area of development.
- **Ground truth.** Evaluating integration quality is difficult without ground truth. Metrics like silhouette score on known cell types are proxies, not definitive measures.

### Emerging Frontiers

- **Spatial multi-omics.** Technologies like MERFISH and spatial CITE-seq add spatial coordinates. Methods like SpatialGLUE and GraphST are leading integration efforts.
- **Foundation models.** Large pretrained models (scGPT, Geneformer, scFoundation) learn general cell state representations that can be fine-tuned for multi-omics integration. Whether they will replace task-specific methods like totalVI remains to be seen.
- **Temporal multi-omics.** Methods for trajectory inference on multi-modal data (e.g., MultiVelo for RNA + ATAC velocity) capture how profiles change during immune responses or disease progression.

---

## Key Takeaways

**Pipeline summary:**

1. Generate or obtain multi-modal data (CITE-seq, Multiome, paired bulk omics)
2. Perform QC on each modality independently
3. Normalize each modality with appropriate methods (log-normalization for RNA, CLR for protein, TF-IDF for ATAC)
4. Select an integration method based on data type and scale (MOFA+ for small/interpretable, totalVI for CITE-seq, scVI for batch correction)
5. Train the model and extract the integrated latent representation
6. Perform downstream analyses (clustering, DE, trajectory) on the integrated space
7. Validate: cross-modality markers should be concordant

**For beginners:**

If you are new to multi-omics, start with totalVI on a CITE-seq dataset. The [scvi-tools tutorials](https://docs.scvi-tools.org/en/stable/tutorials/index.html) provide excellent step-by-step notebooks. CITE-seq data combines the familiarity of scRNA-seq with the added dimension of surface proteins, making it a gentle entry point into multi-modal analysis.

**Resources:**

- [scvi-tools docs](https://docs.scvi-tools.org) | [MOFA+ tutorial](https://biofam.github.io/MOFA2) | [muon docs](https://muon.readthedocs.io)
- Argelaguet et al. (2020), *Genome Biology* -- MOFA+
- Gayoso et al. (2021), *Nature Methods* -- totalVI
- Hao et al. (2021), *Cell* -- Seurat v4 WNN
- Cao & Gao (2022), *Nature Biotechnology* -- GLUE

Multi-omics integration is not just a computational exercise -- it is a lens through which we can see the immune system as it truly operates: a coordinated, multi-layered molecular machine. The methods described here are tools for building that understanding, one dataset at a time.
