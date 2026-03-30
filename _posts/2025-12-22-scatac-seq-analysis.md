---
layout: post
title: "scATAC-seq Analysis: Chromatin Accessibility with ArchR and Signac"
date: 2025-12-22
permalink: /blog/scatac-seq-analysis/
published: true
categories: [tutorial]
tags:
  - bioinformatics
  - single-cell
  - epigenomics
  - R
  - tutorial
---

Single-cell ATAC-seq (scATAC-seq) measures chromatin accessibility at single-cell resolution, revealing the regulatory landscape that governs gene expression. Unlike scRNA-seq, which captures transcriptional output, scATAC-seq captures the **potential** for transcription by identifying open chromatin regions where transcription factors can bind. This distinction is critical for understanding how cell identity is established and maintained through epigenetic regulation.

In immunology, chromatin accessibility profiling has proven transformative. Different immune cell types -- naive T cells, effector T cells, regulatory T cells, macrophages, B cells -- share a common genome but deploy vastly different regulatory programs. scATAC-seq reveals these cell-type-specific regulatory elements: enhancers that activate cytokine genes in Th1 cells, silencers that repress effector programs in Tregs, and poised promoters in naive cells awaiting activation signals.

This tutorial walks through a complete scATAC-seq analysis pipeline using **ArchR** and **Signac** in R, with brief Python alternatives. We will use a PBMC dataset as a running example, focusing on identifying immune cell populations and their regulatory landscapes.

---

## 1. What scATAC-seq Measures

The Assay for Transposase-Accessible Chromatin (ATAC-seq) uses a hyperactive Tn5 transposase to insert sequencing adapters into regions of open chromatin. In the single-cell version, each cell is individually barcoded before tagmentation, producing a library where each fragment is associated with a specific cell.

The resulting data is inherently **sparse and binary**: for any given genomic region in a single cell, chromatin is either accessible (fragment detected) or not. This sparsity -- far more extreme than scRNA-seq -- demands specialized computational methods.

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure1-atac-overview.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. Overview of the scATAC-seq assay. Tn5 transposase inserts adapters into open chromatin regions. After single-cell barcoding and sequencing, fragments are mapped back to the genome to identify accessible regions per cell.
</div>

---

## 2. Data Formats

Before diving into analysis, it is important to understand the key input files:

**Fragments file** (`fragments.tsv.gz`): The primary input for most scATAC-seq tools. Each row represents a single Tn5 insertion fragment with columns: chromosome, start, end, cell barcode, and duplicate count.

```
chr1    10073   10344   ACGTACGTACGTACGT-1    1
chr1    10073   10496   TGCATGCATGCATGCA-1    2
chr1    10109   10344   GCTAGCTAGCTAGCTA-1    1
```

**Peak-by-cell matrix**: A sparse matrix where rows are genomic peaks (open chromatin regions) and columns are cells. Values indicate the number of fragments overlapping each peak per cell.

**Cell barcodes** (`singlecell.csv` or `barcodes.tsv`): A list of valid cell barcodes that passed Cell Ranger ATAC quality filters, often with per-barcode QC metrics.

---

## 3. ArchR Workflow

[ArchR](https://www.archrproject.com/) is a scalable R package for scATAC-seq analysis that uses Arrow files (an HDF5-based format) to handle large datasets with minimal memory footprint. It is well suited for datasets with hundreds of thousands of cells.

### 3.1 Installation and Setup

```r
# Install ArchR
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())

# Load and set parameters
library(ArchR)
set.seed(42)

# Set the number of threads for parallel processing
addArchRThreads(threads = 16)

# Set the reference genome
# Use hg38 for human, mm10 for mouse
addArchRGenome("hg38")
```

### 3.2 Creating an ArchR Project from Fragments

The first step is to convert fragments files into Arrow files, which are ArchR's on-disk data format optimized for random access.

```r
# Define input fragment files
# For a typical 10x Genomics experiment:
inputFiles <- c(
  "PBMC_10k" = "/path/to/pbmc_10k_fragments.tsv.gz"
)

# Create Arrow files from fragments
# This step performs initial QC filtering
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,           # Minimum TSS enrichment score
  minFrags = 1000,       # Minimum number of fragments per cell
  addTileMat = TRUE,     # Add 500-bp tile matrix
  addGeneScoreMat = TRUE # Add gene activity score matrix
)

# Infer doublets using ArchR's projection-based method
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,                # Number of nearest neighbors
  knnMethod = "UMAP",
  LSIMethod = 1
)

# Create the ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "PBMC_ArchR",
  copyArrows = TRUE
)

# Filter doublets
proj <- filterDoublets(proj)

# Check project summary
proj
# Output:
# class: ArchRProject
# outputDirectory: PBMC_ArchR
# Samples: PBMC_10k
# nCells: 8,923
```

### 3.3 Quality Control

Two key QC metrics for scATAC-seq are the **TSS enrichment score** and the **fragment size distribution**.

TSS enrichment measures the enrichment of fragments around transcription start sites relative to flanking regions. High-quality cells show strong enrichment (score > 4), indicating that the assay preferentially captured open chromatin at promoters.

The fragment size distribution should show a clear nucleosomal banding pattern: a prominent sub-nucleosomal peak (< 147 bp) from accessible regions, followed by periodic peaks at mono-, di-, and tri-nucleosomal sizes (~200 bp, ~400 bp, ~600 bp).

```r
# Plot TSS enrichment vs. number of fragments (log10)
# This is the most informative QC plot for scATAC-seq
p1 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)

p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
)

# Fragment size distribution
p3 <- plotFragmentSizes(ArchRProj = proj)

# TSS enrichment profile
p4 <- plotTSSEnrichment(ArchRProj = proj)

# Save plots
plotPDF(p1, p2, p3, p4,
        name = "QC-Plots",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
```

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure2-qc-plots.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. Quality control metrics for scATAC-seq data. (A) TSS enrichment score distribution. (B) Log10 fragment count distribution. (C) Fragment size distribution showing nucleosomal banding pattern. (D) TSS enrichment profile showing enrichment of fragments around transcription start sites.
</div>

### 3.4 Dimensionality Reduction with Iterative LSI

ArchR uses **Iterative Latent Semantic Indexing (LSI)** for dimensionality reduction. LSI applies TF-IDF normalization followed by SVD, conceptually similar to PCA but adapted for the sparse, binary nature of scATAC-seq data. The iterative approach refines the feature set across rounds, selecting the most variable features in each iteration.

```r
# Iterative LSI on the tile matrix (500-bp genome-wide bins)
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  seed = 42
)
```

### 3.5 Clustering and UMAP Visualization

```r
# Add clusters using the Leiden algorithm on the LSI-reduced dimensions
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",      # Uses Seurat's FindClusters under the hood
  name = "Clusters",
  resolution = 0.8,
  seed = 42
)

# Add UMAP embedding
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  seed = 42
)

# Plot UMAP colored by clusters
p_clusters <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

# Plot UMAP colored by sample
p_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)

plotPDF(p_clusters, p_sample,
        name = "UMAP-Clusters",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
```

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure3-umap-clusters.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. UMAP visualization of PBMC scATAC-seq data colored by cluster identity. Distinct clusters correspond to different immune cell populations identified by their chromatin accessibility profiles.
</div>

### 3.6 Gene Activity Scores

Since scATAC-seq does not directly measure gene expression, ArchR computes **gene activity scores** as a proxy. These scores aggregate chromatin accessibility signals across the gene body and promoter region, weighted by distance from the TSS. Genes with accessible promoters and regulatory elements receive higher scores.

```r
# Gene activity scores were already computed during Arrow file creation.
# We can visualize them on the UMAP:

# Immune cell marker genes
markerGenes <- c(
  "CD3D", "CD3E",       # T cells
  "CD4",                 # CD4+ T cells
  "CD8A", "CD8B",        # CD8+ T cells
  "FOXP3", "IL2RA",      # Regulatory T cells
  "MS4A1", "CD79A",      # B cells
  "CD14", "LYZ",         # Monocytes (CD14+)
  "FCGR3A",              # Monocytes (CD16+)
  "GNLY", "NKG7",        # NK cells
  "FCER1A", "CST3"       # Dendritic cells
)

p_gene <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj)
)

# Arrange marker gene plots in a grid
p_gene_grid <- lapply(p_gene, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 4), p_gene_grid))
```

To use imputation (which smooths sparse gene scores for cleaner visualization), add imputation weights first:

```r
# Add imputation weights using MAGIC-based approach
proj <- addImputeWeights(proj)

# Now re-plot with imputed gene scores
p_imputed <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)
```

### 3.7 Peak Calling with MACS2

ArchR calls peaks on pseudo-bulk replicates -- aggregated fragments from cells within each cluster. This approach overcomes the extreme sparsity of individual cells while maintaining cluster-specific peak sets.

```r
# Ensure MACS2 is installed and in PATH
# pip install macs2

# Specify the path to MACS2 if it is not in the system PATH
pathToMacs2 <- findMacs2()

# Call peaks using MACS2 on pseudo-bulk replicates per cluster
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters"
)

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = 2.7e9,    # Effective genome size for hg38
  cutOff = 0.05,         # FDR cutoff
  extsize = 150,         # Extension size for paired-end fragments
  shift = -75
)

# Add peak matrix (peaks x cells)
proj <- addPeakMatrix(proj)

# Examine the peak set
getPeakSet(proj)
# GRanges object with ~150,000-300,000 peaks
```

### 3.8 Identifying Marker Peaks per Cluster

```r
# Find marker peaks for each cluster using a Wilcoxon test
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get a list of significant marker peaks per cluster
markerList <- getMarkers(
  markerPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  returnGR = TRUE
)

# Visualize as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.5",
  transpose = TRUE
)

draw(heatmapPeaks,
     heatmap_legend_side = "bot",
     annotation_legend_side = "bot")
```

### 3.9 Motif Enrichment with chromVAR

**chromVAR** identifies transcription factor motifs enriched in cell-type-specific accessible regions. This is one of the most powerful downstream analyses for scATAC-seq, linking chromatin accessibility patterns to the transcription factors that establish them.

```r
# Add motif annotations from the JASPAR or CIS-BP database
proj <- addMotifAnnotations(
  ArchRProj = proj,
  motifSet = "cisbp",    # or "JASPAR2020"
  name = "Motif"
)

# Run chromVAR to compute per-cell motif deviation scores
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

# Plot motif deviations on UMAP
# Key immune transcription factors:
motifs_of_interest <- c(
  "RUNX3",    # CD8+ T cell differentiation
  "TBX21",    # Th1 commitment (T-bet)
  "GATA3",    # Th2 commitment
  "FOXP3",    # Treg identity
  "SPI1",     # PU.1 -- myeloid lineage
  "PAX5",     # B cell identity
  "TCF7",     # Naive/memory T cells (TCF-1)
  "BATF"      # AP-1 family, effector T cells
)

p_motifs <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = paste0("z:", motifs_of_interest),
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p_motifs,
        name = "Motif-Deviations-UMAP",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
```

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure4-motif-deviations.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. chromVAR motif deviation scores projected onto UMAP. Each panel shows the activity of a key immune transcription factor motif. RUNX3 marks CD8+ T cells, SPI1 (PU.1) marks myeloid cells, PAX5 marks B cells, and FOXP3 marks regulatory T cells.
</div>

We can also perform enrichment analysis on marker peaks to identify motifs enriched in specific clusters:

```r
# Motif enrichment in marker peaks
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

# Plot enrichment as a heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs,
  n = 7,
  transpose = TRUE
)

draw(heatmapEM,
     heatmap_legend_side = "bot",
     annotation_legend_side = "bot")
```

### 3.10 Browser Track Visualization

ArchR can generate genome browser-style tracks showing accessibility across clusters:

```r
# Plot browser tracks for specific loci
p_tracks <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "Clusters",
  geneSymbol = c("CD3D", "MS4A1", "CD14"),
  upstream = 50000,
  downstream = 50000
)

plotPDF(p_tracks,
        name = "BrowserTracks",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 8, height = 6)
```

---

## 4. Signac Workflow (Comparison)

[Signac](https://stuartlab.org/signac/) extends the Seurat framework for chromatin data analysis. If you are already familiar with Seurat for scRNA-seq, Signac provides a natural extension with a consistent API. The core differences from ArchR:

| Feature | ArchR | Signac |
|---------|-------|--------|
| Data storage | HDF5-based Arrow files | In-memory Seurat object |
| Scalability | Excellent (100k+ cells) | Moderate (requires more RAM) |
| Peak calling | Built-in MACS2 wrapper | Manual or via `CallPeaks()` |
| Motif analysis | Built-in chromVAR | Via `RunChromVAR()` wrapper |
| Integration | Built-in label transfer | Via Seurat v5 bridge |
| Learning curve | Standalone API | Familiar if you know Seurat |

### 4.1 Creating a Signac Object

```r
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

# Read 10x Cell Ranger ATAC output
counts <- Read10X_h5("/path/to/filtered_peak_bc_matrix.h5")
metadata <- read.csv("/path/to/singlecell.csv",
                     header = TRUE, row.names = 1)

# Create a ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = "/path/to/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Create Seurat object with the chromatin assay
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations
```

### 4.2 QC and Processing

```r
# Compute QC metrics
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments /
                           pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments /
                        pbmc$peak_region_fragments

# Nucleosome signal and TSS enrichment
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# Filter cells
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

# Normalization and dimensionality reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

# Check that the first LSI component does not correlate with sequencing depth
DepthCor(pbmc)
# If component 1 correlates with depth, exclude it from downstream steps

# UMAP and clustering
pbmc <- RunUMAP(object = pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3, resolution = 0.5)

DimPlot(pbmc, label = TRUE) + NoLegend()
```

### 4.3 Gene Activity and Motif Analysis in Signac

```r
# Gene activity scores
gene.activities <- GeneActivity(pbmc)
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc,
            features = c("CD3D", "MS4A1", "CD14", "NKG7"),
            max.cutoff = "q95")

# chromVAR motif analysis in Signac
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(pbmc) <- "peaks"

# Get JASPAR motif position frequency matrices
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE",
              tax_group = "vertebrates",
              all_versions = FALSE)
)

# Add motif information to the Seurat object
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Run chromVAR
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(pbmc) <- "chromvar"
FeaturePlot(pbmc,
            features = c("MA0139.1"),  # CTCF motif
            min.cutoff = "q10",
            max.cutoff = "q90")
```

---

## 5. Integration with scRNA-seq (Label Transfer)

One of the most valuable analyses is integrating scATAC-seq with a matched or reference scRNA-seq dataset. This enables cell type annotation of ATAC clusters using RNA-based labels and provides a multi-omic view of gene regulation.

### 5.1 Label Transfer in ArchR

```r
# Load a Seurat scRNA-seq reference object
rna_ref <- readRDS("/path/to/pbmc_scrna_reference.rds")

# Perform unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna_ref,
  addToArrow = FALSE,
  groupBy = "celltype",        # Column in RNA metadata with labels
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

# Visualize transferred labels
p_labels <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAP"
)

plotPDF(p_labels,
        name = "LabelTransfer-UMAP",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 7, height = 5)
```

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure5-label-transfer.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. UMAP of scATAC-seq data colored by cell type labels transferred from a reference scRNA-seq dataset. Label transfer uses gene activity scores from ATAC data and gene expression from RNA data to find shared nearest neighbors across modalities.
</div>

### 5.2 Label Transfer in Signac

```r
# Using Seurat v5 for cross-modality integration
transfer.anchors <- FindTransferAnchors(
  reference = rna_ref,
  query = pbmc,
  reduction = "cca",
  query.assay = "RNA"    # Uses the gene activity assay
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna_ref$celltype,
  weight.reduction = pbmc[["lsi"]],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
DimPlot(pbmc,
        group.by = "predicted.id",
        label = TRUE,
        repel = TRUE) + NoLegend()
```

---

## 6. Identifying Cell-Type-Specific Peaks and Regulatory Elements

With cluster identities established (either by clustering or label transfer), we can identify regulatory elements unique to each immune cell type.

### 6.1 Differential Accessibility Analysis

```r
# In ArchR: find marker peaks per cell type
markerPeaks_celltype <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Extract peaks specific to CD8+ T cells
cd8_markers <- getMarkers(
  markerPeaks_celltype,
  cutOff = "FDR <= 0.01 & Log2FC >= 2",
  returnGR = TRUE
)$`CD8 T`

# Annotate peaks with nearest gene
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

peakAnno <- annotatePeak(
  cd8_markers,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  annoDb = "org.Hs.eg.db",
  level = "gene"
)

plotAnnoPie(peakAnno)
```

### 6.2 Linking Peaks to Genes (Co-accessibility)

ArchR can identify **peak-to-gene links** -- distal peaks whose accessibility correlates with the activity of a target gene, indicating potential enhancer-gene regulatory relationships.

```r
# Add peak-to-gene links
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix"  # or "GeneScoreMatrix"
)

# Retrieve links
p2gLinks <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

# Visualize peak-to-gene links as a heatmap
p_p2g <- plotPeak2GeneHeatmap(
  ArchRProj = proj,
  groupBy = "predictedGroup"
)

# Browser tracks with peak-to-gene links for a specific locus
p_links <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "predictedGroup",
  geneSymbol = "IFNG",
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj, corCutOff = 0.45, resolution = 1000,
                            returnLoops = TRUE)
)
```

{% include figure.liquid loading="eager" path="assets/img/blog/scatac-seq/figure6-peak2gene-browser.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 6. Genome browser view of the IFNG locus showing cell-type-specific chromatin accessibility tracks and peak-to-gene links (arcs). Distal regulatory elements with correlated accessibility-expression patterns are connected to the IFNG promoter, revealing the enhancer landscape that drives interferon-gamma production in effector T cells.
</div>

### 6.3 Motif Footprinting

Transcription factor footprinting uses the fine-grained pattern of Tn5 insertions around binding sites to infer factor occupancy. A bound transcription factor protects DNA from Tn5 insertion, creating a "footprint" -- a dip in signal surrounded by accessible flanking DNA.

```r
# Compute motif footprints
motifPositions <- getPositions(proj)

# Select specific motifs for footprinting
seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[c("RUNX3", "SPI1", "PAX5")],
  groupBy = "predictedGroup"
)

# Plot footprints
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Subtract",
  plotName = "Footprints-Subtract",
  addDOC = FALSE,
  smoothWindow = 5
)
```

---

## 7. Python Alternatives

While R dominates the scATAC-seq analysis ecosystem, Python tools are maturing rapidly and offer integration with the broader scverse ecosystem.

### 7.1 episcanpy

[episcanpy](https://episcanpy.readthedocs.io/) extends Scanpy for epigenomic data, providing a familiar API for those already using Scanpy for scRNA-seq.

```python
import episcanpy as epi
import scanpy as sc
import anndata as ad

# Load a count matrix (peaks x cells) into AnnData
adata = epi.ct.peak_mtx(
    matrix_file="/path/to/matrix.mtx",
    features_file="/path/to/peaks.bed",
    barcodes_file="/path/to/barcodes.tsv"
)

# Basic QC
epi.pp.filter_cells(adata, min_features=1000)
epi.pp.filter_features(adata, min_cells=10)

# Binarize the matrix (recommended for ATAC data)
epi.pp.binarize(adata)

# Feature selection and dimensionality reduction
epi.pp.select_var_feature(adata, nb_features=50000, show=True)
epi.pp.lazy(adata)  # TF-IDF normalization

# PCA, neighbors, clustering, UMAP
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.leiden(adata, resolution=0.8)
sc.tl.umap(adata)

# Visualize
sc.pl.umap(adata, color="leiden", frameon=False)
```

### 7.2 pycisTopic

[pycisTopic](https://pycistopic.readthedocs.io/) uses probabilistic topic modeling (LDA) to decompose the peak-by-cell matrix into interpretable "topics" that correspond to regulatory programs.

```python
import pycisTopic

# Create a cisTopic object from a fragments file
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments="/path/to/fragments.tsv.gz",
    path_to_regions="/path/to/consensus_peaks.bed",
    path_to_blacklist="/path/to/hg38-blacklist.v2.bed"
)

# Run topic modeling with multiple numbers of topics
from pycisTopic.lda_models import run_cgs_models
models = run_cgs_models(
    cistopic_obj,
    n_topics=[10, 20, 30, 40, 50],
    n_cpu=8,
    n_iter=500,
    random_state=42
)

# Select the best model
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(models, select_model=30)
cistopic_obj.add_LDA_model(model)

# Dimensionality reduction and clustering
from pycisTopic.clust_vis import run_umap, run_leiden
run_umap(cistopic_obj)
run_leiden(cistopic_obj, resolution=0.8)
```

### 7.3 SnapATAC2

[SnapATAC2](https://kzhang.org/SnapATAC2/) is a newer Python-native tool offering fast preprocessing and analysis with Rust-backed operations:

```python
import snapatac2 as snap

# Import data from fragments file
data = snap.pp.import_data(
    fragment_file="/path/to/fragments.tsv.gz",
    genome=snap.genome.hg38,
    min_num_fragments=1000,
    min_tsse=4
)

# Feature selection using variable bins
snap.pp.select_features(data)

# Spectral embedding (similar to LSI)
snap.tl.spectral(data)

# Clustering and UMAP
snap.tl.umap(data)
snap.pp.knn(data)
snap.tl.leiden(data)

# Visualization
snap.pl.umap(data, color="leiden", interactive=False)
```

---

## 8. Key Takeaways

**Choosing a tool.** ArchR is the most full-featured R package for scATAC-seq, with excellent scalability and a complete analysis pipeline from fragments to motif enrichment. Signac integrates naturally with the Seurat ecosystem, making it ideal when you are already analyzing matched scRNA-seq data. In Python, SnapATAC2 offers the best balance of speed and features, while pycisTopic provides unique topic-modeling-based decomposition.

**Critical QC metrics.** Always check TSS enrichment (> 4), fragment count (> 1000), and the nucleosomal banding pattern in fragment size distributions. Poor libraries show flat TSS profiles and lack nucleosomal periodicity.

**The sparsity problem.** scATAC-seq data is far sparser than scRNA-seq. Strategies to mitigate this include using large genomic bins (500-bp tiles) for initial dimensionality reduction, aggregating cells into pseudo-bulk profiles for peak calling, and imputing gene activity scores for visualization.

**Integration is essential.** Standalone scATAC-seq analysis provides chromatin accessibility patterns, but integration with scRNA-seq -- through label transfer or multi-omic assays like 10x Multiome -- connects these regulatory landscapes to transcriptional output, enabling a mechanistic understanding of gene regulation.

**Motif analysis reveals regulatory logic.** chromVAR deviation scores and motif enrichment analyses link accessible chromatin regions to the transcription factors that bind them. In immunology, this reveals how lineage-defining factors like T-bet, GATA3, FOXP3, PU.1, and PAX5 establish and maintain cell identity through distinct chromatin programs.

**Reproducibility note.** Throughout this tutorial, we set `seed = 42` for reproducible results. In practice, scATAC-seq analyses involve stochastic steps (LSI, UMAP, clustering) that should be evaluated for robustness by varying random seeds and resolution parameters.

---

## References

- Granja, J.M. et al. "ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis." *Nature Genetics* 53, 403--411 (2021). [DOI: 10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6)
- Stuart, T. et al. "Multimodal single-cell chromatin analysis with Signac." *Nature Methods* 21, 789--797 (2024). [DOI: 10.1038/s41592-023-02036-9](https://doi.org/10.1038/s41592-023-02036-9)
- Schep, A.N. et al. "chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data." *Nature Methods* 14, 975--978 (2017). [DOI: 10.1038/nmeth.4401](https://doi.org/10.1038/nmeth.4401)
- Zhang, K. et al. "SnapATAC2: a fast, scalable and versatile tool for analysis of single-cell omics data." *Nature Methods* 21, 217--227 (2024). [DOI: 10.1038/s41592-023-02139-9](https://doi.org/10.1038/s41592-023-02139-9)
- Gonzalez-Blas, C.B. et al. "cisTopic: cis-regulatory topic modeling on single-cell ATAC-seq data." *Nature Methods* 16, 397--400 (2019). [DOI: 10.1038/s41592-019-0367-1](https://doi.org/10.1038/s41592-019-0367-1)
- Buenrostro, J.D. et al. "Single-cell chromatin accessibility reveals principles of regulatory variation." *Nature* 523, 486--490 (2015). [DOI: 10.1038/nature14590](https://doi.org/10.1038/nature14590)
