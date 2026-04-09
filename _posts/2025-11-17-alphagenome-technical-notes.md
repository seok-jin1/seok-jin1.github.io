---
layout: post
title: "AlphaGenome: The AI That Reads DNA's Regulatory Code"
date: 2025-11-17
permalink: /blog/alphagenome-technical-notes/
published: true
categories: [paper-review]
tags:
  - AI
  - genomics
  - deep-learning
  - transformer
  - biology
---

Imagine you have a 3-billion-letter instruction manual for building a human, but 98% of it isn't written in the language of proteins—it's written in a mysterious regulatory code that controls when, where, and how much of each protein gets made. For decades, this "dark matter" of the genome has confounded scientists. How do distant DNA sequences turn genes on or off? How does a single mutation in a seemingly empty region cause disease?

In June 2025, Google DeepMind unveiled **AlphaGenome**, an AI that can read this hidden regulatory code across megabase-scale genomic regions at single-nucleotide resolution.

Formally, we are given a DNA sequence

$$
s = (s_1, s_2, \dots, s_L), \quad s_i \in \{A, C, G, T\},
$$

where $$s$$ is the DNA sequence, each $$s_i$$ is one of four nucleotides, and $$L$$ can be up to 1,000,000 base pairs (1 megabase).

and asked to predict thousands of regulatory properties

$$
\mathbf{Y} = \{y^{(1)}, y^{(2)}, \dots, y^{(K)}\},
$$

where each $$y^{(k)} \in \mathbb{R}^L$$ represents a regulatory track (gene expression, chromatin accessibility, etc.) at single-nucleotide resolution, and $$K$$ can be 5,930 tracks for human genome.

In doing so, AlphaGenome outperformed the best external models on 22 out of 24 genomic prediction tasks and recovered twice as many gene expression variants as previous state-of-the-art methods.

Explore the full [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.06.25.661532v1) and [supplementary material](https://www.biorxiv.org/content/biorxiv/early/2025/06/26/2025.06.25.661532/DC1/embed/media-1.pdf) for deeper technical details.

---

## Introduction

The authors outline seven novel contributions that distinguish the AlphaGenome model.

1. **Hybrid CNN-Transformer architecture.** A new design combines <span style="background-color: #fff3b0;">convolutional neural networks</span> for local pattern detection with <span style="background-color: #fff3b0;">transformers</span> for long-range interactions, enabling efficient processing of megabase-scale sequences.

   Mathematically, the architecture processes sequences through three stages:

   $$
   h^{(0)} = \text{Conv}_{\text{enc}}(\mathbf{X}), \quad h^{(L_{\text{trans}})} = \text{Transformer}(h^{(0)}), \quad \mathbf{Y} = \text{Decoder}(h^{(L_{\text{trans}})})
   $$

   where $$\mathbf{X} \in \{0,1\}^{L \times 4}$$ is the one-hot encoded DNA sequence, $$h^{(0)} \in \mathbb{R}^{L' \times d_{\text{model}}}$$ is the downsampled representation, and $$\mathbf{Y}$$ collects all predicted tracks.

2. **Megabase-scale context at single-base resolution.** By using sequence parallelism across TPU chips, AlphaGenome achieves 1 million base pair context at 1 base pair resolution—5× longer context and 128× finer resolution than previous models like Enformer.

   Previous models faced a fundamental trade-off:

   $$
   \text{context length} \times \text{resolution} \leq C,
   $$

   where $$C$$ is constrained by memory. AlphaGenome breaks this by distributing computation across devices.

3. **Unified multimodal prediction.** A single model jointly predicts across 11 biological modalities: RNA expression, <span style="background-color: #fff3b0;">CAGE</span> (transcription start sites), <span style="background-color: #fff3b0;">PRO-cap</span> (nascent RNA), splice sites and junctions, <span style="background-color: #fff3b0;">DNase</span> and <span style="background-color: #fff3b0;">ATAC-seq</span> (chromatin accessibility), histone modifications, transcription factor binding, and <span style="background-color: #fff3b0;">Hi-C</span> contact maps.

   The multi-task objective:

   $$
   \mathcal{L}_{\text{total}} = \sum_{k=1}^{K} w_k \mathcal{L}_k(\mathbf{y}^{(k)}, \hat{\mathbf{y}}^{(k)}),
   $$

   where $$w_k$$ are learned weights balancing different modalities, and $$\mathcal{L}_k$$ are modality-specific loss functions.

4. **Explicit splice junction modeling.** Unlike previous models that only predict splice donor and acceptor sites independently, AlphaGenome directly predicts which exons connect together and how strongly.

   For a gene with $$n$$ exons, the model predicts junction probabilities:

   $$
   p(\text{junction}_{i \to j} \mid s), \quad 1 \leq i < j \leq n,
   $$

   where $$p(\text{junction}_{i \to j})$$ represents the probability that exon $$i$$ splices to exon $$j$$.

5. **2D pairwise interaction module for Hi-C.** A dedicated pathway models pairwise interactions between genomic positions to predict 3D genome organization.

   The model constructs a pair representation:

   $$
   Z_{ij} = g(h_i, h_j), \quad Z \in \mathbb{R}^{L' \times L' \times d_z},
   $$

   where $$g$$ is a learned function combining features from positions $$i$$ and $$j$$, producing a contact probability matrix $$C \in \mathbb{R}^{L \times L}$$.

6. **Efficient knowledge distillation.** A teacher-student training framework enables deployment on consumer GPUs while maintaining accuracy.

   The student model learns from a teacher ensemble:

   $$
   \mathcal{L}_{\text{distill}} = \text{KL}\left( p_{\text{teacher}}(\mathbf{Y} \mid s) \,\|\, p_{\text{student}}(\mathbf{Y} \mid s) \right),
   $$

   where the teacher is an ensemble of 4 models from cross-validation, and KL is the Kullback-Leibler divergence.

7. **Comprehensive variant effect scoring.** AlphaGenome scores how single-nucleotide changes affect all regulatory layers in approximately 1 second per variant.

   In silico mutagenesis computes:

   $$
   \Delta \mathbf{Y} = f_\theta(s_{\text{alt}}) - f_\theta(s_{\text{ref}}),
   $$

   where $$s_{\text{ref}}$$ and $$s_{\text{alt}}$$ are reference and alternate alleles, and $$\Delta \mathbf{Y}$$ quantifies regulatory impact.

<img src="/assets/img/alphagenome/1.PNG" alt="AlphaGenome model overview" class="zoomable" style="width:90%;max-width:1000px;display:block;margin:20px auto;" />
*Figure 1. AlphaGenome model overview. Input: 1Mb DNA sequence processed through sequence parallelism across interconnected devices. Output: Predictions across 11 modalities—RNA-seq (667 human, 546 mouse tracks), CAGE (173, 188), DNase (305, 67), ATAC (167, 18), histone modifications (1116, 183), TF binding (1617, 127), splice sites (4, 4), and splice site usage (734, 180)—with Pearson correlation values showing performance on held-out test chromosomes.*

---

## How Is AlphaGenome Different?

### From Enformer to AlphaGenome

Classic sequence-to-function models like Enformer (Avsec et al., Nature Methods 2021) faced limitations in both context length and output resolution. They could either look at long sequences with coarse predictions, or make fine-grained predictions over short windows—but not both.

AlphaGenome breaks this trade-off through architectural innovations and distributed computing:

| Feature               | Enformer          | AlphaGenome                |
| :-------------------- | :---------------- | :------------------------- |
| **Context window**    | 200kb             | 1Mb (5×)                   |
| **Output resolution** | 128bp bins        | 1bp (128×)                 |
| **Architecture**      | CNN + Transformer | CNN + Transformer + U-Net  |
| **Modalities**        | Separate heads    | Unified (11 modalities)    |
| **Splice junctions**  | No                | Yes (explicit)             |
| **Hi-C contacts**     | Limited           | Dedicated 2D module        |
| **Training time**     | Days (single TPU) | Distributed across TPU pod |
| **Inference**         | ~2 seconds        | ~1 second                  |
| **Parameters**        | 220M              | 450M                       |

Mathematically, the key difference lies in resolution restoration:

- **Enformer:** Bins sequence into 128bp windows, transformer operates on $$L/128$$ tokens

  $$
  \text{Output} \in \mathbb{R}^{(L/128) \times K}
  $$

- **AlphaGenome:** Sequence parallelism over 131kb chunks, U-Net decoder restores full resolution

  $$
  \text{Output} \in \mathbb{R}^{L \times K}
  $$

### End-to-End Multimodal Learning

Previous approaches trained separate specialized models for each regulatory layer: one for gene expression, another for accessibility, a third for splicing. AlphaGenome instead uses a shared encoder with modality-specific heads, learning cross-regulatory relationships.

The architecture ensures that understanding one regulatory layer (e.g., histone marks) improves predictions for related layers (e.g., gene expression), because enhancers marked by H3K27ac typically drive expression.

---

## A Mathematical Glimpse Inside AlphaGenome

This section summarizes the core computations in slightly more formal terms, while staying high-level enough for a blog post.

### Input Representation

**DNA encoding:**

- Raw sequence:
  $$\text{DNA} = \text{"ATCGATCG..."}$$

- One-hot encoding:
  $$\mathbf{X} \in \{0, 1\}^{L \times 4}$$

  where:
  $$\mathbf{X}_{i,k} = \mathbb{1}[s_i = \text{nucleotide}_k], \quad k \in \{A, C, G, T\}$$

**Positional encoding:**

AlphaGenome uses learned relative positional embeddings. For positions $$i$$ and $$j$$:

$$
\text{relpos}_{ij} = \text{Embedding}(i - j),
$$

where $$\text{relpos}_{ij}$$ is added to attention scores during the transformer stage.

### Stage 1: Convolutional Encoder

The encoder consists of a DNA embedder followed by 6 downsampling blocks (Dowres blocks), progressively compressing the sequence while detecting increasingly complex regulatory motifs:

```
Input:       X ∈ R^(1M × 4)          # 1 megabase × 4 nucleotides
DNA embed:   → R^(500k × 64)         # 2bp resolution, detect TF binding motifs
Dowres1:     → R^(250k × 128)        # Composite regulatory elements
Dowres2:     → R^(125k × 256)        # Enhancer modules
Dowres3:     → R^(62.5k × 512)       # Regulatory domains
Dowres4:     → R^(31.2k × 768)       # Long-range interactions
Dowres5:     → R^(15.6k × 1024)      # Regulatory neighborhoods
Dowres6:     → R^(7.8k × 1536)       # 128bp resolution, d_model dimension
```

Formally:

$$
h^{(0)} = \text{Conv}_{\text{enc}}(\mathbf{X}), \quad h^{(0)} \in \mathbb{R}^{L' \times d_{\text{model}}},
$$

where $$L' = L / 128 \approx 7812$$ for 1Mb input, and $$d_{\text{model}} = 1536$$ is the hidden dimension.

Each convolutional layer uses:

$$
h^{(\ell+1)} = \text{ReLU}\left(\text{BatchNorm}\left(\text{Conv}(h^{(\ell)})\right)\right),
$$

with kernel sizes decreasing from 15bp to 3bp as resolution coarsens.

<img src="/assets/img/alphagenome/2.PNG" alt="AlphaGenome architecture diagram" class="zoomable" style="width:50%;max-width:500px;display:block;margin:20px auto;" />
*Figure 2. AlphaGenome detailed architecture. The model processes 1Mb input through: (1) DNA embedder producing 2bp resolution features, (2) Six downsampling blocks (Dowres blocks 1-6) progressively reducing resolution to 128bp, (3) Transformer tower maintaining 128bp resolution for long-range interactions, (4) Seven upsampling blocks (Upres blocks 1-7) restoring single-base resolution, (5) Output embedder producing modality-specific predictions. Tensor shapes (batch × sequence length × channels) shown at each stage. Skip connections (gray arrows) from encoder to decoder preserve fine-grained patterns.*

### Stage 2: Transformer Tower

Multi-head self-attention captures long-range regulatory interactions:

$$
\text{Attn}(h)_i = \sum_{j=1}^{L'} \alpha_{ij} V_j,
$$

where attention weights are:

$$
\alpha_{ij} = \frac{\exp(Q_i^\top K_j / \sqrt{d_k} + \text{bias}_{ij})}{\sum_{j'=1}^{L'} \exp(Q_i^\top K_{j'} / \sqrt{d_k} + \text{bias}_{ij'})},
$$

and $$\text{bias}_{ij} = \text{relpos}_{ij}$$ incorporates relative position information.

**Key innovation: Sequence parallelism**

To handle $$L' \approx 7812$$ tokens, the sequence is split into chunks distributed across 8 TPU chips. Each chip processes:

$$
h^{(\ell)}_{\text{chip}_m} = \text{LocalAttn}(h^{(\ell)}[(m-1) \cdot C : m \cdot C]),
$$

where $$C = L'/8 \approx 976$$ tokens per chip. Global information sharing occurs through:

$$
h^{(\ell+1)} = \text{GlobalCommunication}\left(\{h^{(\ell)}_{\text{chip}_1}, \ldots, h^{(\ell)}_{\text{chip}_8}\}\right).
$$

The transformer consists of approximately 12-24 layers (architecture details not fully disclosed in preprint):

$$
h^{(\ell+1)} = h^{(\ell)} + \text{FFN}\left(\text{LayerNorm}\left(h^{(\ell)} + \text{Attn}(h^{(\ell)})\right)\right),
$$

where $$\text{FFN}$$ is a feed-forward network with GELU activation.

### Stage 3: 2D Pairwise Module (for Hi-C Contacts)

For predicting 3D genome organization, AlphaGenome constructs a pairwise representation similar to AlphaFold's approach:

$$
Z_{ij} = \text{Linear}\left([h_i \,;\, h_j \,;\, h_i \odot h_j \,;\, |h_i - h_j|]\right),
$$

where $$[\cdot;\cdot]$$ denotes concatenation, $$\odot$$ is element-wise product, and $$Z \in \mathbb{R}^{L' \times L' \times d_z}$$.

**Triangle multiplicative updates** (inspired by AlphaFold) enable the model to reason about genomic triplets $$(i, j, k)$$:

$$
Z'_{ij} = Z_{ij} + \sum_{k=1}^{L'} \sigma(W_{\text{left}} Z_{ik}) \odot \sigma(W_{\text{right}} Z_{kj}),
$$

where $$\sigma$$ is a gating function, and $$W_{\text{left}}, W_{\text{right}}$$ are learned projections.

The final contact map is:

$$
C_{ij} = \text{sigmoid}\left(\text{Projection}(Z_{ij})\right), \quad C \in \mathbb{R}^{L \times L},
$$

representing the probability of physical contact between positions $$i$$ and $$j$$ in 3D space.

### Stage 4: U-Net Decoder with Skip Connections

To restore single-nucleotide resolution, AlphaGenome uses a U-Net decoder that combines:

1. **Upsampling** from compressed representation $$h^{(L_{\text{trans}})} \in \mathbb{R}^{L' \times d}$$
2. **Skip connections** from convolutional encoder stages

$$
\hat{\mathbf{y}}^{(k)} = \text{Conv}_{\text{dec}}^{(k)}\left(h^{(L_{\text{trans}})}, \text{skip}_7, \text{skip}_6, \ldots, \text{skip}_1\right),
$$

where $$\text{skip}_\ell$$ are activations from encoder layer $$\ell$$, and $$\text{Conv}_{\text{dec}}^{(k)}$$ is the modality-specific decoder head.

The decoder uses 7 upsampling blocks (Upres blocks) with transposed convolutions:

```
h^(L_trans):  R^(7.8k × 1536)      # Transformer output, 128bp resolution
Upres1:       R^(15.6k × 1024)     + skip from Dowres5
Upres2:       R^(31.2k × 768)      + skip from Dowres4
Upres3:       R^(62.5k × 512)      + skip from Dowres3
Upres4:       R^(125k × 256)       + skip from Dowres2
Upres5:       R^(250k × 128)       + skip from Dowres1
Upres6:       R^(500k × 64)        + skip from DNA embedder
Upres7:       R^(1M × 64)          # 1bp resolution restored
Output:       R^(1M × K)           # Modality-specific predictions
```

This enables the model to preserve fine-grained patterns (e.g., exact TF binding positions) while incorporating long-range context (e.g., distant enhancers).

### Splice Junction Prediction (Novel Component)

Traditional models predict splice donor and acceptor sites independently:

$$
p(\text{donor})_i, \quad p(\text{acceptor})_j.
$$

AlphaGenome additionally predicts:

1. **Which exons connect**: For every pair of potential donor site $$i$$ and acceptor site $$j$$:

   $$
   p(\text{junction}_{i \to j} \mid s),
   $$

2. **Junction strength**: The expected number of reads supporting junction $$(i, j)$$:

   $$
   \text{strength}_{ij} = \mathbb{E}[\text{read count}_{ij} \mid s].
   $$

The splice loss combines three terms:

$$
\mathcal{L}_{\text{splice}} = \mathcal{L}_{\text{donor}} + \mathcal{L}_{\text{acceptor}} + \lambda_{\text{junc}} \mathcal{L}_{\text{junction}},
$$

where:

$$
\mathcal{L}_{\text{junction}} = -\sum_{(i,j) \in \mathcal{J}_{\text{true}}} \log p_\theta(\text{junction}_{i \to j} \mid s),
$$

and $$\mathcal{J}_{\text{true}}$$ is the set of true junctions from RNA-seq data.

This enables ab initio prediction of alternative splicing patterns directly from genomic sequence.

<img src="/assets/img/alphagenome/3.PNG" alt="Splice junction prediction and DLG1 example" class="zoomable" style="width:90%;max-width:1000px;display:block;margin:20px auto;" />
*Figure 3. AlphaGenome's comprehensive splice predictions. (a) Comparison with specialized splicing models. AlphaGenome uniquely predicts all four aspects: RNA-seq coverage (32bp resolution), splice sites, splice site usage, and explicit splice junctions—whereas models like SpliceAI, Borzoi, and Pangolin predict only subsets. (b) **Example: DLG1 exon skipping.** In tibial artery tissue (GTEx), the variant chr3:197081044:TACTC>T (deletion) causes exon skipping in the DLG1 gene (Discs Large Homolog 1, a synaptic scaffolding protein). AlphaGenome correctly predicts: (1) reduced junction strength from 0.54 → 0.32 (40% decrease), (2) corresponding drop in splice site usage, and (3) decreased RNA-seq coverage over the skipped exon. GTEx data: 20 homozygous REF vs 1 heterozygous ALT individual confirms the prediction.*

### Loss Functions and Training

**Multi-task loss:**

The overall objective combines losses across all $$K$$ tracks:

$$
\mathcal{L}_{\text{total}} = \sum_{k=1}^{K} w_k \mathcal{L}_k + \lambda_{\text{distill}} \mathcal{L}_{\text{distill}},
$$

where $$w_k$$ are per-track weights (learned or set heuristically based on data abundance).

**Modality-specific losses:**

1. **Count-based tracks** (RNA-seq, CAGE, ChIP-seq):

   Negative Poisson log-likelihood:

   $$
   \mathcal{L}_{\text{Poisson}} = -\sum_{i=1}^{L} \left[ y_i \log \hat{y}_i - \hat{y}_i - \log \Gamma(y_i + 1) \right],
   $$

   or negative binomial for overdispersed data.

2. **Binary tracks** (DNase peaks, TF binding sites):

   Binary cross-entropy:

   $$
   \mathcal{L}_{\text{BCE}} = -\sum_{i=1}^{L} \left[ y_i \log \hat{y}_i + (1 - y_i) \log(1 - \hat{y}_i) \right].
   $$

3. **Contact maps** (Hi-C):

   Mean squared error on log-transformed contact frequencies:

   $$
   \mathcal{L}_{\text{Hi-C}} = \sum_{i,j} \left(\log(1 + C_{ij}) - \log(1 + \hat{C}_{ij})\right)^2.
   $$

**Knowledge distillation:**

The student model learns from an ensemble of 4 teacher models (from 4-fold cross-validation):

$$
p_{\text{teacher}}(\mathbf{Y} \mid s) = \frac{1}{4}\sum_{m=1}^{4} p_{\text{teacher}_m}(\mathbf{Y} \mid s),
$$

$$
\mathcal{L}_{\text{distill}} = \sum_{k=1}^{K} \text{KL}\left(p_{\text{teacher}}(y^{(k)} \mid s) \,\|\, p_{\text{student}}(y^{(k)} \mid s)\right).
$$

This enables a smaller student model (450M parameters) to match or exceed the ensemble's performance while being deployable on consumer hardware.

---

## Variant Effect Prediction

A key application of AlphaGenome is predicting how genetic variants affect gene regulation—critical for interpreting disease-associated mutations.

### In Silico Mutagenesis Workflow

For a variant at position $$v$$:

1. **Reference prediction:**
   $$\mathbf{Y}_{\text{ref}} = f_\theta(s_{\text{ref}})$$

2. **Create alternate sequence:**
   $$s_{\text{alt},i} = \begin{cases} s_{\text{ref},i} & \text{if } i \neq v \\ \text{alternate allele} & \text{if } i = v \end{cases}$$

3. **Alternate prediction:**
   $$\mathbf{Y}_{\text{alt}} = f_\theta(s_{\text{alt}})$$

4. **Effect score:**
   $$\Delta \mathbf{Y} = \mathbf{Y}_{\text{alt}} - \mathbf{Y}_{\text{ref}}$$

This produces a $$K$$-dimensional vector quantifying the variant's impact across all regulatory layers.

### Aggregation Strategies by Modality

Different regulatory questions require different aggregation of the $$\Delta \mathbf{Y}$$ signal:

**Gene expression (eQTL prediction):**

Sum effect over the gene's transcription start site (TSS) and gene body:

$$
\Delta E_g = \sum_{i \in [\text{TSS}_g - 10\text{kb}, \text{TSS}_g + 10\text{kb}]} \Delta y_i^{(\text{CAGE})},
$$

where positive $$\Delta E_g$$ indicates increased expression of gene $$g$$.

<img src="/assets/img/alphagenome/4.PNG" alt="RNA-seq variant scoring and APOL4 eQTL example" class="zoomable" style="width:90%;max-width:1000px;display:block;margin:20px auto;" />
*Figure 4. RNA-seq variant effect prediction with concrete eQTL example. (a) Variant scoring methodology: (1) Apply exon mask to focus on gene body, (2) Compute MEAN of masked predictions (REF: 1.8, ALT: 1.5 in this example), (3) Apply log-transform: log(x+1e-3), (4) Calculate effect: ALT - REF = -0.2. (b) **Example: APOL4 expression QTL in colon tissue.** The variant chr22:36201698:A>C creates a new transcription factor binding motif (shown in sequence logo on right). AlphaGenome predicts increased APOL4 (Apolipoprotein L4, involved in lipid metabolism) expression in sigmoid colon tissue. The prediction accurately captures: (1) Reference vs alternate RNA-seq coverage differences across 25kb region, (2) Tissue-specific effect (sigmoid colon), and (3) Direction of effect confirmed by GTEx data (20 homozygous REF vs 2 heterozygous ALT samples). The motif creation mechanism explains the quantitative expression change.*

**Splice site disruption:**

Maximum absolute change at annotated splice positions:

$$
\Delta S = \max_{i \in \{\text{donors} \cup \text{acceptors}\}} |\Delta y_i^{(\text{splice})}|.
$$

**Transcription factor binding:**

Change in peak height:

$$
\Delta B_{\text{TF}} = \max_i \hat{y}_{\text{alt},i}^{(\text{ChIP-TF})} - \max_i \hat{y}_{\text{ref},i}^{(\text{ChIP-TF})}.
$$

**Chromatin accessibility:**

Sum effect in a window around the variant:

$$
\Delta A = \sum_{i \in [v - 500, v + 500]} \Delta y_i^{(\text{ATAC})}.
$$

### Applications

**1. eQTL recovery:**

For each variant-gene pair, AlphaGenome predicts $$\Delta E_g$$. Variants are ranked by $$|\Delta E_g|$$, and the top predictions are compared against experimentally validated eQTLs from GTEx.

**Result:** AlphaGenome recovers **41% of GTEx eQTLs** at FDR < 5%, compared to 19% for the previous best model (Borzoi).

**2. GWAS interpretation:**

GWAS studies identify genomic loci associated with disease, but the causal variant within each locus is often unclear. AlphaGenome scores all variants in GWAS credible sets and identifies the most likely causal variant by regulatory impact.

**Result:** AlphaGenome resolves **49% of GWAS loci** to a single variant, compared to ~30% via traditional colocalization methods.

**3. Disease mutation analysis:**

Example: T-cell acute lymphoblastic leukemia (T-ALL) patients harbor recurrent mutations in the TAL1 enhancer region. AlphaGenome reveals that the mutation creates a _de novo_ MYB transcription factor binding motif, leading to aberrant TAL1 activation.

---

## Performance and Benchmarks

### Sequence Prediction Tasks

AlphaGenome was evaluated against the best external models (Enformer, Borzoi, Sei, Basenji2) on 24 held-out genomic regions.

**Table: Pearson correlation on human test chromosomes (chr8, chr9)**

| Modality               | Best Prior Model | Prior r | AlphaGenome | Improvement |
| :--------------------- | :--------------- | :------ | :---------- | :---------- |
| Gene expression (CAGE) | Borzoi           | 0.46    | 0.54        | +17.4%      |
| TSS activity (CAGE)    | Enformer         | 0.58    | 0.62        | +6.9%       |
| Nascent RNA (PRO-cap)  | Enformer         | 0.51    | 0.56        | +9.8%       |
| DNase accessibility    | Borzoi           | 0.72    | 0.76        | +5.6%       |
| ATAC-seq               | Borzoi           | 0.68    | 0.71        | +4.4%       |
| H3K4me3 (promoters)    | Borzoi           | 0.61    | 0.65        | +6.6%       |
| H3K27ac (enhancers)    | Enformer         | 0.64    | 0.68        | +6.3%       |
| H3K36me3 (gene bodies) | Enformer         | 0.59    | 0.63        | +6.8%       |
| CTCF binding           | Borzoi           | 0.53    | 0.57        | +7.5%       |
| Hi-C contacts          | Orca             | 0.52    | 0.55        | +5.8%       |
| Splice donor usage     | Pangolin         | 0.71    | 0.80        | +12.7%      |
| Splice acceptor usage  | Pangolin         | 0.69    | 0.78        | +13.0%      |

**Overall:** AlphaGenome wins **22 out of 24** tasks.

### Variant Effect Prediction Tasks

**Table: Variant interpretation benchmarks**

| Task              | Dataset      | Metric           | AlphaGenome | 2nd Best Model        |
| :---------------- | :----------- | :--------------- | :---------- | :-------------------- |
| eQTL recovery     | GTEx v8      | Recall @ 5% FDR  | 41%         | 19% (Borzoi)          |
| eQTL direction    | GTEx v8      | Accuracy         | 74%         | 68% (Enformer)        |
| ATAC QTL          | UK Biobank   | Accuracy         | 74%         | 68% (ChromBPNet)      |
| GWAS resolution   | GWAS Catalog | % single variant | 49%         | ~30% (colocalization) |
| Splice disruption | ClinVar      | auROC            | 0.88        | 0.82 (Pangolin)       |
| Promoter variants | MPRA         | Spearman ρ       | 0.61        | 0.54 (Enformer)       |

**Overall:** AlphaGenome wins **24 out of 26** variant tasks.

<img src="/assets/img/alphagenome/5.PNG" alt="Accessibility variant scoring and QTL benchmarks" class="zoomable" style="width:90%;max-width:1000px;display:block;margin:20px auto;" />
*Figure 5. Comprehensive variant effect prediction performance. (a) Accessibility (ATAC/DNase) variant scoring methodology using center mask aggregation: (1) Apply center mask around variant, (2) SUM the masked predictions (REF: 0.4, ALT: 1.3 in example), (3) Log-transform, (4) Compute ALT - REF difference = 0.7. (b) Causality tasks: Average precision for chromatin accessibility QTLs (caQTL), DNase QTLs (dsQTL), and binding QTLs (bQTL) across multiple populations (African, European, Yoruba) and datasets (SPI1). AlphaGenome (distilled, blue) consistently outperforms Borzoi ensemble (yellow) and ChromBPNet (gray). (c) Coefficient tasks: Pearson correlation for predicting effect sizes across same QTL types and cell types (including SMC smooth muscle cells and microglia). AlphaGenome achieves highest or competitive performance on all tasks, demonstrating robust variant interpretation across diverse regulatory modalities and populations.*

### Computational Efficiency

Despite being 2× larger than Enformer (450M vs 220M parameters), AlphaGenome achieves efficient training and inference through distributed computing:

- **Training hardware:** Distributed across TPU pod with sequence parallelism
- **Inference time:** ~1 second per 1Mb sequence
- **Deployment:** Distilled student model can run on a single GPU

This efficiency comes from:

1. Knowledge distillation reducing deployment model size
2. Sequence parallelism enabling distributed training across chips
3. Optimized convolution and attention kernels

---

## The Brain Behind AlphaGenome

AlphaGenome has three "brains" working in concert, much like how your own brain processes information at multiple scales:

### 1. The Pattern Spotter (Convolutional Encoder)

This module scans DNA for regulatory "words" and "phrases":

- **Short words (6-12bp):** TATA box, E-box, CAAT box—promoter elements
- **Phrases (20-50bp):** Transcription factor binding site clusters
- **Paragraphs (100-500bp):** Enhancer modules combining multiple TF sites
- **Chapters (1-10kb):** Super-enhancers controlling developmental genes

Like reading text, it starts with letters (individual nucleotides), builds words (motifs), then sentences (regulatory modules), and finally understands meaning (gene regulation).

### 2. The Relationship Builder (Transformer Tower)

Enhancers can lie 100,000 base pairs away from the genes they control. The transformer asks:

- "Which distant enhancer controls this gene's expression?"
- "If I change this CTCF site, which genes' expression will change?"
- "How do multiple enhancers cooperate to fine-tune expression?"

It's like understanding how clauses in a legal document relate to each other—even when separated by many pages.

### 3. The Precision Artist (U-Net Decoder)

The transformer works with coarse 128bp bins. The U-Net decoder refines predictions to single nucleotides:

- Pinpoints exact TF binding positions (often 8-10bp motifs)
- Identifies precise splice donor/acceptor sites (2bp GT-AG dinucleotides)
- Delineates nucleosome positioning (147bp boundaries)

Like an artist refining a charcoal sketch into a photorealistic painting, it adds fine detail while preserving the coarse composition.

---

## Real-World Impact

### Case Study: TAL1 Oncogenic Mutations in T-ALL

T-cell acute lymphoblastic leukemia (T-ALL) patients harbor recurrent non-coding mutations near the **TAL1** (T-cell acute lymphocytic leukemia 1) gene. Three distinct mutation clusters—located 6-21kb upstream, at the transcription start site (TSS), and 7-18kb downstream—all converge on upregulating TAL1 expression, driving oncogenesis.

AlphaGenome reveals the molecular mechanism of these oncogenic variants through multimodal regulatory analysis (Figure 6):

**The chr1:47239296:C>ACG insertion:**

This variant creates a _de novo_ **MYB transcription factor binding motif**. In silico mutagenesis comparing reference vs alternate sequences shows:

- **Reference sequence**: Scanning 40bp around the variant position → no predicted impact on TAL1 expression
- **Alternate sequence**: The 3bp insertion (C→ACG) creates MYB motif → predicted increase in:
  - TAL1 RNA-seq expression (0.5 log-fold increase in CD34+ common myeloid progenitors)
  - DNase accessibility (chromatin opening at the variant site)
  - Active histone marks: H3K27ac (enhancer activity), H3K4me1/3 (active chromatin), H3K36me3 (transcribed gene body)
  - Decreased repressive marks: H3K27me3 (polycomb silencing), H3K9me3 (heterochromatin)

The model additionally identifies a second **ETS-like motif** nearby that affects TAL1 expression only in the alternate sequence, not the reference—a previously unknown regulatory element.

**Validation and clinical relevance:**

When comparing all oncogenic TAL1 mutations against length-matched random insertions/deletions, AlphaGenome's variant scores correctly separate pathogenic mutations from benign background variants. Unsupervised clustering based on predicted multimodal effects cleanly groups oncogenic variants together, showing they share a common regulatory mechanism despite occurring at different genomic positions.

This mechanistic insight has therapeutic implications: understanding that TAL1 overexpression is driven by MYB motif creation suggests potential therapeutic strategies targeting MYB activity or the newly formed enhancer region.

<img src="/assets/img/alphagenome/6.PNG" alt="TAL1 T-ALL oncogenic mutations multimodal analysis" class="zoomable" style="width:90%;max-width:1000px;display:block;margin:20px auto;" />
*Figure 6. Interpreting T-ALL oncogenic mutations affecting TAL1 through multimodal variant analysis. (a) Three groups of mutations spanning 32.7kb (chr1:47209255-47242023) all target the TAL1 gene: 5' neo-enhancer mutations (6bp, 21bp upstream), TSS region mutations (C>T variants), and downstream cluster (7-18bp, 3bp, 2bp, 1bp insertions). (b) Predicted regulatory changes for chr1:47239296:C>ACG in CD34+ common myeloid progenitors (CMPs). Tracks show ALT-REF differences: RNA-seq predicts increased TAL1 expression 7.5kb away, DNase shows chromatin opening, active histone marks increase (H3K27ac, H3K4me1, H3K4me3, H3K36me3), while repressive marks decrease (H3K9me3, H3K27me3). (c) Variant effect scores for TAL1 RNA-seq expression in CD34+ CMPs comparing oncogenic mutations (orange) vs shuffled length-matched controls (gray). Oncogenic variants show significantly higher predicted impact. (d) Multimodal heatmap clustering variants by their predicted effects across all tracks. Oncogenic mutations cluster together, showing coherent regulatory mechanism. (e) In silico mutagenesis reveals mechanism: Reference sequence (top) shows no TAL1 expression sensitivity within 40bp. Alternate sequence (bottom) shows the ACG insertion creates MYB binding motif (yellow highlight matches UniProbe MYB consensus), which drives increased TAL1 expression, DNase, and H3K27ac. An ETS-like motif nearby (unlabeled) also contributes in the alternate sequence only.*

### 1. Precision Medicine

**Interpreting GWAS variants:**

Genome-wide association studies link genetic variants to diseases, but most disease-associated variants lie in non-coding regions. AlphaGenome answers: "How does this variant cause disease?"

Example: A variant associated with Type 2 diabetes lies 50kb from the nearest gene. AlphaGenome predicts it disrupts a FOXA2 binding site in a pancreatic enhancer, reducing insulin gene expression.

**Clinical variant interpretation:**

When sequencing a patient's genome, doctors find thousands of rare variants. AlphaGenome prioritizes which variants likely affect disease:

$$
\text{Pathogenicity score} = \sum_{k=1}^{K} \beta_k |\Delta y^{(k)}|,
$$

where $$\beta_k$$ weights each modality by disease relevance.

### 2. Drug Target Discovery

**Identifying causal genes at GWAS loci:**

GWAS identifies loci containing 10-50 genes. Which gene is the true target? AlphaGenome ranks genes by predicted expression change:

$$
\text{Target score}_g = \frac{\sum_{v \in \text{locus}} |\Delta E_{g,v}|}{\text{distance}(v, g)},
$$

prioritizing genes with large regulatory effects from nearby variants.

**Screening regulatory sequences:**

Pharmaceutical companies design antisense oligonucleotides (ASOs) or CRISPR therapies targeting specific genes. AlphaGenome predicts:

- Which splice sites to target for exon skipping
- Optimal guide RNA positions for minimal off-target effects
- Enhancer sequences to activate/silence for gene therapy

### 3. Synthetic Biology

**Designing synthetic promoters:**

Want a promoter that drives 10× higher expression in liver than brain? AlphaGenome enables:

1. Start with a weak constitutive promoter
2. Computationally screen 10,000 variants
3. Select variants predicted to increase liver-specific enhancer activity
4. Synthesize and test only top 10 designs

**Result:** Instead of testing 10,000 constructs in the lab, test only 10—saving months and millions of dollars.

**Optimizing codon usage:**

AlphaGenome predicts how synonymous codon changes affect:

- mRNA stability (via secondary structure)
- Splicing (exonic splicing enhancers)
- Translation efficiency (codon optimality)

This enables multi-objective optimization:

$$
\text{Maximize: } \text{Expression} \\
\text{Subject to: } \text{No splice site creation}, \text{GC content} \in [40\%, 60\%]
$$

### 4. Evolutionary Biology

**Comparing human and mouse regulatory code:**

AlphaGenome was trained on both human and mouse genomes (5,930 + 1,128 tracks). Cross-species predictions reveal:

- **Conserved regulatory grammar:** TATA boxes, splice sites work the same
- **Species-specific enhancers:** ~30% of enhancers are human-specific
- **Regulatory turnover:** Expression levels conserved, but enhancer sequences diverge

**Dating regulatory mutations:**

By analyzing ancient DNA (Neanderthal, Denisovan genomes), AlphaGenome identifies regulatory changes unique to modern humans:

- FOXP2 enhancer changes linked to language
- NOTCH2NL duplications driving brain expansion
- Metabolic enhancers adapting to dietary changes

---

## Limitations and Future Directions

### Current Limitations

**1. Very long-range regulation:**

While 1Mb is 5× longer than Enformer, some enhancers lie >1Mb from target genes (e.g., sonic hedgehog limb enhancer is 1.5Mb away). For such cases, AlphaGenome may miss the regulatory connection.

**2. Cell-type specificity:**

Training data is sparse for rare cell types (e.g., retinal photoreceptors, enteroendocrine cells). Predictions for these cells are less reliable due to limited examples.

**3. Dynamic processes:**

AlphaGenome is trained on static snapshots. It doesn't model:

- Developmental trajectories (embryonic → adult)
- Circadian rhythms (day/night gene expression oscillations)
- Response to stimuli (immune activation, stress response)

**4. 3D genome structure:**

Hi-C contact maps are predicted at coarse resolution. AlphaGenome doesn't explicitly model:

- Topologically associating domain (TAD) boundaries
- CTCF-mediated chromatin loops
- Phase-separated condensates (e.g., transcriptional hubs)

**5. Epigenetic memory:**

DNA methylation—a key epigenetic mark—is not predicted. This means AlphaGenome can't model:

- Imprinting (parent-of-origin effects)
- X-chromosome inactivation
- Cellular reprogramming

**6. Complex genetics:**

The in silico mutagenesis framework assumes:

- Single-variant effects (no epistasis)
- Additive contributions
- No environmental interactions

In reality, many diseases involve complex multi-variant haplotypes and gene-environment interactions.

### What's Next?

**AlphaGenome-Multimer:**

Just as AlphaFold-Multimer predicts protein complexes, an AlphaGenome-Multimer could model:

- Cooperative TF binding (e.g., AP-1 dimers)
- Enhanceosome assembly (multiple TFs on one enhancer)
- Chromatin remodeler positioning

**Dynamic AlphaGenome:**

Incorporate time-series data:

- Single-cell RNA-seq trajectories
- Circadian transcriptome datasets
- Developmental atlases (embryo → adult)

This would enable predictions like: "At 6 hours post-fertilization, this enhancer activates; at 24 hours, it switches off."

**Cross-species transfer:**

Current model handles human and mouse. Extending to:

- Model organisms (zebrafish, drosophila, C. elegans)
- Agricultural species (crops, livestock)
- Non-model organisms with limited training data

**Integration with AlphaFold:**

Joint prediction of:

- Gene regulatory networks
- Protein-DNA binding structures
- Regulatory protein complexes

This would answer: "How does SNP rs123 alter CTCF-DNA binding, change chromatin looping, and affect gene X expression?"

**Therapeutic design:**

Use AlphaGenome for:

- CRISPR guide RNA design (minimize off-targets)
- Antisense oligonucleotide optimization (maximize efficacy)
- Gene therapy enhancer engineering (tissue-specific expression)

---

## Key Takeaways

- AlphaGenome reads the genome's regulatory code—the 98% "dark matter" that controls when, where, and how much of each gene is expressed
- It achieves megabase-scale context (1 million bp) at single-nucleotide resolution, breaking previous computational barriers
- Unlike specialized tools, it jointly predicts 11 regulatory layers: expression, splicing, accessibility, 3D structure, and more
- It recovers twice as many gene expression variants as previous methods and resolves nearly half of GWAS loci to causal variants
- Most importantly, it democratizes genome interpretation—researchers can now score regulatory effects of any variant in ~1 second, accelerating precision medicine and synthetic biology

The genome is no longer a static blueprint. With AlphaGenome, it becomes a programmable regulatory computer—and we're just beginning to learn its instruction set.

---

_Reference: Avsec et al., "AlphaGenome: Nucleotide-resolution prediction of regulatory function at megabase scale,"_ bioRxiv (2025). DOI: [10.1101/2025.06.25.661532](https://doi.org/10.1101/2025.06.25.661532)

**Official resources:**

- [AlphaGenome homepage](https://deepmind.google.com/science/alphagenome/)
- [GitHub repository](https://github.com/google-deepmind/alphagenome)
- [DeepMind blog](https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/)

<style>
.post-content img.zoomable {
  cursor: zoom-in;
}

.image-zoom-overlay {
  position: fixed;
  inset: 0;
  background: rgba(0, 0, 0, 0.7);
  display: none;
  align-items: center;
  justify-content: center;
  padding: 24px;
  z-index: 9999;
}

.image-zoom-overlay.visible {
  display: flex;
}

.image-zoom-overlay img {
  max-width: 90vw;
  max-height: 90vh;
  border-radius: 8px;
  box-shadow: 0 12px 30px rgba(0, 0, 0, 0.35);
}
</style>

<script>
document.addEventListener("DOMContentLoaded", function () {
  const images = document.querySelectorAll(".post-content img.zoomable");
  if (!images.length) return;

  const overlay = document.createElement("div");
  overlay.className = "image-zoom-overlay";
  overlay.innerHTML = '<img alt="Expanded diagram" />';
  const overlayImg = overlay.querySelector("img");

  overlay.addEventListener("click", function () {
    overlay.classList.remove("visible");
  });

  document.body.appendChild(overlay);

  images.forEach(function (img) {
    img.addEventListener("click", function () {
      overlayImg.src = img.src;
      overlayImg.alt = img.alt || "Expanded image";
      overlay.classList.add("visible");
    });
  });
});
</script>
