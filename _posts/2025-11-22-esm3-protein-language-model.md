---
layout: post
title: "ESM-3: Simulating Evolution with Protein Language Models"
date: 2025-11-22
permalink: /blog/esm3-protein-language-model/
published: true
categories: [paper-review]
tags:
  - AI
  - biology
  - deep-learning
  - protein-design
  - language-model
---

Imagine you could speak the language of proteins -- not just reading their amino acid sequences, but simultaneously understanding their three-dimensional shapes and biological functions. Now imagine you could use that language to write entirely new proteins that nature has never seen, as if fast-forwarding evolution by hundreds of millions of years. This is what ESM-3 attempts to do: a single AI model that reads and generates proteins across sequence, structure, and function all at once.

Formally, a protein is described by three complementary representations:

$$
\mathbf{s} = (s_1, s_2, \dots, s_L), \quad s_i \in \{1, \dots, 20\},
$$

where $\mathbf{s}$ is the amino acid sequence of length $L$, with each $s_i$ indexing one of 20 standard amino acids,

$$
\mathbf{C} = \{(\mathbf{x}_{i,1}, \mathbf{x}_{i,2}, \mathbf{x}_{i,3}) \in \mathbb{R}^{3 \times 3} \mid i = 1, \dots, L\},
$$

where $\mathbf{C}$ collects the 3D coordinates of the three backbone atoms (N, C$_\alpha$, C) for each residue, and

$$
\mathbf{f} = \{f_k\}_{k=1}^{K},
$$

where $\mathbf{f}$ encodes functional annotations such as Gene Ontology terms, enzyme commission numbers, and InterPro domain labels. ESM-3 is a single model that reasons over all three modalities jointly.

Explore the full [bioRxiv preprint for ESM-3](https://www.biorxiv.org/content/10.1101/2024.07.01.600583v1) and the [Science paper for ESMFold](https://www.science.org/doi/10.1126/science.ade2574) for deeper technical details.

---

## Introduction

The story of protein language models begins with ESM-2 and ESMFold, which demonstrated that <span style="background-color: #fff3b0;">a single protein language model can learn enough about protein biology from sequences alone to predict three-dimensional structure</span> -- without the multiple sequence alignments (MSAs) that AlphaFold2 requires.

**ESMFold** (Lin et al., Science 2023) showed that a large language model (ESM-2, with up to 15 billion parameters) trained on millions of protein sequences develops an internal representation so rich that a relatively lightweight "folding trunk" can convert it into accurate 3D structures. The key insight was that evolutionary information, traditionally extracted from MSAs of related sequences, is implicitly encoded in the language model's learned representations. ESMFold runs an order of magnitude faster than AlphaFold2 because it skips the expensive MSA search step, enabling structure prediction for entire metagenomic databases.

**ESM-3** (Hayes et al., bioRxiv 2024), developed by EvolutionaryScale, extends this vision dramatically. Rather than just predicting structure from sequence, ESM-3 is a <span style="background-color: #fff3b0;">generative multimodal masked language model</span> that operates simultaneously over sequence, structure, and function. It can generate new proteins by iteratively unmasking tokens across all three modalities, enabling programmable protein design with control over shape, function, or both. The largest ESM-3 model has 98 billion parameters and was trained on a combined dataset of 2.78 billion protein sequences, 236 million structures, and 539 million function annotations -- making it one of the largest and most data-rich protein models ever built.

The conceptual progression from ESM-2 to ESMFold to ESM-3 reflects a broader trend in AI: moving from discriminative models (that classify or predict) to generative models (that create). Just as GPT moved from understanding text to generating it, the ESM family has moved from understanding proteins to designing them. This shift has profound implications for biotechnology, drug discovery, and synthetic biology, as it opens the door to designing proteins with properties that go beyond what natural evolution has explored.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure1-overview.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 1: ESM-3 architecture overview.</strong> ESM-3 processes proteins as sequences of discrete tokens across three tracks: amino acid sequence, structure, and function. All three tracks are embedded, fused, and processed by a single transformer backbone. The model can accept any combination of inputs (e.g., partial sequence and desired function) and generate the remaining modalities through iterative masked decoding.
</div>

---

## Novel Contributions

The ESM-3 paper introduces several key innovations that distinguish it from prior protein language models and structure prediction methods.

1. **Multimodal protein tokenization.** ESM-3 represents proteins as discrete tokens across three modalities -- sequence, structure, and function -- enabling a unified language modeling framework. Structure is tokenized using a <span style="background-color: #fff3b0;">Vector Quantized Variational Autoencoder (VQ-VAE)</span> that encodes local atomic neighborhoods into a learned codebook. Function is tokenized by encoding keyword sets (GO terms, EC numbers, InterPro labels) into binary vectors and quantizing them.

   The tokenization maps each modality into discrete tokens:

   $$
   t_i^{\text{seq}} \in \{1, \dots, V_s\}, \quad t_i^{\text{struct}} \in \{1, \dots, V_c\}, \quad t_i^{\text{func}} \in \{1, \dots, V_f\},
   $$

   where $V_s = 20$ (amino acids), $V_c = 4{,}096$ (structure codebook), and $V_f = 8$ per-residue function token slots, each drawn from a learned function codebook.

2. **All-to-all masked generation.** Unlike models that predict structure from sequence (ESMFold, AlphaFold2) or generate sequence from structure (inverse folding), ESM-3 can generate any modality conditioned on any combination of the others. This is achieved through a <span style="background-color: #fff3b0;">generalized masked language modeling objective</span> applied across all tracks simultaneously.

3. **Geometric attention for structural reasoning.** To handle the geometric nature of protein structure, ESM-3 incorporates a geometric attention mechanism that operates on 3D coordinates, allowing the transformer to reason about spatial relationships between residues while maintaining the discrete token framework.

4. **Unified architecture at scale.** A single transformer backbone processes all three modalities together. The authors train models at three scales -- 1.4 billion, 7 billion, and 98 billion parameters -- showing consistent improvements with scale. The 98B model represents one of the largest protein models trained to date.

5. **Function conditioning for programmable design.** By including function tokens, ESM-3 enables <span style="background-color: #fff3b0;">function-directed protein generation</span>: a user can specify desired functional properties (e.g., "fluorescent protein" or "kinase activity") and have the model generate sequences and structures that satisfy those constraints.

6. **Iterative decoding for generation.** Rather than generating all tokens in one pass, ESM-3 uses an iterative masking-unmasking procedure inspired by masked diffusion models. At each step, the model predicts masked positions, samples from the predicted distributions, and re-masks a fraction of tokens for the next step. This allows the model to refine its outputs over multiple rounds, resolving long-range dependencies that would be difficult to capture in a single autoregressive pass.

7. **Massive training data.** ESM-3 is trained on a dataset of 2.78 billion natural protein sequences, 236 million predicted structures (from ESMFold and other sources), and 539 million proteins with function annotations, making it by far the largest multimodal protein training set assembled. The diversity of training data sources -- including UniRef, predicted structures from ESMFold, and functional annotations from InterPro and Gene Ontology -- ensures broad coverage of protein space.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure2-tokenization.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 2: Multimodal tokenization in ESM-3.</strong> The three tracks of ESM-3: sequence tokens are the standard 20 amino acids; structure tokens are learned via a VQ-VAE that encodes local atomic coordinate neighborhoods into discrete codes from a codebook of 4,096 entries; function tokens encode Gene Ontology terms, enzyme commission numbers, and InterPro domain annotations as quantized keyword vectors. Each residue position carries tokens from all three tracks.
</div>

---

## How Is ESM-3 Different?

The following table compares ESM-3 with AlphaFold2, ESMFold, and other representative protein language models across key dimensions.

| Feature                      | AlphaFold2                     | ESMFold                           | ProteinMPNN                       | ESM-3                               |
| ---------------------------- | ------------------------------ | --------------------------------- | --------------------------------- | ----------------------------------- |
| **Primary task**             | Structure prediction           | Structure prediction              | Inverse folding (seq from struct) | Multimodal generation               |
| **Input**                    | Sequence + MSA + templates     | Sequence only                     | Backbone structure                | Any subset of seq/struct/func       |
| **Output**                   | 3D coordinates                 | 3D coordinates                    | Amino acid sequence               | Seq + struct + func tokens          |
| **MSA required**             | Yes                            | No                                | No                                | No                                  |
| **Generative**               | No                             | No                                | Yes (sequence only)               | Yes (all modalities)                |
| **Function awareness**       | No                             | No                                | No                                | Yes                                 |
| **Model size**               | ~93M params                    | Up to 15B (ESM-2) + folding trunk | ~1.7M params                      | 1.4B / 7B / 98B                     |
| **Training data**            | ~170K structures               | ~65M sequences (ESM-2)            | ~19K structures                   | 2.78B seq + 236M struct + 539M func |
| **Speed**                    | Minutes per protein (with MSA) | Seconds per protein               | Milliseconds                      | Seconds (per decoding iteration)    |
| **Structure representation** | Continuous (frames + torsions) | Continuous (frames + torsions)    | Continuous (backbone coords)      | Discrete (VQ-VAE tokens)            |

The most fundamental distinction is that ESM-3 treats protein design as a <span style="background-color: #fff3b0;">language generation problem across multiple modalities</span>, whereas AlphaFold2 and ESMFold are discriminative models focused on structure prediction, and ProteinMPNN is a conditional generative model limited to a single modality.

Another important distinction is the treatment of structure. AlphaFold2 and ESMFold represent structure continuously using backbone frames and torsion angles, while ESM-3 discretizes structure into tokens from a learned codebook. This discretization sacrifices some structural precision but enables a unified transformer architecture that processes all modalities in the same way. The trade-off is between the expressiveness of continuous representations and the simplicity and scalability of discrete token-based processing.

The training data requirements also differ substantially. AlphaFold2 is trained primarily on experimentally determined structures (roughly 170,000 at the time of its publication), while ESM-3 leverages billions of sequences, hundreds of millions of predicted structures, and hundreds of millions of function annotations. This massive data advantage, combined with the multimodal training objective, allows ESM-3 to learn richer representations of protein biology.

---

## A Mathematical Glimpse Inside

This section walks through the core mathematical machinery of ESM-3 and briefly revisits ESMFold as a precursor.

### ESMFold's Folding Trunk (Precursor)

ESMFold demonstrated that a protein language model's representations can be directly converted into 3D structures. The architecture has two stages:

1. **ESM-2 language model:** A transformer trained with masked language modeling on protein sequences. The ESM-2 model family ranges from 8 million to 15 billion parameters, with the largest model using 48 transformer layers with embedding dimension $d = 5{,}120$. Given a sequence $\mathbf{s}$, it produces per-residue representations:

   $$
   \mathbf{h} = \text{ESM-2}(\mathbf{s}), \quad \mathbf{h} \in \mathbb{R}^{L \times d},
   $$

   where $\mathbf{h}$ is the sequence of hidden representations and $d$ is the embedding dimension. The masked language modeling objective during pre-training is:

   $$
   \mathcal{L}_{\text{MLM}} = -\sum_{i \in \mathcal{M}} \log p_\theta(s_i \mid \mathbf{s}_{\setminus \mathcal{M}}),
   $$

   where $\mathcal{M}$ is the set of masked positions, $s_i$ is the true amino acid at position $i$, and $\mathbf{s}_{\setminus \mathcal{M}}$ denotes the unmasked sequence context. Through this objective, the model learns rich per-residue and pairwise representations that capture evolutionary constraints.

2. **Folding trunk with IPA:** The folding trunk takes ESM-2 representations and predicts backbone frames using a structure module with <span style="background-color: #fff3b0;">Invariant Point Attention (IPA)</span>, similar to AlphaFold2's structure module. The trunk first converts single-sequence representations into pair representations via an outer product operation:

   $$
   Z_{ij} = \text{Linear}(\mathbf{h}_i) \otimes \text{Linear}(\mathbf{h}_j), \quad Z \in \mathbb{R}^{L \times L \times d_z},
   $$

   where $Z_{ij}$ is the pair representation between residues $i$ and $j$, constructed from outer products of projected single representations. In AlphaFold2, this pair representation is derived from MSAs through the Evoformer; in ESMFold, it comes entirely from the language model's single-sequence embeddings.

   The structure module then iteratively refines backbone frames:

   $$
   F_i^{(t+1)} = \text{IPA}(F_i^{(t)}, \mathbf{h}_i, Z), \quad F_i = (R_i, \mathbf{t}_i) \in \text{SE}(3),
   $$

   where $F_i^{(t)}$ is the backbone frame of residue $i$ at refinement iteration $t$, composed of a rotation $R_i \in \text{SO}(3)$ and translation $\mathbf{t}_i \in \mathbb{R}^3$. The IPA mechanism ensures that the attention is equivariant under rigid body transformations:

   $$
   \text{IPA}(QF_i + \mathbf{u}, \mathbf{h}_i, Z) = Q \cdot \text{IPA}(F_i, \mathbf{h}_i, Z) + \mathbf{u},
   $$

   for any global rotation $Q \in \text{SO}(3)$ and translation $\mathbf{u} \in \mathbb{R}^3$, meaning the predicted structure is independent of the coordinate frame.

   The critical insight of ESMFold is that the language model representations $\mathbf{h}$ implicitly encode enough co-evolutionary information to replace MSAs, enabling single-sequence structure prediction. According to Lin et al., the attention maps of large ESM-2 models show patterns that closely mirror residue-residue contact maps, suggesting that the model has internalized the co-evolutionary signal that MSAs make explicit.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure7-esmfold.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 3: ESMFold architecture.</strong> ESMFold uses ESM-2 (a protein language model with up to 15 billion parameters) as its backbone, feeding learned representations into a folding trunk that includes Invariant Point Attention (IPA) to predict 3D backbone frames. Unlike AlphaFold2, ESMFold requires no MSA computation, enabling structure prediction in seconds rather than minutes. The folding trunk converts single-sequence embeddings into pair representations and iteratively refines residue frames.
</div>

### Tokenization of Sequence, Structure, and Function

ESM-3 discretizes all three protein modalities into tokens that a standard transformer can process.

**Sequence tokenization** is straightforward: each amino acid maps to one of 20 tokens (plus special tokens for mask, padding, etc.).

**Structure tokenization** uses a VQ-VAE trained on local atomic coordinate frames. For each residue $i$, the encoder extracts the local geometric neighborhood -- the positions of backbone and side-chain atoms relative to the residue's local frame -- and maps it to the nearest entry in a learned codebook:

$$
t_i^{\text{struct}} = \arg\min_{k \in \{1, \dots, V_c\}} \| \mathbf{z}_i^{\text{enc}} - \mathbf{e}_k \|_2,
$$

where $\mathbf{z}_i^{\text{enc}} \in \mathbb{R}^{d_c}$ is the encoder output for residue $i$, $\mathbf{e}_k$ is the $k$-th codebook embedding, and $V_c = 4{,}096$ is the codebook size. The VQ-VAE is trained to reconstruct local coordinates from the discrete codes with high fidelity.

The VQ-VAE training objective combines reconstruction loss, commitment loss, and codebook loss:

$$
\mathcal{L}_{\text{VQ-VAE}} = \underbrace{\| \mathbf{C}_i - \hat{\mathbf{C}}_i \|^2}_{\text{reconstruction}} + \beta \underbrace{\| \mathbf{z}_i^{\text{enc}} - \text{sg}[\mathbf{e}_{k^*}] \|^2}_{\text{commitment}} + \underbrace{\| \text{sg}[\mathbf{z}_i^{\text{enc}}] - \mathbf{e}_{k^*} \|^2}_{\text{codebook}},
$$

where $\hat{\mathbf{C}}_i$ are the reconstructed coordinates, $k^* = t_i^{\text{struct}}$ is the selected codebook index, $\text{sg}[\cdot]$ denotes the stop-gradient operator, and $\beta$ is a hyperparameter balancing the commitment loss. The codebook is updated via exponential moving averages during training.

**Function tokenization** encodes per-residue functional annotations. For each residue, the set of applicable functional keywords (from InterPro, Gene Ontology, and EC annotations) is represented as a binary vector and then quantized into a small number of discrete tokens using locality-sensitive hashing or a learned quantizer.

The overall input to the transformer at each position $i$ is the sum of embeddings from all available tracks:

$$
\mathbf{x}_i = \text{Embed}_{\text{seq}}(t_i^{\text{seq}}) + \text{Embed}_{\text{struct}}(t_i^{\text{struct}}) + \text{Embed}_{\text{func}}(t_i^{\text{func}}) + \text{Embed}_{\text{pos}}(i),
$$

where each $\text{Embed}$ is a learned embedding table for the corresponding track and $\text{Embed}_{\text{pos}}(i)$ adds positional information. When a track is masked, its embedding is replaced by a learned mask embedding $\mathbf{m}^{(m)}$. This additive fusion allows the transformer to process any combination of available modalities through a single forward pass.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure3-generation.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 4: Iterative decoding in ESM-3.</strong> Generation begins with all positions masked across all tracks. At each step, the model predicts probability distributions over tokens at masked positions, samples the most confident predictions, and unmasks them. The remaining positions are re-masked and the process repeats. Over many iterations, the protein is progressively constructed across sequence, structure, and function simultaneously. The decoding schedule controls how many tokens are unmasked per step, analogous to the noise schedule in diffusion models.
</div>

### Masked Language Model Training Across Modalities

ESM-3 is trained with a <span style="background-color: #fff3b0;">masked language modeling (MLM)</span> objective that operates across all three tracks simultaneously. During training, a random fraction of tokens from each track is independently masked, and the model learns to predict the masked tokens given the unmasked context.

For a protein with token sequences $\mathbf{t}^{\text{seq}}$, $\mathbf{t}^{\text{struct}}$, and $\mathbf{t}^{\text{func}}$, let $\mathcal{M}$ denote the set of masked positions across all tracks. The training objective is:

$$
\mathcal{L} = -\sum_{(i, m) \in \mathcal{M}} \log p_\theta(t_i^{(m)} \mid \mathbf{t}_{\setminus \mathcal{M}}),
$$

where $(i, m)$ indexes residue $i$ and modality $m \in \{\text{seq}, \text{struct}, \text{func}\}$, $t_i^{(m)}$ is the true token at that position, and $\mathbf{t}_{\setminus \mathcal{M}}$ denotes all unmasked tokens across all tracks.

The masking rate is sampled uniformly for each training example:

$$
r \sim \text{Uniform}(0, 1), \quad \text{each token masked independently with probability } r.
$$

This uniform masking schedule is critical: at $r \approx 0$ the model learns to fill in a few missing residues (akin to protein engineering), while at $r \approx 1$ it must generate nearly complete proteins from minimal context (akin to de novo design). Training across the full range of masking rates ensures the model can handle both scenarios.

A key design choice is that <span style="background-color: #fff3b0;">tracks are masked independently</span>. This means the model may see the full sequence but a masked structure, or a partial function annotation with no sequence -- forcing it to learn genuine cross-modal relationships rather than simply copying from one track to another.

To understand why independent masking is essential, consider the alternative: if all tracks were always masked or unmasked together at each position, the model could learn to simply copy sequence tokens to predict structure tokens (or vice versa) without learning the underlying physical relationship between them. Independent masking creates training scenarios where the model must, for example:

- Predict sequence from structure alone (inverse folding)
- Predict structure from sequence alone (forward folding)
- Predict function from structure (function annotation)
- Complete partial sequences given functional constraints (constrained design)
- Generate all three tracks from scratch (de novo design)

This diversity of training conditions is what enables ESM-3's versatility at inference time.

The output of the transformer produces per-position logits for each track through separate prediction heads:

$$
p(t_i^{(m)} \mid \mathbf{t}_{\setminus \mathcal{M}}) = \text{softmax}\left( W^{(m)} \mathbf{h}_i^{(\text{final})} + \mathbf{b}^{(m)} \right),
$$

where $W^{(m)}$ and $\mathbf{b}^{(m)}$ are the weight matrix and bias of the prediction head for modality $m$, and $\mathbf{h}_i^{(\text{final})}$ is the final-layer hidden state at position $i$. For structure tokens, an additional structure decoder (the VQ-VAE decoder) can convert the predicted discrete tokens back into 3D coordinates.

### Geometric Attention for Structure Tokens

While the discrete structure tokens capture local geometry, the model also needs to reason about global 3D relationships. ESM-3 incorporates a <span style="background-color: #fff3b0;">geometric attention</span> mechanism that augments the standard transformer attention with spatial information.

Given predicted or input 3D coordinates, the geometric attention computes pairwise distance and orientation features:

$$
d_{ij} = \| \mathbf{x}_i - \mathbf{x}_j \|_2, \quad \mathbf{o}_{ij} = R_i^T (\mathbf{x}_j - \mathbf{x}_i),
$$

where $d_{ij}$ is the Euclidean distance between residues $i$ and $j$, and $\mathbf{o}_{ij}$ is the relative position of residue $j$ in the local frame of residue $i$, with $R_i$ being the local rotation matrix. These geometric features are projected and added as biases to the attention logits:

$$
\alpha_{ij} = \text{softmax}_j \left( \frac{Q_i^T K_j}{\sqrt{d_k}} + b(d_{ij}, \mathbf{o}_{ij}) \right),
$$

where $Q_i$ and $K_j$ are the query and key vectors from the standard transformer attention, and $b(d_{ij}, \mathbf{o}_{ij})$ is a learned function of the geometric features. This allows the model to attend preferentially to spatially proximal residues even when they are distant in sequence.

This geometric bias is particularly important for proteins because residues that are far apart in the linear sequence are often close together in 3D space. For example, in a beta-sheet, two residues separated by 50 or more positions in the sequence may be hydrogen-bonded neighbors in the folded structure. Without geometric attention, the transformer would need to learn these long-range spatial relationships purely from positional encodings and data, which is much less efficient than providing explicit geometric information.

The geometric attention also differs from the Invariant Point Attention (IPA) used in AlphaFold2 and ESMFold. While IPA operates on continuous backbone frames and is SE(3)-equivariant by construction, ESM-3's geometric attention operates as a bias on the standard transformer attention mechanism. This design choice keeps the core architecture as a standard transformer (enabling efficient scaling) while still incorporating geometric reasoning where coordinates are available.

### Generation via Iterative Decoding

At inference time, ESM-3 generates proteins through an iterative decoding process. Given some conditioning information (e.g., desired function keywords, a partial sequence, or a target fold), the model proceeds as follows:

1. Initialize all unconditioned positions as masked tokens across all tracks.
2. Run a forward pass to obtain predicted distributions $p_\theta(t_i^{(m)} \mid \text{context})$ for all masked positions.
3. Select the top-$k$ most confident predictions (highest predicted probability) and unmask them.
4. Repeat from step 2 with the updated (partially unmasked) input until all positions are unmasked.

Formally, at decoding step $n$, the number of tokens to unmask follows a schedule:

$$
k_n = \left\lfloor L_{\text{masked}}^{(n)} \cdot \gamma(n/N) \right\rfloor,
$$

where $L_{\text{masked}}^{(n)}$ is the number of remaining masked tokens at step $n$, $N$ is the total number of decoding steps, and $\gamma: [0,1] \to [0,1]$ is a monotonically increasing schedule function. This process is analogous to the reverse process in discrete diffusion models, where noise is progressively removed.

The iterative nature of decoding is crucial: it allows the model to resolve dependencies between positions. Early steps establish the global fold and key functional motifs, while later steps refine local details and ensure consistency across modalities.

This process can be contrasted with autoregressive generation (used by natural language models like GPT), where tokens are generated left-to-right. Proteins do not have a natural left-to-right ordering for design purposes -- a residue at position 200 may be spatially adjacent to a residue at position 10 and must be designed in concert. The iterative masked decoding allows ESM-3 to handle these long-range dependencies by making globally informed decisions at each step, rather than committing to a fixed generation order.

The connection to discrete diffusion models is worth noting formally. In a discrete diffusion framework, a forward process progressively masks tokens:

$$
q(\mathbf{t}^{(n)} \mid \mathbf{t}^{(0)}) = \prod_i \left[ (1 - \beta_n) \cdot \mathbb{1}[t_i^{(n)} = t_i^{(0)}] + \beta_n \cdot \mathbb{1}[t_i^{(n)} = \texttt{[MASK]}] \right],
$$

where $\beta_n$ is the masking probability at noise level $n$, $\mathbf{t}^{(0)}$ is the original token sequence, and $\mathbf{t}^{(n)}$ is the noised version. The reverse process (generation) then progressively unmasks tokens using the model's predictions. ESM-3's training with uniform masking rates and its iterative decoding procedure can be understood as an instance of this framework.

---

## The GFP Story

Perhaps the most striking result from the ESM-3 paper is the generation of <span style="background-color: #fff3b0;">esmGFP</span>, a novel green fluorescent protein that was designed by the model and experimentally validated to be functional.

### Background

Green fluorescent protein (GFP) is one of the most important tools in biology, earning the 2008 Nobel Prize in Chemistry. GFP enables researchers to visualize proteins and cellular processes under a microscope by fusing it to proteins of interest. The GFP chromophore -- the part of the protein that actually glows -- forms autocatalytically from three specific amino acids (Ser65-Tyr66-Gly67 in wild-type GFP from _Aequorea victoria_). This chromophore formation is exquisitely sensitive to the surrounding protein structure: the three residues must be precisely positioned within an 11-stranded beta-barrel, and even small perturbations to the barrel's hydrogen bonding network or interior packing can abolish fluorescence entirely.

This makes GFP an exceptionally demanding test case for protein design. The protein must simultaneously satisfy:

- **Structural requirements:** A complete beta-barrel fold with correct strand topology
- **Chemical requirements:** Precise positioning of the chromophore triad for autocatalytic cyclization
- **Dynamic requirements:** The barrel must be rigid enough to exclude water from the chromophore environment (solvent exposure quenches fluorescence)

Prior protein engineering efforts on GFP have typically modified only a handful of residues at a time, staying within 90-95% sequence identity to natural GFPs.

### How esmGFP Was Generated

The authors used ESM-3 to generate new fluorescent proteins through a multi-step process. The generation was conditioned on:

- The GFP functional annotation (fluorescence-related GO terms and InterPro annotations)
- Key structural features of the beta-barrel fold (structure tokens encoding the 11-stranded barrel geometry)
- The critical chromophore-forming residues (sequence constraints at the active site positions)

The iterative decoding procedure first established the global fold through structure tokens, then progressively filled in the amino acid sequence while maintaining consistency between all three tracks. Through experimental screening of generated candidates, one design -- esmGFP -- was confirmed to exhibit bright green fluorescence with an excitation/emission profile similar to natural GFPs.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure4-gfp.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 5: The esmGFP story.</strong> ESM-3 generated a novel fluorescent protein (esmGFP) by conditioning on the GFP functional annotation and key structural constraints. The generated protein was experimentally validated to exhibit bright green fluorescence despite having only 58% sequence identity to the closest known fluorescent protein. This level of sequence divergence is comparable to the evolutionary distance between proteins separated by approximately 500 million years of natural evolution, roughly the time since the Cambrian explosion. The figure shows the generated structure, the fluorescence spectrum, and evolutionary distance analysis.
</div>

### Evolutionary Distance Analysis

What makes esmGFP remarkable is its sequence divergence from known GFPs. According to the authors, esmGFP has only <span style="background-color: #fff3b0;">58% sequence identity</span> to the closest known fluorescent protein in nature. To put this in perspective:

- Typical protein engineering campaigns modify a few percent of residues
- Directed evolution experiments rarely reach below 80% identity while maintaining function
- The 58% identity level corresponds to an evolutionary distance equivalent to approximately <span style="background-color: #fff3b0;">500 million years</span> of natural evolution -- roughly the time since the Cambrian explosion, when most modern animal phyla first appeared

To quantify this divergence, consider that sequence identity between homologous proteins decays roughly exponentially with evolutionary time. The relationship can be approximated as:

$$
\text{Identity}(t) \approx \text{Id}_{\text{random}} + (\text{Id}_0 - \text{Id}_{\text{random}}) \cdot e^{-\lambda t},
$$

where $\text{Id}_0$ is the initial identity (100% for the same protein), $\text{Id}_{\text{random}} \approx 5\%$ is the identity expected by chance for random sequences, $\lambda$ is the substitution rate, and $t$ is evolutionary time. For GFP-family proteins, 58% identity places esmGFP at a divergence point corresponding to roughly the Cambrian explosion.

This result suggests that ESM-3 has learned enough about the relationship between sequence, structure, and function to "simulate" evolutionary-scale exploration of protein space, generating functional proteins that are far outside the space of known natural sequences while preserving the delicate structural requirements for fluorescence.

The significance of this result extends beyond GFP itself. It demonstrates that a language model trained on natural protein sequences can generalize to regions of sequence space that nature has not explored -- or at least, that we have not yet observed in sequence databases. The model is not merely interpolating between known sequences; it is extrapolating to genuinely novel proteins that maintain functional integrity.

The name of the paper -- "Simulating 500 million years of evolution with a language model" -- comes directly from this finding.

---

## Benchmarks and Performance

### ESMFold: Structure Prediction

ESMFold demonstrated competitive structure prediction quality while being dramatically faster than MSA-based methods. According to Lin et al., ESMFold achieves accuracy approaching AlphaFold2 on proteins where high-quality MSAs are available, and the gap narrows with increasing language model scale.

The primary evaluation metric is the predicted local distance difference test (pLDDT), which measures per-residue structural accuracy on a 0-100 scale:

$$
\text{pLDDT}_i = 100 \times \mathbb{E}\left[ \frac{1}{|\mathcal{R}_i|} \sum_{j \in \mathcal{R}_i} \mathbb{1}\left[ |d_{ij}^{\text{pred}} - d_{ij}^{\text{true}}| < 0.5\text{\AA} \right] \right],
$$

where $\mathcal{R}_i$ is the set of residues within a distance cutoff of residue $i$, $d_{ij}^{\text{pred}}$ and $d_{ij}^{\text{true}}$ are predicted and true inter-residue distances, and the indicator function checks whether the distance error is below the threshold. Scores above 70 generally indicate reliable predictions, while scores above 90 suggest near-experimental accuracy.

Key performance characteristics of ESMFold:

- **Speed:** Order-of-magnitude faster than AlphaFold2 due to eliminating MSA search. The MSA computation step, which involves searching large sequence databases (UniRef, BFD, MGnify), typically takes minutes to hours per protein and dominates AlphaFold2's total runtime. ESMFold replaces this with a single forward pass through the language model.
- **Accuracy on well-studied proteins:** Near AlphaFold2 quality for proteins with many homologs in sequence databases
- **Accuracy on orphan proteins:** Lower accuracy when few homologous sequences exist, since the language model has less implicit evolutionary information to draw on. This represents the fundamental trade-off: MSA-based methods can explicitly leverage evolutionary co-variation, while language models must have learned these patterns during pre-training.
- **Metagenomic applications:** Enabled structure prediction for hundreds of millions of metagenomic sequences from environmental samples, a task computationally infeasible with MSA-based methods. According to the authors, this allowed structural characterization of proteins from organisms that have never been cultured in a laboratory.
- **Scaling behavior:** The quality of structure predictions improves systematically with the size of the underlying language model. Moving from ESM-2 (650M parameters) to ESM-2 (3B) to ESM-2 (15B) yields progressively better structures, suggesting that larger models capture more subtle evolutionary constraints.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure5-structure-prediction.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 6: Structure prediction benchmarks.</strong> Comparison of ESMFold and AlphaFold2 on protein structure prediction tasks. ESMFold achieves competitive accuracy on well-characterized proteins while running significantly faster due to its single-sequence input design. The accuracy gap between ESMFold and AlphaFold2 depends on MSA depth: for proteins with many homologs, the methods perform similarly, while AlphaFold2 retains an advantage on hard targets with sparse evolutionary information.
</div>

### ESM-3: Scaling and Generation Quality

The ESM-3 paper reports consistent improvements across model scales (1.4B, 7B, 98B parameters), with performance scaling log-linearly with compute across multiple evaluation metrics.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/esm3/figure6-scaling.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 7: Scaling laws in ESM-3.</strong> Performance on protein-related tasks improves log-linearly with compute across the three model scales (1.4B, 7B, 98B parameters). This scaling behavior, reminiscent of scaling laws observed in natural language models, suggests that further scaling may continue to yield improvements in protein understanding and generation capabilities.
</div>

Key findings reported by the authors:

- **Unconditional generation quality:** Generated proteins exhibit natural-like properties including realistic secondary structure composition, hydrophobic packing, and predicted structural confidence scores. The authors evaluate generated proteins using metrics such as predicted TM-score (measuring structural plausibility), sequence perplexity, and biophysical property distributions.

- **Conditional generation:** When conditioned on specific folds or functions, ESM-3 generates proteins that are predicted to adopt the desired structure with high confidence. For example, given structure tokens encoding a TIM barrel fold, the model generates sequences that, when folded by ESMFold or AlphaFold2, adopt the specified fold with high TM-scores.

- **Representation learning:** The model's internal representations capture meaningful biological features, with performance on downstream tasks (fold classification, function prediction, fitness landscape prediction) improving with model scale. The representations learned by ESM-3 outperform those from sequence-only models on several representation benchmarks, suggesting that multimodal training produces richer protein embeddings.

- **Cross-modal reasoning:** The model demonstrates genuine cross-modal capabilities -- for example, generating plausible structures from function descriptions alone, or inferring function from structural motifs. This bidirectional reasoning between modalities is a capability that no single-task model can provide.

- **Fitness landscape prediction:** When used to score mutations, ESM-3's log-likelihoods correlate with experimentally measured fitness effects, suggesting that the model has learned the relationship between sequence variation and functional consequences.

### Comparison with Specialized Models

While ESM-3 is a generalist model, the authors evaluate it against specialized baselines across multiple tasks:

- **Structure prediction:** ESMFold (the precursor) and AlphaFold2 remain strong baselines. ESM-3 is not primarily designed for structure prediction, but its structure token predictions can be decoded into coordinates.
- **Inverse folding:** For designing sequences that fold into given structures, ProteinMPNN serves as a comparison point. ESM-3 can perform inverse folding by conditioning on structure tokens and generating sequence tokens.
- **Unconditional protein generation:** Diffusion-based models like RFdiffusion provide a reference for generating novel protein structures. ESM-3's iterative decoding offers a discrete alternative to the continuous diffusion process.
- **Function prediction:** Models like DeepFRI and InterProScan provide baselines for predicting protein function from sequence or structure.

The advantage of ESM-3 is not necessarily surpassing each specialist at their particular task, but rather providing a single unified model that can handle all of these tasks and, critically, the compositional combinations between them (e.g., "generate a protein with this fold AND this function"). This compositionality is a unique capability that emerges from the multimodal training framework.

---

## Limitations and Future Directions

### 1. Structure Token Fidelity

The VQ-VAE discretization of structure inherently involves information loss. While the codebook of 4,096 entries captures local geometry well, subtle structural details (e.g., precise side-chain rotamer states, backbone torsion angles at sub-angstrom precision) may be lost in the quantization process. This sets a ceiling on the structural accuracy achievable through the token-based representation.

The fundamental trade-off can be expressed as a rate-distortion problem:

$$
D(R) = \min_{\text{encoder, decoder}} \mathbb{E}\left[ \| \mathbf{C}_{\text{true}} - \hat{\mathbf{C}} \|^2 \right], \quad \text{s.t. } R \leq \log_2(V_c) \text{ bits per residue},
$$

where $D(R)$ is the minimum achievable reconstruction distortion at rate $R$, and $V_c = 4{,}096$ gives $R = 12$ bits per residue. Increasing the codebook size would reduce distortion but increase the vocabulary the transformer must handle, potentially degrading generation quality.

### 2. Experimental Validation Gap

While esmGFP is an impressive demonstration, the experimental validation in the paper is limited to a small number of designed proteins. The success rate of ESM-3 designs across diverse protein families and functions remains to be systematically characterized. Fluorescent proteins, while challenging, have well-established experimental assays -- extending to more complex functions (e.g., enzymatic catalysis, protein-protein interactions) will require more extensive wet-lab validation campaigns.

### 3. Computational Cost

The 98B parameter model requires substantial computational resources for both training and inference. While the iterative decoding procedure is more efficient than some diffusion-based approaches, generating a single protein still requires multiple forward passes through a very large model. This may limit accessibility for researchers without access to large-scale computing infrastructure.

### 4. Limited Multimodal Training Data

While the sequence database is vast (2.78 billion proteins), the number of proteins with all three modalities annotated is much smaller. Many proteins have sequence data but no experimental structure or detailed functional annotation. This data imbalance means the model's cross-modal capabilities are primarily trained on predicted (rather than experimental) structures and computationally inferred (rather than experimentally verified) functions.

### 5. Single-Chain Focus

ESM-3 currently operates on individual protein chains. Many biological functions arise from protein complexes, multi-domain interactions, and protein-nucleic acid assemblies. For context, according to some estimates, over 80% of proteins function as part of multi-protein complexes in the cell. Extending the framework to handle multi-chain systems and non-protein biomolecules (DNA, RNA, small molecules) represents a natural but challenging next step. AlphaFold-Multimer has shown that structure prediction can be extended to complexes, but generative modeling of complexes introduces additional challenges around stoichiometry, interface design, and co-evolutionary constraints between chains.

### 6. Controllability and Reliability

While ESM-3 demonstrates function-conditioned generation, the degree of control and the reliability of functional outcomes remain open questions. Generating a protein that is predicted to have a desired function is different from generating one that actually works in a cell. Improving the precision of functional control and developing better computational validators for generated designs are important directions.

### Future Extensions

- **Multimodal fine-tuning with experimental feedback:** Integrating experimental screening results (e.g., from high-throughput functional assays) into the training loop to improve design success rates. This could leverage reinforcement learning from human feedback (RLHF) or direct preference optimization (DPO) approaches adapted from natural language models.
- **Multi-chain and complex generation:** Extending ESM-3 to generate protein complexes, protein-ligand systems, and protein-nucleic acid assemblies. This would require extending the tokenization scheme to handle multiple chains and inter-chain interactions.
- **Integration with laboratory automation:** Connecting the model to automated protein synthesis and testing pipelines for rapid design-build-test cycles. The combination of computational protein design with robotic wet-lab platforms could dramatically accelerate the discovery of functional proteins.
- **Alignment and safety:** Developing methods to ensure generated proteins are safe and beneficial, particularly as the model's capabilities improve. The authors of ESM-3 have noted the importance of responsible development, and EvolutionaryScale has implemented a responsible use framework. As protein generation models become more capable, ensuring that they cannot be misused to design harmful biological agents becomes increasingly important.
- **Higher-resolution structure generation:** Moving beyond discrete structure tokens to directly generate atomic coordinates, potentially through hybrid approaches that combine the advantages of discrete tokenization (for the transformer backbone) with continuous refinement (for final structure output).

---

## Key Takeaways

- **ESMFold proved that protein language models encode enough evolutionary information to predict 3D structure from a single sequence**, eliminating the need for MSA computation and enabling structure prediction at metagenomic scale.
- **ESM-3 extends the language model paradigm from understanding to generation**, treating proteins as sequences of multimodal tokens (sequence, structure, function) that can be generated through iterative masked decoding.
- **The multimodal tokenization framework** -- particularly the VQ-VAE for structure -- enables a standard transformer to reason about protein geometry without specialized structural modules.
- **esmGFP demonstrates that ESM-3 can generate genuinely novel functional proteins**, with the generated fluorescent protein showing only 58% sequence identity to any known GFP, a level of divergence equivalent to roughly 500 million years of natural evolution.
- **Scaling laws hold for protein language models**, with performance improving log-linearly with compute across the 1.4B to 98B parameter range, suggesting that further scaling may unlock additional capabilities.
- **The unified multimodal approach enables compositional protein design** -- specifying desired function, fold, or partial sequence and generating the rest -- opening new possibilities for protein engineering that go beyond what specialized single-task models can achieve.
- **Key limitations remain** in experimental validation breadth, multi-chain modeling, computational cost, and the fidelity of discrete structure representations, pointing to rich opportunities for future work.
- **The broader significance** lies in demonstrating that the language model paradigm -- which has transformed natural language processing, computer vision, and code generation -- can be successfully applied to the fundamental challenge of understanding and designing the molecular machines of life. Proteins are arguably the most complex "language" that exists, and the success of ESM-3 suggests we are only beginning to understand what large-scale models can learn about biology.

---

_References:_

_Hayes, T. et al., "Simulating 500 million years of evolution with a language model," bioRxiv (2024). DOI: [10.1101/2024.07.01.600583](https://doi.org/10.1101/2024.07.01.600583)_

_Lin, Z. et al., "Evolutionary-scale prediction of atomic-level protein structure with a language model," Science 379, 1123--1130 (2023). DOI: [10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)_

_Code and models: [ESM on GitHub](https://github.com/evolutionaryscale/esm)_
