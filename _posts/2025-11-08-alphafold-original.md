---
layout: post
title: "AlphaFold: The AI That Learned How Proteins Fold"
date: 2025-11-08
permalink: /blog/alphafold-original/
published: false
categories: [paper-review]
tags:
  - AI
  - biology
  - deep-learning
  - science
---

Imagine you're given a long string of beads, each bead representing one of the twenty amino acids that make up life's proteins. Now, without touching it, you must predict exactly how that string twists and folds into a three dimensional shape that decides whether it becomes silk, muscle, or an enzyme. For decades this challenge, known as the **protein folding problem**, baffled scientists.

Formally, we are given an amino-acid sequence

$$
x = (x_1, x_2, \dots, x_L), \quad x_i \in \{1,\dots,20\},
$$

where $$x$$ is the amino-acid sequence, each $$x_i$$ indexes one residue in the 20-letter alphabet, and $$L$$ is the total sequence length.

and asked to predict the 3D coordinates of all atoms

$$
\hat{\mathbf{X}} = \{\hat{\mathbf{x}}_{i,a} \in \mathbb{R}^3
\mid i = 1,\dots,L,\ a \in \text{atoms of residue } i\}.
$$

where $\hat{\mathbf{X}}$ collects predicted coordinates, $\hat{\mathbf{x}}_{i,a}$ is the 3D position of atom $a$ in residue $i$, and $\mathbb{R}^3$ denotes ordinary 3D space.

In 2021, Google DeepMind's **AlphaFold** shocked the world by _solving_ much of it. Published in _Nature_, the model predicted protein shapes with almost experimental accuracy, earning headlines like "the greatest breakthrough in biology since the human genome."

Explore the full [Nature article](https://www.nature.com/articles/s41586-021-03819-2) and its [supplementary material (PDF)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03819-2/MediaObjects/41586_2021_3819_MOESM1_ESM.pdf) for deeper technical details.

---

## Introduction

The authors outline seven novel contributions that distinguish the AlphaFold2 model.

1. **Evoformer for joint embeddings.** A new architecture jointly embeds <span style="background-color: #fff3b0;">multiple sequence alignments (MSAs)</span> and pairwise features so AlphaFold can reason about evolutionary couplings and spatial context at the same time.

   Mathematically, AlphaFold maintains two main latent tensors:

   $$
   M \in \mathbb{R}^{N_{\text{msa}} \times L \times d_m}
   \quad\text{(MSA representation)},
   $$

   where $$N_{\text{msa}}$$ is the number of alignment rows, $$L$$ the residue length, and $$d_m$$ the channel width of the MSA embedding.

   $$
   Z \in \mathbb{R}^{L \times L \times d_z}
   \quad\text{(pair representation)}.
   $$

   where $d_z$ controls the dimensionality of the pair features over all residue pairs $L \times L$.

   Evoformer is a deep stack of layers that iteratively update $M$ and $Z$ with cross-talk between them.

2. **Backbone-frame output representation.** By predicting residue-level rigid frames plus atom-level torsion angles, the network adds a loss that directly supervises 3-D geometry, enabling end-to-end structure prediction.

   Each residue $i$ is represented by a rigid body frame

   $$
   F_i = (R_i, \mathbf{t}_i),\quad
   R_i \in \mathrm{SO}(3),\ \mathbf{t}_i \in \mathbb{R}^3,
   $$

   where $$F_i$$ encodes residue $$i$$, $$R_i$$ is its rotation in $$\mathrm{SO}(3)$$, and $$\mathbf{t}_i$$ is its translation vector in $$\mathbb{R}^3$$, and internal atom positions are defined via learned torsion angles relative to that frame.

3. **Invariant Point Attention (IPA).** The equivariant attention design lets AlphaFold consider 3-D distances between residues while ignoring global translations or rotations, so attention focuses on true spatial relationships. If we transform all frames by a global rotation $Q \in \mathrm{SO}(3)$ and translation $\mathbf{u} \in \mathbb{R}^3$,

   $$
   (R_i, \mathbf{t}_i) \mapsto (Q R_i, Q\mathbf{t}_i + \mathbf{u}),
   $$

   where $Q$ is any global rotation and $\mathbf{u}$ a translation applied uniformly to every residue frame. IPA is constructed so that its outputs transform in the same way (equivariance), and its attention weights depend only on relative geometry.

4. **Intermediate loss signals.** Losses applied at several depths encourage iterative refinement of the coordinates, making each pass of the network sharpen the structure. Instead of only supervising the final structure, AlphaFold adds structure losses at multiple recycling iterations.

5. **Masked MSA loss.** Similar to BERT, parts of the MSA are deliberately masked and reconstructed, forcing the model to learn richer sequence statistics.

   If $M$ denotes the MSA tokens and $\mathcal{M}$ the set of masked positions, the model optimizes a cross-entropy loss:

   $$
   \mathcal{L}_{\text{msa-mask}}
   =
   - \sum_{(s,i)\in\mathcal{M}}
     \log p_{\theta}\big(M_{s,i} \mid M_{\text{masked}}\big),
   $$

   where:

   - $$\mathcal{L}_{\text{msa-mask}}$$ is the masked-token loss,
   - $$(s,i)$$ indexes MSA row and residue,
   - $$\mathcal{M}$$ is the masked set,
   - $$p_{\theta}$$ is the model's predicted distribution.

6. **Self-distillation on unlabeled sequences.** A noisy-student self-training loop (inspired by the CVPR 2020 paper _Self-training with Noisy Student improves ImageNet classification_) lets AlphaFold learn from vast unlabeled protein sequences. A "teacher" AlphaFold model creates pseudo-labels (structures and distances), and a "student" model learns to match them under strong augmentations.

7. **Self-estimated accuracy.** AlphaFold outputs per-residue pLDDT scores, giving scientists a built-in estimate of how trustworthy each part of the predicted structure is. Concretely, the network predicts a 50-bin distribution over "local distance difference test" (LDDT) scores; the expected value becomes pLDDT.

---

## How Is AlphaFold2 Different?

### Feature Representation Level

Classic AlphaFold and similar models ingest an $L \times L \times c$ tensor (sequence length $L$, feature channels $c$) by tiling sequence-length features such as MSAs across both axes so they match the pairwise feature grid, then feeding that stack into the network. The pair features (covariation, contact priors, etc.) are naturally $L \times L \times c$, but the sequence-length features are duplicated to fit that same shape.

AlphaFold2 instead embeds the **MSA track** and the **pair track** separately, then lets the **Evoformer** exchange information between them.

- MSA representation:
  $M \in \mathbb{R}^{N_{\text{msa}} \times L \times d_m}$
- Pair representation:
  $Z \in \mathbb{R}^{L \times L \times d_z}$

At a high level, a single Evoformer block performs:

1. **MSA row & column attention** (information flow across both sequence and alignment axes),
2. **MSA–pair cross-talk** via outer-product mean,
3. **Triangle updates** on the pair representation (two-body reasoning on residue triplets).

<img src="/assets/img/alphafold/triangle-update-types.png" alt="Triangle update motifs" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
*Figure 1. Triangle multiplicative updates and triangle self-attention reason over triples $(i,j,k)$ by passing information along outgoing or incoming edges before returning to the $ij$ pair representation.*

Formally, one can summarize an Evoformer block as

$$
M' = \Phi_{\text{msa}}(M, Z), \qquad
Z' = \Phi_{\text{pair}}(Z, M),
$$

where $$\Phi_{\text{msa}}$$ updates the MSA tensor using both $$M$$ and $$Z$$, and $$\Phi_{\text{pair}}$$ updates the pair tensor with information from $$Z$$ and $$M$$. These operators are compositions of attention, outer-product, and feed-forward sublayers with residual connections.

As far as the literature shows, no prior deep-learning protein-folding system independently embedded the MSA and pair representations while allowing structured cross-talk the way Evoformer does.

<img src="/assets/img/alphafold/evoformer-overview.png" alt="Evoformer overview diagram" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
*Figure 2. Each Evoformer block alternates row/column attention on the MSA track with triangle updates and transitions on the pair track, keeping both tensors in lock-step across roughly 48 unshared layers.*

### End-to-End Modeling

AlphaFold, trRosetta (Yang et al., PNAS 2019), and related models first train a network to predict residue–residue distance distributions (distograms) and then solve a downstream optimization problem to find a structure that fits those constraints, essentially a two-step pipeline of distogram prediction followed by energy minimization.

AlphaFold2 represents each residue as a rigid backbone frame $(R_i, \mathbf{t}_i)$ and assumes residue internals depend only on torsion angles once that frame is set. By directly predicting the frame transforms and the torsion angles, AlphaFold2 computes every atom's 3-D position, compares it with experimental structures, and **backpropagates the loss in one sweep**. That makes the whole system an **end-to-end model** rather than a cascade.

---

## A Mathematical Glimpse Inside AlphaFold2

This section summarizes the core computations in slightly more formal terms, while staying high-level enough for a blog post.

### Representations and Input Embedding

- Target sequence (one-hot):
  $$X \in \{0,1\}^{L \times 21}$$
- Residue indices:
  $$\text{idx} \in \mathbb{Z}^{L}$$
- Clustered MSA features:
  $$\text{MSA}_{\text{feat}} \in \mathbb{R}^{N_{\text{msa}} \times L \times d_{\text{in}}}$$
- Template features, extra MSA, etc.

The **InputEmbedder** maps these into initial MSA and pair embeddings:

$$
M^{(0)} = f_{\text{msa}}(\text{MSA}_{\text{feat}}, X),
\quad
Z^{(0)} = f_{\text{pair}}(X, \text{idx}),
$$

where $$f_{\text{msa}}$$ embeds raw MSA features and the target sequence $$X$$, while $$f_{\text{pair}}$$ uses sequence tokens and residue indices $$\text{idx}$$ to seed the pair representation.

<img src="/assets/img/alphafold/input-embedding-pipeline.png" alt="Input embedding pipeline" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
*Figure 3. InputEmbedder combines clustered MSAs, residue indices, templates, and extra MSA rows before handing the aligned MSA/pair tensors to the Evoformer stack.*

including relative positional encodings of the form

$$
d_{ij} = \mathrm{clip}(\text{idx}_i - \text{idx}_j, -K, K),
\quad
\text{relpos}_{ij} = \mathrm{one\_hot}(d_{ij}),
$$

where $$d_{ij}$$ is the clipped residue index difference between $$i$$ and $$j$$, $$K$$ is the clipping threshold, and $$\text{relpos}_{ij}$$ is its one-hot encoding. These are linearly projected and added to $$Z^{(0)}$$.

### Recycling, Templates, and Extra MSA (Conceptual Overview)

Before each Evoformer pass, AlphaFold adds three sources of information:

1. **Recycling:**
   From the previous iteration's predicted frames $\{F_i^{\text{prev}}\}$ and representations, it computes:

   - Distances between Cβ atoms:

   $$
   d^{\text{prev}}_{ij} = \|\mathbf{x}^{\text{prev}}_{i, C_\beta} - \mathbf{x}^{\text{prev}}_{j, C_\beta}\|_2,
   $$

   - Linear projections of previous $M$ and $Z$,

   and adds them as biases to the current $M^{(0)}$ and $Z^{(0)}$.

2. **Templates:**
   Template-based angle features per residue and per pair are embedded and (i) some are treated as extra MSA rows, (ii) some are passed through a small Evoformer-like stack and injected into the pair representation $Z$.

3. **Extra MSA:**
   A large set of additional alignment rows is processed by a compact "ExtraMsaStack" whose output is used to update the pair representation $Z$ without blowing up memory.

Symbolically,

$$
M^{(0)}_{\text{final}} = M^{(0)} + \Delta M_{\text{recycle}} + \Delta M_{\text{template}},
$$

where $$M^{(0)}_{\text{final}}$$ mixes the base embedding with corrections from recycling ($$\Delta M_{\text{recycle}}$$) and template cues ($$\Delta M_{\text{template}}$$).

$$
Z^{(0)}_{\text{final}} = Z^{(0)} + \Delta Z_{\text{recycle}}
                        + \Delta Z_{\text{template}} + \Delta Z_{\text{extra\_msa}},
$$

where the pair seed combines the raw pair features with recycled, template, and extra-MSA contributions.

These are then fed through the main Evoformer stack.

### Evoformer Block (Very High-Level)

A single Evoformer block repeatedly applies:

1. **MSA row attention**

   For each MSA row $s$ and residue $i$,

   $$
     \mathrm{Attn}^{\text{row}}(M)_{s,i}
     = \sum_{j} \alpha_{s,i,j} V_{s,j},
   $$

   where $$\mathrm{Attn}^{\text{row}}$$ aggregates row information for MSA row $$s$$ at residue $$i$$, using attention weights $$\alpha_{s,i,j}$$ over values $$V_{s,j}$$, and

   $$
     \alpha_{s,i,j}
     = \mathrm{softmax}_j\!\big(
         Q_{s,i}^\top K_{s,j} / \sqrt{d_m} + b_{ij}
       \big),
   $$

   where $$Q_{s,i}$$ and $$K_{s,j}$$ are query and key vectors, $$d_m$$ scales the dot product, $$b_{ij}$$ is a bias derived from the pair representation $$Z_{ij}$$ (allowing geometric context to modulate sequence attention), and the softmax normalizes across positions $$j$$.

   <img src="/assets/img/alphafold/row-attention-block.png" alt="Row attention block" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
   *Figure 4. Row-wise attention mixes residues within a single MSA row while injecting pair-derived biases before writing the updates back into $$M$$.*

2. **MSA column attention** (attend along the MSA axis instead of the residue axis).

<img src="/assets/img/alphafold/column-attention-block.png" alt="Column attention block" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
*Figure 5. Column-wise attention processes the stack of MSA rows for a single residue position, ensuring homologous sequences vote on residues $$i$$ consistently.*

3. **Outer-product mean** from MSA to pair:

   $$
   \mathrm{OP}_{ij}
   =
   \frac{1}{N_{\text{msa}}}
   \sum_{s=1}^{N_{\text{msa}}}
   (W_1 M_{s,i}) \otimes (W_2 M_{s,j}),
   $$

   where $$\mathrm{OP}_{ij}$$ captures correlations between residues $$i$$ and $$j$$, $$W_1$$ and $$W_2$$ are learned projections, and $$\otimes$$ denotes an outer product averaged over all $$N_{\text{msa}}$$ rows, which is added to the pair representation $$Z_{ij}$$.

4. **Triangle multiplicative updates** and **triangle attention** on $$Z$$, which model interactions among triplets of residues.
   Multiplicative paths treat the pair matrix as edges of a fully connected graph and pass messages around triangles $$(i,j,k)$$, while triangle attention aggregates those messages with dot-product attention across both "left" and "right" edges before projecting back to the $$ij$$ entry.

   <img src="/assets/img/alphafold/triangle-multiplicative-block.png" alt="Triangle multiplicative block" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
   *Figure 6. Triangle multiplicative updates gate information separately along left and right edges before normalizing and writing the result to the $$ij$$ slot of $$Z$$.*

   <img src="/assets/img/alphafold/triangle-attention-block.png" alt="Triangle attention block" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
   *Figure 7. Triangle attention treats the triplet edges as attention keys/values so that each pair chooses which intermediate residue $$k$$ best explains its geometry.*

   All sublayers live inside residual blocks with layer normalization and feed-forward networks.

### Structure Module and FAPE Loss

After several Evoformer blocks (and recycling iterations), AlphaFold passes the "single" representation (roughly, the first MSA row after processing) and the pair representation into the **Structure Module**. This module uses **Invariant Point Attention** (IPA) to update residue frames $F_i = (R_i, \mathbf{t}_i)$.

<img src="/assets/img/alphafold/structure-module-ipa.png" alt="Structure module IPA" class="zoomable" style="width:82%;max-width:900px;display:block;margin:0 auto;" />
*Figure 8. The Structure Module combines pair biases, single-sequence features, and learned query/key/value points to perform IPA before emitting updated frames and coordinates.*

The key geometric loss is the **Frame Aligned Point Error** (FAPE). For each residue $i$, each atom $a$, and some supervising frame $F_i^\star = (R_i^\star, \mathbf{t}_i^\star)$ from the experimental structure, we look at the atom expressed in the local frame:

- Ground truth in local frame:

$$
  \mathbf{y}_{i,a}
  = R_i^{\star\top}(\mathbf{x}^\star_{i,a} - \mathbf{t}_i^\star),
$$

where $$\mathbf{y}_{i,a}$$ is the ground-truth atom $$a$$ from residue $$i$$ expressed in the supervising local frame $$(R_i^\star, \mathbf{t}_i^\star)$$.

- Prediction in the same frame:

$$
  \hat{\mathbf{y}}_{i,a}
  = R_i^{\star\top}(\hat{\mathbf{x}}_{i,a} - \mathbf{t}_i^\star),
$$

where $$\hat{\mathbf{y}}_{i,a}$$ is the predicted atom positioned in that same local frame for comparison.

The FAPE term is

$$
\mathcal{L}_{\text{FAPE}}
=
\frac{1}{N_{\text{atoms}}}
\sum_{i,a}
\mathrm{clamp}\big(
\|\hat{\mathbf{y}}_{i,a} - \mathbf{y}_{i,a}\|_2,\ d_{\text{cut}}
\big),
$$

where $$\mathcal{L}_{\text{FAPE}}$$ averages the clamped distance error over all atoms, $$N_{\text{atoms}}$$ is the normalization count, and $$d_{\text{cut}}$$ bounds extreme penalties.

which is **invariant** to any global rigid-body motion of the whole protein but still sensitive to local geometry.

The full training objective is a weighted sum of:

- FAPE and other coordinate-based terms,
- torsion-angle losses,
- distogram/angle distribution losses,
- pLDDT calibration losses,
- masked MSA reconstruction loss,
- self-distillation losses (matching teacher predictions).

---

## The Brain Behind AlphaFold

AlphaFold has two main "brains," shown in the classic diagrams of colored blocks flowing from left to right:

1. **Evoformer: The Relationship Builder**

   - Reads a _multiple sequence alignment_ (MSA), thousands of related sequences that reveal which amino acids evolve together.
   - Maintains the two latent tensors $M$ (MSA) and $Z$ (pair).
   - Uses Transformer-style attention and geometric updates (triangle rules, outer-product mean) to learn which parts of a protein likely touch or move together.

2. **Structure Module: The Sculptor**
   - Takes those relationships and builds a 3-D model atom by atom.
   - Uses **Invariant Point Attention**, which "looks" at the structure in 3-D space while staying unaffected by rotations, as if holding the molecule and spinning it in your hand.
   - Predicts residue frames $F_i = (R_i, \mathbf{t}_i)$ and torsion angles, then deterministically computes full-atom coordinates.

The two modules are run multiple times with **recycling**, so each pass refines the coordinates in light of the previous prediction—like an artist going over the same sculpture with finer tools.

---

## How Accurate Is It?

In the international **CASP14** competition, AlphaFold stunned everyone: for many target domains, the backbone error was **around 1 Ångström**, approaching the width of a single atom. That's the level where even crystallography experiments start to disagree with each other.

To help users judge trust in each prediction, AlphaFold reports:

| Score          | What It Means          | Typical Use                  |
| :------------- | :--------------------- | :--------------------------- |
| **pLDDT > 90** | Nearly atomic accuracy | safe for detailed modeling   |
| **70 to 90**   | Domain-level reliable  | good for backbone tracing    |
| **< 70**       | Uncertain or flexible  | may indicate loops or motion |

🖼️ _Visual idea:_ imagine a rainbow-colored protein model, where bright blue regions are rock-solid while orange and red show where the AI isn't sure.

---

## Why It Matters

AlphaFold changed biology overnight. Within months, millions of protein structures from bacteria to humans were predicted and released in the **AlphaFold Protein Structure Database**, a free, searchable atlas for researchers everywhere.

Scientists now use these models to:

- design new enzymes for green chemistry,
- understand disease mutations,
- and even build custom proteins that never existed in nature.

---

## Limits and What's Next

Like any expert, AlphaFold still has blind spots:

- It struggles when few related sequences exist (a "shallow MSA").
- It models single proteins best, while complexes of multiple chains remain trickier.
- It doesn't explicitly handle small molecules, metals, or dynamic motions.

DeepMind and others have since expanded it: **AlphaFold-Multimer** for complexes, **ESMFold** for faster predictions, and new hybrids that blend AI with physics.

---

### Key Takeaways

- AlphaFold taught AI to understand the rules of life's most fundamental building blocks.
- It bridged the gap between biological data and 3-D reality.
- Most importantly, it showed how learning from patterns, not brute force, can decode nature itself.

---

_Reference: Jumper et al.,_ **Nature 596**, 583–589 (2021). DOI: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)

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
