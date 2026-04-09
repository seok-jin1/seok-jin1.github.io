---
layout: post
title: "AlphaFold 3: Predicting the Structure of All Life's Molecules"
date: 2025-11-20
permalink: /blog/alphafold3-nature2024-explained/
published: true
categories: [paper-review]
tags:
  - AI
  - biology
  - deep-learning
  - structural-biology
  - drug-discovery
---

Imagine trying to assemble a jigsaw puzzle where the pieces are not rigid plastic but floppy, vibrating molecular chains that constantly shift shape depending on what they touch. Now imagine that the puzzle includes not just protein pieces but also DNA strands, RNA loops, drug-like small molecules, metal ions, and sugar modifications---all interacting simultaneously. For decades, structural biologists solved these puzzles one painstaking experiment at a time, using X-ray crystallography, cryo-EM, or NMR. AlphaFold 2 revolutionized the protein-only version of this puzzle in 2020, but the full molecular jigsaw---proteins interacting with nucleic acids, ligands, and each other---remained largely unsolved computationally.

In 2024, Google DeepMind's **AlphaFold 3** (AF3) tackled this broader challenge. Published in _Nature_, AF3 introduced a unified deep learning architecture capable of jointly predicting the 3-D structure of complexes containing proteins, DNA, RNA, small molecules, ions, and modified residues. The key architectural innovation is a <span style="background-color: #fff3b0;">diffusion-based structure module</span> that replaces AlphaFold 2's deterministic structure module, enabling the model to generate diverse, high-quality structural predictions by learning to denoise atomic coordinates from random noise.

This post walks through the Nature paper, examining how AF3 combines a simplified trunk architecture (Pairformer), a powerful diffusion module, and a unified tokenization scheme to predict the structure of virtually any biomolecular complex.

Explore the full [Nature article](https://www.nature.com/articles/s41586-024-07487-w) and its supplementary materials for comprehensive technical details. AF3 predictions are accessible through the [AlphaFold Server](https://alphafoldserver.com/).

---

## Novel Contributions

The authors present seven key innovations that distinguish AlphaFold 3 from its predecessor and competing methods:

1. **Unified biomolecular structure prediction.** AF3 is the first model to jointly predict the structure of complexes containing <span style="background-color: #fff3b0;">proteins, DNA, RNA, small molecules, ions, and covalent modifications</span> within a single architecture. Previous methods required separate specialized tools for each interaction type (e.g., protein-protein docking, molecular docking for ligands, RNA structure prediction).

2. **Diffusion-based structure generation.** AF3 replaces AlphaFold 2's deterministic Structure Module with a <span style="background-color: #fff3b0;">diffusion module</span> that operates directly on raw atom coordinates. This module learns to iteratively denoise atomic positions from Gaussian noise, producing full-atom structures in a generative framework. The diffusion approach naturally handles stochasticity and can generate multiple plausible conformations for a given input.

3. **Simplified trunk architecture (Pairformer).** AF3 replaces the Evoformer's coupled MSA-pair representation with a simpler <span style="background-color: #fff3b0;">Pairformer</span> that operates only on single and pair representations during the main trunk processing. MSA information is processed separately in an earlier module and summarized into the pair representation, substantially reducing computational cost.

4. **Unified token-and-atom representation.** AF3 introduces a two-level representation: a token level (one token per standard residue or ligand atom) and an atom level (all heavy atoms). This allows the model to handle standard polymer residues (proteins, nucleic acids) at residue resolution while representing small molecules and modified residues at atomic resolution within the same framework.

5. **Cross-distillation from AlphaFold 2.** To expand the effective training set, AF3 uses <span style="background-color: #fff3b0;">cross-distillation</span>: AlphaFold 2 predictions on sequences without experimental structures are used as additional training targets. This allows AF3 to learn from millions of predicted structures beyond the roughly 180,000 experimental structures in the Protein Data Bank.

6. **Confidence metrics for interfaces.** AF3 extends AlphaFold 2's confidence predictions to biomolecular interfaces, introducing metrics like <span style="background-color: #fff3b0;">interface predicted TM-score (iPTM)</span> that specifically assess the quality of predicted inter-chain contacts. A composite ranking score combines pLDDT, PAE, and iPTM to select the best prediction from multiple diffusion samples.

7. **State-of-the-art across interaction types.** The paper reports that AF3 substantially outperforms specialized methods across protein-ligand, protein-nucleic acid, and protein-protein interaction benchmarks, demonstrating that a single generalist model can surpass dedicated tools.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure1-overview.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 1: AlphaFold 3 model overview.</strong> The figure illustrates the complete AF3 pipeline. <strong>Left:</strong> Inputs include protein sequences, nucleic acid sequences, small molecule representations (as SMILES or CCD codes), ion identities, and covalent modification annotations. <strong>Center:</strong> The architecture proceeds through three stages: (1) an MSA module that processes multiple sequence alignments and summarizes evolutionary information into pair representations, (2) the Pairformer trunk that refines single and pair representations through attention-based updates without maintaining a full MSA track, and (3) the diffusion module that generates 3-D atomic coordinates by iteratively denoising from random noise, conditioned on the trunk representations. <strong>Right:</strong> The output is a full-atom 3-D structure of the biomolecular complex, along with per-residue and pairwise confidence estimates (pLDDT, PAE, iPTM).
</div>

---

## How Is AlphaFold 3 Different?

**_Architecture and scope._** AlphaFold 2 was designed exclusively for single-chain protein structure prediction (later extended to multimers via AlphaFold-Multimer). It used a 48-layer Evoformer that jointly processed MSA and pair representations, followed by a deterministic Structure Module that predicted backbone frames and side-chain torsion angles. AF3 takes a fundamentally different approach:

| Feature                  | AlphaFold 2                                | AlphaFold 3                                               |
| :----------------------- | :----------------------------------------- | :-------------------------------------------------------- |
| **Scope**                | Proteins only                              | Proteins, DNA, RNA, ligands, ions, modifications          |
| **MSA processing**       | Evoformer (joint MSA + pair, 48 layers)    | Separate MSA module, then Pairformer (pair + single only) |
| **Structure generation** | Deterministic (backbone frames + torsions) | Diffusion-based (raw atom coordinates)                    |
| **Representation**       | Residue-level only                         | Token-level + atom-level (unified)                        |
| **Output format**        | Backbone frames + chi angles               | Full-atom Cartesian coordinates                           |
| **Sampling**             | Recycling (deterministic)                  | Stochastic diffusion (multiple samples)                   |
| **Training data**        | PDB proteins                               | PDB + cross-distilled AF2 predictions                     |
| **Confidence**           | pLDDT, PAE                                 | pLDDT, PAE, iPTM, ranking score                           |

**_Why diffusion?_** The deterministic Structure Module in AF2 produces a single prediction per input, making it difficult to capture conformational heterogeneity or assess uncertainty through sampling. AF3's diffusion module generates structures by sampling from a learned distribution, allowing multiple plausible conformations to be produced and ranked. This is particularly important for flexible interfaces (protein-ligand binding, protein-RNA contacts) where a single structure may not adequately represent the ensemble of possible conformations.

**_Why simplify the trunk?_** The Evoformer's MSA track was computationally expensive, requiring $O(N_{\text{seq}} \times L)$ memory and compute at each of 48 layers. By processing MSAs separately and summarizing the information into pair features before the main trunk, AF3 reduces compute while retaining the essential evolutionary information. The paper reports that this simplification does not degrade accuracy and may even improve generalization.

---

## A Mathematical Glimpse Inside

AlphaFold 3's architecture can be decomposed into four major components: input processing, the Pairformer trunk, the diffusion module, and confidence prediction. Let's formalize each.

### Input Representation

AF3 accepts a complex consisting of $K$ chains, where each chain $k$ has a sequence of tokens $\mathbf{t}^{(k)} = (t_1^{(k)}, \ldots, t_{L_k}^{(k)})$. For proteins, each token is an amino acid residue; for nucleic acids, each token is a nucleotide; for small molecules, each atom is a separate token.

The model constructs two initial representations:

1. **Single representation** $s_i \in \mathbb{R}^{d_s}$ for each token $i$, encoding its identity, chain membership, and local features.
2. **Pair representation** $z_{ij} \in \mathbb{R}^{d_z}$ for each pair of tokens $(i, j)$, encoding relative position, chain identity, and evolutionary coupling information.

Additionally, at the atom level, each token $i$ maps to a set of atoms $\mathcal{A}_i$ with local reference frames. For a standard amino acid, $\mathcal{A}_i$ contains all heavy atoms of that residue; for a ligand atom token, $\mathcal{A}_i$ contains just that atom.

### MSA Module

Multiple sequence alignments are processed in a separate module before the main trunk. Given an MSA matrix $M \in \mathbb{R}^{N_{\text{seq}} \times L \times d}$, the MSA module applies row-wise and column-wise attention to extract coevolutionary patterns:

$$
M' = \text{ColumnAttn}\bigl(\text{RowAttn}(M)\bigr)
$$

The processed MSA is then summarized into the pair representation via an outer product mean operation:

$$
z_{ij}^{\text{MSA}} = \frac{1}{N_{\text{seq}}} \sum_{n=1}^{N_{\text{seq}}} m_{ni} \otimes m_{nj}
$$

where $m_{ni}$ is the embedding of sequence $n$ at position $i$, and $\otimes$ denotes outer product. This summary is added to the pair representation $z_{ij}$ before it enters the Pairformer. Crucially, the full MSA matrix $M$ is not carried forward into the trunk, substantially reducing memory.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure3-pairformer.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 2: Pairformer architecture.</strong> The Pairformer replaces AlphaFold 2's Evoformer by operating only on single ($s$) and pair ($z$) representations without maintaining a full MSA track. Each Pairformer block updates the pair representation through triangle multiplicative updates and triangle self-attention (operations borrowed from the Evoformer's pair stack), then updates the single representation via attention biased by the pair representation. MSA information has already been summarized into the pair representation before the Pairformer begins processing. This simplification reduces memory requirements from $O(N_{\text{seq}} \times L \times d)$ to $O(L^2 \times d_z + L \times d_s)$ while preserving the model's ability to capture long-range dependencies and coevolutionary signals.
</div>

### Pairformer (Simplified Evoformer)

The Pairformer operates on the single representation $s \in \mathbb{R}^{L \times d_s}$ and pair representation $z \in \mathbb{R}^{L \times L \times d_z}$. At each layer $\ell$, updates proceed as:

$$
\begin{aligned}
z^{(\ell+1)} &= z^{(\ell)} + \text{TriangleMult}(z^{(\ell)}) + \text{TriangleAttn}(z^{(\ell)}) \\
s^{(\ell+1)} &= s^{(\ell)} + \text{Attention}(s^{(\ell)}, \text{bias}=z^{(\ell+1)})
\end{aligned}
$$

The **Triangle Multiplicative Update** enforces geometric consistency by updating edge $(i,j)$ using information from triangles $(i,k)$ and $(k,j)$:

$$
z_{ij} \leftarrow z_{ij} + \sum_k g(z_{ik}) \cdot h(z_{kj})
$$

where $g$ and $h$ are learned projections. This operation ensures that the pair representation captures transitive spatial relationships: if residue $i$ is close to $k$ and $k$ is close to $j$, the model can infer that $i$ and $j$ are likely nearby.

The **single representation update** uses attention where the pair features serve as additive biases to the attention logits, coupling the single-residue features with pairwise context.

### Diffusion Module

The diffusion module is AF3's central innovation. It operates on the full set of atom coordinates $\mathbf{x} = (x_1, \ldots, x_N) \in \mathbb{R}^{N \times 3}$, where $N$ is the total number of atoms across all chains.

**Forward (noising) process.** Starting from the true atom coordinates $\mathbf{x}_0$, Gaussian noise is progressively added according to a variance schedule $\{\sigma_t\}_{t=1}^T$:

$$
\mathbf{x}_t = \mathbf{x}_0 + \sigma_t \boldsymbol{\epsilon}, \quad \boldsymbol{\epsilon} \sim \mathcal{N}(\mathbf{0}, \mathbf{I})
$$

where $\sigma_t$ increases with $t$, so that $\mathbf{x}_T$ is approximately pure Gaussian noise.

**Reverse (denoising) process.** The model learns a denoising network $\hat{\mathbf{x}}_0 = D_\theta(\mathbf{x}_t, t, c)$ that predicts the clean coordinates $\mathbf{x}_0$ given noisy coordinates $\mathbf{x}_t$, the noise level $t$, and conditioning information $c$ (the single and pair representations from the Pairformer trunk). At inference, the model starts from random noise $\mathbf{x}_T \sim \mathcal{N}(\mathbf{0}, \sigma_T^2 \mathbf{I})$ and iteratively denoises:

$$
\mathbf{x}_{t-1} = D_\theta(\mathbf{x}_t, t, c) + \sigma_{t-1} \boldsymbol{\epsilon}', \quad \boldsymbol{\epsilon}' \sim \mathcal{N}(\mathbf{0}, \mathbf{I})
$$

This iterative refinement progressively resolves the structure from a coarse blob into a detailed atomic arrangement.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure2-diffusion.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 3: Diffusion module for structure generation.</strong> <strong>Left:</strong> During training, ground-truth atom coordinates $\mathbf{x}_0$ are corrupted with Gaussian noise at a randomly sampled noise level $\sigma_t$ to produce noisy coordinates $\mathbf{x}_t$. The denoising network takes $\mathbf{x}_t$, the noise level $t$, and conditioning features from the Pairformer trunk (single and pair representations) as input, and predicts the clean coordinates $\hat{\mathbf{x}}_0$. The training loss penalizes the squared error between $\hat{\mathbf{x}}_0$ and $\mathbf{x}_0$. <strong>Right:</strong> At inference, the model starts from pure noise $\mathbf{x}_T \sim \mathcal{N}(\mathbf{0}, \sigma_T^2 \mathbf{I})$ and applies 200 denoising steps, progressively refining atom positions from a formless cloud into a precise molecular structure. Multiple independent samples can be generated by starting from different noise realizations, enabling ensemble-based confidence estimation.
</div>

**Training loss.** The diffusion module is trained with a weighted mean squared error (MSE) loss between predicted and true coordinates:

$$
\mathcal{L}_{\text{diffusion}} = \mathbb{E}_{t, \boldsymbol{\epsilon}} \left[ w(t) \| D_\theta(\mathbf{x}_t, t, c) - \mathbf{x}_0 \|^2 \right]
$$

where $w(t)$ is a weighting function that emphasizes certain noise levels (the paper reports using a weighting that focuses on intermediate noise levels where the model must resolve both global arrangement and local geometry). The expectation is over uniformly sampled time steps $t$ and noise realizations $\boldsymbol{\epsilon}$.

The total training objective combines the diffusion loss with auxiliary losses:

$$
\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{diffusion}} + \lambda_{\text{conf}} \mathcal{L}_{\text{confidence}} + \lambda_{\text{dist}} \mathcal{L}_{\text{distogram}}
$$

where $\mathcal{L}_{\text{confidence}}$ trains the confidence heads (pLDDT, PAE) and $\mathcal{L}_{\text{distogram}}$ is a pairwise distance prediction loss that provides additional structural supervision.

### Confidence Metrics

AF3 produces several confidence estimates:

- **pLDDT (predicted Local Distance Difference Test):** Per-atom confidence in local structure quality, ranging from 0 to 100. Computed by predicting a distribution over LDDT score bins.

- **PAE (Predicted Aligned Error):** For each pair of tokens $(i, j)$, the predicted error in the position of token $j$ when the predicted and true structures are aligned on token $i$:

$$
\text{PAE}(i,j) = \mathbb{E}\left[\| \mathbf{x}_j^{\text{pred}} - \mathbf{x}_j^{\text{true}} \|_2 \;\Big|\; \text{align on } i \right]
$$

- **iPTM (interface Predicted TM-score):** A predicted TM-score computed only over inter-chain residue pairs, assessing the quality of predicted interfaces:

$$
\text{iPTM} = \frac{1}{|\mathcal{P}_{\text{inter}}|} \sum_{(i,j) \in \mathcal{P}_{\text{inter}}} \frac{1}{1 + \left(\frac{\text{PAE}(i,j)}{d_0}\right)^2}
$$

where $\mathcal{P}_{\text{inter}}$ is the set of inter-chain token pairs and $d_0$ is a length-dependent normalization factor.

- **Ranking score:** A composite metric used to select the best prediction from multiple diffusion samples:

$$
\text{ranking} = 0.2 \cdot \text{pTM} + 0.8 \cdot \text{iPTM} + \text{disorder penalty}
$$

where the disorder penalty down-weights predictions of known disordered regions to avoid inflating confidence scores.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure7-confidence.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 4: Confidence metrics in AlphaFold 3.</strong> Visualization of AF3's confidence outputs for an example biomolecular complex. <strong>Left:</strong> Per-residue pLDDT scores mapped onto the predicted structure, with blue indicating high confidence (>90) and red indicating low confidence (<50). Well-ordered regions such as protein cores and base-paired nucleic acid regions typically show high pLDDT, while flexible loops and disordered tails show low pLDDT. <strong>Center:</strong> PAE matrix showing the predicted aligned error between all pairs of residues. Dark blue regions along the diagonal indicate high confidence in intradomain structure; off-diagonal blue blocks indicate confidently predicted interdomain or interchain contacts. <strong>Right:</strong> iPTM scores for different interfaces within the complex, providing a per-interface quality assessment that is critical for evaluating the reliability of predicted binding modes.
</div>

---

## Real-World Impact

AF3's ability to model diverse biomolecular interactions opens new avenues across multiple fields:

### Drug Discovery: Protein-Ligand Structure Prediction

Traditional drug design relies heavily on experimental protein-ligand co-crystal structures or physics-based docking tools like AutoDock Vina and Glide. AF3 substantially outperforms these tools on protein-ligand structure prediction benchmarks, as the paper reports. For drug-like small molecules, AF3 achieves higher success rates (fraction of predictions below 2 Angstrom RMSD from the experimental pose) compared to the best traditional docking methods and even specialized deep learning approaches.

This capability enables:

- **Virtual screening** of drug candidates against predicted binding poses
- **Lead optimization** by predicting how chemical modifications affect binding geometry
- **Polypharmacology analysis** by modeling a drug's interactions with multiple protein targets

### Nucleic Acid Interactions

AF3 is the first general-purpose model to accurately predict protein-DNA and protein-RNA complex structures. This is critical for understanding:

- **Transcription factor-DNA binding:** How regulatory proteins recognize specific DNA sequences
- **CRISPR-Cas complexes:** Structural basis of genome editing machinery
- **Ribosome and spliceosome assemblies:** RNA-protein machines essential for gene expression
- **RNA aptamer design:** Engineering RNA molecules that bind specific targets

The paper reports that AF3 substantially outperforms existing protein-nucleic acid docking tools on established benchmarks.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure6-nucleic-acids.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 5: Nucleic acid complex predictions.</strong> Examples of AF3 predictions for protein-nucleic acid complexes. <strong>Left:</strong> A protein-DNA complex showing the predicted structure (colored) overlaid with the experimental structure (gray), demonstrating accurate prediction of both the protein fold and DNA conformation, as well as the protein-DNA interface geometry. <strong>Center:</strong> A protein-RNA complex illustrating AF3's ability to capture RNA secondary structure elements (stems, loops, bulges) and their positioning relative to the protein binding surface. <strong>Right:</strong> An RNA-only structure prediction, showing that AF3 can predict RNA tertiary folds including pseudoknots and long-range base-pairing interactions, a task that was previously very challenging for computational methods.
</div>

### Post-Translational Modifications and Covalent Ligands

AF3 natively handles covalent modifications to proteins, including:

- **Glycosylation:** Sugar chains attached to asparagine or serine/threonine residues
- **Phosphorylation:** Addition of phosphate groups that regulate protein activity
- **Bonded ligands:** Covalent inhibitors, cofactors (heme, FAD, NAD+), and prosthetic groups

This is a significant advance over AF2, which treated proteins as unmodified polypeptide chains.

### Protein-Protein Interactions

While AlphaFold-Multimer already predicted protein complex structures, AF3 improves accuracy on antibody-antigen interfaces, a particularly challenging category due to the high sequence variability of antibody complementarity-determining regions (CDRs). The paper reports that AF3 achieves improved performance on recent antibody-antigen targets compared to AlphaFold-Multimer.

---

## Benchmarks and Performance

AF3 was evaluated on a range of biomolecular interaction benchmarks. All evaluations used targets with deposition dates after the training data cutoff to prevent data leakage.

### Protein-Ligand Binding

On the PoseBusters benchmark (a curated set of recent protein-ligand co-crystal structures designed to test physically valid predictions), the paper reports that AF3 achieves substantially higher success rates than traditional docking methods and the specialized deep learning tool DiffDock:

| Method          | Category                                   |
| :-------------- | :----------------------------------------- |
| **AlphaFold 3** | Generalist (deep learning)                 |
| DiffDock        | Specialist (deep learning, ligand docking) |
| Vina (AutoDock) | Traditional (physics-based docking)        |
| Gold            | Traditional (physics-based docking)        |

The paper reports that AF3 achieves the highest fraction of successful ligand poses (RMSD < 2 Angstrom) on PoseBusters while also passing physical validity checks (no steric clashes, correct bond geometry). Notably, many deep learning docking methods that achieve high RMSD success rates fail physical validity checks, a problem that AF3 largely avoids.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure5-ligands.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 6: Protein-ligand prediction performance.</strong> Comparison of AF3 with specialized docking methods on the PoseBusters benchmark. <strong>Left:</strong> Fraction of predictions with ligand RMSD below 2 Angstrom (success rate), showing AF3 achieving the highest success rate among all methods. <strong>Center:</strong> Example predictions showing AF3-predicted ligand poses (colored) overlaid on experimental poses (gray) for diverse drug-like molecules, including cases with multiple rotatable bonds, ring systems, and charged groups. <strong>Right:</strong> Physical validity assessment showing that AF3 predictions pass stereochemistry and clash checks at higher rates than competing deep learning methods, which often generate physically implausible poses despite low RMSD values.
</div>

### Protein-Nucleic Acid Complexes

On recent protein-DNA and protein-RNA complex structures, the paper reports that AF3 substantially outperforms existing methods. The evaluation uses interface LDDT (iLDDT) and DockQ scores to assess prediction quality at the binding interface:

| Interaction Type | AF3 Performance                                |
| :--------------- | :--------------------------------------------- |
| **Protein-DNA**  | Substantially improved over prior methods      |
| **Protein-RNA**  | Substantially improved over prior methods      |
| **RNA-only**     | Competitive with RNA-specific prediction tools |

### Protein-Protein Interfaces

On recent PDB structures (post training cutoff), the paper reports AF3 performance across different protein-protein interaction categories:

| Category                    | AF3 vs. AF2-Multimer   |
| :-------------------------- | :--------------------- |
| **General protein-protein** | Comparable or improved |
| **Antibody-antigen**        | Notably improved       |

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/alphafold3/figure4-benchmarks.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    <strong>Figure 7: Comprehensive benchmark comparisons.</strong> Performance of AF3 across all biomolecular interaction categories compared to the best existing methods for each category. <strong>Top row:</strong> Protein structure prediction (single chain), where AF3 maintains AF2-level accuracy. <strong>Middle row:</strong> Protein-protein, protein-DNA, and protein-RNA interfaces, where AF3 shows substantial improvements over both AF2-Multimer and specialized docking tools. <strong>Bottom row:</strong> Protein-ligand and protein-ion predictions, where AF3 outperforms both physics-based and deep learning docking methods. Each bar represents the fraction of targets achieving a quality threshold (e.g., DockQ > 0.23 for acceptable quality, RMSD < 2 Angstrom for ligands). Error bars indicate 95% confidence intervals from bootstrapping.
</div>

### Ablation Studies

The paper presents ablation experiments demonstrating the importance of key design choices:

- **Diffusion module:** Removing the diffusion module and reverting to a deterministic structure prediction approach degrades performance across all interaction types, with the largest drops on protein-ligand and protein-nucleic acid targets.
- **Cross-distillation:** Removing AF2-distilled training data reduces accuracy, particularly for targets with limited experimental homologs.
- **Pairformer vs. Evoformer:** The simplified Pairformer achieves comparable performance to the full Evoformer while being more computationally efficient.
- **MSA processing:** Removing MSA information entirely substantially degrades performance, confirming that evolutionary covariance remains essential even in the diffusion framework.

---

## Limitations and Future Directions

Despite its broad capabilities, AF3 has several important limitations:

### 1. Hallucination of Plausible but Incorrect Structures

The diffusion module can generate structures that look physically reasonable but are incorrect. Unlike AF2's deterministic predictions, which tended to produce obviously wrong structures when uncertain (e.g., unfolded chains), AF3's diffusion samples can produce confident-looking but wrong conformations. The paper acknowledges this as a key concern and introduces confidence metrics to help identify unreliable predictions, but the risk of <span style="background-color: #fff3b0;">hallucinated structures</span> remains.

### 2. Stereochemistry and Chirality Issues

Diffusion models operating on raw Cartesian coordinates do not inherently enforce chemical validity. While AF3 includes post-processing steps to correct bond lengths and chirality violations, some samples may still contain:

- Incorrect stereochemistry at chiral centers
- Strained bond angles
- Physically implausible conformations for small molecules

The paper notes that generating multiple samples and filtering with confidence metrics mitigates but does not eliminate these issues.

### 3. Limited Conformational Diversity

While the diffusion module can in principle generate diverse conformations, the paper reports that AF3 tends to produce similar structures across samples for well-determined targets. For targets with genuine conformational heterogeneity (e.g., allosteric proteins, intrinsically disordered regions), the diversity of AF3 samples may not fully capture the structural ensemble.

### 4. Training Data Biases

AF3 is trained primarily on crystallographic structures from the PDB, which are biased toward:

- Well-folded, stable proteins (underrepresenting disordered and membrane proteins)
- Specific organism biases (overrepresenting human and model organisms)
- Crystallization-compatible conformations (excluding transient states)

### 5. Computational Cost

AF3 requires substantial computational resources, particularly for large complexes. The diffusion module requires multiple denoising steps (the paper reports using 200 steps), and generating multiple samples for confidence estimation multiplies the cost further.

### Future Extensions

- **Conformational ensembles:** Training diffusion models to explicitly capture conformational distributions rather than single structures
- **Dynamics prediction:** Extending from static structures to molecular dynamics trajectories
- **Covalent drug design:** Leveraging AF3's ability to model covalent modifications for targeted covalent inhibitor design
- **Integration with experimental data:** Using cryo-EM density maps or crosslinking mass spectrometry data as additional conditioning inputs
- **Protein design:** Inverting the structure prediction model to design sequences that fold into desired structures with specific binding properties

---

## Key Takeaways

- **AlphaFold 3 unifies biomolecular structure prediction** across proteins, nucleic acids, small molecules, ions, and covalent modifications within a single architecture, eliminating the need for separate specialized tools.

- **The diffusion-based structure module is the key innovation**, replacing AF2's deterministic structure prediction with a generative approach that learns to denoise atomic coordinates from random noise, enabling stochastic sampling of structural predictions.

- **The Pairformer simplifies the Evoformer** by separating MSA processing from the main trunk, operating only on single and pair representations. This reduces computational cost without sacrificing accuracy.

- **AF3 substantially outperforms specialized methods** on protein-ligand, protein-nucleic acid, and protein-protein benchmarks, demonstrating that a single generalist model can surpass dedicated tools across diverse interaction types.

- **Confidence metrics (pLDDT, PAE, iPTM) are critical for practical use**, as the diffusion module can hallucinate plausible but incorrect structures. The ranking score enables selection of the best prediction from multiple samples.

- **Limitations include hallucination risk, stereochemistry issues, and training data biases**, highlighting that AF3 predictions should be validated experimentally, particularly for novel interaction types not well represented in the PDB.

- **Cross-distillation from AF2 predictions** expands the effective training set beyond experimental structures, and the unified token-atom representation enables seamless handling of diverse molecular entities.

---

_Reference: Abramson et al.,_ **Nature 630**, 493--500 (2024). DOI: [10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w)

_Access predictions: [AlphaFold Server](https://alphafoldserver.com/)_
