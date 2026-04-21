---
layout: post
title: "AlphaFold + ProteinMPNN: A Practical Protein Design Workflow"
date: 2025-12-05
permalink: /blog/alphafold-proteinmpnn-workflow/
published: true
categories: [tutorial]
tags:
  - AI
  - protein-design
  - structural-biology
  - tutorial
  - immunology
description: "An end-to-end antibody-design workflow combining AlphaFold / ColabFold structure prediction with ProteinMPNN inverse folding for CDR redesign, with notes on Chothia numbering and TCR-pMHC engineering."
---

Computational protein design has entered a new era. With **AlphaFold2** providing near-experimental accuracy in structure prediction and **ProteinMPNN** enabling rapid inverse folding (sequence design from structure), researchers now have a powerful closed-loop workflow for engineering novel proteins. This tutorial walks through the complete pipeline---from structure prediction to sequence design to validation---with a focus on immunology applications such as TCR and antibody engineering.

---

## 1. The Design Loop

The core idea behind modern computational protein design is a three-step cycle:

1. **Structure Prediction**: Given a sequence, predict its 3D structure using AlphaFold2.
2. **Sequence Design**: Given a target backbone structure, design new sequences that fold into it using ProteinMPNN.
3. **Validation**: Predict the structure of the designed sequence with AlphaFold2 and check whether it matches the target (self-consistent RMSD).

This loop enables iterative refinement: generate candidates, filter by structural confidence, and converge on designs that are both novel and likely to fold correctly.

### Why This Matters for Immunology

- **Antibody optimization**: Redesign CDR loops to improve affinity or developability while preserving the framework scaffold.
- **TCR engineering**: Design TCR variants with altered peptide-MHC specificity for adoptive cell therapy.
- **De novo binder design**: Create proteins that bind specific epitopes on viral or tumor antigens.

### What We Will Build

By the end of this tutorial, you will have a working pipeline that:

- Predicts a protein complex structure with ColabFold
- Prepares design inputs specifying which residues to redesign
- Generates 100+ sequence variants with ProteinMPNN
- Validates designs with AlphaFold2 using self-consistent RMSD filtering

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure1-design-loop.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 1. The computational protein design loop: Structure Prediction → Sequence Design → Validation. Designs passing the self-consistency filter are candidates for experimental testing.
</div>

---

## 2. Environment Setup

This workflow requires GPU access for AlphaFold2 predictions. ProteinMPNN runs efficiently on CPU but benefits from GPU acceleration. We recommend a Linux machine with at least one NVIDIA GPU (16 GB+ VRAM) or a cloud instance (e.g., Google Colab Pro, AWS p3).

### 2.1 Install ColabFold (LocalColabFold)

[ColabFold](https://github.com/sokrypton/ColabFold) provides a fast, MMseqs2-based frontend for AlphaFold2 that eliminates the need for large genetic databases.

```bash
# Create a dedicated conda environment
conda create -n protein-design python=3.10 -y
conda activate protein-design

# Install ColabFold with AlphaFold2 backend
pip install "colabfold[alphafold]"

# Verify installation
colabfold_batch --help
```

> **Note**: The first run will download AlphaFold2 model parameters (~3.5 GB). MSA generation uses the ColabFold MMseqs2 server by default, requiring internet access.

### 2.2 Install ProteinMPNN

```bash
# Clone the official repository
git clone https://github.com/dauparas/ProteinMPNN.git
cd ProteinMPNN

# Install dependencies
pip install numpy torch

# Verify by listing available scripts
ls *.py
# protein_mpnn_run.py  protein_mpnn_utils.py  ...
```

### 2.3 Install Visualization and Analysis Tools

```bash
pip install py3Dmol biopython matplotlib mdanalysis
```

We will use **py3Dmol** for inline 3D visualization in Jupyter notebooks and **BioPython** for PDB parsing and RMSD calculations.

### 2.4 Directory Structure

Set up a clean project layout:

```bash
mkdir -p protein_design/{inputs,af2_predictions,mpnn_outputs,validation}
```

```
protein_design/
├── inputs/              # Input sequences and structures
├── af2_predictions/     # AlphaFold2 prediction outputs
├── mpnn_outputs/        # ProteinMPNN designed sequences
└── validation/          # AF2 predictions of designed sequences
```

---

## 3. Step 1: Structure Prediction with ColabFold

### 3.1 Prepare Input Sequences

Create a FASTA file with the target sequence. For a multi-chain complex (e.g., antibody heavy + light chain), separate chains with a colon (`:`).

```bash
# Note: the sequences below cover only the VH and VL (variable) domains, which
# together form an Fv / scFv. A true Fab also requires the CH1 (heavy) and
# CL (light) constant domains — we omit them here to keep the tutorial fast.
cat > protein_design/inputs/antibody.fasta << 'EOF'
>anti_HER2_Fv
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS:DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
EOF
```

### 3.2 Run ColabFold Batch Prediction

```bash
colabfold_batch \
  protein_design/inputs/antibody.fasta \
  protein_design/af2_predictions/ \
  --num-recycle 3 \
  --num-models 5 \
  --amber \
  --use-gpu-relax
```

**Key parameters:**

- `--num-recycle 3`: Number of recycling iterations (improves accuracy for complexes).
- `--num-models 5`: Run all five AF2 model variants and rank them.
- `--amber`: Apply AMBER force field relaxation to reduce steric clashes.
- `--use-gpu-relax`: Use GPU for the relaxation step (faster).

> **GPU requirement**: A single Fv prediction (VH + VL) with 5 models takes approximately 10--20 minutes on an A100. On a T4, expect 30--60 minutes.

### 3.3 Interpreting Outputs

ColabFold produces several output files per prediction:

```bash
ls protein_design/af2_predictions/
# anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_*.pdb
# anti_HER2_Fab_scores_rank_001_*.json
# anti_HER2_Fab_coverage.png
# anti_HER2_Fab_pae.png
# anti_HER2_Fab_plddt.png
```

**Key quality metrics:**

- **pLDDT** (predicted Local Distance Difference Test): Per-residue confidence score (0--100). Regions with pLDDT > 90 are modeled with high confidence; 70--90 is acceptable; below 70 suggests disorder or poor prediction.
- **PAE** (Predicted Aligned Error): A matrix of predicted positional errors between all residue pairs. Low inter-chain PAE indicates confident prediction of the complex interface.

### 3.4 Visualize the Prediction

```python
import py3Dmol
from pathlib import Path

def show_structure(pdb_path, color_by="pLDDT"):
    """Visualize a PDB structure colored by pLDDT or chain."""
    pdb_text = Path(pdb_path).read_text()
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_text, "pdb")

    if color_by == "pLDDT":
        view.setStyle({
            "cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "roygb",
                    "min": 50,
                    "max": 100,
                }
            }
        })
    elif color_by == "chain":
        view.setStyle({"chain": "A"}, {"cartoon": {"color": "#1f77b4"}})
        view.setStyle({"chain": "B"}, {"cartoon": {"color": "#ff7f0e"}})

    view.zoomTo()
    return view.show()

# Visualize the top-ranked prediction
show_structure("protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb")
```

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure2-plddt-structure.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 2. AlphaFold2 prediction of an anti-HER2 Fab colored by pLDDT. Blue regions (pLDDT > 90) are high confidence; CDR loops often show moderate confidence (yellow/orange).
</div>

### 3.5 Extract pLDDT Scores

```python
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt

def plot_plddt(pdb_path, title="pLDDT per residue"):
    """Plot per-residue pLDDT from the B-factor column."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pred", pdb_path)

    plddts = []
    residue_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # skip heteroatoms
                    ca = residue["CA"] if "CA" in residue else list(residue.get_atoms())[0]
                    plddts.append(ca.get_bfactor())
                    residue_ids.append(f"{chain.id}{residue.id[1]}")

    plddts = np.array(plddts)
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(range(len(plddts)), plddts, width=1.0,
           color=["#1a9850" if v > 90 else "#fee08b" if v > 70 else "#d73027" for v in plddts])
    ax.set_xlabel("Residue index")
    ax.set_ylabel("pLDDT")
    ax.set_title(title)
    ax.set_ylim(0, 100)
    ax.axhline(y=90, color="gray", linestyle="--", alpha=0.5, label="90")
    ax.axhline(y=70, color="gray", linestyle=":", alpha=0.5, label="70")
    ax.legend()
    plt.tight_layout()
    plt.savefig("plddt_plot.png", dpi=150)
    plt.show()
    print(f"Mean pLDDT: {plddts.mean():.1f}")

plot_plddt("protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb")
```

---

## 4. Step 2: Preparing Structures for Design

Before running ProteinMPNN, we need to define exactly which residues should be redesigned and which should remain fixed.

### 4.1 Parse the PDB and Identify Chains

```python
from Bio.PDB import PDBParser, PDBIO, Select

def get_chain_residues(pdb_path):
    """List all chains and their residue ranges."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    for model in structure:
        for chain in model:
            residues = [r for r in chain if r.id[0] == " "]
            if residues:
                start = residues[0].id[1]
                end = residues[-1].id[1]
                print(f"Chain {chain.id}: residues {start}-{end} ({len(residues)} residues)")

get_chain_residues("protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb")
# Chain A: residues 1-120 (120 residues)  -- Heavy chain
# Chain B: residues 1-107 (107 residues)  -- Light chain
```

### 4.2 Define Designable vs. Fixed Positions

For antibody design, we typically fix the framework regions and redesign the CDR loops. Here we define CDR positions using **Chothia numbering**.

> ⚠️ **Important caveat**: AF2 writes residues out using raw, sequential PDB numbering (1, 2, 3, ...) — **not** Chothia numbering. Before mapping any CDR range onto the predicted structure you must renumber the heavy and light chains with a dedicated tool such as [ANARCI](https://github.com/oxpig/ANARCI) or [AbNumber](https://github.com/prihoda/AbNumber). If you skip this step and apply the Chothia indices below to the raw AF2 PDB directly, you will freeze/design the **wrong** residues. The code below assumes `pdb_path` has already been renumbered into Chothia scheme (e.g., via `anarci -i input.pdb --scheme chothia -o renumbered.pdb`).

```python
import json

# CDR definitions (Chothia numbering; assumes the PDB has been renumbered)
CDR_DEFINITIONS = {
    "H1": (26, 32),   # Heavy chain CDR1
    "H2": (52, 56),   # Heavy chain CDR2
    "H3": (95, 102),  # Heavy chain CDR3
    "L1": (24, 34),   # Light chain CDR1
    "L2": (50, 56),   # Light chain CDR2
    "L3": (89, 97),   # Light chain CDR3
}

def make_fixed_positions(pdb_path, chain_cdr_map, output_path):
    """
    Create a fixed positions JSON for ProteinMPNN.
    Everything is fixed EXCEPT the specified CDR positions.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    fixed_positions = {}
    pdb_name = pdb_path.split("/")[-1].replace(".pdb", "")

    for model in structure:
        chain_fixed = {}
        for chain in model:
            chain_id = chain.id
            residues = [r.id[1] for r in chain if r.id[0] == " "]

            # Determine which CDRs belong to this chain
            designable = set()
            for cdr_name, (start, end) in chain_cdr_map.get(chain_id, {}).items():
                designable.update(range(start, end + 1))

            # Fixed = all residues NOT in designable set
            fixed = [r for r in residues if r not in designable]
            if fixed:
                chain_fixed[chain_id] = fixed

        fixed_positions[pdb_name] = chain_fixed

    with open(output_path, "w") as f:
        json.dump(fixed_positions, f)

    return fixed_positions

# Heavy chain (A): redesign H1, H2, H3
# Light chain (B): redesign L1, L2, L3
chain_cdr_map = {
    "A": {"H1": (26, 32), "H2": (52, 56), "H3": (95, 102)},
    "B": {"L1": (24, 34), "L2": (50, 56), "L3": (89, 97)},
}

fixed_pos = make_fixed_positions(
    "protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb",
    chain_cdr_map,
    "protein_design/inputs/fixed_positions.jsonl"
)
```

### 4.3 Build ProteinMPNN's JSONL Inputs with the Helper Scripts

ProteinMPNN's JSONL formats are finicky — the safest way to generate them is to use the helper scripts bundled with the repository. We first parse all predicted PDBs into a `parsed_pdbs.jsonl`, then layer on the "which chain is designed" and "which residues are fixed" assignments:

```bash
cd ProteinMPNN/helper_scripts

# (a) Parse every PDB in the input folder into ProteinMPNN's JSONL format.
#     Output contains backbone coords + sequence + chain info for each structure.
python parse_multiple_chains.py \
  --input_path ../../protein_design/af2_predictions/ \
  --output_path ../../protein_design/inputs/parsed_pdbs.jsonl

# (b) Tell ProteinMPNN which chains are designed and which are fixed.
#     Here both VH (A) and VL (B) are designed; if you had an antigen chain
#     (e.g. "C"), you would pass --chain_list "A B" and it would stay fixed.
python assign_fixed_chains.py \
  --input_path ../../protein_design/inputs/parsed_pdbs.jsonl \
  --output_path ../../protein_design/inputs/assigned_pdbs.jsonl \
  --chain_list "A B"

# (c) Lock every framework position and leave only the CDR ranges designable.
#     `--position_list` / `--chain_list` use ProteinMPNN's own conventions —
#     pass the CDR residue indices you generated in section 4.2 here.
python make_fixed_positions_dict.py \
  --input_path ../../protein_design/inputs/parsed_pdbs.jsonl \
  --output_path ../../protein_design/inputs/fixed_pdbs.jsonl \
  --chain_list "A B" \
  --position_list "$CDR_H_POSITIONS $CDR_L_POSITIONS" \
  --specify_non_fixed
```

The three files `parsed_pdbs.jsonl`, `assigned_pdbs.jsonl`, `fixed_pdbs.jsonl` are the canonical inputs expected by `protein_mpnn_run.py` in section 5. Hand-crafted JSONL files (as in earlier versions of this tutorial) are easy to get subtly wrong and can silently misalign CDR positions — **always prefer the helper scripts**.

---

## 5. Step 3: Sequence Design with ProteinMPNN

### 5.1 Running ProteinMPNN

ProteinMPNN performs **inverse folding**: given a fixed backbone structure, it predicts amino acid sequences likely to fold into that structure.

```bash
cd ProteinMPNN

python protein_mpnn_run.py \
  --jsonl_path ../protein_design/inputs/parsed_pdbs.jsonl \
  --chain_id_jsonl ../protein_design/inputs/assigned_pdbs.jsonl \
  --fixed_positions_jsonl ../protein_design/inputs/fixed_pdbs.jsonl \
  --out_folder ../protein_design/mpnn_outputs/ \
  --num_seq_per_target 100 \
  --sampling_temp "0.1 0.2 0.3" \
  --seed 42 \
  --batch_size 1
```

> The trio `--jsonl_path parsed_pdbs.jsonl`, `--chain_id_jsonl assigned_pdbs.jsonl`, `--fixed_positions_jsonl fixed_pdbs.jsonl` is the combination the official [ProteinMPNN examples](https://github.com/dauparas/ProteinMPNN/tree/main/examples) use. Mixing `--pdb_path` (single-PDB shortcut) with helper-script JSONLs will error out because the formats don't match.

**Key parameters explained:**

| Parameter                 | Description                                                                   | Recommended                       |
| ------------------------- | ----------------------------------------------------------------------------- | --------------------------------- |
| `--num_seq_per_target`    | Number of sequences to generate per temperature                               | 100 for screening                 |
| `--sampling_temp`         | Temperature(s) for sampling. Lower = more conservative, higher = more diverse | `"0.1 0.2 0.3"`                   |
| `--fixed_positions_jsonl` | Positions to keep fixed (not redesigned)                                      | Required for CDR design           |
| `--batch_size`            | Batch size for inference                                                      | 1 (increase if GPU memory allows) |
| `--seed`                  | Random seed for reproducibility                                               | Any integer                       |

> **CPU vs GPU**: ProteinMPNN runs well on CPU for single structures. For batch design of hundreds of targets, GPU is recommended. On CPU, 100 sequences take approximately 30 seconds per target.

### 5.2 Understanding the Temperature Parameter

Temperature controls the diversity-conservation trade-off in sequence sampling:

- **T = 0.1**: Very conservative. Sequences closely resemble the most probable amino acid at each position. High sequence recovery but low diversity.
- **T = 0.2**: Moderate. Good balance between novelty and foldability. **Recommended starting point.**
- **T = 0.3--0.5**: More exploratory. Greater sequence diversity but some designs may not fold correctly.
- **T > 0.5**: Generally not recommended---too much randomness leads to poorly folding sequences.

### 5.3 Interpreting ProteinMPNN Output

```bash
ls ../protein_design/mpnn_outputs/seqs/
# anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.fa
```

The output is a FASTA file with designed sequences:

```python
def parse_mpnn_output(fasta_path):
    """Parse ProteinMPNN output FASTA and extract scores."""
    sequences = []
    with open(fasta_path) as f:
        header = None
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq:
                    sequences.append({"header": header, "sequence": seq})
                header = line
                seq = ""
            else:
                seq += line
        if header and seq:
            sequences.append({"header": header, "sequence": seq})

    for entry in sequences[:5]:
        # Headers contain: sample number, T (temperature),
        # sample (index), score, global_score, seq_recovery
        print(entry["header"])
        print(f"  Length: {len(entry['sequence'])}")
        print()

    return sequences

seqs = parse_mpnn_output(
    "protein_design/mpnn_outputs/seqs/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.fa"
)
```

Each header line contains:

- **score**: Negative log-likelihood of the sequence given the structure (lower is better).
- **global_score**: Average score across all designed positions.
- **seq_recovery**: Fraction of positions matching the original sequence.

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure3-mpnn-scores.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 3. ProteinMPNN output analysis. Left: distribution of global scores across 300 designs at three temperatures. Right: sequence recovery vs. global score, showing the expected trade-off.
</div>

### 5.4 Analyzing Designed Sequences

```python
import re
import matplotlib.pyplot as plt
import numpy as np

def analyze_mpnn_designs(fasta_path):
    """Analyze score distributions from ProteinMPNN output."""
    scores = []
    recoveries = []
    temps = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">") and "score" in line:
                # Parse score and recovery from header
                score_match = re.search(r"score=([0-9.]+)", line)
                recovery_match = re.search(r"seq_recovery=([0-9.]+)", line)
                temp_match = re.search(r"T=([0-9.]+)", line)

                if score_match:
                    scores.append(float(score_match.group(1)))
                if recovery_match:
                    recoveries.append(float(recovery_match.group(1)))
                if temp_match:
                    temps.append(float(temp_match.group(1)))

    scores = np.array(scores)
    recoveries = np.array(recoveries)
    temps = np.array(temps)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Score distribution by temperature
    for t in sorted(set(temps)):
        mask = temps == t
        axes[0].hist(scores[mask], bins=20, alpha=0.6, label=f"T={t}")
    axes[0].set_xlabel("Global Score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Score Distribution by Temperature")
    axes[0].legend()

    # Recovery vs score
    scatter = axes[1].scatter(scores, recoveries, c=temps, cmap="viridis", alpha=0.6, s=20)
    axes[1].set_xlabel("Global Score")
    axes[1].set_ylabel("Sequence Recovery")
    axes[1].set_title("Recovery vs Score")
    plt.colorbar(scatter, ax=axes[1], label="Temperature")

    plt.tight_layout()
    plt.savefig("mpnn_analysis.png", dpi=150)
    plt.show()

    print(f"Total designs: {len(scores)}")
    print(f"Score range: {scores.min():.3f} - {scores.max():.3f}")
    print(f"Mean recovery: {recoveries.mean():.3f}")

analyze_mpnn_designs(
    "protein_design/mpnn_outputs/seqs/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.fa"
)
```

---

## 6. Step 4: Validation with AlphaFold2

The critical validation step: predict the structure of each designed sequence with AlphaFold2 and check if it matches the original target backbone.

### 6.1 Prepare Designed Sequences for AF2

```python
from pathlib import Path

def prepare_validation_fastas(mpnn_fasta, output_dir, top_n=20):
    """
    Extract top N designs by score and write individual FASTAs
    for ColabFold validation.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse all sequences with scores
    designs = []
    with open(mpnn_fasta) as f:
        header = None
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq:
                    score_match = re.search(r"score=([0-9.]+)", header)
                    score = float(score_match.group(1)) if score_match else 999
                    designs.append({"header": header, "sequence": seq, "score": score})
                header = line
                seq = ""
            else:
                seq += line
        if header and seq:
            score_match = re.search(r"score=([0-9.]+)", header)
            score = float(score_match.group(1)) if score_match else 999
            designs.append({"header": header, "sequence": seq, "score": score})

    # Skip the first entry (original sequence) and sort by score
    designs = sorted(designs[1:], key=lambda x: x["score"])[:top_n]

    # Write individual FASTAs for ColabFold
    # Multi-chain sequences use ':' separator
    for i, design in enumerate(designs):
        fasta_path = output_dir / f"design_{i:03d}.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">design_{i:03d}_score_{design['score']:.3f}\n")
            f.write(design["sequence"] + "\n")
        print(f"Wrote {fasta_path} (score={design['score']:.3f})")

    return designs

import re
top_designs = prepare_validation_fastas(
    "protein_design/mpnn_outputs/seqs/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.fa",
    "protein_design/validation/fastas/",
    top_n=20
)
```

### 6.2 Batch Validate with ColabFold

```bash
# Predict structures for all designed sequences
colabfold_batch \
  protein_design/validation/fastas/ \
  protein_design/validation/predictions/ \
  --num-recycle 3 \
  --num-models 1 \
  --amber \
  --use-gpu-relax
```

> **Tip**: For validation, using `--num-models 1` (instead of 5) saves significant compute while still providing reliable structures. Reserve 5-model predictions for your final top candidates.

### 6.3 Compute Self-Consistent RMSD (scRMSD)

The scRMSD measures how well the designed sequence reproduces the target backbone when folded by AlphaFold2. Low scRMSD indicates the design is self-consistent.

```python
from Bio.PDB import PDBParser, Superimposer
import numpy as np
from pathlib import Path

def compute_ca_rmsd(pdb_path1, pdb_path2, chain_id="A"):
    """
    Compute C-alpha RMSD between two structures after superposition.
    """
    parser = PDBParser(QUIET=True)
    struct1 = parser.get_structure("s1", pdb_path1)
    struct2 = parser.get_structure("s2", pdb_path2)

    # Extract CA atoms
    def get_ca_atoms(structure, chain):
        atoms = []
        for model in structure:
            if chain in model:
                for residue in model[chain]:
                    if residue.id[0] == " " and "CA" in residue:
                        atoms.append(residue["CA"])
        return atoms

    ca1 = get_ca_atoms(struct1, chain_id)
    ca2 = get_ca_atoms(struct2, chain_id)

    # Align on the shorter set
    n = min(len(ca1), len(ca2))
    if n == 0:
        return float("inf")

    sup = Superimposer()
    sup.set_atoms(ca1[:n], ca2[:n])
    return sup.rms


def validate_designs(target_pdb, validation_dir, chain_ids=["A", "B"]):
    """
    Compare all validation predictions against the target structure.
    """
    results = []
    val_path = Path(validation_dir)

    for pdb_file in sorted(val_path.glob("*.pdb")):
        if "rank_001" not in pdb_file.name:
            continue

        rmsds = {}
        for chain in chain_ids:
            rmsd = compute_ca_rmsd(target_pdb, str(pdb_file), chain)
            rmsds[f"rmsd_chain_{chain}"] = rmsd

        avg_rmsd = np.mean(list(rmsds.values()))
        results.append({
            "design": pdb_file.stem,
            "avg_scRMSD": avg_rmsd,
            **rmsds
        })

    # Sort by average scRMSD
    results.sort(key=lambda x: x["avg_scRMSD"])

    print(f"{'Design':<50} {'Avg scRMSD':>10} {'Chain A':>10} {'Chain B':>10}")
    print("-" * 82)
    for r in results:
        status = "PASS" if r["avg_scRMSD"] < 2.0 else "FAIL"
        print(f"{r['design']:<50} {r['avg_scRMSD']:>10.2f} "
              f"{r['rmsd_chain_A']:>10.2f} {r['rmsd_chain_B']:>10.2f}  {status}")

    passed = [r for r in results if r["avg_scRMSD"] < 2.0]
    print(f"\n{len(passed)}/{len(results)} designs passed (scRMSD < 2.0 A)")

    return results

results = validate_designs(
    "protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb",
    "protein_design/validation/predictions/"
)
```

### 6.4 Filter by pLDDT and scRMSD

Apply dual filtering: the designed sequence must both fold correctly (low scRMSD) and fold confidently (high pLDDT).

```python
def filter_designs(results, validation_dir, plddt_threshold=80, rmsd_threshold=2.0):
    """Filter designs by pLDDT and scRMSD thresholds."""
    parser = PDBParser(QUIET=True)
    filtered = []

    for r in results:
        if r["avg_scRMSD"] >= rmsd_threshold:
            continue

        # Find the corresponding PDB
        pdb_files = list(Path(validation_dir).glob(f"*{r['design']}*.pdb"))
        if not pdb_files:
            continue

        # Compute mean pLDDT from B-factors
        structure = parser.get_structure("s", str(pdb_files[0]))
        plddts = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " " and "CA" in residue:
                        plddts.append(residue["CA"].get_bfactor())

        mean_plddt = np.mean(plddts)
        if mean_plddt >= plddt_threshold:
            filtered.append({
                **r,
                "mean_pLDDT": mean_plddt,
            })

    print(f"\nFiltered designs (pLDDT >= {plddt_threshold}, scRMSD < {rmsd_threshold} A):")
    for d in filtered:
        print(f"  {d['design']}: scRMSD={d['avg_scRMSD']:.2f}, pLDDT={d['mean_pLDDT']:.1f}")

    return filtered

final_designs = filter_designs(results, "protein_design/validation/predictions/")
```

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure4-validation-scatter.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 4. Validation results. Each dot is a designed sequence. X-axis: self-consistent RMSD (lower is better). Y-axis: mean pLDDT of the AF2 prediction (higher is better). Designs in the upper-left quadrant (green box) pass both filters.
</div>

---

## 7. Application: Antibody CDR Redesign

Let us walk through a concrete example: redesigning the CDR-H3 loop of an anti-HER2 antibody to explore sequence diversity while maintaining structural compatibility.

### 7.1 Load and Inspect the Antibody

```python
from Bio.PDB import PDBParser

pdb_path = "protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("ab", pdb_path)

# Extract CDR-H3 sequence (residues 95-102 on chain A).
# residue.get_resname() returns 3-letter codes (e.g. "ALA"), so concatenating
# them gives "ALAGLYSER..." rather than a proper one-letter sequence — we
# need to convert each 3-letter code to its one-letter equivalent first.
from Bio.PDB.Polypeptide import protein_letters_3to1

cdr_h3_seq = ""
for residue in structure[0]["A"]:
    if residue.id[0] == " " and 95 <= residue.id[1] <= 102:
        three = residue.get_resname().upper()
        cdr_h3_seq += protein_letters_3to1.get(three, "X")

print(f"Original CDR-H3: {cdr_h3_seq}")
```

### 7.2 Design CDR-H3 Only

For focused CDR redesign, fix everything except CDR-H3:

```python
import json

pdb_name = "anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000"

# Fix all positions EXCEPT CDR-H3 (95-102) on chain A
parser = PDBParser(QUIET=True)
structure = parser.get_structure("ab", pdb_path)

fixed_a = [r.id[1] for r in structure[0]["A"]
           if r.id[0] == " " and not (95 <= r.id[1] <= 102)]
fixed_b = [r.id[1] for r in structure[0]["B"] if r.id[0] == " "]

fixed_positions = {pdb_name: {"A": fixed_a, "B": fixed_b}}
with open("protein_design/inputs/fixed_cdrh3_only.jsonl", "w") as f:
    f.write(json.dumps(fixed_positions) + "\n")
```

```bash
cd ProteinMPNN

python protein_mpnn_run.py \
  --pdb_path ../protein_design/af2_predictions/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb \
  --chain_id_jsonl ../protein_design/inputs/chain_id.jsonl \
  --fixed_positions_jsonl ../protein_design/inputs/fixed_cdrh3_only.jsonl \
  --out_folder ../protein_design/mpnn_outputs/cdrh3_only/ \
  --num_seq_per_target 200 \
  --sampling_temp "0.1 0.15 0.2" \
  --seed 42 \
  --batch_size 1
```

### 7.3 Analyze CDR-H3 Diversity

```python
from collections import Counter

def analyze_cdr_diversity(fasta_path, cdr_start, cdr_end, chain_sep_pos):
    """
    Analyze amino acid diversity at CDR positions across designs.
    chain_sep_pos: position of ':' separator in the multi-chain sequence.
    """
    cdr_sequences = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                # For multi-chain, split by '/' and take chain A
                chains = line.split("/")
                if len(chains) >= 1:
                    chain_a = chains[0]
                    cdr_seq = chain_a[cdr_start:cdr_end]
                    if len(cdr_seq) == (cdr_end - cdr_start):
                        cdr_sequences.append(cdr_seq)

    unique = set(cdr_sequences)
    print(f"Total CDR sequences: {len(cdr_sequences)}")
    print(f"Unique CDR sequences: {len(unique)}")
    print(f"Diversity ratio: {len(unique)/max(len(cdr_sequences),1):.2%}")

    # Position-wise amino acid frequencies
    if cdr_sequences:
        cdr_len = len(cdr_sequences[0])
        print(f"\nPosition-wise amino acid frequencies (CDR length={cdr_len}):")
        for pos in range(cdr_len):
            aa_counts = Counter(seq[pos] for seq in cdr_sequences)
            top3 = aa_counts.most_common(3)
            top3_str = ", ".join(f"{aa}({count})" for aa, count in top3)
            print(f"  Position {cdr_start + pos + 1}: {top3_str}")

analyze_cdr_diversity(
    "protein_design/mpnn_outputs/cdrh3_only/seqs/anti_HER2_Fab_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.fa",
    cdr_start=94, cdr_end=102, chain_sep_pos=120
)
```

### 7.4 Structural Superposition of Top Designs

```python
import py3Dmol

def overlay_designs(target_pdb, design_pdbs, max_show=5):
    """Overlay designed structures on the target for visual comparison."""
    view = py3Dmol.view(width=800, height=600)

    # Add target in gray
    target_text = Path(target_pdb).read_text()
    view.addModel(target_text, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "gray", "opacity": 0.5}})

    # Add designs in different colors
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    for i, dpdb in enumerate(design_pdbs[:max_show]):
        design_text = Path(dpdb).read_text()
        view.addModel(design_text, "pdb")
        view.setStyle({"model": i + 1}, {"cartoon": {"color": colors[i % len(colors)]}})

    view.zoomTo()
    return view.show()
```

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure5-cdr-overlay.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 5. Structural overlay of top 5 CDR-H3 redesigns (colored) on the original antibody (gray). Despite sequence differences, the backbone conformations remain highly consistent, confirming successful inverse folding.
</div>

---

## 8. Application: TCR-pMHC Engineering (Conceptual)

The same workflow applies directly to TCR engineering, an area of growing importance for cancer immunotherapy.

### 8.1 The TCR Design Challenge

T-cell receptors recognize short peptide fragments presented by MHC molecules on cell surfaces. Engineering TCRs with enhanced or altered specificity could improve:

- **TCR-T cell therapy** (distinct from CAR-T): TCR-engineered T cells use a full αβ TCR that recognises peptide-MHC, so boosting pMHC affinity / specificity is exactly the inverse-folding problem addressed here. CAR-T cells, by contrast, use a synthetic scFv against a surface antigen and are **not** pMHC-restricted.
- **Neoantigen targeting**: Design TCRs that recognize patient-specific mutation-derived peptides.
- **Safety optimization**: Reduce cross-reactivity with self-peptides to minimize autoimmune toxicity.

### 8.2 Workflow Adaptation for TCR Design

The workflow maps directly:

1. **Structure Prediction**: Use ColabFold in multimer mode to predict the TCR alpha/beta chains in complex with peptide-MHC.

```bash
cat > tcr_pmhc.fasta << 'EOF'
>TCR_pMHC_complex
TCRA_SEQUENCE:TCRB_SEQUENCE:PEPTIDE_SEQUENCE:MHC_SEQUENCE
EOF

colabfold_batch tcr_pmhc.fasta tcr_predictions/ \
  --num-recycle 6 \
  --num-models 5 \
  --amber
```

> **Note**: Use `--num-recycle 6` for large complexes (>800 residues). The increased recycling helps AlphaFold2 resolve inter-chain contacts more accurately.

2. **Sequence Design**: Fix the MHC and peptide chains entirely. Fix the TCR framework regions. Redesign only the CDR3 loops (alpha and beta), which are the primary determinants of peptide-MHC specificity.

3. **Validation**: Predict structures of redesigned TCRs and filter by:
   - scRMSD < 2.0 angstroms for the entire complex
   - pLDDT > 80 at the binding interface
   - Low inter-chain PAE between TCR CDR3 and peptide (indicates confident interface prediction)

### 8.3 Important Caveats

- **AlphaFold2 confidence does not equal binding affinity.** A high-pLDDT prediction means the structure is likely correct, not that binding is strong.
- **TCR-pMHC complexes are challenging** for AlphaFold2 multimer. Always validate predictions against known crystal structures when available.
- **Experimental validation is essential.** Computational designs should be treated as hypotheses to test with binding assays (SPR, flow cytometry) and functional assays (T cell activation, cytokine release).

---

## 9. Tips and Pitfalls

### MSA Quality Matters for AlphaFold2

- ColabFold uses the MMseqs2 server for fast MSA generation, which works well for most proteins.
- For engineered or heavily mutated sequences, MSA quality may be poor. Consider using `--msa-mode single_sequence` for purely sequence-based prediction (less accurate but avoids MSA bias).
- For antibodies and TCRs, the variable region MSAs can be noisy. Paired MSAs (from OAS or similar databases) can help.

### Temperature Tuning in ProteinMPNN

```python
# Quick experiment: generate 50 sequences at each temperature
# and compare score distributions
temps_to_test = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
```

- **Start at T=0.1** for conservative redesign (high recovery, low risk).
- **Use T=0.2** for moderate exploration (recommended default).
- **Go above T=0.3** only if you specifically need high diversity and plan to filter aggressively.

### When to Use Tied vs. Untied Design

- **Tied design**: Symmetric complexes (homodimers). Both chains get the same designed sequence.
- **Untied design** (default): Asymmetric complexes (Fab heavy + light chain). Each chain gets an independent sequence.

```bash
# For tied design of a homodimer (chains A and B get same sequence)
python protein_mpnn_run.py \
  --pdb_path homodimer.pdb \
  --tied_positions_jsonl tied_positions.jsonl \
  --out_folder output/ \
  --num_seq_per_target 100 \
  --sampling_temp "0.2"
```

### Common Failure Modes

1. **High scRMSD despite good ProteinMPNN score**: The designed sequence folds into a different structure. Usually caused by designing too many positions simultaneously or using high temperature.
2. **Low pLDDT in redesigned regions**: AlphaFold2 is uncertain about the redesigned region. Consider using a more conservative temperature or redesigning fewer positions.
3. **Chain separation in AF2 prediction**: The chains drift apart in the predicted structure. This often indicates the designed interface residues are incompatible. Try fixing key interface contacts.

### Limitations to Keep in Mind

- **AF2 confidence is not binding affinity.** pLDDT and PAE measure structural confidence, not thermodynamic stability or binding energy.
- **ProteinMPNN optimizes for backbone compatibility**, not for function. A sequence that folds correctly may not bind its target.
- **Neither tool models post-translational modifications**, glycosylation, or solvent effects at high fidelity.
- **Experimental validation remains mandatory** for any design intended for therapeutic use.

---

## 10. Key Takeaways

### Pipeline Summary

```
Input Sequence
     │
     ▼
┌─────────────────┐
│  AlphaFold2      │  Structure prediction
│  (ColabFold)     │  → PDB + confidence scores
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Define Design   │  Select designable positions
│  Specifications  │  (CDRs, interface, etc.)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  ProteinMPNN     │  Inverse folding
│                  │  → 100s of candidate sequences
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  AlphaFold2      │  Validation
│  (ColabFold)     │  → scRMSD + pLDDT filtering
└────────┬────────┘
         │
         ▼
  Top Designs for
  Experimental Testing
```

### When to Use This Workflow vs. Alternatives

| Scenario                                  | Recommended Tool                                          |
| ----------------------------------------- | --------------------------------------------------------- |
| Redesign specific regions (CDRs, loops)   | **This workflow** (AF2 + ProteinMPNN)                     |
| De novo backbone generation               | **RFdiffusion** + ProteinMPNN                             |
| Small molecule binding site design        | **RFdiffusion** or **Rosetta**                            |
| Sequence optimization only (no structure) | **ESM-IF** or **ProteinMPNN** with experimental structure |
| Ultra-large scale virtual screening       | **ESMFold** (faster than AF2) + ProteinMPNN               |

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/blog/alphafold-proteinmpnn/figure6-pipeline-comparison.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>
<div class="caption">
    Figure 6. Decision tree for choosing a computational protein design pipeline. The AF2 + ProteinMPNN workflow (highlighted) is ideal for fixed-backbone sequence design tasks.
</div>

### Resources

- **ColabFold**: [github.com/sokrypton/ColabFold](https://github.com/sokrypton/ColabFold)
- **ProteinMPNN**: [github.com/dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
- **ProteinMPNN paper**: Dauparas et al., _Science_ (2022). [DOI: 10.1126/science.add2187](https://doi.org/10.1126/science.add2187)
- **AlphaFold2 paper**: Jumper et al., _Nature_ (2021). [DOI: 10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- **ColabFold paper**: Mirdita et al., _Nature Methods_ (2022). [DOI: 10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)
- **RFdiffusion**: Watson et al., _Nature_ (2023). [DOI: 10.1038/s41586-023-06415-8](https://doi.org/10.1038/s41586-023-06415-8)

---

_This tutorial provides a computational starting point. All designs should be validated experimentally before drawing biological conclusions or pursuing therapeutic applications._
