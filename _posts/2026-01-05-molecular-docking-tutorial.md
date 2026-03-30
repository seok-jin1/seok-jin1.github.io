---
layout: post
title: "Molecular Docking: DiffDock vs AutoDock Vina for Drug Discovery"
date: 2026-01-05
permalink: /blog/molecular-docking-tutorial/
published: true
categories: [tutorial]
tags:
  - CADD
  - drug-discovery
  - python
  - deep-learning
  - tutorial
---

Molecular docking is a cornerstone of computer-aided drug discovery (CADD). Given a protein target and a small molecule, docking predicts **where** and **how strongly** the ligand binds. This prediction guides medicinal chemists toward promising drug candidates before expensive wet-lab experiments begin.

In this tutorial we compare two paradigms side by side:

- **AutoDock Vina** -- the most widely cited physics-based docking engine, which uses a scoring function derived from empirical force fields.
- **DiffDock** -- a deep-learning method that frames docking as a generative diffusion process over the SE(3) manifold of ligand poses.

By the end you will have working Python code for both approaches, understand their trade-offs, and know how to plug them into a virtual screening pipeline.

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure1-docking-overview.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. Molecular docking overview. A ligand is placed into a protein binding pocket, and a scoring function evaluates the predicted binding affinity.
</div>

---

## 1. Why Molecular Docking Matters

Drug discovery is expensive: the average cost to bring a single drug to market exceeds **$2 billion** and takes more than a decade. Molecular docking accelerates the earliest stages -- hit identification and lead optimization -- by computationally screening millions of compounds against a target protein. A good docking workflow can:

1. **Rank** a compound library by predicted binding affinity.
2. **Reveal** binding poses that explain structure-activity relationships (SAR).
3. **Filter** out molecules unlikely to bind, saving wet-lab resources.

The two dominant approaches are:

| Approach | Representative Tools | Scoring | Speed |
|----------|---------------------|---------|-------|
| Physics-based | AutoDock Vina, Glide, GOLD | Empirical / force-field | Minutes per ligand |
| AI-based | DiffDock, EquiBind, TankBind | Learned from crystal structures | Seconds per ligand |

---

## 2. Traditional Docking with AutoDock Vina

### 2.1 Installation

```bash
# Create a dedicated conda environment
conda create -n docking python=3.10 -y
conda activate docking

# Install core packages
pip install vina meeko rdkit-pypi numpy scipy
pip install prody biopython matplotlib

# Verify Vina
python -c "import vina; print(vina.__version__)"
```

### 2.2 Preparing the Receptor (PDB to PDBQT)

AutoDock Vina requires receptor files in **PDBQT** format, which adds partial charges and atom types to standard PDB coordinates. We use `MolKit` via the `meeko` ecosystem or the `prepare_receptor` utility.

```python
"""
prepare_receptor.py
Convert a PDB file to PDBQT format for AutoDock Vina.
"""
import subprocess
from pathlib import Path
from prody import parsePDB, writePDB

# ── 1. Download the receptor from the PDB ──────────────────────────
pdb_id = "6LU7"  # SARS-CoV-2 main protease (Mpro)
receptor_pdb = Path(f"{pdb_id}_receptor.pdb")

structure = parsePDB(pdb_id)

# Keep only protein atoms (remove water, ligands, ions)
protein = structure.select("protein")
writePDB(str(receptor_pdb), protein)
print(f"Saved protein-only PDB: {receptor_pdb}")

# ── 2. Add hydrogens with reduce (AmberTools) ──────────────────────
receptor_h = receptor_pdb.with_suffix(".h.pdb")
subprocess.run(
    ["reduce", "-BUILD", str(receptor_pdb)],
    stdout=open(receptor_h, "w"),
    check=True,
)
print(f"Added hydrogens: {receptor_h}")

# ── 3. Convert to PDBQT using ADFR Suite's prepare_receptor ───────
receptor_pdbqt = receptor_pdb.with_suffix(".pdbqt")
subprocess.run(
    [
        "prepare_receptor",
        "-r", str(receptor_h),
        "-o", str(receptor_pdbqt),
        "-A", "hydrogens",
    ],
    check=True,
)
print(f"Receptor PDBQT ready: {receptor_pdbqt}")
```

> **Tip**: If `prepare_receptor` is not available, install the [ADFR Suite](https://ccsb.scripps.edu/adfr/) or use the `meeko` Python API as shown below.

### 2.3 Preparing the Ligand (SMILES to PDBQT via Meeko + RDKit)

```python
"""
prepare_ligand.py
Convert a SMILES string to a 3D PDBQT file for docking.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from meeko import MoleculePreparation, PDBQTWriterLegacy

# ── 1. Define the ligand ────────────────────────────────────────────
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
name = "aspirin"

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# ── 2. Generate 3D coordinates ──────────────────────────────────────
params = AllChem.ETKDGv3()
params.randomSeed = 42
AllChem.EmbedMolecule(mol, params)
AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

# ── 3. Save as SDF for reference ────────────────────────────────────
writer = Chem.SDWriter(f"{name}.sdf")
writer.write(mol)
writer.close()

# ── 4. Convert to PDBQT using Meeko ────────────────────────────────
preparator = MoleculePreparation()
mol_setups = preparator.prepare(mol)

for setup in mol_setups:
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
    if is_ok:
        with open(f"{name}.pdbqt", "w") as f:
            f.write(pdbqt_string)
        print(f"Ligand PDBQT saved: {name}.pdbqt")
    else:
        print(f"Error: {error_msg}")
```

### 2.4 Defining the Search Box

The search box tells Vina where to look for binding poses. Coordinates are usually centered on the known binding site or determined from a co-crystallized ligand.

```python
"""
define_search_box.py
Determine the search box center and dimensions from a reference ligand.
"""
import numpy as np
from prody import parsePDB

# Parse the original PDB with the co-crystallized ligand
structure = parsePDB("6LU7")

# Select the co-crystallized inhibitor (residue name N3)
ligand = structure.select("resname N3")
coords = ligand.getCoords()

# Compute center and size
center = coords.mean(axis=0)
padding = 10.0  # Angstroms of padding around the ligand
box_size = coords.ptp(axis=0) + 2 * padding

print(f"Center (x, y, z): {center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}")
print(f"Box size (x, y, z): {box_size[0]:.1f}, {box_size[1]:.1f}, {box_size[2]:.1f}")

# Typical output for 6LU7:
# Center (x, y, z): -10.44, 12.44, 68.92
# Box size (x, y, z): 30.0, 30.0, 30.0
```

### 2.5 Running AutoDock Vina

```python
"""
run_vina.py
Execute molecular docking with AutoDock Vina Python bindings.
"""
from vina import Vina

# ── Initialize Vina ────────────────────────────────────────────────
v = Vina(sf_name="vina")  # Scoring function: 'vina' or 'vinardo'

# ── Set receptor ───────────────────────────────────────────────────
v.set_receptor(rigid_pdbqt="6LU7_receptor.pdbqt")

# ── Set ligand ─────────────────────────────────────────────────────
v.set_ligand_from_file("aspirin.pdbqt")

# ── Define search space ────────────────────────────────────────────
v.compute_vina_maps(
    center=[-10.44, 12.44, 68.92],
    box_size=[30, 30, 30],
)

# ── Run docking ────────────────────────────────────────────────────
v.dock(
    exhaustiveness=32,   # Higher = more thorough (default 8)
    n_poses=10,          # Number of output poses
)

# ── Retrieve results ───────────────────────────────────────────────
energies = v.energies()
print("\nDocking Results:")
print(f"{'Pose':<6} {'Affinity (kcal/mol)':<22} {'RMSD l.b.':<12} {'RMSD u.b.'}")
for i, row in enumerate(energies):
    print(f"{i+1:<6} {row[0]:<22.1f} {row[1]:<12.3f} {row[2]:.3f}")

# ── Save output poses ─────────────────────────────────────────────
v.write_poses("aspirin_vina_out.pdbqt", n_poses=10, overwrite=True)
print("\nPoses saved to aspirin_vina_out.pdbqt")
```

**Interpreting Vina Scores**: Vina reports binding free energy in **kcal/mol**. More negative values indicate stronger predicted binding. A rough guide:

| Affinity (kcal/mol) | Interpretation |
|---------------------|----------------|
| < -10 | Very strong binding |
| -8 to -10 | Strong binding |
| -6 to -8 | Moderate binding |
| > -6 | Weak binding |

> **Caveat**: Vina scores correlate only moderately with experimental binding affinities. They are best used for **ranking** compounds rather than predicting absolute $$K_d$$ values.

### 2.6 Visualizing Docking Poses

```python
"""
visualize_vina_poses.py
Visualize docking results using py3Dmol in a Jupyter notebook.
"""
import py3Dmol
from pathlib import Path

# ── Load receptor and docked ligand ─────────────────────────────────
receptor_pdb = Path("6LU7_receptor.pdb").read_text()
ligand_pdbqt = Path("aspirin_vina_out.pdbqt").read_text()

view = py3Dmol.view(width=800, height=600)

# Add receptor as cartoon
view.addModel(receptor_pdb, "pdb")
view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})

# Add docked ligand (first model only)
view.addModel(ligand_pdbqt, "pdbqt")
view.setStyle(
    {"model": 1},
    {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}},
)

# Add surface around the binding site
view.addSurface(
    py3Dmol.VDW,
    {"opacity": 0.15, "color": "white"},
    {"model": 0, "within": {"distance": 6, "sel": {"model": 1}}},
)

view.zoomTo({"model": 1})
view.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure2-vina-pose.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. AutoDock Vina docking result for aspirin in the SARS-CoV-2 Mpro active site. The protein is shown in cartoon representation, the ligand in green sticks, and the binding pocket surface in transparent white.
</div>

---

## 3. AI-Based Docking with DiffDock

### 3.1 How DiffDock Works

DiffDock treats molecular docking as a **generative modeling** problem. Instead of sampling poses through stochastic search with a hand-crafted scoring function, DiffDock learns a **diffusion process over SE(3)** -- the space of 3D rotations and translations -- from thousands of experimentally determined protein-ligand complexes.

The key ideas:

1. **Forward diffusion**: Starting from the crystal pose, progressively add noise to the ligand's position, orientation, and torsion angles until the pose is random.
2. **Reverse diffusion**: A neural network learns to reverse this process -- given a noisy pose and the protein, predict the denoising step.
3. **Sampling**: At inference time, start from random poses and iteratively denoise to generate plausible binding modes.
4. **Confidence model**: A separate model scores each generated pose, predicting how close it is to the true binding mode.

The diffusion operates on three degrees of freedom:

$$
\mathbf{p}_t = (\mathbf{r}_t, \mathbf{R}_t, \boldsymbol{\tau}_t) \in \mathbb{R}^3 \times SO(3) \times \mathbb{T}^k
$$

where $$\mathbf{r}_t$$ is the translational component, $$\mathbf{R}_t$$ is the rotational component in SO(3), and $$\boldsymbol{\tau}_t$$ represents $$k$$ torsion angles on the torus $$\mathbb{T}^k$$.

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure3-diffdock-architecture.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. DiffDock architecture. A score model learns to reverse the diffusion over translational, rotational, and torsional degrees of freedom. A separate confidence model ranks the generated poses.
</div>

### 3.2 Installation

```bash
# Clone DiffDock repository
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock

# Create environment
conda create -n diffdock python=3.10 -y
conda activate diffdock

# Install PyTorch with CUDA support
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install PyTorch Geometric and dependencies
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch-geometric

# Install DiffDock dependencies
pip install e3nn biopandas scipy networkx spyrmsd
pip install -r requirements.txt

# Download pre-trained model weights
python -c "from utils.download import download_and_extract; download_and_extract()"
```

### 3.3 Running DiffDock

DiffDock accepts a CSV file specifying protein-ligand pairs to dock:

```python
"""
run_diffdock.py
Run DiffDock for molecular docking.
"""
import csv
import subprocess
from pathlib import Path

# ── 1. Prepare input CSV ───────────────────────────────────────────
input_csv = Path("diffdock_input.csv")
with open(input_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
    writer.writerow([
        "mpro_aspirin",          # Complex name
        "6LU7_receptor.pdb",     # Receptor PDB file
        "CC(=O)Oc1ccccc1C(=O)O", # Ligand SMILES
        "",                       # Leave empty when providing PDB
    ])

# ── 2. Run DiffDock inference ──────────────────────────────────────
diffdock_dir = Path("DiffDock")
output_dir = Path("diffdock_results")

cmd = [
    "python", str(diffdock_dir / "inference.py"),
    "--config", str(diffdock_dir / "default_inference_args.yaml"),
    "--protein_ligand_csv", str(input_csv),
    "--out_dir", str(output_dir),
    "--inference_steps", "20",
    "--samples_per_complex", "10",
    "--batch_size", "10",
    "--actual_steps", "18",
    "--no_final_step_noise",
]

print("Running DiffDock inference...")
result = subprocess.run(cmd, capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print(f"Error: {result.stderr}")
```

### 3.4 Parsing DiffDock Results

```python
"""
parse_diffdock_results.py
Parse and analyze DiffDock output.
"""
import re
from pathlib import Path

import numpy as np
from rdkit import Chem

results_dir = Path("diffdock_results/mpro_aspirin")

# ── Collect confidence scores and poses ─────────────────────────────
poses = []
for sdf_file in sorted(results_dir.glob("rank*_confidence*.sdf")):
    # Extract rank and confidence from filename
    match = re.search(r"rank(\d+)_confidence(-?[\d.]+)", sdf_file.name)
    if match:
        rank = int(match.group(1))
        confidence = float(match.group(2))

        mol = Chem.SDMolSupplier(str(sdf_file), removeHs=False)[0]
        if mol is not None:
            poses.append({
                "rank": rank,
                "confidence": confidence,
                "file": sdf_file.name,
                "mol": mol,
                "centroid": np.mean(mol.GetConformer().GetPositions(), axis=0),
            })

# ── Display results ─────────────────────────────────────────────────
poses.sort(key=lambda x: x["rank"])
print(f"\nDiffDock Results ({len(poses)} poses):")
print(f"{'Rank':<6} {'Confidence':<14} {'Centroid (x,y,z)':<30} {'File'}")
print("-" * 70)
for p in poses:
    cx, cy, cz = p["centroid"]
    print(f"{p['rank']:<6} {p['confidence']:<14.3f} ({cx:.1f}, {cy:.1f}, {cz:.1f})     {p['file']}")

# ── Confidence interpretation ───────────────────────────────────────
top_pose = poses[0]
print(f"\nTop pose confidence: {top_pose['confidence']:.3f}")
if top_pose["confidence"] > 0:
    print("  -> High confidence: likely within 2A RMSD of true pose")
elif top_pose["confidence"] > -1.5:
    print("  -> Moderate confidence: may be reasonable")
else:
    print("  -> Low confidence: treat with caution")
```

### 3.5 Comparing Vina and DiffDock Poses

```python
"""
compare_poses.py
Compare poses from Vina and DiffDock by computing RMSD.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from spyrmsd import rmsd as spyrmsd

# ── Load DiffDock top pose ──────────────────────────────────────────
diffdock_mol = Chem.SDMolSupplier("diffdock_results/mpro_aspirin/rank1_confidence0.42.sdf")[0]

# ── Load Vina top pose (convert PDBQT -> SDF via RDKit) ────────────
vina_mol = Chem.MolFromPDBFile("aspirin_vina_pose1.pdb", removeHs=False)

# ── Compute RMSD between the two poses ─────────────────────────────
if diffdock_mol and vina_mol:
    # Align by maximum common substructure
    mcs = Chem.MolFromSmarts(
        Chem.MolToSmarts(
            AllChem.GetBestRMS(diffdock_mol, vina_mol)
        )
    )
    rmsd_value = rdMolAlign.GetBestRMS(diffdock_mol, vina_mol)
    print(f"RMSD between DiffDock (rank1) and Vina (pose1): {rmsd_value:.2f} A")

    # Symmetry-corrected RMSD using spyrmsd
    coords_dd = diffdock_mol.GetConformer().GetPositions()
    coords_vina = vina_mol.GetConformer().GetPositions()
    symmrmsd = spyrmsd.symmrmsd(
        coords_dd, [coords_vina],
        diffdock_mol.GetNumAtoms(),
        diffdock_mol.GetNumAtoms(),
    )
    print(f"Symmetry-corrected RMSD: {symmrmsd[0]:.2f} A")
```

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure4-vina-vs-diffdock.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. Overlay of top-ranked poses from AutoDock Vina (green) and DiffDock (magenta) in the Mpro binding site. Despite different methodologies, both methods often identify similar binding modes for well-defined pockets.
</div>

---

## 4. Preparing Structures from AlphaFold Predictions

When no experimental structure is available, **AlphaFold2** predictions provide a viable alternative. However, predicted structures require extra care.

```python
"""
prepare_alphafold_structure.py
Download and prepare an AlphaFold-predicted structure for docking.
"""
import requests
import numpy as np
from pathlib import Path
from prody import parsePDB, writePDB

# ── 1. Download AlphaFold structure by UniProt ID ───────────────────
uniprot_id = "P0DTD1"  # SARS-CoV-2 polyprotein (contains Mpro)
af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

response = requests.get(af_url)
af_pdb = Path(f"AF-{uniprot_id}.pdb")
af_pdb.write_text(response.text)
print(f"Downloaded: {af_pdb}")

# ── 2. Assess prediction confidence (pLDDT) ────────────────────────
structure = parsePDB(str(af_pdb))
bfactors = structure.select("protein and name CA").getBetas()

print(f"\npLDDT statistics:")
print(f"  Mean:   {np.mean(bfactors):.1f}")
print(f"  Median: {np.median(bfactors):.1f}")
print(f"  Min:    {np.min(bfactors):.1f}")
print(f"  Max:    {np.max(bfactors):.1f}")

# ── 3. Filter low-confidence regions ───────────────────────────────
# Residues with pLDDT < 70 are unreliable for docking
high_conf = structure.select("protein and beta > 70")
writePDB(f"AF-{uniprot_id}_high_conf.pdb", high_conf)
print(f"\nRetained {high_conf.numResidues()} / {structure.select('protein').numResidues()} residues")

# ── 4. Quick energy minimization with OpenMM ────────────────────────
"""
AlphaFold structures lack proper hydrogen placement and may have
minor clashes. A brief minimization resolves these issues.
"""
from openmm.app import PDBFile, ForceField, Simulation, Modeller
from openmm import LangevinMiddleIntegrator, unit

pdb = PDBFile(f"AF-{uniprot_id}_high_conf.pdb")
forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4)

system = forcefield.createSystem(modeller.topology)
integrator = LangevinMiddleIntegrator(
    300 * unit.kelvin, 1.0 / unit.picoseconds, 0.002 * unit.picoseconds
)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Minimize
simulation.minimizeEnergy(maxIterations=500)

# Save minimized structure
positions = simulation.context.getState(getPositions=True).getPositions()
with open(f"AF-{uniprot_id}_minimized.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)
print("Energy minimization complete")
```

> **Best practices for docking into AlphaFold structures**:
> - Only dock into regions with **pLDDT > 70** (confident predictions).
> - Always energy-minimize before docking.
> - Be cautious with **disordered loops** -- AlphaFold may place them incorrectly.
> - Validate results against known actives if possible.

---

## 5. Virtual Screening Workflow

Now let us scale up: screen a library of compounds against our target.

```python
"""
virtual_screening.py
Screen a compound library using AutoDock Vina.
"""
import csv
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

# ── 1. Load and filter compound library ─────────────────────────────
def load_library(sdf_file: str, max_compounds: int = 1000):
    """Load compounds and apply basic drug-likeness filters."""
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    compounds = []

    for mol in supplier:
        if mol is None:
            continue

        # Lipinski's Rule of Five
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
            compounds.append(mol)

        if len(compounds) >= max_compounds:
            break

    print(f"Loaded {len(compounds)} drug-like compounds")
    return compounds


# ── 2. Prepare ligand PDBQT from RDKit mol ──────────────────────────
def mol_to_pdbqt(mol):
    """Convert an RDKit mol to PDBQT string."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    for setup in mol_setups:
        pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            return pdbqt_string
    return None


# ── 3. Dock a single compound ───────────────────────────────────────
def dock_single(args):
    """Dock one ligand and return the best score."""
    idx, smiles, receptor_pdbqt, center, box_size = args

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return idx, smiles, None, "Invalid SMILES"

    pdbqt_string = mol_to_pdbqt(mol)
    if pdbqt_string is None:
        return idx, smiles, None, "PDBQT preparation failed"

    try:
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rigid_pdbqt=receptor_pdbqt)
        v.set_ligand_from_string(pdbqt_string)
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=8, n_poses=1)
        score = v.energies()[0][0]
        return idx, smiles, score, "OK"
    except Exception as e:
        return idx, smiles, None, str(e)


# ── 4. Run virtual screening ───────────────────────────────────────
def run_screening(
    smiles_list: list[str],
    receptor_pdbqt: str,
    center: list[float],
    box_size: list[float],
    n_workers: int = 4,
    output_csv: str = "screening_results.csv",
):
    """Screen a list of SMILES against a receptor."""
    args_list = [
        (i, smi, receptor_pdbqt, center, box_size)
        for i, smi in enumerate(smiles_list)
    ]

    results = []
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(dock_single, args): args for args in args_list}
        for future in as_completed(futures):
            idx, smiles, score, status = future.result()
            results.append((idx, smiles, score, status))

            if len(results) % 50 == 0:
                elapsed = time.time() - start_time
                rate = len(results) / elapsed
                print(f"  Docked {len(results)}/{len(smiles_list)} "
                      f"({rate:.1f} compounds/sec)")

    # Sort by score (best first)
    results.sort(key=lambda x: x[2] if x[2] is not None else 0)

    # Save results
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "smiles", "vina_score_kcal_mol", "status"])
        for rank, (idx, smiles, score, status) in enumerate(results, 1):
            writer.writerow([rank, smiles, f"{score:.2f}" if score else "N/A", status])

    elapsed = time.time() - start_time
    print(f"\nScreening complete: {len(results)} compounds in {elapsed:.0f}s")
    print(f"Results saved to {output_csv}")

    return results


# ── Example usage ───────────────────────────────────────────────────
if __name__ == "__main__":
    # Example SMILES library (replace with your actual library)
    smiles_library = [
        "CC(=O)Oc1ccccc1C(=O)O",           # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",       # Ibuprofen
        "OC(=O)c1ccccc1O",                   # Salicylic acid
        "CC(=O)Nc1ccc(O)cc1",               # Acetaminophen
        "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",  # Pyrene (negative control)
        # ... add more SMILES
    ]

    results = run_screening(
        smiles_list=smiles_library,
        receptor_pdbqt="6LU7_receptor.pdbqt",
        center=[-10.44, 12.44, 68.92],
        box_size=[30, 30, 30],
        n_workers=4,
    )

    # Print top hits
    print("\nTop 5 Hits:")
    for idx, smi, score, status in results[:5]:
        print(f"  {score:>8.1f} kcal/mol  {smi}")
```

### 5.1 Batch DiffDock Screening

```python
"""
batch_diffdock.py
Screen a compound library using DiffDock.
"""
import csv
from pathlib import Path
import subprocess

# ── 1. Generate input CSV for DiffDock ──────────────────────────────
smiles_library = {
    "aspirin":       "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen":     "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "salicylic_acid": "OC(=O)c1ccccc1O",
    "acetaminophen": "CC(=O)Nc1ccc(O)cc1",
}

input_csv = Path("diffdock_batch_input.csv")
with open(input_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
    for name, smiles in smiles_library.items():
        writer.writerow([name, "6LU7_receptor.pdb", smiles, ""])

# ── 2. Run DiffDock batch inference ─────────────────────────────────
cmd = [
    "python", "DiffDock/inference.py",
    "--config", "DiffDock/default_inference_args.yaml",
    "--protein_ligand_csv", str(input_csv),
    "--out_dir", "diffdock_screening_results",
    "--inference_steps", "20",
    "--samples_per_complex", "5",
    "--batch_size", "10",
    "--no_final_step_noise",
]

subprocess.run(cmd, check=True)

# ── 3. Collect top confidence scores ────────────────────────────────
import re

results = []
results_dir = Path("diffdock_screening_results")
for complex_dir in sorted(results_dir.iterdir()):
    if not complex_dir.is_dir():
        continue

    best_confidence = -float("inf")
    best_file = None
    for sdf in complex_dir.glob("rank1_confidence*.sdf"):
        match = re.search(r"confidence(-?[\d.]+)", sdf.name)
        if match:
            conf = float(match.group(1))
            if conf > best_confidence:
                best_confidence = conf
                best_file = sdf.name

    results.append({
        "name": complex_dir.name,
        "confidence": best_confidence,
        "file": best_file,
    })

results.sort(key=lambda x: x["confidence"], reverse=True)
print("\nDiffDock Screening Results:")
print(f"{'Compound':<20} {'Confidence':<14} {'File'}")
for r in results:
    print(f"{r['name']:<20} {r['confidence']:<14.3f} {r['file']}")
```

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure5-screening-workflow.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Virtual screening pipeline. Compounds are filtered for drug-likeness, docked with Vina or DiffDock, ranked by score, and top hits proceed to interaction analysis.
</div>

---

## 6. Post-Docking Analysis: Interaction Fingerprints

Scores alone do not tell the full story. We need to understand **which** protein-ligand interactions drive binding. Two excellent tools for this are **ProLIF** (Protein-Ligand Interaction Fingerprints) and **PLIP** (Protein-Ligand Interaction Profiler).

### 6.1 Interaction Fingerprints with ProLIF

```python
"""
interaction_analysis.py
Compute protein-ligand interaction fingerprints using ProLIF.
"""
import prolif as plf
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt

# ── 1. Load protein and docked ligand ───────────────────────────────
protein = mda.Universe("6LU7_receptor.pdb")
ligand = mda.Universe("aspirin_vina_pose1.sdf")

# ── 2. Compute interaction fingerprint ─────────────────────────────
fp = plf.Fingerprint(
    interactions=[
        "HBDonor",
        "HBAcceptor",
        "PiStacking",
        "PiCation",
        "CationPi",
        "Hydrophobic",
        "VdWContact",
        "SaltBridge",
    ]
)

# Generate fingerprint for the docked pose
fp.run(
    ligand.select_atoms("all"),
    protein.select_atoms("all"),
    residues=None,  # Analyze all residues
)

# ── 3. Convert to DataFrame for analysis ───────────────────────────
df = fp.to_dataframe()
print("\nInteraction Fingerprint:")
print(df)

# ── 4. Identify key interactions ────────────────────────────────────
bv = fp.to_bitvectors()[0]  # First (only) frame
ifp = fp.ifp[0]  # Interaction details

print("\nKey interactions:")
for (lig_res, prot_res), interactions in ifp.items():
    for interaction_type, metadata in interactions.items():
        if metadata:
            print(f"  {interaction_type}: {lig_res} -- {prot_res}")

# ── 5. Visualize interaction fingerprint ────────────────────────────
fig, ax = plt.subplots(figsize=(12, 4))
fp_df = fp.to_dataframe(return_atoms=False)
fp_df_plot = fp_df.droplevel("ligand", axis=1)

# Plot as heatmap
im = ax.imshow(
    fp_df_plot.values,
    cmap="YlOrRd",
    aspect="auto",
    interpolation="none",
)
ax.set_yticks([0])
ax.set_yticklabels(["Aspirin"])
ax.set_xticks(range(len(fp_df_plot.columns)))
ax.set_xticklabels(
    [f"{res}\n{itype}" for res, itype in fp_df_plot.columns],
    rotation=90,
    fontsize=7,
)
ax.set_title("Protein-Ligand Interaction Fingerprint")
plt.colorbar(im, ax=ax, label="Interaction Present")
plt.tight_layout()
plt.savefig("interaction_fingerprint.png", dpi=300, bbox_inches="tight")
plt.show()
```

### 6.2 Pose Clustering

When multiple docking runs or multiple poses are generated, clustering helps identify distinct binding modes.

```python
"""
pose_clustering.py
Cluster docked poses by RMSD to identify distinct binding modes.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# ── 1. Load all poses ──────────────────────────────────────────────
pose_files = sorted(Path("diffdock_results/mpro_aspirin").glob("rank*_confidence*.sdf"))
mols = []
for f in pose_files:
    mol = Chem.SDMolSupplier(str(f), removeHs=False)[0]
    if mol is not None:
        mols.append(mol)

print(f"Loaded {len(mols)} poses")

# ── 2. Compute pairwise RMSD matrix ────────────────────────────────
n = len(mols)
rmsd_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):
        rmsd = rdMolAlign.GetBestRMS(mols[i], mols[j])
        rmsd_matrix[i, j] = rmsd
        rmsd_matrix[j, i] = rmsd

# ── 3. Hierarchical clustering ─────────────────────────────────────
condensed = squareform(rmsd_matrix)
Z = linkage(condensed, method="average")
clusters = fcluster(Z, t=2.0, criterion="distance")  # 2 A cutoff

print(f"\nClusters (2.0 A RMSD cutoff):")
for c in sorted(set(clusters)):
    members = np.where(clusters == c)[0]
    print(f"  Cluster {c}: {len(members)} poses (indices: {members.tolist()})")

# ── 4. Dendrogram ──────────────────────────────────────────────────
from scipy.cluster.hierarchy import dendrogram

fig, ax = plt.subplots(figsize=(10, 5))
dendrogram(Z, ax=ax, labels=[f"Pose {i+1}" for i in range(n)])
ax.set_ylabel("RMSD (A)")
ax.set_title("Hierarchical Clustering of Docked Poses")
ax.axhline(y=2.0, color="red", linestyle="--", label="Cutoff = 2.0 A")
ax.legend()
plt.tight_layout()
plt.savefig("pose_clustering.png", dpi=300, bbox_inches="tight")
plt.show()
```

{% include figure.liquid loading="eager" path="assets/img/blog/molecular-docking/figure6-interaction-analysis.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 6. Post-docking analysis. Left: Protein-ligand interaction fingerprint showing hydrogen bonds, hydrophobic contacts, and pi-stacking. Right: Hierarchical clustering dendrogram of docked poses reveals two distinct binding modes.
</div>

---

## 7. Comparison: Vina vs DiffDock vs Glide

| Feature | AutoDock Vina | DiffDock | Glide (Schrodinger) |
|---------|:------------:|:--------:|:-------------------:|
| **License** | Open source (Apache 2.0) | Open source (MIT) | Commercial |
| **Scoring** | Empirical (hybrid) | Learned (confidence) | Empirical (GlideScore) |
| **Speed (per ligand)** | 1-5 min | 5-30 sec (GPU) | 1-10 min |
| **GPU acceleration** | No | Yes (required) | Optional |
| **Receptor flexibility** | Limited (side chains) | Implicit (learned) | Side-chain sampling |
| **Success rate (RMSD < 2A)** | ~50-60% | ~35-40% (v1.1) | ~70-80% |
| **Blind docking** | Poor | Good | Poor |
| **Ensemble docking** | Manual | Built-in | Manual |
| **Virtual screening** | Excellent | Good (batched) | Excellent |
| **Water molecules** | Not modeled | Not modeled | Explicit waters |
| **Metal coordination** | Limited | Not supported | Supported |
| **Covalent docking** | No | No | Yes (CovDock) |
| **Best for** | Large-scale screening | Novel pockets, blind docking | Accuracy-critical projects |

> **Key insight**: DiffDock excels at **blind docking** (no prior knowledge of the binding site) and produces diverse poses, but its success rate on standard benchmarks currently trails physics-based methods. The best strategy is often to **use both** and consensus-score the results.

---

## 8. Limitations and Pitfalls

Molecular docking is powerful but comes with important caveats.

### 8.1 Scoring Function Accuracy

No scoring function reliably predicts absolute binding affinities. Pearson correlations between docking scores and experimental $$K_d$$ values typically range from **0.3 to 0.5**. Use docking scores for **ranking** (enrichment) rather than affinity prediction. For more accurate free energy estimates, consider:

- **MM-GBSA / MM-PBSA** rescoring
- **Free energy perturbation (FEP)** calculations
- **Machine learning rescoring** (e.g., RF-Score, OnionNet)

### 8.2 Receptor Flexibility

Proteins are not rigid. Induced-fit effects can drastically change binding pocket shape upon ligand binding. Strategies to address this:

```python
"""
ensemble_docking.py
Dock into multiple receptor conformations (ensemble docking).
"""
from vina import Vina
import numpy as np

# List of receptor conformations from MD simulation or NMR ensemble
receptor_ensemble = [
    "receptor_conf1.pdbqt",
    "receptor_conf2.pdbqt",
    "receptor_conf3.pdbqt",
    "receptor_conf4.pdbqt",
    "receptor_conf5.pdbqt",
]

all_scores = []
for conf_file in receptor_ensemble:
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(rigid_pdbqt=conf_file)
    v.set_ligand_from_file("aspirin.pdbqt")
    v.compute_vina_maps(center=[-10.44, 12.44, 68.92], box_size=[30, 30, 30])
    v.dock(exhaustiveness=16, n_poses=1)

    score = v.energies()[0][0]
    all_scores.append(score)
    print(f"  {conf_file}: {score:.1f} kcal/mol")

# Ensemble scoring: take the best score across conformations
best_score = min(all_scores)
mean_score = np.mean(all_scores)
print(f"\nEnsemble best:  {best_score:.1f} kcal/mol")
print(f"Ensemble mean:  {mean_score:.1f} kcal/mol")
```

### 8.3 Water and Metal Handling

Crystallographic water molecules can mediate protein-ligand interactions. Most docking programs ignore them, which can lead to:

- **False negatives**: Compounds that bind via water-mediated hydrogen bonds score poorly.
- **False positives**: Hydrophobic compounds score well in pockets that are actually hydrated.

For metal-containing binding sites (e.g., zinc proteases, metalloenzymes), specialized parameters or constraints are needed.

### 8.4 Common Mistakes to Avoid

1. **Forgetting to protonate at physiological pH** -- always check ionization states.
2. **Ignoring tautomers** -- enumerate tautomers with RDKit before docking.
3. **Using too small a search box** -- you may miss the binding site.
4. **Over-interpreting scores** -- a 0.5 kcal/mol difference is within noise.
5. **Not validating** -- always run known actives/decoys (e.g., DUD-E benchmarks) to calibrate.

```python
"""
enumerate_tautomers.py
Properly enumerate tautomers and protomers before docking.
"""
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

smiles = "c1[nH]c2c(=O)[nH]c(nc2n1)N"  # Guanine

mol = Chem.MolFromSmiles(smiles)
enumerator = rdMolStandardize.TautomerEnumerator()
tautomers = enumerator.Enumerate(mol)

print(f"Found {len(tautomers)} tautomers for {smiles}:")
for i, taut in enumerate(tautomers):
    print(f"  {i+1}. {Chem.MolToSmiles(taut)}")

# Best practice: dock all reasonable tautomers and keep the best score
```

---

## 9. Key Takeaways

1. **AutoDock Vina** remains the workhorse for large-scale virtual screening. It is fast, well-validated, and easy to use through its Python bindings. Start here if you have a well-defined binding site.

2. **DiffDock** shines in scenarios where the binding site is unknown (blind docking) or when you want to quickly explore diverse binding modes. Its generative approach produces conformationally diverse poses that physics-based methods may miss.

3. **Use both methods together**. Consensus scoring -- keeping compounds that rank well in both Vina and DiffDock -- often outperforms either method alone.

4. **Post-docking analysis is essential**. Raw scores are noisy. Always inspect binding poses visually, compute interaction fingerprints, and cluster poses before selecting compounds for experimental testing.

5. **Validate your workflow** with known actives and decoys before trusting it for prospective screening. Tools like the DUD-E benchmark set make this straightforward.

6. **AlphaFold structures** can substitute for experimental structures, but require careful preparation -- filter by pLDDT, minimize, and validate against known binders when possible.

7. **Know the limitations**: no docking method handles receptor flexibility, water molecules, or entropic contributions perfectly. For high-stakes decisions, follow up with more rigorous methods like FEP or experimental assays.

---

**Useful Resources**:

- [AutoDock Vina Documentation](https://autodock-vina.readthedocs.io/)
- [DiffDock Paper (Corso et al., 2023)](https://arxiv.org/abs/2210.01776)
- [Meeko Documentation](https://meeko.readthedocs.io/)
- [ProLIF Documentation](https://prolif.readthedocs.io/)
- [PDBe Binding Site Validation](https://www.ebi.ac.uk/pdbe/)
- [DUD-E Benchmarking](http://dude.docking.org/)
