---
layout: post
title: "Molecular Dynamics with OpenMM: From Setup to Analysis"
date: 2025-12-15
permalink: /blog/openmm-molecular-dynamics/
published: true
categories: [tutorial]
tags:
  - molecular-dynamics
  - computational-biology
  - python
  - tutorial
  - CADD
---

Proteins are not static sculptures. They breathe, twist, and fluctuate on timescales ranging from femtoseconds to milliseconds, and these motions are often the key to understanding function. Molecular dynamics (MD) simulation lets us watch these motions unfold in atomistic detail by numerically integrating Newton's equations of motion for every atom in the system.

Why does this matter? Three concrete examples:

- **Drug binding**: A ligand must navigate a dynamic binding pocket whose shape changes over time. MD reveals transient sub-pockets and binding pathways invisible to static crystal structures.
- **Intrinsically disordered regions (IDRs)**: Roughly 30% of eukaryotic proteins contain disordered segments that lack a single folded state. MD is one of the few tools that can characterize their conformational ensembles.
- **Protein stability**: Point mutations can shift the balance between folded and unfolded states. Free energy perturbation calculations built on MD quantify these shifts.

In this tutorial, we will walk through a complete MD workflow using **OpenMM** for simulation and **MDAnalysis** for trajectory analysis. Every code block is designed to be runnable end-to-end on a protein of your choice.

---

## 1. Environment Setup

Install the required packages. OpenMM is best installed via conda because it includes GPU-accelerated kernels (CUDA/OpenCL). MDAnalysis and nglview can be installed with pip.

```bash
# Create a dedicated environment (recommended)
conda create -n md-tutorial python=3.11 -y
conda activate md-tutorial

# OpenMM via conda-forge
conda install -c conda-forge openmm -y

# Analysis and visualization
pip install mdanalysis nglview matplotlib numpy
```

Verify that OpenMM detects your hardware:

```python
import openmm
print(openmm.__version__)
print(openmm.Platform.getNumPlatforms())
for i in range(openmm.Platform.getNumPlatforms()):
    print(f"  {openmm.Platform.getPlatform(i).getName()}")
```

You should see platforms like `Reference`, `CPU`, `CUDA`, or `OpenCL`. If you have an NVIDIA GPU and see `CUDA`, you are in good shape for production-scale simulations.

---

## 2. Preparing a Protein System

We will use lysozyme (PDB: 1AKI) as our model protein. It is small (129 residues), well-studied, and ideal for learning.

### 2.1 Loading and Cleaning the PDB

```python
from openmm.app import PDBFile, Modeller, ForceField
from openmm.app import PME, HBonds, NoCutoff
from openmm import unit
import urllib.request
import os

# Download PDB file
pdb_id = "1AKI"
pdb_filename = f"{pdb_id}.pdb"
if not os.path.exists(pdb_filename):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, pdb_filename)
    print(f"Downloaded {pdb_filename}")

# Load the structure
pdb = PDBFile(pdb_filename)
print(f"Loaded {pdb_id}: {pdb.topology.getNumAtoms()} atoms, "
      f"{pdb.topology.getNumResidues()} residues, "
      f"{pdb.topology.getNumChains()} chain(s)")
```

### 2.2 Adding Hydrogens and Solvent

Crystal structures typically lack hydrogen atoms. We also need to solvate the protein in a water box and add counterions to neutralize the system.

```python
# Define force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create modeller from the loaded PDB
modeller = Modeller(pdb.topology, pdb.positions)

# Add missing hydrogens at pH 7.0
modeller.addHydrogens(forcefield, pH=7.0)
print(f"After adding hydrogens: {modeller.topology.getNumAtoms()} atoms")

# Add solvent: cubic box with 1.0 nm padding
# Also add 0.15 M NaCl to mimic physiological ionic strength
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0 * unit.nanometers,
    ionicStrength=0.15 * unit.molar,
    positiveIon='Na+',
    negativeIon='Cl-'
)
print(f"After solvation: {modeller.topology.getNumAtoms()} atoms")
```

The padding parameter controls how much water surrounds the protein on each side. A value of 1.0 nm is standard for most applications. The ionic strength of 0.15 M NaCl approximates physiological salt concentration.

{% include figure.liquid loading="eager" path="assets/img/blog/openmm-md/figure1-system-setup.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. System preparation workflow: the raw crystal structure is protonated, placed in a cubic water box with periodic boundary conditions, and neutralized with counterions.
</div>

### 2.3 Saving the Prepared System

```python
# Save the solvated system for reproducibility
with open('system_solvated.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("Saved solvated system to system_solvated.pdb")
```

---

## 3. Creating the Simulation System

### 3.1 Force Field and System Parameters

The force field defines how atoms interact. We use **AMBER14** with the **TIP3P-FB** water model, a modern combination with good accuracy for protein simulations.

```python
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import Simulation, DCDReporter, StateDataReporter
import sys

# Create the OpenMM System
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,          # Particle Mesh Ewald for long-range electrostatics
    nonbondedCutoff=1.0 * unit.nanometers,  # Short-range cutoff
    constraints=HBonds,           # Constrain bonds involving hydrogen
    hydrogenMass=1.5 * unit.amu   # Hydrogen mass repartitioning for 4 fs timestep
)
print(f"System created: {system.getNumParticles()} particles, "
      f"{system.getNumForces()} force terms")
```

Key choices explained:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `nonbondedMethod` | PME | Accurate long-range electrostatics for periodic systems |
| `nonbondedCutoff` | 1.0 nm | Standard cutoff; balances accuracy and speed |
| `constraints` | HBonds | Allows 2 fs or larger timestep by freezing fast H-bond vibrations |
| `hydrogenMass` | 1.5 amu | Hydrogen mass repartitioning enables 4 fs timestep without loss of accuracy |

### 3.2 Integrator and Barostat

We use a **Langevin integrator** for temperature control and a **Monte Carlo barostat** for pressure control.

```python
# Langevin integrator for NVT / temperature coupling
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
timestep = 4.0 * unit.femtoseconds  # Possible due to hydrogen mass repartitioning

integrator = LangevinMiddleIntegrator(temperature, friction, timestep)

# Monte Carlo barostat for NPT / pressure coupling
pressure = 1.0 * unit.atmospheres
barostat_frequency = 25  # Attempt volume change every 25 steps
barostat = MonteCarloBarostat(pressure, temperature, barostat_frequency)
system.addForce(barostat)
```

{% include figure.liquid loading="eager" path="assets/img/blog/openmm-md/figure2-simulation-workflow.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. Overview of the simulation workflow: energy minimization removes steric clashes, NVT equilibration stabilizes temperature, NPT equilibration stabilizes density, and the production run generates the trajectory for analysis.
</div>

### 3.3 Building the Simulation Object

```python
# Choose the fastest available platform
platform = openmm.Platform.getPlatformByName('CUDA')  # or 'CPU' if no GPU
properties = {'CudaPrecision': 'mixed'}  # Mixed precision: fast and accurate

# For CPU fallback:
# platform = openmm.Platform.getPlatformByName('CPU')
# properties = {}

simulation = Simulation(
    modeller.topology,
    system,
    integrator,
    platform,
    properties
)
simulation.context.setPositions(modeller.positions)
print(f"Simulation built on platform: {simulation.context.getPlatform().getName()}")
```

---

## 4. Energy Minimization

The initial structure from the PDB likely contains steric clashes, especially after adding hydrogens and solvent. Energy minimization resolves these by moving atoms downhill on the potential energy surface.

```python
import time

# Report initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"Initial potential energy: {initial_energy}")

# Run energy minimization
print("Running energy minimization...")
t0 = time.time()
simulation.minimizeEnergy(
    maxIterations=1000,
    tolerance=10.0 * unit.kilojoules_per_mole / unit.nanometer
)
t1 = time.time()

# Report final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"Final potential energy:   {final_energy}")
print(f"Energy reduction:         {initial_energy - final_energy}")
print(f"Minimization took:        {t1 - t0:.1f} seconds")

# Save minimized coordinates
positions_min = state.getPositions()
with open('minimized.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, positions_min, f)
```

You should see the potential energy drop from a large positive value (indicating clashes) to a negative value. If the initial energy is extremely positive (e.g., > 10^6 kJ/mol), something went wrong during system setup -- check for missing atoms or overlapping waters.

---

## 5. Equilibration

Equilibration is performed in two stages: NVT (constant volume, constant temperature) followed by NPT (constant pressure, constant temperature).

### 5.1 NVT Equilibration

The goal is to bring the system to the target temperature while keeping the volume fixed.

```python
# Remove the barostat for NVT phase
# We need to remove and re-add it for the NPT phase
for i in range(system.getNumForces()):
    force = system.getForce(i)
    if isinstance(force, MonteCarloBarostat):
        system.removeForce(i)
        break

# Reinitialize context after modifying the system
simulation.context.reinitialize(preserveState=True)

# Set initial velocities from Maxwell-Boltzmann distribution
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

# Add reporters for NVT equilibration
simulation.reporters.clear()
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        1000,                          # Report every 1000 steps (4 ps)
        step=True,
        temperature=True,
        potentialEnergy=True,
        speed=True,
        remainingTime=True,
        totalSteps=25000
    )
)

# Run 100 ps of NVT equilibration (25,000 steps x 4 fs)
print("Running NVT equilibration (100 ps)...")
simulation.step(25000)
print("NVT equilibration complete.")
```

### 5.2 NPT Equilibration

Now we add the barostat back and allow the box volume to fluctuate, equilibrating the density.

```python
# Re-add the barostat for NPT
barostat = MonteCarloBarostat(
    1.0 * unit.atmospheres,
    300 * unit.kelvin,
    25
)
system.addForce(barostat)
simulation.context.reinitialize(preserveState=True)

# Update reporters for NPT
simulation.reporters.clear()
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        1000,
        step=True,
        temperature=True,
        potentialEnergy=True,
        density=True,
        speed=True,
        remainingTime=True,
        totalSteps=50000
    )
)

# Run 200 ps of NPT equilibration (50,000 steps x 4 fs)
print("Running NPT equilibration (200 ps)...")
simulation.step(50000)
print("NPT equilibration complete.")

# Save equilibrated state
simulation.saveCheckpoint('equilibrated.chk')
print("Saved equilibrated checkpoint.")
```

{% include figure.liquid loading="eager" path="assets/img/blog/openmm-md/figure3-equilibration.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. Equilibration monitoring: (a) temperature stabilizes around 300 K during NVT, (b) density converges to ~1.0 g/cm<sup>3</sup> during NPT, indicating a well-equilibrated system.
</div>

During NPT equilibration, watch for density convergence to approximately 1.0 g/cm^3 for aqueous systems. If the density drifts or oscillates wildly, the system may have issues.

---

## 6. Production Run

The production run generates the trajectory we will analyze. We attach two reporters: `DCDReporter` for coordinates (binary, compact) and `StateDataReporter` for thermodynamic properties.

```python
# Clear previous reporters
simulation.reporters.clear()

# Production parameters
production_steps = 2500000   # 10 ns at 4 fs/step
save_frequency = 2500        # Save frame every 10 ps
log_frequency = 25000        # Log thermodynamics every 100 ps

# DCD trajectory reporter
simulation.reporters.append(
    DCDReporter('trajectory.dcd', save_frequency)
)

# Thermodynamic log
simulation.reporters.append(
    StateDataReporter(
        'production_log.csv',
        log_frequency,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        volume=True,
        speed=True,
        separator=','
    )
)

# Console output for monitoring
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        log_frequency,
        step=True,
        time=True,
        temperature=True,
        speed=True,
        remainingTime=True,
        totalSteps=production_steps
    )
)

# Run production
print(f"Starting production run: {production_steps * 4 / 1e6:.1f} ns")
print(f"Trajectory frames: {production_steps // save_frequency}")
t0 = time.time()
simulation.step(production_steps)
t1 = time.time()

elapsed = t1 - t0
ns_per_day = (production_steps * 4e-6) / (elapsed / 86400)
print(f"\nProduction complete in {elapsed:.0f} seconds")
print(f"Performance: {ns_per_day:.1f} ns/day")

# Save final state
simulation.saveCheckpoint('production_final.chk')
simulation.saveState('production_final.xml')

# Save final coordinates
state = simulation.context.getState(getPositions=True)
with open('production_final.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, state.getPositions(), f)
```

**Typical performance benchmarks** (lysozyme in water, ~25,000 atoms):

| Platform | Approximate Speed |
|----------|------------------|
| CPU (8 cores) | 5--15 ns/day |
| GTX 1080 Ti | 100--200 ns/day |
| RTX 3090 | 200--400 ns/day |
| A100 | 400--800 ns/day |

For a 10 ns simulation on a modern GPU, expect roughly 30 minutes to a few hours.

---

## 7. Trajectory Analysis with MDAnalysis

Now comes the rewarding part: extracting biophysical insights from the trajectory.

### 7.1 Loading the Trajectory

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import matplotlib.pyplot as plt

# Load topology (PDB) and trajectory (DCD)
u = mda.Universe('system_solvated.pdb', 'trajectory.dcd')
print(f"Loaded trajectory: {u.trajectory.n_frames} frames, "
      f"{u.atoms.n_atoms} atoms")
print(f"Time range: {u.trajectory[0].time} to {u.trajectory[-1].time} ps")

# Select protein atoms
protein = u.select_atoms('protein')
backbone = u.select_atoms('protein and backbone')
ca_atoms = u.select_atoms('protein and name CA')
print(f"Protein: {protein.n_atoms} atoms, "
      f"Backbone: {backbone.n_atoms} atoms, "
      f"C-alpha: {ca_atoms.n_atoms} atoms")
```

### 7.2 RMSD: Structural Drift

Root Mean Square Deviation (RMSD) measures how far the protein has moved from its starting structure. A plateau indicates the simulation has reached a stable conformational basin.

```python
# Align trajectory to the first frame using backbone atoms
align.AlignTraj(u, u, select='protein and backbone', in_memory=True).run()

# Calculate backbone RMSD
rmsd_analysis = rms.RMSD(u, u, select='backbone', ref_frame=0)
rmsd_analysis.run()

# Extract results: columns are [frame, time (ps), RMSD (A)]
rmsd_data = rmsd_analysis.results.rmsd
time_ns = rmsd_data[:, 1] / 1000.0  # Convert ps to ns
rmsd_values = rmsd_data[:, 2]       # RMSD in Angstroms

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(time_ns, rmsd_values, linewidth=0.8, color='#2c7bb6')
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Backbone RMSD ($\\AA$)', fontsize=12)
ax.set_title('Backbone RMSD over Simulation Time', fontsize=14)
ax.axhline(y=np.mean(rmsd_values[len(rmsd_values)//2:]),
           color='#d7191c', linestyle='--', linewidth=1,
           label=f'Mean (last half): {np.mean(rmsd_values[len(rmsd_values)//2:]):.2f} $\\AA$')
ax.legend()
ax.set_xlim(0, time_ns[-1])
plt.tight_layout()
plt.savefig('rmsd.png', dpi=150, bbox_inches='tight')
plt.show()
print(f"Mean RMSD (last 5 ns): {np.mean(rmsd_values[len(rmsd_values)//2:]):.2f} A")
```

For a well-behaved simulation of a stable protein like lysozyme, expect RMSD to plateau at 1--2 Angstroms within the first few nanoseconds.

{% include figure.liquid loading="eager" path="assets/img/blog/openmm-md/figure4-rmsd-rmsf.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. (a) Backbone RMSD over time: the plateau after ~2 ns indicates structural convergence. (b) Per-residue RMSF: loop regions (highlighted) show the highest fluctuations, consistent with their known flexibility.
</div>

### 7.3 RMSF: Per-Residue Flexibility

Root Mean Square Fluctuation (RMSF) quantifies the average mobility of each residue. High RMSF values pinpoint flexible loops, termini, and disordered regions.

```python
# Calculate per-residue RMSF using C-alpha atoms
# First, align to average structure
average_pos = align.AverageStructure(u, u, select='protein and name CA',
                                      ref_frame=0).run()
ref = average_pos.results.universe

align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()

# Compute RMSF
rmsf_analysis = rms.RMSF(ca_atoms).run()
rmsf_values = rmsf_analysis.results.rmsf
residue_ids = ca_atoms.resids

# Plot
fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(residue_ids, rmsf_values, linewidth=1.0, color='#2c7bb6')
ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color='#abd9e9')
ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('RMSF ($\\AA$)', fontsize=12)
ax.set_title('Per-Residue C$\\alpha$ RMSF', fontsize=14)

# Highlight flexible regions (RMSF > mean + 1 std)
threshold = np.mean(rmsf_values) + np.std(rmsf_values)
flexible = residue_ids[rmsf_values > threshold]
ax.axhline(y=threshold, color='#d7191c', linestyle='--',
           label=f'Threshold: {threshold:.2f} $\\AA$')
if len(flexible) > 0:
    ax.scatter(flexible, rmsf_values[rmsf_values > threshold],
               color='#d7191c', s=20, zorder=5, label='Flexible residues')
ax.legend()
plt.tight_layout()
plt.savefig('rmsf.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Most flexible residues: {flexible}")
print(f"Mean RMSF: {np.mean(rmsf_values):.2f} A")
```

### 7.4 Radius of Gyration: Compactness

The radius of gyration (R_g) tracks the overall compactness of the protein. Changes in R_g can indicate unfolding, compaction, or large-scale conformational transitions.

```python
# Calculate radius of gyration for each frame
rg_values = []
time_points = []
for ts in u.trajectory:
    rg = protein.radius_of_gyration()
    rg_values.append(rg)
    time_points.append(ts.time / 1000.0)  # ps to ns

rg_values = np.array(rg_values)
time_points = np.array(time_points)

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(time_points, rg_values, linewidth=0.8, color='#1a9641')
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Radius of Gyration ($\\AA$)', fontsize=12)
ax.set_title('Radius of Gyration over Time', fontsize=14)
ax.axhline(y=np.mean(rg_values), color='#d7191c', linestyle='--',
           label=f'Mean: {np.mean(rg_values):.2f} $\\AA$')
ax.legend()
ax.set_xlim(0, time_points[-1])
plt.tight_layout()
plt.savefig('radius_of_gyration.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Mean Rg: {np.mean(rg_values):.2f} +/- {np.std(rg_values):.2f} A")
```

### 7.5 Hydrogen Bond Analysis

Hydrogen bonds are critical for protein secondary structure stability. We can track intra-protein hydrogen bonds over time.

```python
# Hydrogen bond analysis
hbond_analysis = HydrogenBondAnalysis(
    universe=u,
    donors_sel='protein',
    hydrogens_sel='protein',
    acceptors_sel='protein',
    d_a_cutoff=3.0,          # Donor-acceptor distance cutoff (Angstroms)
    d_h_a_angle_cutoff=150   # D-H...A angle cutoff (degrees)
)
hbond_analysis.run()

# Count hydrogen bonds per frame
hbond_counts = hbond_analysis.count_by_time()
hbond_times = np.array([t / 1000.0 for t in hbond_analysis.times])  # ps to ns

# Plot
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(hbond_times, hbond_counts, linewidth=0.8, color='#7b3294', alpha=0.7)

# Running average for clarity
window = min(50, len(hbond_counts) // 5)
if window > 1:
    running_avg = np.convolve(hbond_counts, np.ones(window)/window, mode='valid')
    avg_times = hbond_times[:len(running_avg)]
    ax.plot(avg_times, running_avg, linewidth=2, color='#d7191c',
            label=f'Running avg (n={window})')

ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Number of H-bonds', fontsize=12)
ax.set_title('Intra-Protein Hydrogen Bonds', fontsize=14)
ax.legend()
ax.set_xlim(0, hbond_times[-1])
plt.tight_layout()
plt.savefig('hydrogen_bonds.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Mean H-bonds: {np.mean(hbond_counts):.1f} +/- {np.std(hbond_counts):.1f}")
```

---

## 8. Visualizing Trajectories

### 8.1 Interactive Visualization with nglview

In a Jupyter notebook, nglview provides an interactive 3D viewer:

```python
import nglview as nv

# Create view from MDAnalysis universe
view = nv.show_mdanalysis(protein)

# Styling
view.clear_representations()
view.add_cartoon(selection='protein', color='sstruc')
view.add_ball_and_stick(selection='protein and (resname LYS or resname ASP) and sidechain')
view.add_surface(selection='protein', opacity=0.15, color='white')

# Set background
view.background = 'white'

# Display
view
```

### 8.2 Static Visualization with MDAnalysis

For publication-quality figures without a notebook, you can render snapshots:

```python
# Extract representative frames for visualization
# Save frames at 0%, 25%, 50%, 75%, 100% of simulation
n_frames = u.trajectory.n_frames
frame_indices = [0,
                 n_frames // 4,
                 n_frames // 2,
                 3 * n_frames // 4,
                 n_frames - 1]

for idx in frame_indices:
    u.trajectory[idx]
    frame_time = u.trajectory[idx].time / 1000.0
    protein.write(f'frame_{frame_time:.1f}ns.pdb')
    print(f"Saved frame at {frame_time:.1f} ns")
```

These PDB snapshots can be loaded into PyMOL, VMD, or ChimeraX for high-quality rendering.

{% include figure.liquid loading="eager" path="assets/img/blog/openmm-md/figure5-trajectory-snapshots.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Trajectory snapshots at different time points showing the protein's conformational evolution. Loop regions exhibit the most visible structural changes, while the core beta-sheet and alpha-helical regions remain stable.
</div>

---

## 9. Putting It All Together: Summary Analysis Script

Here is a compact analysis script that generates all key metrics in one pass:

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
import matplotlib.pyplot as plt

def analyze_trajectory(topology, trajectory, output_prefix='analysis'):
    """Complete analysis pipeline for an MD trajectory."""
    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms('protein')
    backbone = u.select_atoms('protein and backbone')
    ca_atoms = u.select_atoms('protein and name CA')

    # Align to first frame
    align.AlignTraj(u, u, select='protein and backbone',
                    in_memory=True).run()

    # --- Collect per-frame data ---
    times, rmsd_vals, rg_vals = [], [], []
    ref_pos = backbone.positions.copy()

    for ts in u.trajectory:
        t = ts.time / 1000.0
        times.append(t)

        # RMSD
        current_rmsd = rms.rmsd(backbone.positions, ref_pos, superposition=False)
        rmsd_vals.append(current_rmsd)

        # Radius of gyration
        rg_vals.append(protein.radius_of_gyration())

    times = np.array(times)
    rmsd_vals = np.array(rmsd_vals)
    rg_vals = np.array(rg_vals)

    # --- RMSF ---
    rmsf_analysis = rms.RMSF(ca_atoms).run()
    rmsf_vals = rmsf_analysis.results.rmsf

    # --- Multi-panel figure ---
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=False)

    # RMSD
    axes[0].plot(times, rmsd_vals, lw=0.8, color='#2c7bb6')
    axes[0].set_ylabel('RMSD ($\\AA$)')
    axes[0].set_title('Backbone RMSD')
    axes[0].set_xlabel('Time (ns)')

    # Rg
    axes[1].plot(times, rg_vals, lw=0.8, color='#1a9641')
    axes[1].set_ylabel('Rg ($\\AA$)')
    axes[1].set_title('Radius of Gyration')
    axes[1].set_xlabel('Time (ns)')

    # RMSF
    axes[2].bar(ca_atoms.resids, rmsf_vals, width=1.0,
                color='#fdae61', edgecolor='#d7191c', linewidth=0.3)
    axes[2].set_ylabel('RMSF ($\\AA$)')
    axes[2].set_xlabel('Residue Number')
    axes[2].set_title('Per-Residue RMSF')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_summary.png', dpi=200, bbox_inches='tight')
    plt.show()

    # --- Print summary statistics ---
    half = len(rmsd_vals) // 2
    print(f"\n{'='*50}")
    print(f"TRAJECTORY ANALYSIS SUMMARY")
    print(f"{'='*50}")
    print(f"Frames:           {u.trajectory.n_frames}")
    print(f"Simulation time:  {times[-1]:.1f} ns")
    print(f"RMSD (last half): {np.mean(rmsd_vals[half:]):.2f} +/- "
          f"{np.std(rmsd_vals[half:]):.2f} A")
    print(f"Rg (mean):        {np.mean(rg_vals):.2f} +/- "
          f"{np.std(rg_vals):.2f} A")
    print(f"RMSF (mean):      {np.mean(rmsf_vals):.2f} A")
    print(f"RMSF (max):       {np.max(rmsf_vals):.2f} A "
          f"(residue {ca_atoms.resids[np.argmax(rmsf_vals)]})")
    print(f"{'='*50}")

    return times, rmsd_vals, rg_vals, rmsf_vals

# Usage
times, rmsd, rg, rmsf = analyze_trajectory('system_solvated.pdb', 'trajectory.dcd')
```

---

## 10. Application Context

### IDR Dynamics

Intrinsically disordered regions (IDRs) do not adopt a single folded structure, making MD an essential tool for characterizing their conformational ensembles. Key considerations for IDR simulations:

- Use longer simulation times (hundreds of nanoseconds to microseconds) since IDRs sample many states.
- Force field choice is critical: AMBER ff19SB with OPC water, or CHARMM36m, have been validated for disordered proteins.
- Enhanced sampling methods (replica exchange MD, metadynamics) may be necessary to adequately sample the conformational landscape.

```python
# Example: monitoring end-to-end distance for an IDR
# Useful for characterizing the degree of compaction
idr_nterm = u.select_atoms('protein and name CA and resid 1')
idr_cterm = u.select_atoms('protein and name CA and resid 129')

e2e_distances = []
for ts in u.trajectory:
    dist = np.linalg.norm(idr_nterm.positions[0] - idr_cterm.positions[0])
    e2e_distances.append(dist)

print(f"End-to-end distance: {np.mean(e2e_distances):.1f} +/- "
      f"{np.std(e2e_distances):.1f} A")
```

### Ligand Binding Studies

For drug discovery applications, MD simulations can reveal:

- **Binding pose stability**: Run MD starting from a docked pose and check if the ligand remains bound.
- **Binding free energies**: Use alchemical free energy methods (FEP/TI) built on MD.
- **Cryptic binding sites**: Transient pockets that open during the simulation.

```python
# Example: monitoring ligand-protein distance during simulation
# (assuming a ligand is present in the system)
# ligand = u.select_atoms('resname LIG')
# binding_site = u.select_atoms('protein and around 5.0 resname LIG')
#
# min_distances = []
# for ts in u.trajectory:
#     dists = mda.lib.distances.distance_array(
#         ligand.positions, binding_site.positions
#     )
#     min_distances.append(np.min(dists))
```

### Protein Stability

Comparing MD simulations of wild-type and mutant proteins reveals how mutations affect dynamics:

- Increased RMSD or RMSF near the mutation site suggests destabilization.
- Loss of hydrogen bonds or salt bridges indicates weakened interactions.
- Changes in R_g can indicate partial unfolding.

---

## 11. Key Takeaways

**What we covered:**

1. **System preparation**: Loading a PDB, adding hydrogens at physiological pH, solvating in a water box with counterions.
2. **Simulation setup**: Choosing force fields (AMBER14), integrators (Langevin), and barostats for temperature and pressure control.
3. **Running MD**: Energy minimization, NVT/NPT equilibration, and production runs with appropriate reporters.
4. **Analysis**: RMSD, RMSF, radius of gyration, and hydrogen bonds using MDAnalysis.
5. **Visualization**: nglview for interactive exploration and snapshot export for publication figures.

**Computational resources:**

| Simulation Goal | System Size | Time | GPU Hours (RTX 3090) |
|----------------|-------------|------|---------------------|
| Quick test | ~25K atoms | 10 ns | ~1 hour |
| Standard analysis | ~25K atoms | 100 ns | ~8 hours |
| Thorough sampling | ~25K atoms | 1 us | ~3 days |
| Protein-ligand FEP | ~50K atoms | 5 ns x 12 windows | ~6 hours |

**Common pitfalls and solutions:**

- **Simulation blows up (NaN energies)**: Usually caused by bad initial geometry. Always minimize before dynamics. Reduce the timestep if problems persist.
- **RMSD never plateaus**: Either the simulation is too short, or the protein is undergoing a genuine conformational change. Extend the simulation and check for artifacts.
- **Density far from 1.0 g/cm^3**: The water model or barostat may be misconfigured. Verify force field files and barostat parameters.
- **Periodic image artifacts**: Protein may diffuse across the periodic boundary. Use MDAnalysis `transformations.unwrap` or `make_whole` to fix coordinates before analysis.

```python
# Fix periodic boundary artifacts before analysis
import MDAnalysis.transformations as trans

u = mda.Universe('system_solvated.pdb', 'trajectory.dcd')
protein = u.select_atoms('protein')

# Define transformations pipeline
workflow = [
    trans.unwrap(u.atoms),
    trans.center_in_box(protein, wrap=True),
    trans.wrap(u.atoms)
]
u.trajectory.add_transformations(*workflow)
```

**Further reading:**

- OpenMM documentation: [openmm.org](http://openmm.org)
- MDAnalysis User Guide: [userguide.mdanalysis.org](https://userguide.mdanalysis.org)
- Best practices for MD: [Living Journal of Computational Molecular Science (LIVECOMP)](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5957)
- AMBER force field parameters: [ambermd.org](https://ambermd.org)

This tutorial provides a foundation. Real research simulations typically require longer timescales, enhanced sampling techniques, and careful validation against experimental data (NMR order parameters, SAXS profiles, crystallographic B-factors). Start with the workflow here, verify it on a known system, and then adapt it to your research question.
