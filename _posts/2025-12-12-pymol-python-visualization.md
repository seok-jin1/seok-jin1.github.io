---
layout: post
title: "PyMOL + Python: Publication-Quality Protein Structure Visualization"
date: 2025-12-12
permalink: /blog/pymol-python-visualization/
published: true
categories: [tutorial]
tags:
  - structural-biology
  - visualization
  - python
  - tutorial
---

Structural biology papers live and die by their figures. **PyMOL** remains the gold standard for molecular visualization, and its Python API turns it into a fully scriptable rendering engine. This tutorial walks through practical recipes for generating publication-quality figures entirely from Python scripts -- no GUI clicking required.

## Why PyMOL?

- **Scriptability**: Every operation maps to a Python function call, making figures 100% reproducible.
- **Ray tracing**: Built-in ray tracer produces smooth, anti-aliased images at arbitrary resolution.
- **Community**: Thousands of published papers use PyMOL, and reviewers expect its visual style.
- **Flexibility**: From simple cartoons to complex multi-panel assemblies, one tool handles it all.

## 1. Setup: Installing PyMOL and pymol2

```bash
# Conda (recommended)
conda create -n pymol-env python=3.10
conda activate pymol-env
conda install -c conda-forge pymol-open-source

# Verify pymol2 API
python -c "import pymol2; print('pymol2 ready')"

# For headless rendering on a server
sudo apt-get install xvfb
xvfb-run -a python my_pymol_script.py
```

Every PyMOL Python script follows this skeleton:

```python
import pymol2

with pymol2.PyMOL() as p:
    cmd = p.cmd
    cmd.load("structure.pdb", "my_protein")
    cmd.show("cartoon", "my_protein")
    cmd.color("marine", "my_protein")
    cmd.ray(2400, 2400)
    cmd.png("output.png", dpi=300)
```

## 2. Basic Commands: Loading, Selecting, Coloring, Representations

```python
# --- Loading ---
cmd.load("6lzg.pdb", "spike_ace2")   # Local file
cmd.fetch("6LZG", name="spike_ace2") # Fetch from PDB

# --- Selections ---
cmd.select("chain_A", "chain A")
cmd.select("binding_loop", "resi 480-505 and chain A")
cmd.select("interface", "chain A within 4.0 of chain B")
cmd.select("polar_if", "interface and (resn ASP+GLU+LYS+ARG+HIS)")

# --- Representations ---
cmd.hide("everything", "all")
cmd.show("cartoon", "all")            # Ribbon diagram
cmd.show("surface", "chain_A")        # Molecular surface
cmd.show("sticks", "binding_loop")    # Sidechain sticks
cmd.show("spheres", "element Zn")     # Space-filling for ions

# --- Coloring ---
cmd.color("marine", "chain_A")
cmd.color("salmon", "chain_B")
cmd.set_color("custom_blue", [0.2, 0.4, 0.8])
cmd.spectrum("b", "blue_white_red", "all")  # Color by B-factor
```

## 3. Python Scripting: Key `cmd` Functions

The core API functions you will use in every script:

```python
cmd.load(filename, object_name)     # Load a structure
cmd.select(name, selection_expr)    # Create a named selection
cmd.color(color, selection)         # Apply color
cmd.show(representation, selection) # Show representation
cmd.hide(representation, selection) # Hide representation
cmd.set(setting, value)             # Change a global setting
cmd.orient(selection)               # Auto-orient the camera
cmd.zoom(selection, buffer)         # Zoom to fit selection
cmd.rotate(axis, angle)             # Rotate the view
cmd.get_view()                      # Get current view matrix (18-tuple)
cmd.set_view(view_tuple)            # Restore a saved view
cmd.ray(width, height)              # Ray trace at given resolution
cmd.png(filename, dpi=300)          # Save image
cmd.save(filename)                  # Save session (.pse) or structure
```

## 4. Recipe 1: Color by pLDDT (AlphaFold Confidence)

AlphaFold stores per-residue confidence (pLDDT) in the B-factor column. This recipe creates the standard confidence color scheme.

{% include figure.liquid loading="eager" path="assets/img/blog/pymol-python/figure1-plddt-coloring.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 1. AlphaFold model colored by pLDDT confidence. Dark blue: very high (>90), light blue: high (70-90), yellow: low (50-70), orange: very low (<50).
</div>

```python
import pymol2

def render_plddt(pdb_path, output_path):
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(pdb_path, "af_model")

        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)
        cmd.set("antialias", 2)
        cmd.set("ray_shadows", 0)

        cmd.hide("everything", "all")
        cmd.show("cartoon", "af_model")

        # Standard AlphaFold pLDDT palette
        cmd.set_color("plddt_vhigh", [0.051, 0.341, 0.827])  # >90
        cmd.set_color("plddt_high",  [0.416, 0.796, 0.945])  # 70-90
        cmd.set_color("plddt_low",   [1.000, 0.859, 0.071])  # 50-70
        cmd.set_color("plddt_vlow",  [1.000, 0.494, 0.271])  # <50

        cmd.color("plddt_vlow",  "b < 50")
        cmd.color("plddt_low",   "b >= 50 and b < 70")
        cmd.color("plddt_high",  "b >= 70 and b < 90")
        cmd.color("plddt_vhigh", "b >= 90")

        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_oval_length", 1.2)
        cmd.orient("af_model")
        cmd.zoom("af_model", buffer=5)
        cmd.ray(2400, 2400)
        cmd.png(output_path, dpi=300)

render_plddt("AF-P0DTC2-F1-model_v4.pdb", "plddt_colored.png")
```

### Batch Processing Multiple Models

If you have many AlphaFold predictions to render:

```python
import glob

pdb_files = glob.glob("alphafold_models/*.pdb")
for pdb in pdb_files:
    name = pdb.split("/")[-1].replace(".pdb", "")
    render_plddt(pdb, f"figures/{name}_plddt.png")
```

## 5. Recipe 2: Highlight a Binding Interface

Showing two interacting chains with interface residues as sticks and hydrogen bonds as dashed lines.

{% include figure.liquid loading="eager" path="assets/img/blog/pymol-python/figure2-binding-interface.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 2. Spike RBD (blue) bound to ACE2 (coral). Interface residues within 4 angstroms are shown as sticks with polar contacts as dashed lines.
</div>

```python
import pymol2

def render_interface(pdb_path, chain_a, chain_b, output_path, cutoff=4.0):
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(pdb_path, "complex")
        cmd.remove("solvent")
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)

        # Define interface selections
        cmd.select("if_A", f"byres (chain {chain_a} within {cutoff} of chain {chain_b})")
        cmd.select("if_B", f"byres (chain {chain_b} within {cutoff} of chain {chain_a})")

        # Representations
        cmd.hide("everything", "all")
        cmd.show("cartoon", "all")
        cmd.show("sticks", "if_A and sidechain")
        cmd.show("sticks", "if_B and sidechain")

        # Colors
        cmd.color("marine", f"chain {chain_a}")
        cmd.color("salmon", f"chain {chain_b}")
        cmd.color("slate", "if_A and sidechain")
        cmd.color("lightorange", "if_B and sidechain")

        # Polar contacts (H-bonds)
        cmd.distance("hbonds", "if_A", "if_B", cutoff=3.5, mode=2)
        cmd.set("dash_color", "gray40", "hbonds")
        cmd.set("dash_width", 2.0, "hbonds")
        cmd.hide("labels", "hbonds")

        cmd.set("cartoon_transparency", 0.15)
        cmd.set("stick_radius", 0.15)
        cmd.orient("if_A or if_B")
        cmd.zoom("complex", buffer=8)
        cmd.ray(2400, 1800)
        cmd.png(output_path, dpi=300)

render_interface("6lzg.pdb", "A", "B", "binding_interface.png")
```

To add residue labels at key positions:

```python
for chain, resi in [("A", 31), ("A", 353), ("B", 501)]:
    cmd.select(f"lbl_{chain}{resi}", f"chain {chain} and resi {resi} and name CA")
    cmd.label(f"lbl_{chain}{resi}", "'%s%s' % (resn, resi)")
cmd.set("label_size", 14)
cmd.set("label_font_id", 7)
```

## 6. Recipe 3: Electrostatic Surface Potential (APBS)

Electrostatic surfaces reveal charge distributions critical for understanding molecular recognition.

{% include figure.liquid loading="eager" path="assets/img/blog/pymol-python/figure3-electrostatics.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 3. Electrostatic surface potential of ACE2. Red = negative, blue = positive, white = neutral. The binding groove shows complementary negative charge to the spike RBD.
</div>

```bash
# Prerequisites
conda install -c conda-forge apbs
pip install pdb2pqr
```

```python
import pymol2, subprocess, os

def render_electrostatics(pdb_path, output_path, apbs_range=5.0):
    base = os.path.splitext(pdb_path)[0]

    # Step 1: PDB -> PQR (adds charges and radii)
    subprocess.run([
        "pdb2pqr", "--ff", "AMBER", "--apbs-input", f"{base}.in",
        "--keep-chain", pdb_path, f"{base}.pqr"
    ], check=True)

    # Step 2: Run APBS calculation
    subprocess.run(["apbs", f"{base}.in"], check=True)

    # Step 3: Visualize in PyMOL
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(pdb_path, "protein")
        cmd.remove("solvent")
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)

        cmd.load(f"{base}.dx", "potential_map")
        cmd.hide("everything", "all")
        cmd.show("surface", "protein")

        cmd.ramp_new("e_ramp", "potential_map",
                     [-apbs_range, 0, apbs_range],
                     ["red", "white", "blue"])
        cmd.set("surface_color", "e_ramp", "protein")
        cmd.set("surface_quality", 2)

        cmd.orient("protein")
        cmd.zoom("protein", buffer=5)
        cmd.ray(2400, 2400)
        cmd.png(output_path, dpi=300)

render_electrostatics("6lzg_chainA.pdb", "electrostatics.png")
```

## 7. Recipe 4: Multi-Panel Figure Assembly

Journal figures often require multiple views. Script consistent viewpoints with `cmd.get_view()` / `cmd.set_view()`.

{% include figure.liquid loading="eager" path="assets/img/blog/pymol-python/figure4-multi-panel.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 4. Multi-panel figure: (A) overview, (B) interface close-up, (C) 90-degree rotation, (D) surface representation. All panels share identical lighting and style.
</div>

```python
import pymol2, os

def render_panels(pdb_path, out_dir, w=1800, h=1800):
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(pdb_path, "complex")
        cmd.remove("solvent")

        # Shared style
        for setting, val in [("ray_opaque_background", 0), ("antialias", 2),
            ("ray_shadows", 1), ("depth_cue", 0), ("fog", 0),
            ("cartoon_fancy_helices", 1)]:
            cmd.set(setting, val)
        cmd.bg_color("white")
        cmd.hide("everything")
        cmd.show("cartoon", "all")
        cmd.color("marine", "chain A")
        cmd.color("salmon", "chain B")

        # Panel A: Overview
        cmd.orient("complex")
        view_a = cmd.get_view()
        cmd.ray(w, h)
        cmd.png(f"{out_dir}/panel_A.png", dpi=300)

        # Panel B: Interface close-up
        cmd.select("iface", "byres (chain A within 4.0 of chain B)")
        cmd.show("sticks", "iface and sidechain")
        cmd.orient("iface")
        cmd.ray(w, h)
        cmd.png(f"{out_dir}/panel_B.png", dpi=300)
        cmd.hide("sticks", "all")

        # Panel C: Rotated view
        cmd.set_view(view_a)
        cmd.rotate("y", 90)
        cmd.ray(w, h)
        cmd.png(f"{out_dir}/panel_C.png", dpi=300)

        # Panel D: Surface
        cmd.set_view(view_a)
        cmd.hide("cartoon", "all")
        cmd.show("surface", "all")
        cmd.set("transparency", 0.1)
        cmd.ray(w, h)
        cmd.png(f"{out_dir}/panel_D.png", dpi=300)

os.makedirs("panels", exist_ok=True)
render_panels("6lzg.pdb", "panels")
```

Assemble panels programmatically with Pillow:

```python
from PIL import Image, ImageDraw, ImageFont

def assemble(panel_paths, output, cols=2, pad=20):
    panels = [Image.open(p) for p in panel_paths]
    pw, ph = panels[0].size
    rows = (len(panels) + cols - 1) // cols
    canvas = Image.new("RGBA",
        (cols * pw + (cols+1)*pad, rows * ph + (rows+1)*pad),
        (255, 255, 255, 255))

    draw = ImageDraw.Draw(canvas)
    for i, panel in enumerate(panels):
        r, c = divmod(i, cols)
        x, y = pad + c*(pw+pad), pad + r*(ph+pad)
        canvas.paste(panel, (x, y))
        draw.text((x+10, y+5), "ABCD"[i], fill="black")

    canvas.save(output, dpi=(300, 300))

assemble(
    [f"panels/panel_{x}.png" for x in "ABCD"],
    "figure1_assembled.png"
)
```

Save view matrices to JSON for version control:

```python
import json
views = {"overview": list(cmd.get_view())}
with open("views.json", "w") as f:
    json.dump(views, f)
```

## 8. Tips for Journal-Quality Figures

{% include figure.liquid loading="eager" path="assets/img/blog/pymol-python/figure5-quality-comparison.png" class="img-fluid rounded z-depth-1" zoomable=true %}
<div class="caption">
    Figure 5. Left: default PyMOL (no ray tracing). Right: publication settings with ray tracing, anti-aliasing, and optimized lighting.
</div>

### Resolution and Ray Tracing

```python
# Single-column (3.5"): 3.5 * 300 = 1050px minimum; use 2x
cmd.ray(2100, 2100)
cmd.png("figure.png", dpi=300)

# Production rendering settings
cmd.set("ray_trace_mode", 1)
cmd.set("ray_shadows", 1)
cmd.set("antialias", 2)
cmd.set("spec_reflect", 0.2)
cmd.set("spec_power", 250)
cmd.set("ambient", 0.3)
cmd.set("direct", 0.7)
cmd.set("depth_cue", 0)
cmd.set("fog", 0)
```

### Transparent Backgrounds

```python
cmd.set("ray_opaque_background", 0)
cmd.png("transparent.png", dpi=300)  # PNG supports alpha; JPEG does not
```

### Colorblind-Friendly Palette

```python
# Wong (2011) palette
PALETTE = {
    "blue": [0.000, 0.447, 0.698], "orange": [0.902, 0.624, 0.000],
    "green": [0.000, 0.620, 0.451], "pink": [0.800, 0.475, 0.655],
    "sky": [0.337, 0.706, 0.914], "red": [0.835, 0.369, 0.000],
}
for name, rgb in PALETTE.items():
    cmd.set_color(f"cb_{name}", rgb)
```

### Common Pitfalls

```python
# Slow ray tracing? Use draft settings during development
cmd.set("surface_quality", 0)
cmd.ray(800, 800)

# Colors differ between screen and ray trace?
cmd.set("ray_trace_mode", 1)  # Matches OpenGL preview

# Script works in GUI but fails headless?
# Use: xvfb-run -a python script.py
# Or:  export QT_QPA_PLATFORM=offscreen
```

### Label Font Sizes

```python
# For a 3.5"-wide figure rendered at 2100px (600 px/inch):
cmd.set("label_size", 24)
cmd.set("label_font_id", 7)       # Bold sans-serif
cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_transparency", 0.4)
```

### Exporting Sessions

```python
cmd.save("figure_session.pse")  # Editable session for collaborators
cmd.save("figure.wrl")          # VRML for vector rendering
```

## 9. Key Takeaways

- **Script everything.** A Python script is a reproducible figure protocol. When reviewers ask for changes, re-run a script instead of clicking through a GUI.
- **Use `pymol2.PyMOL()` context managers** for clean, self-contained rendering sessions.
- **Store `cmd.get_view()` matrices** in version control. Views are the hardest thing to recreate.
- **Ray trace at 2x your target resolution** and let the journal downsample.
- **Use transparent backgrounds** (`ray_opaque_background`, 0) for easy compositing.
- **Adopt a colorblind-friendly palette** from the start -- retrofitting is painful.
- **Save `.pse` session files** alongside final figures so collaborators can modify the exact scene.

Every structural biology figure should be one `python render_figure.py` away from regeneration.
