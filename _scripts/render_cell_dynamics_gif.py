from pathlib import Path
import math

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "assets" / "img" / "cellular_dynamics_perturbation"
GIF = OUT_DIR / "set_mmd_flow_dynamics.gif"
FRAMES = OUT_DIR / "dynamics_frames"

WIDTH = 560
HEIGHT = 360
FPS_MS = 70
N_FRAMES = 86
N_CELLS = 170


def ease(t):
    return 0.5 - 0.5 * math.cos(math.pi * t)


def make_cells(seed=7):
    rng = np.random.default_rng(seed)
    n1 = int(N_CELLS * 0.68)
    n2 = N_CELLS - n1

    cluster1 = rng.normal([0.0, 0.0], [0.42, 0.27], size=(n1, 2))
    cluster2 = rng.normal([0.75, 0.42], [0.22, 0.18], size=(n2, 2))
    ctrl = np.vstack([cluster1, cluster2])

    additive = ctrl + np.array([1.35, 0.58])

    rel = ctrl - np.array([0.18, 0.08])
    theta = np.deg2rad(28)
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    warped = rel @ rot.T
    scale = np.array([1.22, 0.62])
    nonlinear = warped * scale + np.array([1.34, 0.78])
    branch = 0.45 / (1 + np.exp(-4.2 * (ctrl[:, 0] - 0.34)))
    nonlinear[:, 1] += branch
    nonlinear[:, 0] += 0.22 * np.sin(ctrl[:, 1] * 3.0)
    return ctrl, additive, nonlinear


def draw_panel(ax, pts, target, title, subtitle, color):
    ax.set_xlim(-1.15, 2.65)
    ax.set_ylim(-0.95, 2.15)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#fbfbfb")
    for spine in ax.spines.values():
        spine.set_color("#dddddd")

    target_mean = target.mean(axis=0)
    target_cov = np.cov(target.T)
    vals, vecs = np.linalg.eigh(target_cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    for scale, alpha in [(2.8, 0.08), (1.8, 0.15)]:
        ell = Ellipse(
            target_mean,
            width=scale * np.sqrt(vals[0]),
            height=scale * np.sqrt(vals[1]),
            angle=angle,
            facecolor=color,
            edgecolor=color,
            lw=1.5,
            alpha=alpha,
        )
        ax.add_patch(ell)

    ax.scatter(target[:, 0], target[:, 1], s=10, c=color, alpha=0.13, linewidths=0)
    ax.scatter(pts[:, 0], pts[:, 1], s=13, c=color, alpha=0.86, edgecolors="white", linewidths=0.25)
    ax.text(0.03, 0.95, title, transform=ax.transAxes, ha="left", va="top", fontsize=11, weight="bold", color="#222222")
    ax.text(0.03, 0.86, subtitle, transform=ax.transAxes, ha="left", va="top", fontsize=8.4, color="#555555")


def render_frame(path, t, ctrl, additive, nonlinear):
    phase = ease(min(max(t, 0), 1))
    add_pts = ctrl + phase * (additive - ctrl)
    flow_pts = ctrl + phase * (nonlinear - ctrl)

    fig = plt.figure(figsize=(WIDTH / 100, HEIGHT / 100), dpi=100)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, 2, wspace=0.08, left=0.04, right=0.98, top=0.86, bottom=0.11)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])

    draw_panel(ax_left, add_pts, additive, "Additive baseline", "same shift for every cell", "#6f7d8f")
    draw_panel(ax_right, flow_pts, nonlinear, "Set-MMD Flow", "cell-aware distribution shift", "#2474d8")

    for ax, target in [(ax_left, additive), (ax_right, nonlinear)]:
        if 0.18 < t < 0.86:
            idx = np.linspace(0, N_CELLS - 1, 16, dtype=int)
            start = ctrl[idx]
            end = target[idx]
            cur = start + phase * (end - start)
            delta = 0.18 * (end - start)
            ax.quiver(cur[:, 0], cur[:, 1], delta[:, 0], delta[:, 1], angles="xy", scale_units="xy", scale=1, width=0.004, alpha=0.42, color="#222222")

    fig.text(0.5, 0.955, "Population-aware perturbation prediction", ha="center", va="top", fontsize=13, weight="bold", color="#111111")
    fig.text(0.5, 0.055, "Control cells move toward predicted perturbed distributions", ha="center", va="bottom", fontsize=8.5, color="#555555")

    fig.savefig(path, dpi=100)
    plt.close(fig)


def build_gif():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if FRAMES.exists():
        for frame in FRAMES.glob("*.png"):
            frame.unlink()
    else:
        FRAMES.mkdir(parents=True)

    ctrl, additive, nonlinear = make_cells()
    for i in range(N_FRAMES):
        t = i / (N_FRAMES - 1)
        if i > N_FRAMES - 12:
            t = 1.0
        render_frame(FRAMES / f"frame_{i:03d}.png", t, ctrl, additive, nonlinear)

    frames = [Image.open(p).convert("RGB") for p in sorted(FRAMES.glob("*.png"))]
    frames[0].save(
        GIF,
        save_all=True,
        append_images=frames[1:],
        optimize=False,
        duration=FPS_MS,
        loop=0,
        disposal=2,
    )

    for frame in FRAMES.glob("*.png"):
        frame.unlink()
    FRAMES.rmdir()
    print(f"Wrote {GIF}")


if __name__ == "__main__":
    build_gif()
