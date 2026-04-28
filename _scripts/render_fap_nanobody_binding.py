from pathlib import Path
import math
import shutil

from pymol import cmd


ROOT = Path(__file__).resolve().parents[1]
STRUCTURE = Path("/home/laugh/FAP_nanobody/fold_anti_fap_nb_and_hfap_2/fold_anti_fap_nb_and_hfap_2_model_0.cif")
OUT_DIR = ROOT / "assets" / "img" / "nanobody_humanization_library"
FRAMES = OUT_DIR / "binding_frames"
GIF = OUT_DIR / "fap_nanobody_binding.gif"

WIDTH = 560
HEIGHT = 360
APPROACH_FRAMES = 42
ZOOM_FRAMES = 18
ROTATE_FRAMES = 28
START_DISTANCE = 22.0

NB_INTERFACE = "31+52+53+54+55+56+57+58+98+99+100+101+102+103+104+105+106+107+108+110+111+114"
CDR1 = "31"
CDR2 = "52+53+54+55+56+57+58"
CDR3 = "98+99+100+101+102+103+104+105+106+107+108+110+111+114"
HOTSPOTS = "101+102+103+104+105+106"
FAP_INTERFACE = "575+591+593+594+595+597+598+601+605+632+636+637+638+675+676+679+680+681+682+684+685+687+688"


def center(selection):
    model = cmd.get_model(selection)
    n = len(model.atom)
    if n == 0:
        raise RuntimeError(f"No atoms found for selection: {selection}")
    xyz = [0.0, 0.0, 0.0]
    for atom in model.atom:
        xyz[0] += atom.coord[0]
        xyz[1] += atom.coord[1]
        xyz[2] += atom.coord[2]
    return [v / n for v in xyz]


def norm(vec):
    length = math.sqrt(sum(v * v for v in vec))
    if length == 0:
        return [1.0, 0.0, 0.0]
    return [v / length for v in vec]


def ease(t):
    return 1 - (1 - t) ** 3


def interpolate_view(start, end, t):
    return tuple(start[i] + (end[i] - start[i]) * t for i in range(len(start)))


def prepare_scene(show_interface=False):
    cmd.hide("everything")
    cmd.show("cartoon", "fap")
    cmd.show("cartoon", "nb")

    cmd.color("gray80", "fap")
    cmd.color("marine", "nb")

    if show_interface:
        cmd.select("nb_interface", f"nb and chain B and resi {NB_INTERFACE}")
        cmd.select("cdr1", f"nb and chain B and resi {CDR1}")
        cmd.select("cdr2", f"nb and chain B and resi {CDR2}")
        cmd.select("cdr3", f"nb and chain B and resi {CDR3}")
        cmd.select("hotspots", f"nb and chain B and resi {HOTSPOTS}")
        cmd.select("fap_interface", f"fap and chain A and resi {FAP_INTERFACE}")

        cmd.color("cyan", "fap_interface")
        cmd.color("blue", "cdr1")
        cmd.color("yellow", "cdr2")
        cmd.color("red", "cdr3")
        cmd.color("orange", "hotspots")
        cmd.show("sticks", "nb_interface or fap_interface")

        cmd.distance("contact1", "(nb and chain B and resi 111 and name OD2)", "(fap and chain A and resi 688 and name OH)")
        cmd.distance("contact2", "(nb and chain B and resi 110 and name NH2)", "(fap and chain A and resi 687 and name OE2)")
        cmd.distance("contact3", "(nb and chain B and resi 52 and name OG1)", "(fap and chain A and resi 598 and name OE1)")
        cmd.hide("labels", "contact*")
        cmd.color("yellow", "contact*")
        cmd.set("dash_width", 2.2)
        cmd.set("dash_gap", 0.25)

    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("surface_quality", 1)
    cmd.set("stick_radius", 0.16)
    cmd.set("ambient", 0.55)
    cmd.set("specular", 0.25)
    cmd.set("shininess", 18)
    cmd.set("depth_cue", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")


def render_frame(idx, offset, axis, view, show_interface=False, spin_degrees=0):
    cmd.delete("nb")
    cmd.create("nb", "nb_bound")
    if offset:
        cmd.translate([axis[0] * offset, axis[1] * offset, axis[2] * offset], "nb", camera=0)
    prepare_scene(show_interface=show_interface)
    cmd.set_view(view)
    if spin_degrees:
        cmd.turn("y", spin_degrees)
    cmd.png(str(FRAMES / f"frame_{idx:03d}.png"), width=WIDTH, height=HEIGHT, ray=1)


def build_gif():
    from PIL import Image, ImageChops

    frame_paths = sorted(FRAMES.glob("frame_*.png"))
    raw = [Image.open(path).convert("RGB") for path in frame_paths]

    union = None
    white = Image.new("RGB", raw[0].size, "white")
    for img in raw:
        diff = ImageChops.difference(img, white).convert("L")
        mask = diff.point(lambda px: 255 if px > 15 else 0)
        bbox = mask.getbbox()
        if bbox:
            union = bbox if union is None else (
                min(union[0], bbox[0]),
                min(union[1], bbox[1]),
                max(union[2], bbox[2]),
                max(union[3], bbox[3]),
            )

    if union is None:
        union = (0, 0, raw[0].width, raw[0].height)

    margin = 22
    x0 = max(union[0] - margin, 0)
    y0 = max(union[1] - margin, 0)
    x1 = min(union[2] + margin, raw[0].width)
    y1 = min(union[3] + margin, raw[0].height)

    aspect = WIDTH / HEIGHT
    crop_w = x1 - x0
    crop_h = y1 - y0
    if crop_w / crop_h > aspect:
        target_h = crop_w / aspect
        pad = (target_h - crop_h) / 2
        y0 = max(int(y0 - pad), 0)
        y1 = min(int(y1 + pad), raw[0].height)
    else:
        target_w = crop_h * aspect
        pad = (target_w - crop_w) / 2
        x0 = max(int(x0 - pad), 0)
        x1 = min(int(x1 + pad), raw[0].width)

    frames = []
    for img in raw:
        cropped = img.crop((x0, y0, x1, y1)).resize((WIDTH, HEIGHT), Image.Resampling.LANCZOS)
        frames.append(cropped.convert("P", palette=Image.Palette.ADAPTIVE, colors=128))

    frames[0].save(
        GIF,
        save_all=True,
        append_images=frames[1:],
        optimize=True,
        duration=70,
        loop=0,
        disposal=2,
    )


def main():
    if not STRUCTURE.exists():
        raise FileNotFoundError(STRUCTURE)

    if FRAMES.exists():
        shutil.rmtree(FRAMES)
    FRAMES.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cmd.reinitialize()
    cmd.set("retain_order", 1)
    cmd.load(str(STRUCTURE), "complex")
    cmd.create("fap", "complex and chain A")
    cmd.create("nb_bound", "complex and chain B")
    cmd.delete("complex")
    cmd.remove("solvent")

    fap_center = center("fap")
    nb_center = center("nb_bound")
    axis = norm([nb_center[i] - fap_center[i] for i in range(3)])

    cmd.create("nb", "nb_bound")
    cmd.translate([axis[0] * START_DISTANCE, axis[1] * START_DISTANCE, axis[2] * START_DISTANCE], "nb", camera=0)
    prepare_scene(show_interface=False)
    cmd.orient("fap or nb")
    cmd.turn("x", -18)
    cmd.turn("y", 28)
    cmd.zoom("fap or nb", 3)
    view = cmd.get_view()
    cmd.delete("nb")

    cmd.create("nb", "nb_bound")
    prepare_scene(show_interface=True)
    cmd.set_view(view)
    cmd.zoom("nb_interface or fap_interface", 4.5)
    close_view = cmd.get_view()
    cmd.delete("nb")

    frame_idx = 0
    for i in range(APPROACH_FRAMES):
        t = i / (APPROACH_FRAMES - 1)
        if t < 0.82:
            local_t = t / 0.82
            offset = START_DISTANCE - ease(local_t) * (START_DISTANCE - 3.2)
        else:
            local_t = (t - 0.82) / 0.18
            offset = 3.2 * (1 - local_t)
        render_frame(frame_idx, offset, axis, view, show_interface=i > APPROACH_FRAMES - 10)
        frame_idx += 1

    for i in range(ZOOM_FRAMES):
        t = ease(i / (ZOOM_FRAMES - 1))
        render_frame(frame_idx, 0.0, axis, interpolate_view(view, close_view, t), show_interface=True)
        frame_idx += 1

    for i in range(ROTATE_FRAMES):
        render_frame(frame_idx, 0.0, axis, close_view, show_interface=True, spin_degrees=i * 2.0)
        frame_idx += 1

    build_gif()
    print(f"Wrote {GIF}")
    cmd.quit()


main()
