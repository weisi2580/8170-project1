#!/usr/bin/env python3
"""Render 3D overlays of reference vs. custom-MSA model vs. default-MSA model
with PyMOL (pymol-open-source, importable from the venv).

For each target the top-ranked model of each condition is superposed onto
the experimental reference by C-alpha pairs within the CASP evaluation unit,
and rendered twice:

  figures/<ID>/overlay.png         reference grey, custom blue, default orange
  figures/<ID>/overlay_plddt.png   the two models side by side, coloured by
                                   pLDDT (AF3 bands: >90, 70-90, 50-70, <50)

Usage:
    python scripts/render_overlays.py [--target T1183]
"""
import argparse
import csv
from pathlib import Path

from PIL import Image, ImageChops
from pymol import cmd

ROOT = Path(__file__).resolve().parents[1]
REFS = {
    "T1183": "data/RCSB_chain/8IFX_B.cif",
    "T1112": "data/RCSB/8ORK.cif",
    "T1122": "data/RCSB_chain/8BBT_A.cif",
}
# PyMOL object names ("custom" is a reserved selection keyword in PyMOL)
OBJ = {"custom": "m_custom", "default": "m_default"}
COLORS = {"reference": [0.62, 0.62, 0.60], "m_custom": [0.165, 0.471, 0.839], "m_default": [0.922, 0.408, 0.204]}
# AlphaFold's pLDDT band colours
PLDDT_BANDS = [(90, 101, [0.0, 0.325, 0.839]), (70, 90, [0.396, 0.796, 0.953]),
               (50, 70, [1.0, 0.859, 0.075]), (0, 50, [1.0, 0.490, 0.271])]


def png(path: Path, width, height, margin=40):
    """Ray-trace to PNG, then crop PyMOL's bounding-sphere whitespace."""
    cmd.png(str(path), width=width, height=height, dpi=200, ray=1)
    img = Image.open(path).convert("RGB")
    bbox = ImageChops.difference(img, Image.new("RGB", img.size, (255, 255, 255))).getbbox()
    if bbox:
        l, t, r, b = bbox
        img.crop((max(l - margin, 0), max(t - margin, 0),
                  min(r + margin, img.width), min(b + margin, img.height))).save(path)


def eval_unit(casp_id):
    with open(ROOT / "data" / "metadata.csv") as fh:
        row = next(r for r in csv.DictReader(fh) if r["casp_id"] == casp_id)
    return row["eval_unit"]


def setup_scene():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0)
    cmd.set("ray_trace_mode", 0)
    cmd.set("cartoon_loop_radius", 0.35)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("specular", 0.2)
    for name, rgb in {**COLORS, **{f"plddt{i}": b[2] for i, b in enumerate(PLDDT_BANDS)}}.items():
        cmd.set_color(f"c_{name}", rgb)


def render(casp_id):
    unit = eval_unit(casp_id)
    cmd.reinitialize()
    setup_scene()
    cmd.load(str(ROOT / REFS[casp_id]), "reference")
    cmd.remove("reference and not polymer.protein")
    cmd.remove("not alt ''+A")
    cmd.alter("all", "alt=''")
    ref_resi = {a.resi for a in cmd.get_model(f"reference and name CA and resi {unit}").atom}
    for cond in ("custom", "default"):
        obj = OBJ[cond]
        cmd.load(str(ROOT / f"af3_{cond}" / casp_id / f"{casp_id}_model.cif"), obj)
        # residue numbers match (CASP FASTA numbering), so fit CA pairs directly
        common = sorted(ref_resi & {a.resi for a in cmd.get_model(f"{obj} and name CA").atom}, key=int)
        sel = "resi " + "+".join(common)
        rms = cmd.pair_fit(f"{obj} and name CA and {sel}", f"reference and name CA and {sel}")
        print(f"[{casp_id}/{cond}] CA RMSD after least-squares fit on {len(common)} residues: {rms:.2f} A")
    cmd.hide("everything")
    cmd.show("cartoon")
    for obj in COLORS:
        cmd.color(f"c_{obj}", obj)
    cmd.set("cartoon_transparency", 0.35, "reference")
    cmd.orient("reference")
    cmd.zoom("all", 0, complete=1)
    out_dir = ROOT / "figures" / casp_id
    out_dir.mkdir(parents=True, exist_ok=True)
    png(out_dir / "overlay.png", 1800, 1350)
    print(f"wrote figures/{casp_id}/overlay.png")

    # pLDDT view: models side by side (default translated along screen x)
    view = cmd.get_view()
    cmd.disable("reference")
    for i, (lo, hi, _) in enumerate(PLDDT_BANDS):
        cmd.color(f"c_plddt{i}", f"(m_custom or m_default) and b > {lo - 0.001} and b < {hi}")
    cmd.set("cartoon_transparency", 0)
    ext = cmd.get_extent("m_custom")
    width = max(hi - lo for lo, hi in zip(*ext)) * 1.1
    cmd.translate([width, 0, 0], "m_default", camera=1)
    cmd.orient("m_custom or m_default")
    cmd.zoom("m_custom or m_default", 0, complete=1)
    png(out_dir / "overlay_plddt.png", 2400, 1200)
    cmd.set_view(view)
    print(f"wrote figures/{casp_id}/overlay_plddt.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=list(REFS))
    args = parser.parse_args()
    for casp_id in [args.target] if args.target else list(REFS):
        render(casp_id)


if __name__ == "__main__":
    main()
