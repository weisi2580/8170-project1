# PyMOL overlay script (Step 7): reference vs. custom-MSA model vs. default-MSA model,
# colored by pLDDT. Requires PyMOL (open-source), and AF3 model files already downloaded.
#
# Usage (from repo root):
#   pymol -cq scripts/overlay.pml -- T1183
#
# Writes figures/<ID>/overlay.png

python
import sys, os
args = sys.argv[1:] if hasattr(sys, "argv") else []
target = args[0] if args else os.environ.get("OVERLAY_TARGET", "T1183")
python end

python
ref_map = {"T1183": ("8IFX", "B"), "T1112": ("8ORK", None), "T1122": ("8BBT", "A")}
pdb_id, chain = ref_map[target]
ref_path = f"data/RCSB_chain/{pdb_id}_{chain}.cif" if chain else f"data/RCSB/{pdb_id}.cif"
custom_glob = f"af3_custom/{target}"
default_glob = f"af3_default/{target}"

import glob
def find_model(d):
    hits = [f for f in glob.glob(f"{d}/*model*.cif") + glob.glob(f"{d}/*model*.pdb")]
    return hits[0] if hits else None

custom_model = find_model(custom_glob)
default_model = find_model(default_glob)
python end

python
cmd.load(ref_path, "reference")
if custom_model:
    cmd.load(custom_model, "custom_msa")
if default_model:
    cmd.load(default_model, "default_msa")
python end

hide everything
show cartoon
color grey70, reference
color skyblue, custom_msa
color salmon, default_msa

align custom_msa, reference
align default_msa, reference

# color predictions by pLDDT (b-factor column, AF3 convention)
spectrum b, red_yellow_green, custom_msa
spectrum b, red_yellow_green, default_msa

bg_color white
set ray_opaque_background, 0
orient reference

python
os.makedirs(f"figures/{target}", exist_ok=True)
cmd.png(f"figures/{target}/overlay.png", width=1600, height=1200, dpi=300, ray=1)
python end
