"""shapes_layer.py: put the proof gallery's shape families on the laws landscape.

Each settlement (numbered dot) on the landscape is linked to the library theorems that
state its law, by name and statement (the LINKS table below: that link is the only
choice made by hand). Each theorem's shape family comes from the gallery's rules
(../proof-gallery/gallery.json), never from this file. The region comes from the dot's
position: below the parabola d = k²/4 two real walls, on it a double wall, above it none.
Output: shapes.json (used by index.html).
"""
import json, pathlib
from collections import Counter

HERE = pathlib.Path(__file__).resolve().parent
G = json.load(open(HERE.parent / "proof-gallery" / "gallery.json"))
by_short = {}
for p in G["proofs"]:
    by_short.setdefault(p["name"].split(".")[-1], []).append(p)
fam = {g["id"]: g for g in G["groups"]}

# dot number, (x, y) as on the map, title, [(module suffix or "", theorem short name)]
LINKS = [
 (1, 1, 0, "The circle operation", [("Base", "multAdd_lift"), ("Base", "multAdd_comm"), ("Base", "multAdd_assoc"),
    ("DarkSector", "redshift_compose"), ("DarkSector", "redshift_power"), ("DataCenter", "pue_corner"),
    ("DataCenter", "additive_books_corner"), ("DataCenter", "foldMA_split"), ("Finance", "fisher"),
    ("InformationTheory", "uniform_tsallis_is_fused_power"), ("InformationTheory", "campbell"),
    ("Landscape", "multadd_point")]),
 (2, 0, 0, "Addition", [("Base", "wAdd_weight_zero"), ("InformationTheory", "hazards_add"),
    ("InformationTheory", "shannon_additive"), ("InformationTheory", "shannon_weight"),
    ("InformationTheory", "hartley"), ("InformationTheory", "lengths_add"),
    ("NavierStokes", "burgers_corner_defect"), ("NavierStokes", "burgers_corner_ok_iff"),
    ("Thermodynamics", "renyi_adds"), ("DialSweep", "shear_is_add")]),
 (3, -1, 0, "Union of failures", [("Base", "wAdd_either"), ("InformationTheory", "union_rule"),
    ("InformationTheory", "survival"), ("InformationTheory", "union_is_flipped_madd"),
    ("InformationTheory", "union_grows"), ("InformationTheory", "series_system"),
    ("InformationTheory", "gini_simpson"), ("InformationTheory", "kraft_slack"),
    ("QuantumNoise", "amplitude_damping"), ("Consensus", "either"), ("Consensus", "caught")]),
 (4, -2, 0, "Parity", [("InformationTheory", "bsc_cascade"), ("InformationTheory", "bsc_cascade_diag"),
    ("InformationTheory", "piling_up"), ("InformationTheory", "n_flips"), ("InformationTheory", "parity3"),
    ("QuantumNoise", "bit_flip"), ("QuantumNoise", "stacked_bit_flips"),
    ("Landscape", "parity_reaches_floor"), ("Thermodynamics", "ising_series")]),
 (5, -1.5, 0, "3-letter channel", [("InformationTheory", "qary_cascade"),
    ("InformationTheory", "ladder_never_cancel"), ("InformationTheory", "ladder_binary"),
    ("Landscape", "channel_floor")]),
 (6, -4/3, 0, "4-letter channel", [("InformationTheory", "ladder_four"), ("QuantumNoise", "depolarizing")]),
 (7, -9/8, 0, "Depolarizing qutrit", []),
 (8, 2, 0, "2× leverage", []),
 (9, -0.5, 0, "Partial damage", []),
 (10, 0, -1, "Velocity addition", [("Landscape", "velocity_point"), ("InformationTheory", "velocity_corners"),
    ("InformationTheory", "odds_to_velocity"), ("InformationTheory", "vadd_assoc"),
    ("NumberPlane", "split_point"), ("Statistics", "amh_velocity"), ("Statistics", "pool_velocity"), ("DialSweep", "boost_walls")]),
 (11, 0, 1, "Slope addition", [("", "tangent_addition"), ("", "complex_point"), ("DialSweep", "rotation_no_wall")]),
 (12, 0, -0.5, "Sine addition", [("", "angle_addition")]),
 (13, 2, 1, "Euler characteristic", [("Landscape", "chi_y_line"), ("Landscape", "chi_y_touch"),
    ("Elliptic", "ochanine_nodal_num")]),
 (14, 2, 0.5, "Rung 2", [("Rungs", "fexp_lift"), ("Rungs", "fexp_iter"), ("Rungs", "fexp_tower"),
    ("Rungs", "fexp_multAdd_right"), ("Rungs", "fexp_sadd"), ("Rungs", "fexp_carries_multAdd")]),
 (15, -1, 0.5, "Rung −1", [("FractionalRungs", "fuse_max"), ("FractionalRungs", "fuse_add")]),
 (16, 0.5, -1/16, "Half rung", [("Landscape", "lens_lift"), ("Landscape", "lens_factor_u"),
    ("Landscape", "lens_factor_v")]),
]
# laws that live on a stretch of the map rather than at one dot
STRETCHES = [
 ("KdV corner, dimension dial", "the straight line between k = 0 and k = 1", "line",
  [("Landscape", "dial_is_multadd"), ("Landscape", "kdv_phase_shift"), ("FractionalRungs", "kdv_weight"),
   ("FractionalRungs", "dial_or"), ("FractionalRungs", "dial_three"), ("Base", "wAdd_scale")]),
 ("The border (double wall)", "the parabola d = k²/4", "border",
  [("Landscape", "parabola_additive"), ("Landscape", "parabola_den"), ("Landscape", "discriminant"), ("DialSweep", "dial_discriminant")]),
 ("Every law at once", "the whole map: F(a, b) = N/D", "all",
  [("Landscape", "num_comm"), ("Landscape", "den_comm"), ("Landscape", "assoc_num"), ("Landscape", "assoc_den"),
   ("Landscape", "inverse_num"), ("Landscape", "units_num"), ("Landscape", "units_den"),
   ("DialSweep", "dial_product"), ("DialSweep", "dial_lens"), ("DialSweep", "dial_det")]),
]

def find(mod, short):
    c = by_short.get(short, [])
    c = [p for p in c if (not mod) or p["module"].endswith(mod)]
    return c[0] if c else None

def region(x, y):
    par = x * x / 4
    if abs(y - par) < 1e-9: return "border"
    return "two" if y < par else "none"

def pack(p):
    return dict(n=p["name"], s=p["statement"], c=p["cluster"], b=p["body"], lead=p["lead"],
                mv=p["moves"], km=p["kmoves"], fr=p["frame"], src=p["lead_src"])

out_dots, missing = [], []
for n, x, y, t, links in LINKS:
    ps = []
    for mod, short in links:
        p = find(mod, short)
        if p: ps.append(pack(p))
        else: missing.append((n, mod, short))
    cnt = Counter(p["c"] for p in ps)
    out_dots.append(dict(n=n, x=x, y=y, t=t, region=region(x, y), proofs=ps,
                         families=cnt.most_common(), main=(cnt.most_common(1)[0][0] if cnt else None)))
out_str = []
for t, where, reg, links in STRETCHES:
    ps = [pack(p) for p in (find(m, s) for m, s in links) if p]
    cnt = Counter(p["c"] for p in ps)
    out_str.append(dict(t=t, where=where, region=reg, proofs=ps, families=cnt.most_common(),
                        main=(cnt.most_common(1)[0][0] if cnt else None)))
used = sorted({f for d in out_dots + out_str for f, _ in d["families"]})
fams = {f: dict(id=f, name=fam[f]["name"], body=fam[f]["body"], action=fam[f]["action"],
               sentence=fam[f]["sentence"], count=fam[f]["count"], height=fam[f]["height"]) for f in used}
# laws grouped by shape: family -> settlements and stretches whose proofs fall in it
grouped = {}
for d in out_dots:
    for f, k in d["families"]:
        grouped.setdefault(f, []).append(dict(kind="dot", n=d["n"], t=d["t"], k=k, region=d["region"]))
for i, d in enumerate(out_str):
    for f, k in d["families"]:
        grouped.setdefault(f, []).append(dict(kind="stretch", n=i, t=d["t"], k=k, region=d["region"]))
json.dump(dict(dots=out_dots, stretches=out_str, families=fams, grouped=grouped,
               stats=dict(linked=sum(len(d["proofs"]) for d in out_dots + out_str),
                          dots_with_proofs=sum(1 for d in out_dots if d["proofs"]),
                          families=len(used))),
          open(HERE / "shapes.json", "w"), ensure_ascii=False, indent=1)
print("missing:", missing)
for d in out_dots + out_str:
    print(d.get("n", "-"), d["t"], d["region"], [(fams[f]["name"], k) for f, k in d["families"]])
