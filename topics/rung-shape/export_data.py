"""Grids for the page's interactive views (writes page/data.json)."""
import sys, os, json, math, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fractional-rungs"))
from regular_iteration import Es
R = lambda v: [None if not np.isfinite(t) else round(float(t), 5) for t in v]
def sheet_b(s, a, c=1.0):
    return Es(s, float(Es(-s, np.array([c]))[0]) - Es(-s, np.asarray(a, float)))
out = {}
S = np.round(np.arange(-3, 4.0001, 0.125), 4); A = np.round(np.linspace(0, 1, 51), 4)
out["sheet"] = {"s": S.tolist(), "a": A.tolist(), "b": [R(np.clip(sheet_b(s, A), 0, 1)) for s in S]}
af = np.linspace(0, 1, 121)
out["rungs"] = {str(s): {"a": R(af), "b": R(np.clip(sheet_b(s, af), 0, 1))} for s in [-2, -1, 0, 1, 2, 3]}
S2 = np.round(np.arange(-2, 2.0001, 0.1), 3); X = np.round(np.linspace(-0.9, 1.5, 49), 4)
fz = []
for s in S2:
    z = Es(s, X); z[(z > 5) | ~np.isfinite(z)] = np.nan; fz.append(R(z))
out["flow"] = {"s": S2.tolist(), "x": X.tolist(), "z": fz}
out["flow_rungs"] = {str(s): R(np.where(np.abs(Es(s, X)) <= 5, Es(s, X), np.nan)) for s in [-2, -1, 0, 1, 2]}
S3 = np.round(np.arange(-1, 2.0001, 0.1), 3); D = np.round(np.linspace(0, 4, 41), 3)
dz = []
for s in S3:
    c = float(Es(-s, np.array([1.0]))[0]); z = Es(s, D*c); z[(z > 16) | ~np.isfinite(z)] = np.nan; dz.append(R(z))
out["dial2"] = {"s": S3.tolist(), "d": D.tolist(), "z": dz}
res = json.load(open("results.json"))
out["area"] = res["area"]; out["uniq"] = res["uniqueness"]; out["loop"] = res["loop"]; out["growth"] = res["growth_half"]
# growth curves for the page: n -> rung s run n times on 1, s = 0, .25, .5, .75, 1
ns = np.linspace(1, 40, 79)
out["growthcurves"] = {str(s): R(np.log10(1 + Es(s, ns*float(Es(-s, np.array([1.0]))[0])))) for s in [0, 0.25, 0.5, 0.75, 1]}
out["growthcurves"]["n"] = R(ns)
json.dump(out, open("page/data.json", "w"), separators=(",", ":"))
print("bytes", os.path.getsize("page/data.json"))
# --- follow-up: the sheet through the flow's lens (u = E^-s(a), v = E^-s(b)) and the flow's "= 1" curve
out = json.load(open("page/data.json"))
S = np.array(out["sheet"]["s"]); A = np.array(out["sheet"]["a"])
lu, lv = [], []
for i, s in enumerate(S):
    b = np.array([np.nan if t is None else t for t in out["sheet"]["b"][i]])
    lu.append(R(Es(-s, A))); lv.append(R(Es(-s, b)))
out["lens"] = {"u": lu, "v": lv, "w": R(np.array([float(Es(-s, np.array([1.0]))[0]) for s in S]))}
json.dump(out, open("page/data.json", "w"), separators=(",", ":"))
print("bytes with lens", os.path.getsize("page/data.json"))
