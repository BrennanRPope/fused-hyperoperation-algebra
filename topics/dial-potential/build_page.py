"""Build index.html (the write-up page) from results.json and running.npz.  Run after search.py."""
import json, math, os, html
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
run = np.load(os.path.join(HERE, "running.npz"))

GROUP = {"flat": "higgs", "cw": "higgs", "veltman": "higgs", "unif": "higgs", "angle": "formula",
         "tree_lam": "formula", "sumrule": "formula", "yt": "formula", "alpha": "formula",
         "gravity": "gravity", "bd": "gravity", "lightdial": "gravity", "hier": "gravity",
         "quint": "dark", "quint_shape": "dark"}

pts = []
for i, rnd in enumerate(R["rounds"]):
    for g in rnd:
        nb = g["net_bits"]
        if nb is None or (isinstance(nb, float) and math.isnan(nb)):
            continue
        tested = bool(g["rows"]) or "cosmo_bits" in g["extra"]
        if not tested:
            continue
        nbv = float(nb) if isinstance(nb, (int, float)) and math.isfinite(nb) else -1e9
        pts.append({"r": i + 1, "b": round(max(nbv, -30), 2), "out": bool(g["ruled_out"]),
                    "g": GROUP.get(g["family"], "formula"), "n": g["name"], "id": g["id"]})

mu = run["mu"]; lam = run["y"][4]
sel = np.unique(np.linspace(0, len(mu) - 1, 160).astype(int))
lam_curve = [[round(math.log10(mu[j]), 3), round(float(lam[j]), 5)] for j in sel]
prof = R["slope_profile"]
cmin = min(prof["chi2"])
slope_curve = [[round(l, 3), round(c - cmin, 3)] for l, c in zip(prof["lam"], prof["chi2"])]

r1 = R["rounds"][0]


def esc(s):
    return html.escape(str(s))


def verdict(g):
    if g["ruled_out"]:
        return "ruled out"
    if not g["rows"] and "cosmo_bits" not in g["extra"]:
        return "nothing to test"
    return "alive"


def fmt_bits(b):
    if b is None:
        return ""
    b = float(b)
    if not math.isfinite(b) or b < -1e4:
        return "far below −1000"
    return f"{b:+.1f}"


rows1 = []
for g in sorted(r1, key=lambda g: g["id"]):
    rows1.append(f"<tr><td>{esc(g['id'])}</td><td>{esc(g['name'])}<div class=sub>{esc(g['shape'])}; {esc(g['coupling'])}</div></td>"
                 f"<td>{esc(g['grade'])}</td><td class=num>{fmt_bits(g['net_bits'])}</td><td>{verdict(g)}</td></tr>")

data = {"pts": pts, "lam": lam_curve, "slope": slope_curve, "best": prof["best"], "band": prof["one_sigma"],
        "geo": math.sqrt(2 / 3)}

tpl = open(os.path.join(HERE, "page_template.html")).read()
out = tpl.replace("/*DATA*/", "const D=" + json.dumps(data) + ";").replace("<!--ROUND1-->", "\n".join(rows1))
for k, v in R.get("page_numbers", {}).items():
    out = out.replace("{{" + k + "}}", str(v))
open(os.path.join(HERE, "index.html"), "w").write(out)
print("wrote index.html", len(out), "bytes;", len(pts), "points")
