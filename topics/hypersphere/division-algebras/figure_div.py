"""DRAFT figures D1, D2. Run after checks_div.py and checks_shape.py."""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "figures")
KNOWN, OURS, BAD, OK = "#d9822b", "#2b6cb0", "#c0392b", "#3a9a5b"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
R = json.load(open(os.path.join(HERE, "results_div.json")))
S = json.load(open(os.path.join(HERE, "results_shape.json")))

# ---------- D1 ----------
fig, axs = plt.subplots(1, 3, figsize=(19, 6.2), gridspec_kw={"width_ratios": [1.45, 1, 0.9]})
ax = axs[0]
rows = [("1", "R  real"), ("2", "C  complex"), ("4", "H  quaternions"), ("8", "O  octonions"), ("16", "S  sedenions")]
cols = [("norm multiplies |xy|=|x||y|", "size-keeping\nproduct\n= polar split\nof ∘ₖ"),
        ("commutative", "order of\nfactors\ndoesn't matter"),
        ("associative", "grouping\ndoesn't matter\n(∘ₖ too)"),
        ("alternative", "alternative\n(weak grouping)\n= superstring\n3-ψ rule"),
        ("zd", "no zero\ndivisors")]
for i, (n, name) in enumerate(rows):
    for j, (key, lab) in enumerate(cols):
        v = (n != "16") if key == "zd" else R["props"][n][key]
        ax.add_patch(plt.Rectangle((j, -i), 0.94, 0.9, color=OK if v else BAD, alpha=0.85))
        ax.text(j + 0.47, -i + 0.45, "yes" if v else "no", ha="center", va="center", color="white", fontweight="bold")
    ax.text(-0.1, -i + 0.45, f"dim {n}:  {name}", ha="right", va="center", fontsize=11)
for j, (key, lab) in enumerate(cols):
    ax.text(j + 0.47, 1.05, lab, ha="center", va="bottom", fontsize=9.5)
ax.set_xlim(-2.6, 5.1); ax.set_ylim(-5.0, 2.6); ax.axis("off")
ax.text(-2.6, -4.85, "From C on, each doubling loses at least one property. The circle operation inherits each one exactly,\n"
        "through its lift 1 + k(a∘b) = (1+ka)(1+kb), which holds in all five.", fontsize=9.5)
ax.set_title("A. What each doubling keeps (checked exactly)")

ax = axs[1]
ns = [1, 2, 4, 8]
for i, n in enumerate(ns):
    ax.bar(i, 2, color=KNOWN); ax.bar(i, n, bottom=2, color=OURS)
    ax.text(i, 1, "time +\n1 space", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
    ax.text(i, 2 + n/2, f"{n}", ha="center", va="center", color="white", fontsize=11, fontweight="bold")
    ax.text(i, n + 2.2, f"{n+2}", ha="center", fontsize=14, fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels(["over R", "over C", "over H", "over O"])
ax.set_ylabel("dimensions of spacetime"); ax.set_ylim(0, 12)
ax.text(-0.45, 11.3, "orange: relativity's plane, size t² − z²\nblue: the algebra K, size |x|²", fontsize=9.5, va="top")
ax.set_title("B. Superstring spacetimes = 2 + (1, 2, 4, 8)\nsize = (t+z)(t−z) − |x|²  (det of a 2×2 matrix over K)")

ax = axs[2]
res = [R["residual"][str(n)][0] for n in [1, 2, 4, 8, 16]]
ax.bar(range(5), res, color=[OK]*4 + [BAD]); ax.set_yscale("log"); ax.set_ylim(1e-18, 3)
ax.set_xticks(range(5)); ax.set_xticklabels(["R", "C", "H", "O", "S"])
for i, v in enumerate(res):
    ax.text(i, v*2, "0 (round-off)" if v < 1e-10 else f"{v:.2f}", ha="center", fontsize=9)
ax.set_ylabel("leftover of the 3-ψ rule (mean, unit spinors)")
ax.set_title("C. The superstring's spinor rule\nholds through octonions, fails at 16")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD1-one-root-1-2-4-8-and-3-4-6-10.png"), dpi=130); plt.close(fig)

# ---------- D2 ----------
def norm(v, p):
    v = np.abs(np.asarray(v, float))
    return v.max() if p == math.inf else (v**p).sum()**(1/p)
fig, axs = plt.subplots(1, 3, figsize=(18, 5.8), gridspec_kw={"width_ratios": [1, 1.1, 1.25]})
ax = axs[0]; ax.set_aspect("equal")
th = np.linspace(0, 2*np.pi, 800)
for p, c, lab in [(1, BAD, "p = 1 diamond (cross-polytope)"), (2, OURS, "p = 2 sphere"), (4, "#999999", "p = 4"), (math.inf, KNOWN, "p = ∞ cube")]:
    pts = np.array([[math.cos(t), math.sin(t)] for t in th]); pts = np.array([q/norm(q, p) for q in pts])
    ax.plot(pts[:, 0], pts[:, 1], color=c, lw=2.4 if p == 2 else 1.6, label=lab)
ax.legend(loc="lower center", fontsize=8.5, bbox_to_anchor=(0.5, -0.32)); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15)
ax.set_title("A. The ℓp shapes: points of size 1")
ax = axs[1]
keys = ["1", "1.5", "2", "3", "4", "inf"]; labs = ["1", "1.5", "2", "3", "4", "∞"]
rot = [S[k]["S2 size change under a full rotation"] for k in keys]
drift = [S[k]["S4 angular momentum drift in a Kepler run"] for k in keys]
prod = [S[k]["S1 product keeps size, plane (max err)"] for k in keys]
x = np.arange(len(keys)); w = 0.27
ax.bar(x - w, rot, w, color="#7aa6d6", label="size change when you turn")
ax.bar(x, drift, w, color=KNOWN, label="angular momentum lost")
ax.bar(x + w, [min(v, 1.0) for v in prod], w, color="#888888", label="product size error (capped at 1)")
ax.set_xticks(x); ax.set_xticklabels([f"p = {l}" for l in labs]); ax.set_ylim(0, 1.45)
ax.text(2, 0.06, "all 0", ha="center", color=OURS, fontweight="bold")
ax.legend(fontsize=9, loc="upper center", ncol=2); ax.set_ylabel("failure size (0 = keeps it)")
ax.set_title("B. Only p = 2 keeps rotations, angular\nmomentum and a size-keeping product")
ax = axs[2]; ax.set_aspect("equal")
for p, c, lab in [(1, BAD, "pull from the p = 1 shape"), (2, OURS, "pull from the sphere (Newton)")]:
    def rhs(t, y, p=p):
        X, Y, vx, vy = y; h = 1e-6; V = lambda a, b: -1/norm([a, b], p)
        return [vx, vy, -(V(X+h, Y)-V(X-h, Y))/(2*h), -(V(X, Y+h)-V(X, Y-h))/(2*h)]
    sol = solve_ivp(rhs, [0, 40], [1, 0.2, 0, 1.0], rtol=1e-10, atol=1e-12, max_step=0.01)
    ax.plot(sol.y[0], sol.y[1], color=c, lw=1.2 if p == 1 else 2.2, label=lab)
ax.plot(0, 0, "k*", ms=10); ax.legend(fontsize=9, loc="lower center", bbox_to_anchor=(0.5, -0.3))
ax.set_title("C. Same start, same time: the sphere's orbit closes;\nthe diamond's wanders and its angular momentum drifts")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD2-gravity-needs-the-sphere.png"), dpi=130); plt.close(fig)
print("ok")
