"""DRAFT figure D4: misaligned charges. Run after checks_misalign.py."""
import os, json, math, itertools
import numpy as np, sympy as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from misalign import W_stu
from precess import u, L, E
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "..", "figures")
KNOWN, OURS, GREY = "#d9822b", "#2b6cb0", "#bbbbbb"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
D = json.load(open(os.path.join(HERE, "results_misalign.json")))["data"]
fig, axs = plt.subplots(1, 3, figsize=(20, 6.4), gridspec_kw={"width_ratios": [1.1, 1.0, 1.05]})

ax = axs[0]
sc = D["scan"]; ang = [s[0] for s in sc]; K = [s[1] for s in sc]; T = [s[2] if s[2] is not None else np.nan for s in sc]
ax.fill_between(ang, 0, [k/2*3.0 for k in K], color=GREY, alpha=0.35, label="pull K (scaled; 0 = no pull)")
ax.plot(ang, T, color=OURS, lw=2.4, label="turn per lap T")
for a0, lab, dy in [(45, "same direction\nas the black hole:\nno pull, no orbit", 1.9), (225, "exactly opposite:\nT = 1 (= n − 1)", 1.35), (135, "one charge flipped:\nT = 1", 1.6), (315, "", 0)]:
    if lab: ax.annotate(lab, xy=(a0, T[a0] if not math.isnan(T[a0]) else 0.05), xytext=(a0-40, dy), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
ax.set_xlim(0, 360); ax.set_ylim(0, 2.3); ax.set_xticks(range(0, 361, 45))
ax.set_xlabel("direction of the visitor's charge in the black hole's two-charge plane (degrees)")
ax.set_ylabel("turn per lap ÷ (π K²/L²)"); ax.legend(loc="upper right", fontsize=9.5)
ax.set_title("A. Around a two-charge black hole, every\ndirection that pulls gives the same turn, 1")

ax = axs[1]; ax.set_aspect("equal")
def orbit(ks, ws, ss, Lv, u0, phimax):
    Wx = W_stu(ks, ws, ss); F = sp.lambdify((u, E, L), Wx/L**2, "numpy"); dF = sp.lambdify((u, E, L), sp.diff(Wx/L**2, u), "numpy")
    from scipy.optimize import brentq
    E0 = brentq(lambda e: F(u0, e, Lv), 0.5, 0.99999)   # u0 is a turning point
    sol = solve_ivp(lambda p, y: [y[1], 0.5*dF(y[0], E0, Lv)], [0, phimax], [u0, 0.0], rtol=1e-11, atol=1e-13, max_step=0.01)
    return sol.t, 1/sol.y[0]
one = sp.Integer(1)
for ks, col, lab in [([one, 0, 0, 0], OURS, "around a one-charge black hole: closes"), ([one, one, one, one], KNOWN, "around a four-charge black hole: turns")]:
    ph, rr = orbit(ks, [one, 0, 0, 0], [-1, 1, 1, 1], 6.0, 1/40.0, 4*2*math.pi)
    ax.plot(rr*np.cos(ph), rr*np.sin(ph), color=col, lw=1.6 if col == KNOWN else 2.4, label=lab)
ax.plot(0, 0, "k*", ms=10); ax.legend(loc="lower center", fontsize=9.5, bbox_to_anchor=(0.5, -0.27))
ax.set_title("B. The same one-charge visitor, flipped,\nfour laps each (exact orbit equation)")

ax = axs[2]; ax.axis("off"); ax.set_xlim(0, 4.2); ax.set_ylim(-0.9, 2.6)
pts = [(0, 0), (1, 0), (0.5, 0.87), (0.5, 0.3)]
for c, n in enumerate([1, 2, 3, 4]):
    ox = c*1.05; P = [(ox+0.1+0.75*x, 1.0+0.75*y) for x, y in pts[:n]]
    for i, j in itertools.combinations(range(n), 2):
        ax.plot([P[i][0], P[j][0]], [P[i][1], P[j][1]], color=OURS, lw=2.5)
    for p in P: ax.plot(*p, "o", color=KNOWN, ms=11)
    pairs = n*(n-1)//2
    ax.text(ox+0.48, 0.6, f"{n} charge{'s' if n>1 else ''}\n{pairs} pair{'s' if pairs != 1 else ''}\nT = {n-1}", ha="center", va="top", fontsize=10.5)
ax.text(0.05, -0.35, "Each pair of the black hole's charges turns the orbit;\na single charge alone never does (exactly Kepler).\nT·K² = 2 Σ over pairs i<j of kᵢkⱼ × (how the visitor lines up with i and j)", fontsize=9.8, va="top")
ax.text(0.05, 2.45, "orange dots: the black hole's charges   blue lines: pairs", fontsize=9.5)
ax.set_title("C. The rule: turning comes from pairs of charges")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD4-turning-comes-from-pairs-of-charges.png"), dpi=130); print("ok")
