"""DRAFT figure D8: two comparable black holes whose charges point different ways, at second order, and the
corrected static term. Run after checks_mixed2.py and checks_static.py."""
import os, sys, json; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = os.path.join(HERE, "..", "..", "..", "figures")
KNOWN, OURS, GREY, GREEN, PURPLE, RED = "#d9822b", "#2b6cb0", "#bbbbbb", "#2f855a", "#805ad5", "#c53030"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
d = json.load(open(os.path.join(HERE, "results_mixed2.json")))["data"]
pairs = [k for k in d if k != "simulated"]
Xs = sp.Symbol('X'); X = np.linspace(0, 1, 201)
def lab(k):
    A, B = k.split("|"); f = lambda s: "(" + ",".join(x.strip() for x in s.strip("() ").split(",")) + ")"
    return f"{f(A)} with {f(B)}"
fig, axs = plt.subplots(1, 4, figsize=(25, 6.6), gridspec_kw={"width_ratios": [1.1, 0.95, 1.0, 0.95]})
cl = [OURS, GREEN, PURPLE, KNOWN, RED, "k", "#319795"]
# A: T2 across the mass ratio for each mixed pair, with both small-visitor ends
ax = axs[0]
for i, k in enumerate(pairs):
    f = sp.lambdify(Xs, sp.sympify(d[k]["T2"]), "numpy"); y = f(X) * np.ones_like(X)
    ax.plot(X, y, color=cl[i], lw=2.2, label=lab(k))
    ax.plot([0, 1], [float(sp.sympify(d[k]["probe_A_round_B"][1])), float(sp.sympify(d[k]["probe_B_round_A"][1]))], "o", color=cl[i], ms=7, mfc="white", mew=2)
ax.set_xlim(0, 1); ax.set_xlabel("share of the total mass carried by the first body, X")
ax.set_ylabel("second-order turn coefficient T₂"); ax.legend(fontsize=8.6, loc="center left", bbox_to_anchor=(0.0, 0.55))
ax.set_title("A. Seven pairs with charges pointing different ways:\ncircles = exact small-visitor values (section 11)")
# B: energy part of second order follows mass shares + nu T on every pair
ax = axs[1]
for i, k in enumerate(pairs):
    T2e = sp.sympify(d[k]["T2e"]); T = sp.sympify(d[k]["T"])
    sh = Xs*sp.sympify(d[k]["probe_B_round_A"][2]) + (1 - Xs)*sp.sympify(d[k]["probe_A_round_B"][2])
    xs = np.linspace(0.05, 0.95, 7)
    lhs = [float((T2e - sh).subs(Xs, x)) for x in xs]; rhs = [float((Xs*(1 - Xs)*T).subs(Xs, x)) for x in xs]
    ax.plot(rhs, lhs, "o", color=cl[i], ms=6)
mx = 1.2; ax.plot([0, mx], [0, mx], color=GREY, lw=1.5, zorder=0); ax.set_xlim(0, mx); ax.set_ylim(0, mx)
ax.set_xlabel("ν × first-order turn T   (ν = X(1 − X))"); ax.set_ylabel("T₂ₑ minus the mass-share average of the two visitors")
ax.set_title("B. The energy part follows one rule on every pair:\nT₂ₑ = shares + νT (what the mass-energy relation requires)")
# C: kappa, the pair constant in the energy-free part
ax = axs[2]
names = [lab(k) for k in pairs] + [f"{n} equal charge{'s' if n > 1 else ''} each, opposite" for n in (1, 2, 3, 4)]
vals = [float(sp.sympify(d[k]["kappa"])) for k in pairs] + [1.5*(4 - n) for n in (1, 2, 3, 4)]
txt = [d[k]["kappa"] for k in pairs] + [str(sp.Rational(3, 2)*(4 - n)) for n in (1, 2, 3, 4)]
cols = cl + [GREY]*4; yy = np.arange(len(vals))[::-1]
ax.barh(yy, vals, color=cols, height=0.65)
for y, v, t in zip(yy, vals, txt): ax.text(v + 0.3, y, t, va="center", fontsize=9.5)
ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=8.6); ax.set_xlim(0, 30)
ax.set_xlabel("κ in T₂ = shares + ν(κ − 3T/2)")
ax.set_title("C. The energy-free part needs one new number\nper pair, κ; no formula found for it yet")
# D: the corrected static term
ax = axs[3]; a2 = np.linspace(0, 3, 301)
ax.plot(a2, 12*(a2 - 1)*(1 + a2)**2, color=KNOWN, lw=2.6, label="as printed (Gupta 2026, Eq. 7.1)")
ax.plot(a2, 0*a2, color=OURS, lw=2.6, label="corrected (our static calculation)")
for n, aa in [(4, 0), (3, 1/3), (2, 1), (1, 3)]:
    ax.plot(aa, 12*(aa - 1)*(1 + aa)**2, "o", color=KNOWN, ms=7); ax.annotate(f"{n} charge{'s' if n > 1 else ''}", (aa, 12*(aa - 1)*(1 + aa)**2), textcoords="offset points", xytext={4: (6, -16), 3: (8, 8), 2: (6, -16), 1: (-60, -4)}[n], fontsize=9.5)
ax.axhline(0, color=GREY, lw=0.8); ax.set_xlim(0, 3.1)
ax.set_xlabel("dilaton coupling a²"); ax.set_ylabel("static second-order energy, units G³m₁²m₂²/r³")
ax.legend(loc="upper left", fontsize=9.5)
ax.set_title("D. Two extremal black holes with like charges\nfeel no force: the printed equations break this")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD8-misaligned-charges-at-second-order.png"), dpi=110, bbox_inches="tight"); print("saved")
