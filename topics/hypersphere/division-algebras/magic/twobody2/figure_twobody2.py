"""DRAFT figure D7 (corrected): two comparable black holes at second order, from Gupta's Lagrangian with its
static G^2 m1 m2 q1 q2 term corrected (pn2.CORRECTION). Run after checks_twobody2.py."""
import os, sys, json; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, sympy as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from orbit2pn import build, measure
OUT = os.path.join(HERE, "..", "..", "..", "figures")
KNOWN, OURS, GREY, GREEN, PURPLE = "#d9822b", "#2b6cb0", "#bbbbbb", "#2f855a", "#805ad5"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
d = json.load(open(os.path.join(HERE, "results_twobody2.json")))["data"]; emd = d["EMD_extremal_opposite"]
cols = {"1": "k", "2": GREEN, "3": OURS, "4": PURPLE}
fig, axs = plt.subplots(1, 4, figsize=(24, 6.4), gridspec_kw={"width_ratios": [1.05, 1.05, 0.85, 0.85]})
X = np.linspace(0, 1, 201); nu = X*(1 - X)
ax = axs[0]
for n in ["4", "3", "2", "1"]:
    a0, a1 = [float(sp.sympify(x)) for x in emd[n]["T2"]]
    ax.plot(X, a0 + a1*nu, color=cols[n], lw=2.6, label=f"{n} charge{'s' if n != '1' else ''} each: {emd[n]['T2'][0]} at the ends, {sp.sympify(emd[n]['T2'][0]) + sp.sympify(emd[n]['T2'][1])/4} at equal masses")
ax.axvline(0.5, color=GREY, ls="--"); ax.set_xlim(0, 1); ax.set_ylim(-0.6, 13.5)
ax.set_xlabel("share of the total mass carried by the first body"); ax.set_ylabel("second-order turn coefficient")
ax.legend(loc="center left", fontsize=9.5, bbox_to_anchor=(0.0, 0.62)); ax.annotate("one charge each: zero only\nfor a small visitor (the ends)", xy=(0.3, 0.95), xytext=(0.04, 2.6), fontsize=10, arrowprops=dict(arrowstyle="->"))
ax.set_title("A. Two black holes with opposite charges, at second order\n(corrected equations): one-charge pairs turn too")
ax = axs[1]
for n in ["1", "2", "3", "4"]:
    a1 = float(sp.sympify(emd[n]["T2"][1]))
    ax.plot(X, a1*nu, color=cols[n], lw=2.6, label=f"{n} charge{'s' if n != '1' else ''} each: {emd[n]['T2'][1]} ν")
gr = d["GR"]; ax.plot(X, float(sp.sympify(gr["T2"][1]))*nu/6, color=KNOWN, lw=2.4, ls="--", label="Einstein's theory ÷ 6 (Damour–Schäfer)")
ax.axhline(0, color=GREY, lw=0.8); ax.set_xlim(0, 1); ax.set_ylim(-1.4, 1.4)
ax.set_xlabel("share of the total mass carried by the first body")
ax.set_ylabel("extra second-order turn from the mass ratio")
ax.annotate("(9/2 − 3T) ν, with T = n − 1\nthe first-order turn from pairs", xy=(0.5, 1.125), xytext=(0.52, 0.55), fontsize=10.5, arrowprops=dict(arrowstyle="->"))
ax.legend(loc="lower left", fontsize=9.5); ax.set_title("B. The mass-ratio part: 9/2 ν from the corrected\nstatic term, less 3T ν from pairs")
def paths(av, Jk=4.0, laps=10):  # Jk = J/k
    mdl = build('EMD', a=av); _, J, E, sol = measure(mdl, Jk, 0.3, laps)
    t = np.linspace(0, sol.t[-1], 6000); y = sol.sol(t); rr, ph = y[0], y[2]
    x, yy = rr*np.cos(ph), rr*np.sin(ph)
    i = [j for j in range(1, len(rr) - 1) if rr[j] < rr[j - 1] and rr[j] <= rr[j + 1]]
    global PERI; PERI = [(x[j]/2, yy[j]/2) for j in i]
    return x/2, yy/2, -x/2, -yy/2
for ax, av, ttl, kw in [(axs[2], sp.sqrt(3), "C. Equal masses, one charge each: the orbit turns at\nsecond order (J/k = 2.5, 6 laps, exaggerated)", dict(Jk=2.5, laps=6)), (axs[3], 0, "D. Equal masses, four charges each:\nthe orbits turn (T = 3 at first order)", {})]:
    xa, ya, xb, yb = paths(av, **kw); sc = np.max(np.hypot(xa, ya))
    ax.plot(xa/sc, ya/sc, color=OURS, lw=1.2, label="first body"); ax.plot(xb/sc, yb/sc, color=KNOWN, lw=1.2, label="second body")
    if av != 0:
        for j, (px, py) in enumerate(PERI):
            ax.plot([0, px/sc*1.1], [0, py/sc*1.1], color="k", lw=1.0, alpha=0.35 + 0.65*j/max(1, len(PERI) - 1))
        ax.annotate("closest approach,\nlap by lap", xy=(PERI[-1][0]/sc*1.1, PERI[-1][1]/sc*1.1), xytext=(-0.95, 0.85), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
    ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(loc="lower center", fontsize=9.5, bbox_to_anchor=(0.5, -0.16), ncol=2); ax.set_title(ttl)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD7-two-black-holes-at-second-order.png"), dpi=110, bbox_inches="tight"); print("saved")
