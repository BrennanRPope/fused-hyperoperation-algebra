"""DRAFT figure D7: two comparable black holes at second order. Run after checks_twobody2.py."""
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
ax.legend(loc="center left", fontsize=9.5, bbox_to_anchor=(0.0, 0.62)); ax.annotate("one charge each:\nstill exactly zero", xy=(0.3, 0), xytext=(0.06, 1.6), fontsize=10, arrowprops=dict(arrowstyle="->"))
ax.set_title("A. Two black holes with opposite charges, at second order:\nequal masses turn a little more than a small visitor")
ax = axs[1]
for n in ["2", "3", "4"]:
    T = float(sp.sympify(emd[n]["T"][0])); a1 = float(sp.sympify(emd[n]["T2"][1]))
    ax.plot(X, a1*nu/T, color=cols[n], lw=6 - int(n), alpha=0.85, label=f"{n} charges each")
gr = d["GR"]; ax.plot(X, float(sp.sympify(gr["T2"][1]))*nu/6, color=KNOWN, lw=2.4, ls="--", label="Einstein's theory, uncharged (Damour–Schäfer)")
ax.axhline(0, color=GREY, lw=0.8); ax.set_xlim(0, 1); ax.set_ylim(-0.75, 0.85)
ax.set_xlabel("share of the total mass carried by the first body")
ax.set_ylabel("extra second-order turn from the mass ratio ÷ first-order turn")
ax.annotate("all three lie on one curve:\n(3/2) × ν × (first-order turn)", xy=(0.5, 0.375), xytext=(0.2, 0.15), fontsize=10.5, arrowprops=dict(arrowstyle="->"))
ax.annotate("Einstein: −15 ν ÷ 6, a different rule", xy=(0.62, -0.58), xytext=(0.66, -0.36), fontsize=10, color=KNOWN, arrowprops=dict(arrowstyle="->", color=KNOWN))
ax.legend(loc="upper right", fontsize=9.5); ax.set_title("B. The mass ratio only enters multiplied by the\nturn that pairs of charges already give")
def paths(av, Jk=4.0, laps=10):
    mdl = build('EMD', a=av); _, J, E, sol = measure(mdl, Jk, 0.3, laps)
    t = np.linspace(0, sol.t[-1], 6000); y = sol.sol(t); rr, ph = y[0], y[2]
    x, yy = rr*np.cos(ph), rr*np.sin(ph); return x/2, yy/2, -x/2, -yy/2
for ax, av, ttl in [(axs[2], sp.sqrt(3), "C. Equal masses, one charge each:\nthe orbits still close at second order"), (axs[3], 0, "D. Equal masses, four charges each:\nthe orbits turn (T = 3 at first order)")]:
    xa, ya, xb, yb = paths(av); sc = np.max(np.hypot(xa, ya))
    ax.plot(xa/sc, ya/sc, color=OURS, lw=1.2, label="first body"); ax.plot(xb/sc, yb/sc, color=KNOWN, lw=1.2, label="second body")
    ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(loc="lower center", fontsize=9.5, bbox_to_anchor=(0.5, -0.16), ncol=2); ax.set_title(ttl)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD7-two-black-holes-at-second-order.png"), dpi=110, bbox_inches="tight"); print("saved")
