# DRAFT. Figures for the bridges study. Colour key on every picture:
# green = proved or derived from stated assumptions, blue = observed in a lab, orange = conjectured,
# grey = contested, magenta = ours (Affine Composition Theory), dashed line = a link that is our reading.
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
C = dict(proved="#2e8b57", observed="#2b6cb0", conj="#dd8a1e", contest="#8a8a8a", ours="#c2187a")
def box(ax, x, y, text, kind, w=2.25, h=0.78, fs=8.2):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.05", fc=C[kind], ec="black", lw=0.8, alpha=0.97, zorder=2))
    ax.text(x, y, text, ha="center", va="center", color="white", fontsize=fs, weight="bold", zorder=3)
def line(ax, p, q, dashed=False, col="black"):
    ax.plot([p[0], q[0]], [p[1], q[1]], ls="--" if dashed else "-", color=col, lw=1.4, zorder=1.5)

# ---- bq1: the bridge map ----
fig, ax = plt.subplots(figsize=(13, 8.2)); ax.set_xlim(0, 13); ax.set_ylim(0, 9.1); ax.axis("off")
ax.text(6.5, 8.85, "bq1. Every quantum-gravity bridge sits on one of two legs, and the legs meet at i", ha="center", fontsize=13, weight="bold")
# legs
ax.add_patch(FancyBboxPatch((0.3, 5.0), 12.4, 2.8, boxstyle="round,pad=0.05", fc="#eef7ee", ec="#2e8b57", lw=1.2))
ax.text(0.5, 7.55, "ENTROPY LEG:  an additive entropy (quantum)  =  area / 4Għ (gravity)", fontsize=10.5, weight="bold", color="#1d5e3a")
ax.add_patch(FancyBboxPatch((0.3, 0.35), 12.4, 2.5, boxstyle="round,pad=0.05", fc="#eef3fb", ec="#2b6cb0", lw=1.2))
ax.text(0.5, 0.55, "PHASE LEG:  a quantum phase  =  action / ħ, with gravity's potential or clock inside it", fontsize=10.5, weight="bold", color="#1f4e85")
E = {"BH": (1.6, 6.6, "Bekenstein-Hawking\nS = A/4Għ (1973-75)", "proved"),
     "UN": (4.0, 6.6, "Unruh / Bisognano-Wichmann\nT = ħa/2πc (1975-76)", "proved"),
     "JA": (6.4, 6.6, "Jacobson 1995\nEinstein eq. from δQ = TδS", "proved"),
     "SV": (1.6, 5.55, "Strominger-Vafa 1996\ncounting of S", "proved"),
     "VE": (6.4, 5.55, "Verlinde 2011\nentropic force", "contest"),
     "RT": (8.9, 6.6, "Ryu-Takayanagi 2006\nS = min area/4G", "proved"),
     "DO": (11.3, 6.6, "Dong 2016: Rényi =\nbrane with tension", "proved"),
     "IS": (8.9, 5.55, "islands / Page curve\n(2019)", "proved"),
     "ER": (11.3, 5.55, "ER = EPR (2013)", "conj")}
P = {"CO": (1.6, 1.75, "COW 1975: gravity\nshifts neutron phase", "observed"),
     "AB": (4.0, 1.75, "gravitational\nAharonov-Bohm (2022)", "observed"),
     "PW": (6.4, 1.75, "Page-Wootters 1983:\ntime from a clock", "conj"),
     "BM": (8.8, 1.75, "gravity-made entanglement\n(BMV 2017, proposed)", "conj"),
     "AN": (11.2, 1.75, "analogue Hawking\nradiation (2016-19)", "observed")}
for d in (E, P):
    for key, (x, y, t, kd) in d.items(): box(ax, x, y, t, kd)
# AdS/CFT frame
ax.add_patch(FancyBboxPatch((7.65, 5.05), 4.9, 2.2, boxstyle="round,pad=0.04", fc="none", ec=C["conj"], lw=2, ls=":"))
ax.text(12.45, 7.32, "inside AdS/CFT (Maldacena 1997, conjectured, strong evidence)", ha="right", fontsize=8, color="#9a5c0b")
# the i hub
box(ax, 5.0, 3.95, "i : a boost is a rotation in imaginary time\n(temperature = κ/2π;  phase = e^{iS/ħ})", "proved", w=4.2, h=0.9, fs=9)
for p in [(4.0, 6.2), (1.9, 6.2), (1.6, 2.15), (6.4, 2.15), (11.3, 2.15)]: line(ax, (5.0, 3.95), p)
# ours
box(ax, 10.6, 3.95, "OURS: linear mixing forces k = 0\n(law a+b+k·ab; Tsallis k = 1-q,\nRényi lift 1+kS = Tr ρⁿ)", "ours", w=3.6, h=1.05, fs=8.6)
line(ax, (11.0, 4.48), (11.3, 6.2), dashed=True, col=C["ours"])
ax.text(11.45, 4.62, "k = 0 is n = 1:\nno brane tension", fontsize=7.6, color=C["ours"])
ax.text(7.2, 4.55, "the additive entropy that Jacobson\nand Ryu-Takayanagi assume", fontsize=7.6, color=C["ours"])
box(ax, 1.45, 3.95, "OURS: spinning river\n(Kerr from a shift by i)", "ours", w=2.3, h=0.8, fs=8.2)
line(ax, (2.6, 3.95), (2.9, 3.95), dashed=True, col=C["ours"])
# legend
for i, (kd, lab) in enumerate([("proved", "proved / derived"), ("observed", "observed"), ("conj", "conjectured or proposed"), ("contest", "contested"), ("ours", "ours (a reading, not a unification)")]):
    ax.add_patch(FancyBboxPatch((0.4+i*2.55, 8.2), 0.28, 0.22, boxstyle="round,pad=0.02", fc=C[kd], ec="black", lw=0.5))
    ax.text(0.75+i*2.55, 8.31, lab, va="center", fontsize=8)
fig.savefig("figures/bq1_bridge_map.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- bq2: the one dial ----
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.9))
nn = np.linspace(0.35, 3.0, 300)
ax = axs[0]
ax.axvspan(0.35, 1, color="#fde7f1"); ax.axvspan(1, 3.0, color="#fde7f1")
ax.plot(nn, (nn-1)/nn, color=C["proved"], lw=2.4, label="brane tension 4G·Tₙ = (n−1)/n  (Dong 2016)")
ax.axhline(0, color="black", lw=0.7); ax.axvline(1, color=C["ours"], lw=2.6)
ax.text(1.04, -1.55, "n = 1  (k = 0)\nRyu-Takayanagi,\nBekenstein-Hawking,\nJacobson's input:\nno tension,\nentropy additive", fontsize=8.3, color="#5a0b38")
ax.text(1.9, 0.12, "k < 0: positive tension,\nthe brane bends the\ngeometry into a cone", fontsize=8.3)
ax.text(0.38, 0.45, "k > 0: negative tension", fontsize=8.3)
ax.set_xlabel("Rényi order n  (= Tsallis q)"); ax.set_ylabel("4G × tension  (= deficit angle / 2π)")
sec = ax.secondary_xaxis("top", functions=(lambda x: 1-x, lambda x: 1-x)); sec.set_xlabel("the law's dial k = 1 − n")
ax.set_title("bq2a. The law's dial is a brane tension; k = 0 is the tensionless point", fontsize=10)
ax.legend(fontsize=8, loc="lower right"); ax.set_ylim(-1.9, 0.75)
ax = axs[1]
for p, col in [(0.5, "#444"), (0.8, "#2b6cb0"), (0.95, "#2e8b57")]:
    lift = p**nn + (1-p)**nn
    ax.plot(nn, lift, color=col, lw=2, label=f"qubit, eigenvalues ({p}, {1-p:.2f})")
ax.axvline(1, color=C["ours"], lw=2.6); ax.axhline(1, color="black", lw=0.7)
ax.text(1.05, 1.32, "lift = 1 at k = 0 for every state:\nthe law's lift is the\nreplica partition function\nZₙ / Z₁ⁿ = Tr ρⁿ", fontsize=8.3, color="#5a0b38")
ax.set_xlabel("Rényi order n"); ax.set_ylabel("lift 1 + k·S_T  =  Tr ρⁿ")
ax.set_title("bq2b. The lift 1 + kS is the replica trace (check B3)", fontsize=10); ax.legend(fontsize=8)
fig.suptitle("bq2. One dial runs through the entropy leg: Tsallis q, Rényi n, brane tension, and the law's k", fontsize=11.5, weight="bold")
fig.tight_layout(); fig.savefig("figures/bq2_one_dial.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- bq3: what correlation does to the law's k (setup for candidate C1) ----
import json; dev = json.load(open("results.json"))["B7_drift"]
cs = np.linspace(0, 0.95, 60)
def ts(p, q): p = p[p > 0]; return (1-(p**q).sum())/(q-1)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for q, col in [(2.0, C["proved"]), (3.0, C["observed"]), (0.5, C["conj"])]:
    ke = []
    for c in cs:
        P = (1-c)*np.full((2, 2), .25) + c*np.diag([.5, .5])
        sa = ts(P.sum(1), q); sb = ts(P.sum(0), q); ke.append((ts(P.ravel(), q)-sa-sb)/(sa*sb) - (1-q))
    ax.plot(cs, ke, color=col, lw=2, label=f"q = {q}")
ax.scatter([0, .1, .3, .6], dev, color="black", zorder=3, label="check B7 (q = 2)")
ax.axhline(0, color="black", lw=0.7)
ax.set_xlabel("correlation between the two pieces"); ax.set_ylabel("k_eff − (1 − q)")
ax.set_title("bq3. Correlation shows up as an extra term on the law's k\n(the reading candidate C1 would test on replica wormholes)", fontsize=9.5)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("figures/bq3_cross_term.png", dpi=150); plt.close(fig)
print("figures written")
