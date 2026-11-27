"""DRAFT round 34 picture: f41_waves_on_pieces.png, from results.json."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.18)
fig.suptitle("Round 34 (DRAFT): gravitational waves on spacetime glued from flat pieces", fontsize=16)

# A: a 2D slice of the gluing with a + ripple on the edge lengths
ax = fig.add_subplot(gs[0, 0])
n = 9; segs, strain = [], []
kx = 2 * np.pi / 6
for i in range(n):
    for j in range(n):
        for (di, dj) in ((1, 0), (0, 1), (1, 1)):
            if i + di < n and j + dj < n:
                p, q = np.array([i, j]), np.array([i + di, j + dj]); u = (q - p) / np.linalg.norm(q - p)
                mid = (p + q) / 2
                s = np.cos(kx * mid[0]) * u[1]**2          # + wave along x: h_yy = cos(kx), edge strain = h_yy u_y^2
                segs.append([p, q]); strain.append(s)
lc = LineCollection(segs, cmap="coolwarm", linewidths=3); lc.set_array(np.array(strain)); lc.set_clim(-1, 1)
ax.add_collection(lc)
ax.set_xlim(-0.5, n + 0.6); ax.set_ylim(-0.5, n + 1.6); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.annotate("", xy=(n - 1, n + 0.6), xytext=(0, n + 0.6), arrowprops=dict(arrowstyle="->", color=GREY, lw=2))
ax.text(0, n + 0.9, "wave moving →  edges stretch (red), shrink (blue)", color=GREY, fontsize=10.5)
ax.set_title("A. What was built: flat pieces glued on a lattice, with a ripple in their edge lengths\n"
             "(a 2D slice; in 4D each cell is 24 flat pieces, 15 edges per corner, 50 kinds of joint)")
ax.text(-0.3, -1.4, "curvature lives only on the joints, as the angle a loop comes back turned by (round 11)", fontsize=11, color=BLUE)

# B: long-wavelength match to Einstein's action
ax = fig.add_subplot(gs[0, 1])
rows = R["G2"]["rows"]; cols = [BLUE, ORANGE, PURPLE, GREY]
for i in range(0, len(rows), 4):
    r = rows[i:i + 4]
    ax.loglog([x["k"] for x in r], [x["rel"] for x in r], "o-", color=cols[i // 4], lw=2, label=f"direction {tuple(r[0]['dir'])}")
ks = np.array([0.005, 0.04]); ax.loglog(ks, 0.3 * ks, "k--", lw=1, label="falls in proportion to k")
ax.set_xlabel("wavenumber k (in units of 1 / gluing length)"); ax.set_ylabel("mismatch with Einstein's action, |E + F/4| / |F/4|")
ax.set_title("B. At long wavelengths the pieces' action is Einstein's (G2)\nRegge's joint size × deficit → ½∫√g R exactly for ripples (sympy ratio −1/8)")
ax.legend(fontsize=10, loc="lower right")

# C: the two new flat directions at real frequency
ax = fig.add_subplot(gs[1, 0])
w = np.array(R["G3"]["scan_w"]); S = np.array(R["G3"]["scan"])     # 7 smallest, descending
ax.semilogy(w, S[:, 0], color=BLUE, lw=2.5, label="new flat direction 1 (+ polarization)")
ax.semilogy(w, S[:, 1], color=ORANGE, lw=2.5, ls="--", label="new flat direction 2 (× polarization)")
for j in range(2, 7): ax.semilogy(w, np.maximum(S[:, j], 1e-12), color=GREY, lw=1, alpha=0.7)
ax.axvline(R["G3"]["cone_w"], color=PURPLE, lw=1.5, ls=":")
ax.text(R["G3"]["cone_w"] + 0.0006, 3e-3, "lattice light cone", color=PURPLE, fontsize=11)
ax.text(0.986, 1.5e-11, "always flat: 4 vertex moves and the lattice's flex", color=GREY, fontsize=11)
ax.set_ylim(1e-12, 1e-1)
ax.set_xlabel("frequency / wavenumber  (ω / k, light speed = 1), at ka = 0.1 along an axis")
ax.set_ylabel("how far each direction is from flat (singular value)")
ax.set_title("C. At real frequency two more flat directions open, both on one cone (G3)\nthey are the + and × shapes of a gravitational wave; no third polarization")
ax.legend(fontsize=10, loc="center left")

# D: the cone, by direction
ax = fig.add_subplot(gs[1, 1])
names = ["axis (1,0,0)", "face (1,1,0)", "body (1,1,1)", "generic (2,-1,1/2)"]
vec = {"axis (1,0,0)": [1, 0, 0], "face (1,1,0)": [1, 1, 0], "body (1,1,1)": [1, 1, 1], "generic (2,-1,1/2)": [2, -1, .5]}
kk = np.linspace(0.001, 1.05, 200)
for c, nm in zip(cols, names):
    d = np.array(vec[nm]) / np.linalg.norm(vec[nm])
    cone = [2 * np.arcsinh(np.sqrt(np.sum(np.sin(x * d / 2)**2))) / x for x in kk]
    ax.plot(kk, cone, color=c, lw=1.5, alpha=0.7)
    r = [x for x in R["G4"]["rows"] if x["dir"] == nm]
    ax.plot([x["k"] for x in r], [x["speed"] for x in r], "o", color=c, ms=9, label=f"{nm}: found (dots), plain lattice wave (line)")
ax.axhline(1, color=GREY, ls="--", lw=1); ax.text(0.02, 1.003, "light speed", color=GREY, fontsize=11)
ax.set_xlabel("wavenumber × gluing length, ka"); ax.set_ylabel("wave speed ω / k")
ax.set_ylim(0.9, 1.015)
ax.set_title("D. The waves move exactly like a plain wave on the same lattice (G4)\nlight speed at long wavelengths, lagging by (ka)²/12 to (ka)²/18")
ax.text(0.40, 1.0035, f"GW170817 (lag ≤ 3×10⁻¹⁵ at 100 Hz) only needs a ≲ {100*R['G4']['a_max_m']:.0f} cm;\nat the Planck length the lag is ~10⁻⁸³", color=RED, fontsize=11)
ax.legend(fontsize=9.5, loc="lower left", bbox_to_anchor=(0, 0.12))
fig.savefig(os.path.join(HERE, "figures", "f41_waves_on_pieces.png"), dpi=110, bbox_inches="tight")
