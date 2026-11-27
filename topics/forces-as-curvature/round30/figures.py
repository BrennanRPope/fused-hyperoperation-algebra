"""DRAFT round 30 picture: f37_lump_tide.png, from results.json."""
import json, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.collections import LineCollection
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8"

fig = plt.figure(figsize=(20, 13))
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.18)
fig.suptitle("Round 30 (DRAFT): a lump's own tidal pull, built from the law's local rule plus fitting", fontsize=16)

# A
ax = fig.add_subplot(gs[0, 0]); L = R["L2"]
r = np.array(L["r"]); qb = np.array(L["qbar"]); ql = np.array(L["qloc"]); tide = np.array(L["tide"])
along = qb - tide
ax.plot(r, ql, color=BLUE, lw=2.5, label="local dial q_loc = q + 8πρ/3 (the law's rule: set by the matter here)")
ax.plot(r, qb, color=ORANGE, lw=2.5, label="bending across the radius, q̄ = q_loc averaged over the ball inside")
ax.plot(r, along, color=RED, lw=2.5, ls="--", label="bending along the radius")
ax.fill_between(r, along, qb, color=ORANGE, alpha=0.15)
ax.axhline(L["q0"], color=GREY, lw=1, ls=":")
ax.text(5.9, L["q0"] + 0.006, "the law's background q", color=GREY, fontsize=10, ha="right")
i = np.searchsorted(r, 2.6)
ax.annotate("shaded gap = the tide,\n3/2 (q̄ − q_loc)\n(outside: 3M/r³)", (r[i], (qb[i] + along[i]) / 2), (3.4, 0.2), fontsize=11,
            color=ORANGE, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.set_xlim(0, 6); ax.set_xlabel("areal radius r"); ax.set_ylabel("curvature of a plane (dial)")
ax.set_title("A. One smooth lump: three dials at each place (L1, L2)")
ax.legend(loc="upper right", fontsize=10)

# B: chain
ax = fig.add_subplot(gs[0, 1]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("B. Where the tide comes from (L2)")
boxes = [(2.45, 8.4, "matter here sets the mean bending\n(the law's local rule, rounds 23 and 25)", BLUE),
         (7.55, 8.4, "every geometry fits together\n(Bianchi identity: no physics, an identity)", GREY),
         (5.0, 5.4, "together: (r³ q̄)′ = 3 r² q_loc\nacross-bending q̄ = q_loc averaged over the ball inside", ORANGE),
         (5.0, 2.4, "tide = across − along = 3/2 (q̄ − q_loc)\noutside a lump: 3M/r³", RED)]
for xx, yy, txt, col in boxes:
    hw = 2.25 if yy > 8 else 3.6
    ax.add_patch(FancyBboxPatch((xx - hw, yy - 0.8), 2 * hw, 1.6, boxstyle="round,pad=0.1", fc="white", ec=col, lw=2))
    ax.text(xx, yy, txt, ha="center", va="center", fontsize=11, color=col)
for (x0, y0), (x1, y1) in [((2.45, 7.5), (4.2, 6.3)), ((7.55, 7.5), (5.8, 6.3)), ((5.0, 4.5), (5.0, 3.3))]:
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="->", color=GREY, lw=2))
ax.text(5.0, 0.6, "Nothing else is needed for a lump at rest, once the dial may differ by direction.", ha="center", fontsize=11, style="italic")

# C: two lumps
ax = fig.add_subplot(gs[1, 0]); L4 = R["L4"]; F = np.array(L4["field"])
mag = np.abs(F[:, 4]); lw = np.clip(np.log10(mag / mag.min()) / np.log10(mag.max() / mag.min()), 0.05, 1)
segs = []
for (xx, zz, vx, vz, w, wmax), s in zip(F, lw):
    l = 0.16 + 0.12 * s
    segs.append([(xx - l * vx, zz - l * vz), (xx + l * vx, zz + l * vz)])
lc = LineCollection(segs, colors=plt.cm.Oranges(0.35 + 0.65 * lw), linewidths=1 + 2.2 * lw)
ax.add_collection(lc)
for zc, m in [(L4["d"] / 2, L4["m1"]), (-L4["d"] / 2, L4["m2"])]:
    ax.add_patch(Circle((0, zc), 0.45 * math.sqrt(m), color=BLUE, alpha=0.85))
    ax.text(0.65, zc, f"lump, mass {m:g}", color=BLUE, fontsize=11, va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=2))
ax.set_xlim(-5.3, 5.3); ax.set_ylim(-5.3, 5.3); ax.set_aspect("equal")
ax.set_xlabel("x"); ax.set_ylabel("z")
ax.set_title("C. Two lumps at rest: direction of strongest stretch (L4)\nlocal dial is the vacuum value everywhere shown; darker = stronger tide")

# D: count
ax = fig.add_subplot(gs[1, 1])
parts = [(3, GREY, "only the choice of coordinates"), (1, BLUE, "the law's local rule\n(mean bending set by matter)"),
         (2, PURPLE, "free: gravitational-wave data\n(the law doesn't fix these)")]
left = 0
for w, col, lab in parts:
    ax.barh(0, w, left=left, color=col, alpha=0.8, height=0.5)
    ax.text(left + w / 2, 0.45 if w != 1 else -0.45, lab, ha="center", va="bottom" if w != 1 else "top", color=col, fontsize=11, weight="bold")
    ax.text(left + w / 2, 0, str(w), ha="center", va="center", color="white", fontsize=16, weight="bold")
    left += w
ax.set_xlim(0, 6); ax.set_ylim(-1.45, 1.1); ax.set_yticks([])
ax.set_xlabel("the 6 functions that give a slice its shape, at each point")
ax.set_title("D. What is left over without round symmetry (L3)")
zz = np.linspace(0.2, 5.8, 300)
ax.plot(zz, -1.0 + 0.12 * np.cos(2 * np.pi * zz / 1.4), color=PURPLE, lw=2)
ax.text(3.0, -1.38, "a ripple in those 2 keeps every local dial unchanged yet stretches by direction: a wave's tide (round 26)",
        ha="center", fontsize=10, color=PURPLE)
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "f37_lump_tide.png"), dpi=120, bbox_inches="tight")
