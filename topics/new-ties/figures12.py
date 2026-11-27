"""Round 12 pictures (DRAFT): fig38 the step function, fig39 one cluster up close, fig40 the cluster shrinking."""
from fractions import Fraction as F
import json, math, random, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cover12 as C

cov = json.load(open("cover12.json"))
od = json.load(open("original_descent.json"))
bad = [p[2] for p in od["O3_descent"]["per_pair"]]

def ev(w, z):
    for x in reversed(w):
        s, c = C.LET[x]; z = s * z + c
    return z
def rand_word(rng, start, L):
    w = start; st = C.nxt("none", start)
    while len(w) < L:
        o = [c for c in C.REDUCED[st] if c != "E"]; x = rng.choice(o); w += x; st = C.nxt(st, x)
    return w
rng = random.Random(38); ratios = []
for _ in range(20000):
    X = rand_word(rng, "N", 24); Y = rand_word(rng, rng.choice("NK"), 24)
    ratios.append(float(ev(X, F(1, 18)) / ev(Y, F(1, 18))))

def panel(ax, name, lo=None, hi=None, rug=True):
    w = cov["windows"][name]; a, b = map(float, map(F, w["window"]))
    lo = a if lo is None else lo; hi = b if hi is None else hi
    for p, q, n in w["proved_pieces"]:
        p, q = float(F(p)), float(F(q))
        if q < lo or p > hi: continue
        ax.plot([max(p, lo), min(q, hi)], [math.log10(max(n, 1))] * 2, color="#1f77b4", lw=1.6, solid_capstyle="butt")
    for p, q in w["open"]:
        p, q = float(F(p)), float(F(q))
        if q < lo or p > hi: continue
        ax.axvspan(max(p, lo), min(q, hi), color="#d62728", alpha=0.25, lw=0)
    if rug:
        rr = [r for r in ratios if lo <= r <= hi]
        ax.plot(rr, [-0.35] * len(rr), "|", color="black", ms=6, alpha=0.4)
        bb = [r for r in bad if lo <= r <= hi]
        ax.plot(bb, [4.15] * len(bb), "|", color="#2ca02c", ms=6, alpha=0.5)
    ax.set_xlim(lo, hi); ax.set_ylim(-0.7, 4.5)
    ax.set_ylabel("search size for the range\n(log10 of nodes)")

# fig38
fig, axs = plt.subplots(2, 1, figsize=(11, 7.2))
panel(axs[0], "N"); panel(axs[1], "K")
axs[0].set_title("Figure 38. The step function: each blue step is a range of stretches proved free of a second collision\n"
                 "by ONE search tree (height = how big that tree was). Red = not proved yet.", fontsize=10.5, loc="left")
axs[1].set_xlabel("stretch ρ = λ^a / μ^b")
axs[0].text(0.003, 0.97, "N window [0.9, 1.111]", transform=axs[0].transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="none", alpha=0.8))
axs[1].text(0.003, 0.97, "K window [0.474, 0.585]", transform=axs[1].transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="none", alpha=0.8))
axs[1].text(0.995, 0.03, "black ticks (bottom): points x/y where two pieces of the fractal line up (the set R)\n"
            "green ticks (top): the 1,406 bad pairs with a ≤ 10,000 (Round 11)", transform=axs[1].transAxes,
            ha="right", va="bottom", fontsize=8.5, bbox=dict(fc="white", ec="0.7", alpha=0.9))
fig.tight_layout(); fig.savefig("figures/fig38-step-function.png", dpi=140); plt.close(fig)

# fig39
fig, ax = plt.subplots(figsize=(11, 4.2))
panel(ax, "N", 0.944, 0.963)
ax.set_xlabel("stretch ρ")
ax.set_title("Figure 39. One cluster up close (0.944 to 0.963). The red gaps sit exactly where the black ticks are:\n"
             "points where two pieces of the fractal line up. Between them the steps close.", fontsize=10.5, loc="left")
ax.axvline(18 / 19, color="0.3", ls=":", lw=1)
ax.text(18 / 19, 3.7, " 18/19: N(N(0)) = (18/19)·N(1/9)", fontsize=8.5, color="0.2")
fig.tight_layout(); fig.savefig("figures/fig39-cluster-zoom.png", dpi=140); plt.close(fig)

# fig40
import glob, numpy as np
pts = []
for fn in glob.glob("cluster12_[0-9]*.txt"):
    f_ = open(fn).readline().split()
    if f_ and f_[0] == "k": pts.append((2.0 ** -int(f_[1]), float(f_[5]) / float(f_[7])))
pts.sort()
fig, ax = plt.subplots(figsize=(6.6, 4.2))
xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
ax.plot(xs, 100 * ys, "o-", color="#d62728")
for x_, y_ in pts: ax.annotate("%.0f%%" % (100 * y_), (x_, 100 * y_), textcoords="offset points", xytext=(4, 7), fontsize=9)
slope = None
if len(pts) >= 3:
    fit = ys < 1   # the coarsest depth has nothing proved yet
    slope = np.polyfit(np.log(xs[fit]), np.log(ys[fit]), 1)[0]
    ax.text(0.03, 0.06, "unproved share ≈ (width)^%.3f\n(so the line-up set R behaves like dimension ≈ %.2f here)" % (slope, 1 - slope),
            transform=ax.transAxes, fontsize=9, bbox=dict(fc="white", ec="0.7"))
ax.set_xscale("log"); ax.set_ylim(0, 100)
ax.set_xlabel("narrowest range allowed (halving stops here)"); ax.set_ylabel("unproved share of 0.944–0.963 (%)")
ax.set_title("Figure 40. Finer ranges prove more of a cluster, but slowly.\nA 4× bigger search per range proves nothing extra.", fontsize=10.5, loc="left")
fig.tight_layout(); fig.savefig("figures/fig40-cluster-shrinks.png", dpi=140); plt.close(fig)
print("figures written", pts, "slope", slope)
