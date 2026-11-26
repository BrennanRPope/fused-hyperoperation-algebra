"""Pictures for fractional-steps/README.md. Reads grid_output.json, far_output_*.json, far2_t*.json."""
import glob
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
G, R, GR = "#2b8a3e", "#c92a2a", "#868e96"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=140)
    plt.close(fig)
    print("wrote", name)


grid = json.load(open(os.path.join(HERE, "grid_output.json")))

# ---- fig 1: half a generation
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
pois = [math.exp(-1) / math.factorial(k) for k in range(8)]
axs[0].bar(range(8), pois, color=G)
axs[0].set_title("(a) One whole generation: Poisson(1) children\nall chances positive", fontsize=10)
axs[0].set_xlabel("number of descendants"); axs[0].set_ylabel("chance")
half = grid["0.5"]["coeffs"][:8]
axs[1].bar(range(8), half, color=[G if v >= 0 else R for v in half])
axs[1].axhline(0, color="k", lw=0.6)
axs[1].annotate(f"P(4) = {half[4]:.5f}\nP(5) = {half[5]:.5f}", (4, half[4]), xytext=(4.6, 0.25),
                color=R, fontsize=9, arrowprops=dict(arrowstyle="->", color=R))
axs[1].set_title("(b) Half a generation, the only candidate\nP(4) and P(5) come out negative", fontsize=10)
axs[1].set_xlabel("number of descendants")
fig.suptitle("Fig 1. The half step is not a real chance distribution", fontsize=11)
save(fig, "fig1-half-generation.png")

# ---- fig 2: generator coefficients
import flow  # noqa: E402
mp.mp.dps = 40
uc = [float(c.real) for c in flow.taylor(flow.u, M=256, rho=0.95)[:70]]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
ax = axs[0]
ks = np.arange(len(uc))
cols = [GR if k == 1 else (G if v >= 0 else R) for k, v in zip(ks, uc)]
ax.bar(ks, uc, color=cols)
ax.set_yscale("symlog", linthresh=1e-16)
ax.set_yticks([-1, -1e-6, -1e-12, 0, 1e-12, 1e-6, 1])
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("k (number of children produced in one event)")
ax.set_ylabel("rate u_k")
ax.set_title("(a) The only possible generator u(s) = 1/alpha'(s)\nred = negative rate (k = 3 is the first; grey k = 1 may be negative)", fontsize=10)
ax = axs[1]
sk = np.sqrt(ks[2:])
ax.scatter(sk, np.sign(uc[2:]) * np.abs(np.array(uc[2:])) ** 0.08, c=[G if v >= 0 else R for v in uc[2:]], s=14)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("square root of k")
ax.set_ylabel("sign x |u_k|^0.08 (squashed)")
ax.set_title("(b) The signs swing like a wave in sqrt(k),\nnegative rates keep coming as far as we checked", fontsize=10)
fig.suptitle("Fig 2. Why no continuous clock: a continuous process would need every red bar to be >= 0", fontsize=11)
save(fig, "fig2-generator.png")

# ---- fig 3: which step sizes work (deepest, most precise search wins for each t)
best = {}  # t -> (rank, first_neg or None, top_of_search)
def note(t, v, top, rank):
    t = round(float(t), 6)
    c = [float(x) for x in v["coeffs"]]
    if rank == 0:      # grid.py already dropped values at its noise level
        neg = v["neg"]
    else:
        # far.py (45 digits): noise near 1e-40; far2.py at 50 digits: read the noise off the run's own far tail
        if rank == 1:
            floor = 1e-38
        elif v.get("dps", 50) >= 90:   # the whole-step control at 90 digits stays smooth and positive past 1e-60
            floor = 1e-60
        else:
            floor = 10 * max(abs(x) for x in c[3 * len(c) // 4:])
        neg = [k for k in v["neg"] if abs(c[k]) > floor]
        top = max([k for k in range(len(c)) if abs(c[k]) > floor], default=top)  # deepest k above noise
    if t not in best or rank > best[t][0]:
        best[t] = (rank, neg[0] if neg else None, top)
for t, v in grid.items():
    if t != "alpha_prime":
        note(t, v, 55, 0)
for f in glob.glob(os.path.join(HERE, "far_output_*.json")):
    for t, v in json.load(open(f)).items():
        note(t, v, 479, 1)
for f in glob.glob(os.path.join(HERE, "far2_t*.json")):
    v = json.load(open(f))
    note(v["t"], v, len(v["coeffs"]) - 1, 2 + int(v.get("dps", 50)))
fig, ax = plt.subplots(figsize=(10, 4.8))
fit_t, fit_k = [], []
for t, (rank, k, top) in sorted(best.items()):
    whole = abs(t - round(t)) < 1e-9
    if k is not None:
        ax.plot(t, k, "o", color=R, ms=4)
        if t >= 1:
            fit_t.append(t); fit_k.append(k)
    elif whole:
        ax.plot(t, top, "^", color=G, ms=8)
    else:
        ax.plot(t, top, "^", color="#f08c00", ms=8)
cfit = sum(k * t * t for t, k in zip(fit_t, fit_k)) / sum(t ** 4 for t in fit_t)
tt = np.linspace(1, max(best) + 0.5, 200)
ax.plot(tt, cfit * tt ** 2, "--", color="k", lw=0.8, label=f"rule of thumb for steps of 1 or more: first failure near k = {cfit:.1f} t^2")
ax.set_yscale("log")
ax.set_xlabel("step size t (generations)")
ax.set_ylabel("first descendant count k with a negative chance")
ax.plot([], [], "o", color=R, label="fractional step: negative chance found at this k")
ax.plot([], [], "^", color=G, label="whole step: clean up to the top of the search")
ax.plot([], [], "^", color="#f08c00", label="fractional step: search not deep enough yet to see it")
ax.legend(fontsize=8, loc="lower right")
ax.set_title("Fig 3. Every fractional step we could settle fails; the first failure moves out like t squared", fontsize=10)
save(fig, "fig3-which-steps.png")
print("fit constant", cfit)

# ---- fig 4: the Abel ruler's slope is positive (why whole steps stay fine far out)
ap = grid["alpha_prime"][:40]
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(range(len(ap)), ap, "o", color=G, ms=4, label="Taylor coefficients of alpha'(s)")
kk = np.arange(len(ap))
ax.plot(kk, 2 * kk + 5 / 3, "k--", lw=0.8, label="2k + 5/3 from 2/(1-s)^2 - (1/3)/(1-s)")
ax.set_xlabel("k"); ax.legend(fontsize=8)
ax.set_title("Fig 4. Check: the ruler's slope has all positive coefficients,\nwith the double pole and the 1/3 residue the height flow already knew", fontsize=10)
save(fig, "fig4-ruler-check.png")
