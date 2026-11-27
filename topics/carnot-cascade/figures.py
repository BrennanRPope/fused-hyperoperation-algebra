"""Pictures for carnot-cascade/README.md. Self-contained (numpy + matplotlib)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
C1, C2, C3, C4, C5, GR = "#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#868e96"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name), dpi=140); plt.close(fig); print("wrote", name)

eC = np.linspace(0, 0.999, 400); tau = 1 - eC
# ------------------------------------------------------------------ Fig 1: the known efficiencies
fig, ax = plt.subplots(figsize=(7.8, 4.8))
ax.fill_between(eC, eC/2, eC/(2 - eC), color=C1, alpha=0.10, lw=0)
ax.plot(eC, eC, color="k", lw=1.2, label="Carnot  η_C")
ax.plot(eC, eC/(2 - eC), color=C1, lw=2, label="upper bound  η_C/(2 − η_C)")
ax.plot(eC, eC/(2 - eC/2), color=C1, lw=1.2, ls="--", label="heat-type, symmetric  η_C/(2 − η_C/2)")
ax.plot(eC, eC/2, color=C1, lw=2, ls=":", label="lower bound  η_C/2")
ax.plot(eC, 1 - np.sqrt(tau), color=C2, lw=2.5, label="Curzon–Ahlborn  1 − √(1 − η_C)")
with np.errstate(invalid="ignore"): fr = (1 - tau)**2/((1 - tau) - tau*np.log(tau))
ax.plot(eC, fr, color=C3, lw=1.5, label="Feynman ratchet (Tu 2008)")
ax.plot(eC, 1 - np.sqrt(tau*(1 + tau)/2), color=C4, lw=1.5, label="ecological (Angulo-Brown 1991)")
ax.annotate("low-dissipation band\nat maximum power", (0.86, 0.58), (0.62, 0.12), color=C1, fontsize=8,
            arrowprops=dict(arrowstyle="->", color=C1, lw=0.8))
ax.set_xlabel("Carnot efficiency η_C = 1 − T_c/T_h"); ax.set_ylabel("efficiency")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.set_title("The known finite-time efficiencies. Only Curzon–Ahlborn is a knob image of Carnot:\n"
             "1 − (1 − η_C)^λ at λ = ½, the only kind that two engines in series reproduce", fontsize=9)
save(fig, "fig1-efficiencies.png")

# ------------------------------------------------------------------ Fig 2: where they sit on the dial
fig, ax = plt.subplots(figsize=(7.8, 5.0))
kk = np.linspace(-3, 3, 300)
ax.plot(kk, kk**2/4, color=GR, lw=1, label="border parabola  k² = 4d (one double wall)")
gs = np.linspace(0, 1, 100); kf = 2*gs - 2; df = -gs*(2 - gs)
ax.plot(kk, kk**2/4 - 1, color=C1, lw=0.8, ls=":", alpha=0.6)
ax.plot(kf, df, color=C1, lw=3, label="low-dissipation family η_C/(2 − γη_C), γ from 0 to 1:\n d = k²/4 − 1, the border moved down by 1")
pts = [(-2, 0, "γ = 0: lower bound η_C/2\n= parity, x + y − 2xy", C1, (-2.9, 0.45)),
       (0, -1, "γ = 1: upper bound η_C/(2 − η_C)\n= velocity addition", C1, (0.15, -1.45)),
       (-1, -0.75, "γ = ½: symmetric heat-type", C1, (-2.4, -1.05)),
       (-1, 0, "Carnot (full speed) and\nCurzon–Ahlborn (half speed):\nthe union line k = −1", C2, (-0.9, 0.35))]
for k_, d_, lab, c, xy in pts:
    ax.plot(k_, d_, "o", color=c, ms=9, mec="white", mew=1.5)
    ax.annotate(lab, (k_, d_), xy, fontsize=7.5, arrowprops=dict(arrowstyle="-", color=GR, lw=0.7))
ax.axhline(0, color=GR, lw=0.6); ax.axvline(0, color=GR, lw=0.6)
ax.set_xlim(-3, 2); ax.set_ylim(-1.7, 1.3)
ax.set_xlabel("k  (how full the corner cell is)"); ax.set_ylabel("d  (how much the next layer bends)")
ax.legend(fontsize=7.5, frameon=False, loc="lower left")
ax.set_title("Read in β = ln(T_h/T_c), every low-dissipation bound is a plane law,\nand the whole family lies on one parabola parallel to the border", fontsize=9)
save(fig, "fig2-dial.png")

# ------------------------------------------------------------------ Fig 3: fine cascades flow to Curzon-Ahlborn
etaC = 0.8; B = -np.log(1 - etaC); ns = np.unique(np.round(np.logspace(0, 3, 40)).astype(int))
fam = [("lower bound η_C/2", lambda t: (1 - t)/2, C1, ":"),
       ("upper bound η_C/(2 − η_C)", lambda t: (1 - t)/(1 + t), C1, "-"),
       ("heat-type symmetric", lambda t: (1 - t)/(2 - (1 - t)/2), C1, "--"),
       ("Feynman ratchet", lambda t: (1 - t)**2/((1 - t) - t*np.log(t)), C3, "-"),
       ("ecological (not max power)", lambda t: 1 - np.sqrt(t*(1 + t)/2), C4, "-")]
fig, ax = plt.subplots(figsize=(7.8, 4.6))
for lab, fn, c, ls in fam:
    ax.plot(ns, [1 - (1 - fn(np.exp(-B/n)))**n for n in ns], color=c, ls=ls, lw=2, label=lab)
ax.axhline(1 - np.sqrt(1 - etaC), color=C2, lw=2.5, label=f"Curzon–Ahlborn {1 - np.sqrt(1 - etaC):.3f} (unchanged by splitting)")
ax.axhline(1 - (1 - etaC)**0.75, color=C4, lw=1, ls=":", label="1 − (1 − η_C)^¾")
ax.set_xscale("log"); ax.set_xlabel("number of stages in the cascade, each run at its own maximum power")
ax.set_ylabel("overall efficiency (η_C = 0.8)")
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
ax.set_title("Split the engine into ever more max-power stages: every slope-½ rule flows to Curzon–Ahlborn,\n"
             "the fixed point of cascading", fontsize=9)
save(fig, "fig3-cascade-flow.png")

# ------------------------------------------------------------------ Fig 4: Riccati test bars
J = json.load(open(os.path.join(HERE, "checks.json")))["_D_residuals"]
names = list(J.keys()); vals = [max(J[n], 1e-45) for n in names]
fig, ax = plt.subplots(figsize=(7.8, 4.2))
cols = [C1 if v < 1e-30 else C5 for v in vals]
y = np.arange(len(names))[::-1]
ax.barh(y, np.log10(vals) + 45, left=-45, color=cols, height=0.6)
ax.axvline(-4, color=GR, lw=0.8, ls=":"); ax.text(-3.8, len(names) - 0.6, "off the dial →", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7.5)
ax.set_xlabel("log₁₀ of the miss when d η/d β is fitted by a quadratic in η (zero miss = a plane law)")
ax.set_xlim(-45, 0)
ax.set_title("Which efficiencies are plane laws? Blue: exact (miss at the 40-digit floor). Pink: not.", fontsize=9)
save(fig, "fig4-riccati-test.png")
