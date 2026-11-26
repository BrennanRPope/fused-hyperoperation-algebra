"""Pictures for liebig-tropical/README.md. Reads checks.json."""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
J = json.load(open(os.path.join(HERE, "checks.json")))
G, R, GR, B, K = "#2b8a3e", "#c92a2a", "#868e96", "#1c7ed6", "k"
CM = plt.get_cmap("viridis")


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name), dpi=140); plt.close(fig); print("wrote", name)


def clay(ta, tb, pv):
    if pv == 0:
        return ta * tb
    if np.isinf(pv):
        return np.minimum(ta, tb)
    sa, sb = -np.log(ta) * pv, -np.log(tb) * pv
    m = np.maximum(sa, sb)
    return np.exp(-(m + np.log(np.exp(sa - m) + np.exp(sb - m) - np.exp(-m))) / pv)


def su(ta, tb):
    ua, ub = 1 / ta - 1, 1 / tb - 1
    return 1 / (1 + ua + ub - ua * ub / (ua + ub + 1e-300))


# ------------------------------------------------------------- Fig 1: one knob, three classic laws
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.5))
ax = axs[0]
u = np.logspace(-2, 2, 400); th = 1 / (1 + u)
for kk in np.linspace(0, 1, 6):
    ax.loglog(u, 1 / (1 + 2 * u + kk * u * u), color=GR, lw=0.8, alpha=0.6)
ax.loglog([], [], color=GR, lw=0.8, label="old one dial, k from 1 to 0 (grey band)")
ps = [0, 1, 3, 10, np.inf]
names = {0: "p = 0: multiplicative", 1: "p = 1: additive", 3: "p = 3", 10: "p = 10", np.inf: "p = ∞: Liebig"}
for i, pv in enumerate(ps):
    col = CM(i / 4)
    ax.loglog(u, clay(th, th, pv), color=col, lw=2.2 if pv in (0, 1, np.inf) else 1.4,
              ls="--" if np.isinf(pv) else "-", label=names[pv])
ax.set_xlabel("scarcity u = K/S (both nutrients equally scarce)")
ax.set_ylabel("growth / maximum")
ax.legend(fontsize=8, loc="lower left")
ax.set_title("(a) The knob p carries the additive law to Liebig;\nthe old dial (grey) stops at additive", fontsize=9)
ax = axs[1]
ax.plot(J["dial_miss"]["k"], J["dial_miss"]["miss"], color=GR, lw=2, label="old one dial: best is k = 0, 17% off")
pp = np.logspace(0, 2.5, 60)
uu = np.logspace(-2, 2, 2000); tt = 1 / (1 + uu)
gap = [np.max(np.abs(clay(tt, tt, pv) - tt)) for pv in pp]
ax2 = ax.twiny()
ax2.semilogx(pp, gap, color=B, lw=2, label="knob family: gap closes as p grows")
ax2.set_xlabel("knob p (top axis, log scale)", color=B)
ax.set_xlabel("old dial k (bottom axis)", color=GR)
ax.set_ylabel("worst gap to Liebig (equal scarcity)")
ax.set_ylim(0, 0.75)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
ax.set_title("(b) The old closure was about the dial, not about Liebig", fontsize=9)
fig.suptitle("Fig 1. Turning the rate knob (a symmetry of the multiplicative law) slides the additive law all the way to Liebig", fontsize=10)
save(fig, "fig1-one-knob.png")

# ------------------------------------------------------------- Fig 2: isoclines, the ecologists' picture
fig, axs = plt.subplots(1, 4, figsize=(13, 4.0))
T = np.linspace(0.005, 1, 300); TA, TB = np.meshgrid(T, T)
laws = [("multiplicative (p = 0)", TA * TB), ("additive (p = 1)", clay(TA, TB, 1)),
        ("synthesizing unit (≈ p = 3.5)", su(TA, TB)), ("Liebig (p = ∞)", np.minimum(TA, TB))]
for ax, (name, Z) in zip(axs, laws):
    cs = ax.contour(TA, TB, Z, levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], cmap="viridis")
    ax.clabel(cs, fontsize=6, fmt="%.1f")
    ax.set_title(name, fontsize=9); ax.set_aspect("equal")
    ax.set_xlabel("θ_A (nutrient A alone)")
axs[0].set_ylabel("θ_B (nutrient B alone)")
fig.suptitle("Fig 2. Growth isoclines: the knob bends them from hyperbolas (multiplicative) to Liebig's right angles", fontsize=10)
save(fig, "fig2-isoclines.png")

# ------------------------------------------------------------- Fig 3: one line of dependence
fig, ax = plt.subplots(figsize=(10.5, 3.6))
ax.set_xlim(-0.3, 1.25); ax.set_ylim(-1.3, 1.5); ax.axis("off")
ax.annotate("", (1.03, 0), (-0.27, 0), arrowprops=dict(arrowstyle="->", lw=1.2))
ax.text(1.04, 0.08, "Kendall's τ: how tied\nthe two shortages are", fontsize=8, ha="left", va="bottom")
kk = np.linspace(2, 0, 50); taus_d = []
for kv in kk:
    a_ = 1 - kv
    # closed form for AMH tau
    if abs(a_) < 1e-9:
        taus_d.append(0.0)
    else:
        taus_d.append(1 - 2 * ((1 - a_) ** 2 * math.log(1 - a_) + a_) / (3 * a_ ** 2) if a_ < 1 else 1 / 3)
ax.plot(taus_d, [0.5] * len(kk), color=GR, lw=6, solid_capstyle="butt")
ax.text(-0.27, 0.75, "old one dial (Ali–Mikhail–Haq), k from 2 to 0", color=GR, fontsize=8.5)
pt = np.linspace(0, 1, 50)
ax.plot(pt, [-0.5] * 50, color=B, lw=6, solid_capstyle="butt")
ax.text(-0.27, -0.95, "knob family (Clayton), p from 0 to ∞:  τ = p/(p + 2)", color=B, fontsize=8.5)
for tv, lab, yy in [(0, "multiplicative\n(independent)", 1.0), (1 / 3, "additive\n(both tracks meet)", 1.0),
                    (J["su"]["tau"], "synthesizing\nunit", -1.25), (1, "Liebig\n(fully tied)", 1.0)]:
    ax.plot([tv, tv], [-0.6, 0.6], color=K, lw=0.8, ls=":")
    ax.plot(tv, 0, "o", color=K, ms=7)
    ax.text(tv, yy if yy > 0 else -0.15, lab, ha="center", fontsize=8, va="bottom" if yy > 0 else "top")
ax.text(min(taus_d) - 0.0, 0.18, f"{min(taus_d):.2f}", fontsize=7, color=GR, ha="center")
ax.set_title("Fig 3. Read as copulas, the two families lie on one line of dependence. The old dial cannot pass τ = 1/3; Liebig sits at τ = 1", fontsize=10)
save(fig, "fig3-dependence-line.png")

# ------------------------------------------------------------- Fig 4: Maslov on the height ruler
fig, ax = plt.subplots(figsize=(7.5, 4.3))
sA = np.linspace(0, 3, 300); sB = 1.5
for i, pv in enumerate([1, 2, 5, 20]):
    h = 1 / pv
    st = h * np.log(np.exp(sA / h) + np.exp(sB / h) - 1)
    ax.plot(sA, st, color=CM(i / 4), lw=1.8, label=f"temperature h = 1/p = {h:g}")
ax.plot(sA, np.maximum(sA, sB), "k--", lw=2, label="Liebig: max (tropical stop)")
ax.set_xlabel("height s_A = ln(1/θ_A)   (s_B fixed at 1.5)")
ax.set_ylabel("height of the combined shortage, ln(1/θ)")
ax.legend(fontsize=8)
ax.set_title("Fig 4. On the height ruler the knob is a temperature: a softened max\nthat cools to Liebig, never more than h·ln 2 above it", fontsize=10)
save(fig, "fig4-maslov.png")
