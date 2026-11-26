"""Pictures for where-k-moves/README.md. Reads checks.json (written by checks.py)."""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
G, R, GR, B, P = "#2b8a3e", "#c92a2a", "#868e96", "#1c7ed6", "#7048e8"
J = json.load(open(os.path.join(HERE, "checks.json")))


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=140)
    plt.close(fig)
    print("wrote", name)


# ------------------------------------------------------------- Fig 1: frozen bath, chances per event
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.3))
ax = axs[0]
ex = J["A_defect_examples"]
labels = ["off + one on level\n(cure model)", "two on levels\n(0.5 and 1)", "off + two on levels", "three on levels\n(0.2, 0.6, 1)"]
vals = [e["defect"] for e in ex]
cols = [G if abs(v) < 1e-12 else R for v in vals]
ax.bar(range(4), vals, color=cols)
for i, (v, e) in enumerate(zip(vals, ex)):
    ax.text(i, v + 0.02, f"gap {v:.3f}\npair k = {e['k']:.2f}", ha="center", fontsize=8)
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("closure gap  (E[h³]E[h] − E[h²]²) / E[h]⁴")
ax.set_ylim(0, max(vals) * 1.35)
ax.set_title("(a) Forget a frozen bath: two events always glue with a flat k,\n"
             "but three events glue the same way only for the cure model", fontsize=9)

ax = axs[1]
hv, wv = np.array([0.2, 0.6, 1.0]), np.ones(3) / 3
M1, M2 = (wv * hv).sum(), (wv * hv ** 2).sum()
kpair = -M2 / M1 ** 2
a = np.linspace(0, 1, 200)
p = a * M1
true3 = np.array([(wv * (1 - (1 - ai * hv) ** 3)).sum() for ai in a])
two = lambda u, v: u + v + kpair * u * v
glued = two(two(p, p), p)
ax.plot(p, true3, color="k", lw=2, label="true chance of at least one of 3 events")
ax.plot(p, glued, color=R, ls="--", label=f"glued with the pair rule, k = {kpair:.2f}")
ax.fill_between(p, glued, true3, color=R, alpha=0.15)
ax.set_xlabel("chance of one event, p")
ax.set_ylabel("chance of at least one")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("(b) Bath with three on levels (0.2, 0.6, 1): the pair rule\n"
             "undercounts three events (shaded gap)", fontsize=9)
fig.suptitle("Fig 1. Frozen hidden state, chances per event: closed only with one 'on' level (Cauchy–Schwarz)", fontsize=10)
save(fig, "fig1-frozen-chances.png")

# ------------------------------------------------------------- Fig 2: where the forgotten-state laws land
fig, ax = plt.subplots(figsize=(9, 6))
dd = np.linspace(0, 4.4, 300)
ax.fill_between(dd, -2 * np.sqrt(dd), 2 * np.sqrt(dd), color="#f1f3f5", zorder=0)
ax.plot(dd, 2 * np.sqrt(dd), "k", lw=1.6); ax.plot(dd, -2 * np.sqrt(dd), "k", lw=1.6)
ax.axhline(0, color=GR, lw=0.5); ax.axvline(0, color=GR, lw=0.5)
# cure model: d = 0, k = -1/pi
for pi_ in (1, 0.5, 1 / 3):
    ax.plot(0, -1 / pi_, "s", color=G, ms=8, zorder=5)
ax.plot([0, 0], [-3.4, -1], color=G, lw=3, alpha=0.5)
ax.annotate("cure model (frozen: off + one on level)\nk = −1/π, d = 0;  π = 1 is the plain 'or' law",
            (0, -3), (-3.9, -3.4), fontsize=8, color=G, arrowprops=dict(arrowstyle="->", color=G))
# geometric frailty line k + d = -1, theta in [0, 1]
th = np.linspace(0, 1, 50)
ax.plot(th, -(1 + th), color=B, lw=3)
ax.annotate("geometric ladder of frozen levels 1, 2, 3, ...\n(Ali–Mikhail–Haq): k = −(1+θ), d = θ",
            (0.5, -1.5), (1.3, -3.6), fontsize=8, color=B, arrowprops=dict(arrowstyle="->", color=B))
ax.plot(1, -2, "*", color=B, ms=15, zorder=6)
ax.annotate("θ = 1: exponential (gamma 1) frailty,\ntouches the parabola", (1, -2), (1.7, -2.2), fontsize=8, color=B,
            arrowprops=dict(arrowstyle="->", color=B))
# map functions: d = -2k - 4
kk = np.linspace(-2, 0, 50)
ax.plot(-2 * kk - 4, kk, color=P, lw=3)
ax.plot(0, -2, "o", color=P, ms=9, zorder=6)
ax.plot(-4, 0, "o", color=P, ms=9, zorder=6)
ax.annotate("Haldane map (no interference) = cure model with π = ½\nk = −2, d = 0", (0, -2), (-3.9, -1.6), fontsize=8, color=P,
            arrowprops=dict(arrowstyle="->", color=P))
ax.annotate("Kosambi map (interference)\nk = 0, d = −4", (-4, 0), (-3.9, 1.0), fontsize=8, color=P,
            arrowprops=dict(arrowstyle="->", color=P))
ax.text(-2.5, -0.55, "Felsenstein's family:\nd = −2k − 4", color=P, fontsize=8, rotation=-26)
ax.text(1.6, 2.6, "inside the parabola: none of\nthese cases lands here", fontsize=8.5, ha="left")
ax.text(2.25, 0.6, "OFF the plane dial (no point):\n- frozen bath, 2+ on levels, chances per event\n"
        "- gamma frailty θ ≠ 1 (Clayton)\n- chi-square interference m = 1\n- Haldane–Wu exclusion statistics",
        fontsize=8, bbox=dict(fc="white", ec=R, lw=0.8), color=R)
ax.set_xlim(-4.3, 4.4); ax.set_ylim(-3.8, 3.4)
ax.set_xlabel("bend d"); ax.set_ylabel("dial k")
ax.set_title("Fig 2. Where the forgotten-hidden-state laws land on the plane dial F = (x + y + kxy)/(1 − dxy)", fontsize=10)
save(fig, "fig2-where-they-land.png")

# ------------------------------------------------------------- Fig 3: interference: drift is not the test
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
C5 = J["C5"]
ax = axs[0]
ax.plot(C5["x_chi"], C5["k_chi"], color=R, lw=2, label="chi-square m = 1 (running two-phase state)")
ax.plot(C5["x_kos"], C5["k_kos"], color=P, lw=2, ls="--", label="Kosambi (on the plane dial, d = −4)")
ax.axhline(-2, color="k", lw=0.8, ls=":"); ax.text(0.02, -1.93, "Haldane map: flat k = −2", fontsize=8)
ax.set_xlabel("recombination fraction r of each half")
ax.set_ylabel("one-dial k read from two halves")
ax.legend(fontsize=8, loc="center left")
ax.set_title("(a) Both drift from 0 to −2: drift in the one-dial k\ncannot tell 'bent' from 'off the family'", fontsize=9)
ax = axs[1]
Dg = np.linspace(0, 3, 600)
for name, r, q, col, ls in [
        ("Haldane: 1 − 2r (straight)", 0.5 * (1 - np.exp(-2 * Dg)), np.exp(-2 * Dg), "k", ":"),
        ("Kosambi: 1 − 4r² (a parabola)", 0.5 * np.tanh(2 * Dg), 1 / np.cosh(2 * Dg) ** 2, P, "--"),
        ("chi-square m = 1", 0.5 * (1 - np.exp(-4 * Dg) * (1 + 2 * Dg)), np.exp(-4 * Dg) * (1 + 4 * Dg), R, "-")]:
    ax.plot(r, q, color=col, ls=ls, lw=2, label=name)
coef = J["C4"]["best_quadratic"]
rr = np.linspace(0, 0.5, 200)
ax.plot(rr, np.polyval(coef, rr), color=R, lw=0.8, alpha=0.6, label=f"best quadratic for chi-square (misses by {J['C4']['max_miss']:.3f})")
ax.set_xlabel("recombination fraction r")
ax.set_ylabel("growth rate r′ at that r")
ax.legend(fontsize=7.5, loc="upper right")
ax.set_title("(b) The real test: on the plane dial exactly when r′ is a\nquadratic in r (a Riccati equation). Chi-square is not.", fontsize=9)
fig.suptitle("Fig 3. Running hidden state (crossover interference): always a law, bent onto the dial only for Kosambi–Felsenstein maps", fontsize=10)
save(fig, "fig3-interference.png")

# ------------------------------------------------------------- Fig 4: Haldane exclusion statistics
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
D = J["D"]
gs = list(D.keys())
N = 12
for i, g in enumerate(gs):
    def coeffs(gv, n):
        return [1.0] + [float(mp.binomial(m * (1 - gv), m - 1) / m) for m in range(1, n + 1)]
    cs = coeffs(mp.mpf(g), N)
    for n in range(1, N + 1):
        c = cs[n]
        col = G if c > 1e-14 else (R if c < -1e-14 else "#dee2e6")
        ax.add_patch(plt.Rectangle((n - 0.45, i - 0.4), 0.9, 0.8, color=col))
ax.set_xlim(0.4, N + 0.6); ax.set_ylim(-0.6, len(gs) - 0.4)
ax.set_yticks(range(len(gs)))
ax.set_yticklabels([f"g = {float(g):.2f}" + (" (Bose)" if float(g) == 0 else " (Fermi)" if float(g) == 1 else "") for g in gs], fontsize=8)
ax.set_xticks(range(1, N + 1))
ax.set_xlabel("number of copies n in one state")
ax.set_title("(a) Counts in Wu's one-state weight: green positive, red negative,\n"
             "grey zero. Every g strictly between 0 and 1 hits a negative count", fontsize=9)
ax = axs[1]
w = np.logspace(-4, 4, 400)
for g, col in ((0.25, "#74c0fc"), (0.5, B), (0.75, "#1864ab")):
    nocc = 1 / (w + g)
    kh = w + g - w ** g * (1 + w) ** (1 - g)
    ax.semilogx(nocc, kh, color=col, label=f"g = {g}: k from {2*g-1:+.2f} (sparse) to {g:.2f} (crowded)")
ax.set_xlabel("occupation n (crowding)")
ax.set_ylabel("one-dial k = 1/n − 1/x")
ax.legend(fontsize=7.5)
ax.set_title("(b) Its k drifts with crowding (from repetition-dial/)", fontsize=9)
fig.suptitle("Fig 4. Haldane–Wu crowding is not a frozen mixture of honest states: any such mixture has only positive counts", fontsize=10)
save(fig, "fig4-haldane.png")
