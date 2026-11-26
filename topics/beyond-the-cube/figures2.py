"""Figures for the second pass (section 3b of README). Run after checks2.py (reads sim2.npz)."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = np.load("sim2.npz")
OUT = "figures/"
BLUE, ORANGE, GREEN, GRAY, RED, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#a0aec0", "#c53030", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
alphas = S["alphas"]; n = len(alphas)
SIDECOL = ["#2b6cb0", "#dd6b20", "#2b6cb0", "#dd6b20", "#dd6b20", "#dd6b20"]


def side_of(a):
    return (np.searchsorted(alphas, a, side="right") - 1) % n


def geo(ax, t, p, **kw):
    d = (p - t) % (2 * np.pi)
    if d > np.pi:
        t, p = p, t; d = 2 * np.pi - d
    mid = t + d / 2
    if d > np.pi - 1e-3:
        ax.plot([math.cos(t), math.cos(p)], [math.sin(t), math.sin(p)], **kw); return
    c = np.array([math.cos(mid), math.sin(mid)]) / math.cos(d / 2)
    r = math.tan(d / 2)
    a0 = math.atan2(math.sin(t) - c[1], math.cos(t) - c[0]); a1 = math.atan2(math.sin(p) - c[1], math.cos(p) - c[0])
    da = (a1 - a0 + np.pi) % (2 * np.pi) - np.pi
    s = np.linspace(a0, a0 + da, 60)
    ax.plot(c[0] + r * np.cos(s), c[1] + r * np.sin(s), **kw)


# ---------------------------------------------------------------- fig9 one universe
rng = np.random.default_rng(7)
def one_universe(eps=0.02):
    ce = 1 / math.tan(eps / 2); N = rng.poisson(math.pi * ce)
    w = rng.uniform(0, 1, N); d = 2 * np.arctan2(1.0, ce - 2 * ce * w)
    t = rng.uniform(0, 2 * np.pi, N); return t, (t + d) % (2 * np.pi)
while True:
    TH, PH = one_universe()
    s1, s2 = side_of(TH), side_of(PH)
    blue = np.any(((s1 == 0) & (s2 == 2)) | ((s1 == 2) & (s2 == 0)))
    orng = np.any(((s1 == 1) & np.isin(s2, [3, 4, 5])) | ((s2 == 1) & np.isin(s1, [3, 4, 5])))
    if blue and orng:
        break
fig, axs = plt.subplots(1, 2, figsize=(13, 6.6))
for ax, title, mode in [(axs[0], "One random universe of geodesics (~314 of them)", "all"),
                        (axs[1], "Chord 02 in this universe: side 0 to side 2 (blue channel)\nand side 1 to sides 3, 4, 5 (orange channel): here both are linked", "chord")]:
    tt = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(tt), np.sin(tt), color="k", lw=0.8)
    for k in range(n):
        a0 = alphas[k]; a1 = alphas[(k + 1) % n] + (2 * np.pi if k == n - 1 else 0)
        s = np.linspace(a0, a1, 80)
        col = SIDECOL[k] if mode == "chord" else "#4a5568"
        ax.plot(1.04 * np.cos(s), 1.04 * np.sin(s), color=col, lw=6, solid_capstyle="butt")
        m = (a0 + a1) / 2
        ax.text(1.17 * math.cos(m), 1.17 * math.sin(m), f"side {k}", ha="center", va="center", fontsize=10)
    for t, p in zip(TH, PH):
        s1, s2 = side_of(t), side_of(p)
        if mode == "all":
            geo(ax, t, p, color=GRAY, lw=0.5, alpha=0.6)
        else:
            pair = {s1, s2}
            if pair == {0, 2}:
                geo(ax, t, p, color=BLUE, lw=2.2)
            elif 1 in pair and (pair - {1}) <= {3, 4, 5} and len(pair) == 2:
                geo(ax, t, p, color=ORANGE, lw=2.2)
            else:
                geo(ax, t, p, color="#e2e8f0", lw=0.4)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=11)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35)
fig.suptitle("Fig 9. The random-geodesic soup (Crofton intensity, rate 1): a chord's two crossing channels are 'linked' "
             "if some geodesic joins them", weight="bold", fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUT + "fig9-geodesic-soup.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig10 simulated vs exact u
fig, ax = plt.subplots(figsize=(6.2, 6))
ax.plot([0, 1], [0, 1], color=GRAY, ls="--")
ax.scatter(S["eu"], S["su"], color=GREEN, s=50, zorder=3)
for (i, j), e, s in zip(S["chords"], S["eu"], S["su"]):
    ax.annotate(f"u{i}{j}", (e, s), (5, -12), textcoords="offset points", fontsize=9)
ax.set_xlabel("exact u (cross-ratio of the two sides)"); ax.set_ylabel("simulated P(no geodesic links the two sides)")
ax.set_title(f"Fig 10. Every coordinate of the hexagon's associahedron\nis a 'no link' probability ({int(S['U'])} universes)",
             fontsize=11, weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig10-u-is-a-probability.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig11 both vs neither
fig, ax = plt.subplots(figsize=(11, 4.8))
xs = np.arange(len(S["chords"])); w = 0.27
ax.bar(xs - w, S["both"], w, color=GREEN, label="both channels linked")
ax.bar(xs, S["neither"], w, color=BLUE, label="neither linked")
ax.bar(xs + w, S["one"], w, color=GRAY, label="exactly one linked")
ax.set_xticks(xs, [f"chord {i}{j}" for i, j in S["chords"]], fontsize=9.5)
ax.set_ylabel("simulated probability")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.13))
ax.set_title("Fig 11. For every chord, green and blue are equal: 'both linked' is exactly as likely as 'neither'.\n"
             "That equality is the chord's u-equation, read as probability.", fontsize=11, weight="bold", pad=34)
fig.tight_layout(); fig.savefig(OUT + "fig11-both-equals-neither.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig12 rate
R = S["rate"]; u = float(S["eu"][0])
lam = np.linspace(0.3, 2.5, 200)
fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.plot(lam, u**lam + (1 - u)**lam, color=PURPLE, lw=2, label="exact: u^λ + (1 − u)^λ")
ax.scatter(R[:, 0], R[:, 1], color=ORANGE, zorder=3, s=45, label="simulated, 12,000 universes each")
ax.axhline(1, color=GRAY, ls="--"); ax.axvline(1, color=GRAY, ls=":")
ax.annotate("the equation holds\nonly at rate 1", (1, 1), (1.35, 1.2), arrowprops=dict(arrowstyle="->"), fontsize=10)
ax.set_xlabel("clock rate λ (geodesics per unit of Crofton measure)")
ax.set_ylabel("P(no link i–j) + P(no link between blocks)")
ax.legend(frameon=False, fontsize=9)
ax.set_title("Fig 12. The associahedron fixes the clock's rate", fontsize=11, weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig12-rate-is-pinned.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig13 certainty
C = S["cert"]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.9))
ax = axs[0]
ax.plot(C[:, 0], C[:, 1], "o-", color=BLUE, label="side 0 linked to side 2 (exact)")
ax.plot(C[:, 0], C[:, 2], "o-", color=ORANGE, label="side 1 linked to 3, 4 or 5 (exact)")
ax.scatter(C[:, 0], C[:, 3], marker="x", color="k", zorder=3, label="simulated")
ax.scatter(C[:, 0], C[:, 4], marker="x", color="k", zorder=3)
ax.set_xscale("log"); ax.invert_xaxis(); ax.minorticks_off(); ax.set_xticks(C[:, 0], [f"{v:g}" for v in C[:, 0]])
ax.set_xlabel("length of side 1 (shrinking to the right)"); ax.set_ylabel("probability")
ax.legend(frameon=False, fontsize=9)
ax.set_title("As side 1 shrinks, sides 0 and 2 touch:\none channel becomes certain, the crossing one impossible", fontsize=10)
ax = axs[1]
ang = np.pi / 2 - 2 * np.pi * np.arange(6) / 6
P = np.c_[np.cos(ang), np.sin(ang)]
for k, (sub, chs) in enumerate([((-1.3, 0), [(0, 2), (0, 3), (0, 4)]), ((1.3, 0), [(1, 3), (3, 5), (5, 1)])]):
    Q = P * 0.9 + np.array(sub)
    ax.fill(*np.vstack([Q, Q[:1]]).T, color="#f7fafc", ec="k", lw=1)
    for a, b in chs:
        ax.plot(*zip(Q[a], Q[b]), color=GREEN, lw=2.5)
ax.plot(*zip(P[0] * 0.9 + (1.3, 0), P[3] * 0.9 + (1.3, 0)), color=RED, lw=2, ls="--")
ax.text(1.3, -1.25, "red crosses green: cannot also be certain", ha="center", fontsize=9.5, color=RED)
ax.text(-1.3, -1.25, "a triangulation: three channels\nthat can all be certain at once", ha="center", fontsize=9.5, color=GREEN)
ax.set_aspect("equal"); ax.axis("off"); ax.set_ylim(-1.6, 1.2)
ax.set_title("So the associahedron's corners (non-crossing chord sets)\nare the largest sets of channels that can be certain together", fontsize=10)
fig.suptitle("Fig 13. Facets and corners of the associahedron, read as certainty", weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig13-facets-as-certainty.png", dpi=140); plt.close(fig)
print("figures2 written")
