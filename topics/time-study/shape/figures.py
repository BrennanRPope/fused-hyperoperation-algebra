"""DRAFT shape figures. Colour key: green = known, blue = check, magenta = his chain inputs/reading, orange = the missing input."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G, B, O, M = "#2e7d32", "#1f5fbf", "#e07b00", "#c2185b"; out = __file__.replace("figures.py", "figures/")
res = json.load(open(__file__.replace("figures.py", "results.json")))
f = res["funnel"]
fig, ax = plt.subplots(figsize=(8, 4))
labs = ["all curvature", "moving at c\n(chain's speed limit)", "balanced for\none observer", "balanced for\nevery observer"]
vals = [f["all curvature"], f["moving at c"], f["balanced for one observer"], f["balanced for every observer"]]
bars = ax.bar(labs, vals, color=["grey", M, O, O])
for b_, v in zip(bars, vals): ax.text(b_.get_x() + b_.get_width()/2, v + 0.3, str(v), ha="center")
ax.set_ylabel("independent curvature pieces")
ax.set_title("sh1. From all curvature to the four-leaf (D1, D4, D5). The chain's speed limit does most of the cut;\n"
             "the last step needs 'balanced for every observer', which is Einstein's vacuum equation in plain words.", fontsize=9)
fig.savefig(out + "sh1_funnel.png", dpi=130, bbox_inches="tight"); plt.close(fig)

th = np.linspace(0, 2*np.pi, 400)
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), subplot_kw=dict(projection="polar"))
for ax, m, col, title in [(axs[0], 0, O, "spin 0: same all round\n(rides along: 2 pieces)"),
                           (axs[1], 1, O, "spin 1: two lobes\n(rides along: 2 pieces)"),
                           (axs[2], 2, M, "spin 2: the four-leaf\n(carries nothing: 2 pieces)")]:
    ax.plot(th, np.abs(np.cos(m*th)) if m else np.ones_like(th), color=col, lw=2.5); ax.set_yticks([]); ax.set_xticks([])
    ax.set_title(title, fontsize=9)
fig.suptitle("sh2. The six shapes a pattern moving at c can have, seen face-on (D2, D6). Orange: shapes that come with something riding\n"
             "along (in GR, energy moving at c, like light). Magenta: the four-leaf, the only shape that carries nothing.", fontsize=9, y=1.08)
fig.savefig(out + "sh2_six_shapes.png", dpi=130, bbox_inches="tight"); plt.close(fig)

g1 = np.arange(-6, 7); Xg, Yg = np.meshgrid(g1, g1, indexing="ij")
fig, axs = plt.subplots(1, 2, figsize=(9, 4))
for ax, F, title in [(axs[0], Xg**2 - Yg**2, "four-leaf x² − y²: averaging changes nothing"),
                     (axs[1], Xg**2 + Yg**2, "trace x² + y²: averaging moves every place")]:
    Fz = F.astype(float); A = Fz.copy()
    A[1:-1, 1:-1] = (Fz[1:-1, 1:-1] + Fz[2:, 1:-1] + Fz[:-2, 1:-1] + Fz[1:-1, 2:] + Fz[1:-1, :-2] + 2*Fz[1:-1, 1:-1])/7
    im = ax.imshow((A - Fz)[1:-1, 1:-1], cmap="coolwarm", vmin=-0.7, vmax=0.7); ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=axs, shrink=0.8, label="change after one neighbour-averaging step")
fig.suptitle("sh3. 'Balanced' in chain words: each place equals the average of its neighbours, the rest state of his neighbour chain (D7).\n"
             "One observer's averaging only removes the trace; every observer's is needed for the four-leaf alone.", fontsize=9, y=1.06)
fig.savefig(out + "sh3_balanced.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("3 figures")
