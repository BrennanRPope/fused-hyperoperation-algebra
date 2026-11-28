"""DRAFT twist figures. Colour key: green = known identity, blue = check, magenta = his reading, orange = not used."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G, B, O, M = "#2e7d32", "#1f5fbf", "#e07b00", "#c2185b"; out = __file__.replace("figures.py", "figures/")
res = json.load(open(__file__.replace("figures.py", "results.json")))

# tw1: what goes in, what comes out
fig, ax = plt.subplots(figsize=(11, 4.4)); ax.axis("off"); ax.set_xlim(0, 11); ax.set_ylim(0, 5)
box = lambda xy, txt, col, w=2.6: (ax.add_patch(plt.Rectangle(xy, w, 1.0, fc="w", ec=col, lw=2)),
                                   ax.text(xy[0] + w/2, xy[1] + 0.5, txt, ha="center", va="center", fontsize=8.5))
box((0.2, 3.6), "rate pattern\n(the four-leaf, E)", M); box((0.2, 2.2), "wave direction\n(along z)", M)
box((0.2, 0.8), "speed limit c\n(chain's cone, ../tilt S1-S3)", M)
box((4.1, 2.2), "steps compose along paths\n(associative chain)\n= Bianchi: dB/dt = curl E (B1)", G, 3.0)
box((8.2, 2.2), "twist B = i·E\nexactly (B4)", B)
for yy in (4.1, 2.7, 1.3): ax.annotate("", (4.1, 2.7), (2.8, yy), arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("", (8.2, 2.7), (7.1, 2.7), arrowprops=dict(arrowstyle="->", color="k", lw=2))
box((4.1, 0.4), "Einstein's field equations\n(not used; cross-check only, B6)", O, 3.0)
ax.set_title("tw1. The tilt's twist from the chain alone. Magenta: the inputs (his picture). Green: a known identity. "
             "Blue: what comes out. Orange: not used.", fontsize=9)
fig.savefig(out + "tw1_inputs_outputs.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# tw2: the given pattern and the forced twist
g = np.linspace(-1, 1, 201); Xg, Yg = np.meshgrid(g, g)
fig, axs = plt.subplots(1, 2, figsize=(9, 4.2))
axs[0].imshow(Xg**2 - Yg**2, extent=[-1, 1, -1, 1], cmap="PiYG", origin="lower"); axs[0].set_title("given: rate pattern E (plus wave)", fontsize=9)
axs[1].imshow(-2*Xg*Yg, extent=[-1, 1, -1, 1], cmap="PiYG", origin="lower"); axs[1].set_title("forced: twist B, the same leaves turned 45°", fontsize=9)
for a in axs: a.set_xticks([]); a.set_yticks([])
fig.suptitle("tw2. Face-on view, wave coming at you. Nothing about B was put in: the chain's consistency and the speed c fix it (B3, B4).", fontsize=9)
fig.savefig(out + "tw2_given_forced.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# tw3: why the speed has to be c
sv = np.array(res["speeds"])
fig, ax = plt.subplots(figsize=(7.5, 4))
ax.plot(sv, 1/sv, color=B, lw=2, label="forced twist size ÷ rate size = 1/s (B5)")
ax.plot(sv, np.abs(sv**2 - 1)/sv**2, color=O, lw=2, ls="--", label="miss in Einstein's partner relation, |s²−1|/s² (B6)")
ax.axvline(1, color=M, lw=1.5); ax.text(1.03, 2.6, "s = c", color=M); ax.axhline(1, color="grey", lw=0.6)
ax.set_xlabel("speed s of the moving pattern (c = 1)"); ax.set_ylim(0, 3.4); ax.legend(fontsize=8)
ax.set_title("tw3. Only at the chain's speed limit is the twist exactly i·E,\nand only there does it also agree with Einstein's equations.", fontsize=9)
fig.savefig(out + "tw3_speed.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("3 figures")
