"""DRAFT time-study figures. Colour key: green = derived/known GR, blue = numeric check, orange = proposed test, magenta = our reading."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G, B, O, M = "#2e7d32", "#1f5fbf", "#e07b00", "#c2185b"; out = __file__.replace("figures.py", "figures/")

# ts1: one shared axis, winding rate per place, for a still star (works) and for a passing wave (rate must swing in time)
fig = plt.figure(figsize=(11, 4.6))
for n, (title, rates) in enumerate([("A. Still star: rate set by place only (works)", lambda T: [2.2, 1.6, 1.0]),
                                    ("B. Passing wave: rate swings with time", lambda T: [1.6 + 0.6*np.sin(2*np.pi*T/4), 1.6, 1.6 - 0.6*np.sin(2*np.pi*T/4)])]):
    ax = fig.add_subplot(1, 2, n + 1, projection="3d"); T = np.linspace(0, 12, 1500)
    for j, xpos in enumerate([-2, 0, 2]):
        r = np.array([rates(tt)[j] for tt in T]); ph = np.cumsum(r)*(T[1] - T[0])*2*np.pi/3
        ax.plot(xpos + 0.35*np.cos(ph), 0.35*np.sin(ph), T, color=[M, B, G][j], lw=1.3)
        ax.plot([xpos, xpos], [0, 0], [0, 12], color="grey", ls="--", lw=0.8)
    ax.plot([-2.6, 2.6], [0, 0], [0, 0], color="k", lw=0.8)
    ax.set_title(title, fontsize=10); ax.set_xlabel("place"); ax.set_zlabel("shared time axis"); ax.set_yticks([])
fig.suptitle("ts1. Each place's clock winds around one shared time axis. Grey dashed: the axis at each place. "
             "In B the winding speeds up and slows down as the wave goes by.", fontsize=9)
fig.savefig(out + "ts1_shared_axis.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ts2: the winding-rate field of a plus wave over half a cycle: a four-leaf pattern that flips (spin 2)
g = np.linspace(-1, 1, 201); X, Y = np.meshgrid(g, g)
fig, axs = plt.subplots(1, 5, figsize=(14, 3.2))
for ax, ph in zip(axs[:4], [0, 0.25, 0.5, 0.75]):
    c = np.cos(2*np.pi*ph); Phi = 0.25*c*(X**2 - Y**2)          # Phi = -1/4 hdd (x^2-y^2), hdd = -w^2 A cos
    im = ax.imshow(Phi, extent=[-1, 1, -1, 1], cmap="PiYG", vmin=-0.25, vmax=0.25)
    ax.set_title(f"phase {ph:.2f} of a wave cycle", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
ax = axs[4]; Phi = 0.5*X*Y; ax.imshow(Phi, extent=[-1, 1, -1, 1], cmap="PiYG", vmin=-0.25, vmax=0.25)
ax.set_title("cross wave: same leaves turned 45°", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=axs, shrink=0.8, label="pink: clock slower than centre   green: faster")
fig.suptitle("ts2. Winding rate per place in a detector as a plus wave passes along the viewing line (check T3, T6). "
             "Four leaves; they swap every half cycle. A 90° turn flips them, a 180° turn returns them: spin 2.", fontsize=9)
fig.savefig(out + "ts2_four_leaf.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ts3: how far 'shared axis + winding field' alone goes: sinc(wL), with LIGO and LISA marked
wl = np.logspace(-3, 1, 400)
fig, ax = plt.subplots(figsize=(7.5, 4))
ax.semilogx(wl, np.sin(wl)/wl, color=G, lw=2, label="exact arm response ÷ winding-only response = sin(ωL)/(ωL)")
for lab, v, col in [("LIGO 100 Hz", 2*np.pi*100*4e3/3e8, B), ("LIGO 1 kHz", 2*np.pi*1e3*4e3/3e8, B),
                    ("LISA 10 mHz", 2*np.pi*0.01*2.5e9/3e8, O), ("LISA 0.1 Hz", 2*np.pi*0.1*2.5e9/3e8, O)]:
    ax.axvline(v, color=col, ls=":"); ax.text(v, 0.05, lab, rotation=90, fontsize=8, color=col)
ax.axhline(1, color="grey", lw=0.6); ax.set_xlabel("ωL: arm length in wavelengths × 2π"); ax.set_ylabel("ratio")
ax.set_title("ts3. Winding rate alone is enough while the detector is small next to the wave (check T5).\n"
             "Past ωL ≈ 1 the axis must also tilt (the shift) and space must stretch.", fontsize=9)
ax.legend(fontsize=8, loc="lower left"); fig.savefig(out + "ts3_how_far.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ts4: the angles on time as a map
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis("off")
hub = (0.5, 0.5); ax.add_patch(plt.Circle(hub, 0.09, color=M, alpha=0.85)); ax.text(*hub, "time as a\ndirection\n(helix)", ha="center", va="center", color="w", fontsize=9)
nodes = [("one shared axis\nvs local clocks", G, "tested here: axis survives,\nrate must swing in time"),
         ("signature: the i\n(d = -1 boosts, +1 turns)", G, "test: compose across\na d sign change"),
         ("arrow and entropy\n(loop phase = 2 x entropy)", M, "test: does k > 0 make\ncombined entropy grow?"),
         ("energy not conserved\n(no time symmetry)", O, "test: wave's energy loss\nfrom a drifting rate"),
         ("probability ≠ 1\n(non-unitary)", O, "test: helix that spirals\nin or out, PT pairs"),
         ("more machinery:\n4 numbers per point", G, "lapse (1) + shift (3);\ntest: two winding rates")]
for n, (lab, col, test) in enumerate(nodes):
    a = np.pi/2 + 2*np.pi*n/len(nodes); p = (0.5 + 0.33*np.cos(a), 0.5 + 0.36*np.sin(a))
    ax.plot([hub[0], p[0]], [hub[1], p[1]], color="grey", lw=1)
    ax.text(*p, lab, ha="center", va="center", fontsize=8.5, bbox=dict(fc="w", ec=col, lw=2))
    ax.text(p[0], p[1] - 0.105, test, ha="center", va="top", fontsize=7.5, color=O if n else B)
ax.set_title("ts4. Angles on time. Box colour: green = known physics, magenta = our reading, orange = proposed relaxation. "
             "Text under each box: one test from the law.", fontsize=9)
fig.savefig(out + "ts4_angles_map.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("4 figures")
