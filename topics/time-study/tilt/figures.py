"""DRAFT tilt/chain/arrow figures. Colour key: green = known GR or derived, blue = numeric check, orange = open, magenta = our reading."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G, B, O, M = "#2e7d32", "#1f5fbf", "#e07b00", "#c2185b"; out = __file__.replace("figures.py", "figures/")
res = json.load(open(__file__.replace("figures.py", "results.json")))

# tl1: rate pattern (what the law's one number per place gives) vs the tilt the wave also needs
g = np.linspace(-1, 1, 201); X, Y = np.meshgrid(g, g); q = np.linspace(-1, 1, 11); QX, QZ = np.meshgrid(q, q)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
axs[0].imshow(0.25*(X**2 - Y**2), extent=[-1, 1, -1, 1], cmap="PiYG", origin="lower")
axs[0].set_title("A. Winding rate (one number per place)\nface-on, plus wave: the four leaves", fontsize=9)
gx, gz = QX*QZ/3, 0*QX                      # plus wave, y = 0 slice: g_0x = x z/3, g_0z = -x^2/3
gz = -QX**2/3
axs[1].quiver(QX, QZ, gx, gz, color=G); axs[1].set_title("B. Tilt of the time axis (3 numbers per place)\nside-on slice y = 0, wave moving up (z)", fontsize=9)
axs[1].set_xlabel("x"); axs[1].set_ylabel("z (wave direction)")
th = np.linspace(0, 2*np.pi, 200)
for ph, col, lab in [(0, M, "winding pattern E"), (np.pi/4, G, "twist pattern B = i·E")]:
    r = np.abs(np.cos(2*(th - ph))); axs[2].plot(r*np.cos(th), r*np.sin(th), color=col, lw=2, label=lab)
axs[2].set_aspect("equal"); axs[2].legend(fontsize=8, loc="lower right"); axs[2].set_xticks([]); axs[2].set_yticks([])
axs[2].set_title("C. For a wave the twist is the winding\npattern turned 45°, i.e. times i (check A8)", fontsize=9)
fig.suptitle("tl1. A rate alone can never tilt the axis (A1). The wave needs a tilt (A2, A3), and its twist is fixed by the rate pattern and the wave's direction.", fontsize=9, y=1.04)
fig.savefig(out + "tl1_rate_vs_tilt.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# tl2: the chain picture: scalar steps close a loop, direction steps leave a twist
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
for i in range(5):
    ax.add_patch(plt.Rectangle((0.3 + 1.9*i, 3.6), 1.3, 1.0, fc="w", ec=M, lw=2)); ax.text(0.95 + 1.9*i, 4.1, f"f at\npoint {i+1}", ha="center", va="center", fontsize=8)
    if i < 4: ax.annotate("", (2.2 + 1.9*i, 4.1), (1.6 + 1.9*i, 4.1), arrowprops=dict(arrowstyle="->", color="k"))
ax.text(5, 5.2, "Brennan's hunch: a pure function at each point feeds the next", ha="center", fontsize=9, color=M)
ax.text(5, 2.6, "one number per step: steps commute, loops close, only a rate (A4)\n"
                "a direction per step: order matters, a loop leaves a turn = the tilt's twist (A5, A6)\n"
                "rung 1 groups freely, so the chain has one answer (A7)", ha="center", va="top", fontsize=8.5)
v = np.linspace(0.01, 0.2, 5); axs[1].plot(v, v**2, color=G, lw=2, label="leading order v²")
axs[1].plot(v, np.array(res["loop_over_v2"])*v**2, "o", color=B, label="numeric: tilt x, tilt y, undo x, undo y")
axs[1].set_xlabel("tilt size v"); axs[1].set_ylabel("leftover turn (radians)"); axs[1].legend(fontsize=8)
axs[1].set_title("Leftover turn after a square loop of tilts (Thomas–Wigner)", fontsize=9)
fig.suptitle("tl2. Time as a composition chain. Magenta: the hunch (ours). Green: known (Thomas 1926, Wigner 1939; time-ordered products, Dyson 1949).", fontsize=9)
fig.savefig(out + "tl2_chain.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# tl3: where the chain picks a direction
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
k = 0.7; law = lambda p, q: p + q + k*p*q; steps = np.random.default_rng(3).normal(size=8)*0.4; ab = np.insert(steps, 4, -1/k)
st = np.linspace(-1, 1, 1000)
for seq, col, lab in [(steps, B, "ordinary steps"), (ab, M, "step 5 = absorbing point −1/k")]:
    xs = st.copy(); cnt = [1000]
    for s_ in seq: xs = law(s_, xs); cnt.append(len(np.unique(np.round(xs, 9))))
    axs[0].step(range(len(cnt)), cnt, where="post", color=col, lw=2, label=lab)
axs[0].set_yscale("log"); axs[0].set_xlabel("step in the chain"); axs[0].set_ylabel("distinct states left"); axs[0].legend(fontsize=8)
axs[0].set_title("A. One-way information loss at −1/k (A9)", fontsize=9)
r3 = res["rung3_bracketings"]; axs[1].scatter([1]*14, [0]*14, color=B, s=30, label="rung 1: all 14 groupings agree")
axs[1].scatter([2]*len(r3), np.array(r3) - np.mean(r3), color=O, s=30, label="rung 3: groupings spread")
axs[1].set_xticks([1, 2]); axs[1].set_xticklabels(["rung 1", "rung 3"]); axs[1].set_ylabel("answer − average")
axs[1].legend(fontsize=8); axs[1].set_title("B. Order matters on rung 3, but no information is lost (A7, A10)", fontsize=9)
fig.suptitle("tl3. The arrow. Only the absorbing point makes the chain one-way: you can run it forward, not back. Rung 3 says direction matters, not which way.", fontsize=9)
fig.savefig(out + "tl3_arrow.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("3 figures")

# tl4: speed limit of the chain
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
H = np.log(np.array(res["cone_hist"])); H[H <= 0] = np.nan
axs[0].imshow(np.log10(H), aspect="auto", origin="lower", cmap="viridis", extent=[-40.5, 40.5, -0.5, 30.5])
axs[0].plot([0, 30], [0, 30], "w--", lw=1); axs[0].plot([0, -30], [0, 30], "w--", lw=1)
axs[0].set_xlabel("place"); axs[0].set_ylabel("step"); axs[0].set_title("A. A nudge spreads one place per step: a light cone (S1)", fontsize=9)
vel = res["vel_chain"]; axs[1].plot(range(len(vel)), vel, "o-", color=B, ms=3, label="steps of 0.5c composed")
axs[1].axhline(1, color=M, lw=2, label="c = the law's absorbing point −1/k (S2)")
axs[1].set_xlabel("number of steps"); axs[1].set_ylabel("speed (c = 1)"); axs[1].legend(fontsize=8, loc="lower right")
axs[1].set_title("B. Composed speeds crowd up to c and never reach it (S3)", fontsize=9)
fig.suptitle("tl4. The speed limit. White dashed: the cone edge. Approaching c is the speed limit; landing on it is the arrow's one-way step.", fontsize=9)
fig.savefig(out + "tl4_speed_limit.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("tl4")
