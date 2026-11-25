import json, math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 150})
INK, BLUE, RED, GREEN, GRAY, AMBER = "#1f2937", "#2563eb", "#dc2626", "#16a34a", "#9ca3af", "#d97706"
R = json.load(open("s3_results.json"))
fig, axs = plt.subplots(1, 2, figsize=(11, 3.9))
x = np.linspace(-1.999, 1.999, 400)
xn = np.array(R["noncm_x"]); xc = np.array(R["cm_x"])
ax = axs[0]
ax.hist(xn, bins=24, range=(-2, 2), density=True, color=GREEN, alpha=.55, label="4 ordinary curves")
ax.plot(x, np.sqrt(4 - x**2)/(2*np.pi), color=INK, lw=2, label="3-sphere: semicircle")
ax.set_title("Ordinary curves: shadow of the 3-sphere", fontsize=11); ax.set_xlabel("a_p / √p  (= 2 cos θ)"); ax.legend(fontsize=9, frameon=False); ax.set_yticks([])
ax = axs[1]
nz = xc[np.abs(xc) > 1e-9]
cnt, ed = np.histogram(nz, bins=24, range=(-2, 2)); ax.bar((ed[1:]+ed[:-1])/2, cnt/(len(xc)*(ed[1]-ed[0])), width=ed[1]-ed[0], color=AMBER, alpha=.6, label="2 CM curves (a_p ≠ 0)")
ax.plot(x, 0.5/(np.pi*np.sqrt(4 - x**2)), color=INK, lw=2, label="circle: arcsine (half the primes)")
ax.axvline(0, color=RED, lw=2.5, label=f"a_p = 0: {np.mean(np.abs(xc) < 1e-9):.0%} of primes")
ax.set_ylim(0, 0.9); ax.set_title("CM curves: shadow of a circle", fontsize=11); ax.set_xlabel("a_p / √p"); ax.legend(fontsize=9, frameon=False); ax.set_yticks([])
plt.tight_layout(); plt.savefig("figures/f7_shadows.png"); plt.close()
# cube layers for 2n = 6
fig, ax = plt.subplots(figsize=(7, 3.4))
L = [1, 6, 15, 20, 15, 6, 1]; freq = [6, 4, 2, 0, -2, -4, -6]
cols = [GRAY]*7; cols[3] = BLUE; cols[2] = cols[4] = RED
ax.bar(range(7), L, color=cols)
for i, v in enumerate(L): ax.text(i, v + .4, str(v), ha="center", fontsize=10)
ax.set_xticks(range(7)); ax.set_xticklabels([f"{j} ones\nturns {f:+d}θ" for j, f in zip(range(7), freq)], fontsize=8.5)
ax.set_yticks([]); ax.set_title("Corners of the 6-cube by layer. Circle keeps 20; 3-sphere keeps 20 − 15 = 5 (Catalan)", fontsize=10)
plt.tight_layout(); plt.savefig("figures/f8_cube_layers.png"); plt.close()
print("ok")
