"""Round 9 figures (README section 13): fig27 the concrete break, fig28 the counts, fig29 the drop over the
collision surface, fig30 what is proved about the size at the break."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"; PURPLE = "#7b4fa0"
B = json.load(open("breaks.json"))
lam, mu, nu, kap = 19/200, 1/20, 1/19, 1/10

# Fig 27: two different routes, one landing place
fig, ax = plt.subplots(figsize=(11, 5.6))
hull = (0.0, 2.0)
def img(c, l, iv): return tuple(sorted((c + l * (iv[0] - c), c + l * (iv[1] - c))))
rows = {"start": 3.0, "mid": 1.8, "end": 0.6}
ax.plot(hull, [rows["start"]] * 2, color="k", lw=6, solid_capstyle="butt")
ax.text(1.0, rows["start"] + 0.25, "the whole fractal sits in [0, 2]", ha="center")
a1 = img(1, nu, hull); a2 = img(0, lam, a1)
b1 = img(2, kap, hull); b2 = img(0, mu, b1)
ax.plot(a1, [rows["mid"]] * 2, color=BLUE, lw=8, solid_capstyle="butt")
ax.plot(b1, [rows["mid"]] * 2, color=ORANGE, lw=8, solid_capstyle="butt")
ax.text(np.mean(a1), rows["mid"] - 0.32, "shrink by 1/19 toward 1", color=BLUE, ha="center")
ax.text(np.mean(b1), rows["mid"] - 0.32, "shrink by 1/10 toward 2", color=ORANGE, ha="center")
ax.plot(a2, [rows["end"] + 0.06] * 2, color=BLUE, lw=8, solid_capstyle="butt")
ax.plot(b2, [rows["end"] - 0.06] * 2, color=ORANGE, lw=8, solid_capstyle="butt")
ax.annotate("", xy=(np.mean(a1) - 0.15, rows["mid"] + 0.1), xytext=(0.7, rows["start"] - 0.08), arrowprops=dict(arrowstyle="->", color=BLUE))
ax.annotate("", xy=(np.mean(b1), rows["mid"] + 0.1), xytext=(1.5, rows["start"] - 0.08), arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.annotate("", xy=(np.mean(a2), rows["end"] + 0.12), xytext=(np.mean(a1) - 0.1, rows["mid"] - 0.4), arrowprops=dict(arrowstyle="->", color=BLUE))
ax.annotate("", xy=(np.mean(b2), rows["end"] - 0.0), xytext=(np.mean(b1) - 0.2, rows["mid"] - 0.4), arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.text(0.16, rows["end"] - 0.45, "then by 19/200 toward 0 (blue) or by 1/20 toward 0 (orange):\n"
        "both land on exactly [0.090, 0.100], the same map z ↦ (18 + z)/200", fontsize=10)
for c in (0, 1, 2): ax.plot([c], [rows["start"] - 0.25], "k^"); ax.text(c, rows["start"] - 0.5, f"centre {c}", ha="center", fontsize=9)
ax.set_xlim(-0.15, 2.25); ax.set_ylim(-0.2, 3.6); ax.axis("off")
ax.set_title("Figure 27. A break: two different words give the same map (exact fractions; Lean example_collision)", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig27-the-break.png", dpi=150); plt.close(fig)

# Fig 28: counts
c = B["B2_counts"]; L = c["lengths"]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]
ax.semilogy(L, c["trace"], "o-", color=GREY, label="different words (what the formula counts)")
ax.semilogy(L, c["distinct"], "s-", color=RED, label="different maps (hashed mod two primes)")
ax.semilogy(L, c["family_model"], "x", color="k", ms=9, label="prediction: only the one break and its family")
ax.set_xlabel("word length"); ax.set_ylabel("how many"); ax.legend(fontsize=9)
ax.set_title("Counts: the prediction matches every length to 9", fontsize=10)
ax = axs[1]
ax.plot(L[1:], [c["trace"][i] / c["trace"][i-1] for i in range(1, len(L))], "o-", color=GREY, label="words")
ax.plot(L[1:], [c["distinct"][i] / c["distinct"][i-1] for i in range(1, len(L))], "s-", color=RED, label="maps")
ax.set_xlabel("word length"); ax.set_ylabel("growth from one length to the next"); ax.legend(fontsize=9)
ax.set_title("The maps grow more slowly: fewer pieces, smaller size", fontsize=10)
fig.suptitle("Figure 28. The break removes maps at every length from 2 on", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig28-counts.png", dpi=150); plt.close(fig)

# Fig 29: drop over the collision surface
g = np.array(B["B4_surface"])
fig, ax = plt.subplots(figsize=(7.4, 6))
sc = ax.scatter(g[:, 0], g[:, 1], c=g[:, 2] - g[:, 3], cmap="magma_r", s=34, marker="s")
plt.colorbar(sc, label="proved drop: D* minus the upper bound")
ax.plot([mu], [kap], "o", ms=12, mfc="none", mec=BLUE, mew=2); ax.text(mu - 0.003, kap + 0.005, "Figure 27's example", color=BLUE, fontsize=9, ha="right")
ax.set_xlabel("μ (class 1's second map)"); ax.set_ylabel("κ (class 3's map)")
ax.set_title("Figure 29. Every point breaks the same way (λ = μ(2 − κ), ν = κ/(2 − κ)),\n"
             "and the size drops below the formula by at least the colour", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig29-drop-surface.png", dpi=150); plt.close(fig)

# Fig 30: what is proved at the example
n = B["B1_numbers"]
fig, ax = plt.subplots(figsize=(10, 3.6))
D = n["Dstar"][0]; F = n["s_family_upper"][1]; Y = n["s_Y_upper"][1]; lo = max(v[0] for v in n["lower_bounds"].values())
ax.barh([0], [F - lo], left=[lo], color=GREEN, height=0.5)
ax.text((lo + F) / 2, 0, "true size is\nin here (proved)", ha="center", va="center", color="white", fontsize=10)
ax.axvline(D, color=RED, lw=2.5); ax.text(D + 0.0015, 0.45, f"formula D* = {D:.4f}", color=RED, fontsize=10)
ax.axvline(Y, color=GREY, ls=":"); ax.text(Y - 0.0015, -0.55, f"bound from one duplicate\npiece alone: {Y:.4f}", color="#555", ha="right", fontsize=9)
ax.text(F + 0.0015, 0.0, f"upper bound {F:.4f}", fontsize=10, va="center")
ax.text(lo - 0.0015, 0.0, f"lower bound {lo:.4f}", fontsize=10, ha="right", va="center")
ax.set_xlim(0.40, 0.535); ax.set_ylim(-0.75, 0.7); ax.set_yticks([]); ax.set_xlabel("dimension")
ax.set_title("Figure 30. At the break the true size is proved at least 0.033 below the formula", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig30-size-at-break.png", dpi=150); plt.close(fig)
print("figures 27-30 written")
