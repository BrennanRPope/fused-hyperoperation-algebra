"""Round 11 figures (README section 15): fig35 every bad pair settled (stretch and size of the exact search),
fig36 the gap at the original break, round by round, fig37 why no single argument covers all pairs yet."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"
D = json.load(open("original_descent.json")); B = json.load(open("original_bracket.json"))
P = D["O3_descent"]["per_pair"]           # (a, b, rho, isbad, status, nodes, depth), bad pairs only
a = np.array([p[0] for p in P]); rho = np.array([p[2] for p in P]); nodes = np.array([p[5] for p in P])
inN = (rho > 0.8)
AM = D["AMAX"]

# Fig 35
fig, axs = plt.subplots(2, 1, figsize=(11.5, 7.2), sharex=True)
ax = axs[0]
ax.axhspan(0.9, 10/9, color=BLUE, alpha=0.12); ax.axhspan(9/19, 100/171, color=ORANGE, alpha=0.12)
ax.scatter(a[inN], rho[inN], s=5, color=BLUE); ax.scatter(a[~inN], rho[~inN], s=5, color=ORANGE)
ax.text(AM * 1.01, 1.0, "N window\n[0.9, 1.11]", color=BLUE, va="center", fontsize=9)
ax.text(AM * 1.01, 0.53, "K window\n[0.474, 0.585]", color=ORANGE, va="center", fontsize=9)
ax.set_ylabel("stretch ρ = λ^a / μ^b"); ax.set_xlim(0, AM * 1.13)
ax.set_title(f"Fig 35. Every bad pair with up to {AM:,} copies of λ, settled exactly (DRAFT)")
ax.text(AM * 0.5, 0.74, f"{len(P):,} bad pairs (a, b): Round 10's push says nothing here. First one: a = 13, b = 10.",
        ha="center", va="center", fontsize=9)
ax = axs[1]
ax.scatter(a, nodes, s=5, color=GREEN); ax.set_yscale("log"); ax.set_ylim(10, max(nodes) * 12)
ax.set_ylabel("nodes in the exact search"); ax.set_xlabel("a = extra copies of λ in the run")
ax.text(0.01, 0.95, f"every search ended with every branch cut off: no second collision for any of them "
        f"(largest search {max(nodes):,} nodes)", transform=ax.transAxes, va="top", fontsize=9, color=GREEN)
plt.tight_layout(); plt.savefig("figures/fig35-bad-pairs-settled.png", dpi=130); plt.close()

# Fig 36
fig, ax = plt.subplots(figsize=(11, 3.8))
labels = ["Round 9\n(collision-free subsystem)", "Round 10\n(runs of at most 12 λ)", f"Round 11\n(runs of at most {AM:,} λ)"]
gaps = [0.4704834 - 0.4367041, B["round10"]["gap"], None]
lg = [np.log10(gaps[0]), np.log10(gaps[1]), B["log10_gap"]]
ax.barh([2, 1, 0], [-v for v in lg], color=[GREY, BLUE, GREEN], height=0.5)
for y, v, txt in zip([2, 1, 0], lg, ["0.034", "1.3 × 10⁻⁷", f"10^{B['log10_gap']:.0f}"]):
    ax.text(-v * 1.12, y, "gap " + txt, va="center", fontsize=10)
ax.set_yticks([2, 1, 0]); ax.set_yticklabels(labels); ax.set_xscale("log")
ax.set_xlabel("digits pinned down (−log₁₀ of the gap between the proved lower and upper bounds)")
ax.set_xlim(0.8, max(-x for x in lg) * 4)
ax.set_title("Fig 36. The original break's size, pinned tighter each round, still not exact (DRAFT)")
plt.tight_layout(); plt.savefig("figures/fig36-gap-by-round.png", dpi=130); plt.close()

# Fig 37
sig = np.where(inN, rho, 19 * rho / 10)
dev = np.abs(sig - 1)
order = np.argsort(a); run = np.minimum.accumulate(dev[order])
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.scatter(a, dev, s=4, color=GREY, label="each bad pair: |σ − 1|")
ax.step(a[order], run, where="post", color=RED, lw=2, label="closest to 1 so far")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("a = extra copies of λ"); ax.set_ylabel("|σ − 1|")
ax.legend(loc="upper left", fontsize=9)
ax.set_title("Fig 37. Why one argument does not cover every pair yet (DRAFT)")
ax.text(0.06, 0.27, "After one more cancel, each bad pair asks: can h + 1 = σ (g + 1)?\n"
        "σ is never 1 (a 19-adic count), but it comes as close to 1 as you like\n"
        "(log λ / log μ is irrational). Near σ = 1 a zone argument has nothing to\n"
        "push with, so each σ needs its own exact search.", transform=ax.transAxes, ha="left", va="top", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9))
plt.tight_layout(); plt.savefig("figures/fig37-why-not-uniform.png", dpi=130); plt.close()
print("figures 35-37 written")
