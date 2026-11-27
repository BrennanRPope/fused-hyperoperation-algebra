"""Round 10 figures (README section 14): fig31 why no other collision happens (the proof picture), fig32 the exact
size against the old bracket, fig33 where on the collision surface the exact formula is proved, fig34 why the
original 19/200 break is pinned but not exact."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"; PURPLE = "#7b4fa0"
E = json.load(open("exact.json"))
lam, mu, nu, kap = 231/3200, 3/80, 3/77, 3/40

# Fig 31: the proof picture
fig, ax = plt.subplots(figsize=(11.5, 5.6))
J1 = (0, 2 * lam); JN = (1 - nu, 1 + nu); JK = (2 - 2 * kap, 2)
def bar(iv, y, col, lw=10):
    ax.plot(iv, [y, y], color=col, lw=lw, solid_capstyle="butt")
bar((0, 2), 3.0, "k", 5); ax.text(1.0, 3.2, "every map sends [0, 2] into itself", ha="center")
bar(J1, 2.1, GREY); bar(JN, 2.1, BLUE); bar(JK, 2.1, ORANGE)
ax.text(np.mean(J1), 2.3, "zone A:\nlast step about 0", ha="center", color="#555555", fontsize=9)
ax.text(np.mean(JN), 2.3, "zone N:\nlast step about 1", ha="center", color=BLUE, fontsize=9)
ax.text(np.mean(JK), 2.3, "zone K:\nlast step about 2", ha="center", color=ORANGE, fontsize=9)
ax.text(2.3, 2.1, "three disjoint landing zones", va="center", fontsize=10)
ax.text(0.0, 1.45, "a second collision would leave, after cancelling,  h' = ρ · g'  with h' landing in zone N.\n"
        "The 3-adic step forces ρ = (77/40)^a with a ≥ 1, so ρ ≥ 1.925:", fontsize=10)
rho = 2 - kap; y = 0.7
bar(JN, y, BLUE, 9); ax.text(np.mean(JN), y - 0.3, "zone N\n(where h' lands)", ha="center", color=BLUE, fontsize=9)
bar((rho * JN[0], rho * JN[1]), y, BLUE, 5); ax.text(rho * np.mean(JN), y + 0.18, "1.925 × zone N", ha="center", color=BLUE, fontsize=9)
bar((rho * JK[0], rho * JK[1]), y, ORANGE, 5); ax.text(rho * np.mean(JK), y + 0.18, "1.925 × zone K", ha="center", color=ORANGE, fontsize=9)
ax.annotate("", xy=(rho * JN[0] - 0.03, y), xytext=(JN[1] + 0.03, y), arrowprops=dict(arrowstyle="->", color=GREEN))
ax.text(1.42, y + 0.18, "pushed clear: no overlap,\nso no second collision", color=GREEN, ha="center", fontsize=9)
ax.text(4.1, y, "larger a\npushes further", fontsize=9, va="center", color="#555555")
ax.set_xlim(-0.1, 4.7); ax.set_ylim(0.0, 3.5); ax.set_yticks([]); ax.set_xlabel("position on the line")
ax.set_title("Fig 31. Why the 3-adic break has no second collision (DRAFT)")
for s_ in ["top", "right", "left"]: ax.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig("figures/fig31-no-second-collision.png", dpi=130); plt.close()

# Fig 32: exact against the bracket
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
o = E["E5_original"]; e1 = E["E1"]
ax = axs[0]
ax.barh([1], [o["upper"][1] - 0.4367041], left=[0.4367041], color=GREY, height=0.35)
ax.barh([0], [o["upper"][1] - o["new_lower"][0]], left=[o["new_lower"][0]], color=GREEN, height=0.35)
ax.plot([o["upper"][1]] * 2, [-0.4, 1.4], color=BLUE, ls="--", lw=1)
ax.text(0.4367041, 1.27, "Round 9: between 0.4367 and 0.4705", fontsize=10)
ax.text(0.437, 0.22, "Round 10: between 0.47048331 and 0.47048344\n(width 0.00000014, a thin line at this scale)", fontsize=10, color=GREEN)
ax.text(o["upper"][1] - 0.0005, -0.33, "family bound 0.4704834", color=BLUE, fontsize=9, ha="right")
ax.set_yticks([0, 1]); ax.set_yticklabels(["now", "before"]); ax.set_xlim(0.432, 0.474)
ax.set_xlabel("dimension"); ax.set_title("the original break (19/200, 1/20, 1/19, 1/10)")
ax = axs[1]
ax.plot([e1["exact_dim"][0]], [0], "o", color=GREEN, ms=11)
ax.plot([e1["Dstar"][0]], [0], "o", color=GREY, ms=11)
ax.annotate("exact size 0.4232693\n(proved, no counting)", (e1["exact_dim"][0], 0), xytext=(0.395, 0.45), color=GREEN,
            arrowprops=dict(arrowstyle="->", color=GREEN))
ax.annotate("formula's guess D* = 0.4531865", (e1["Dstar"][0], 0), xytext=(0.43, -0.5), color="#555555",
            arrowprops=dict(arrowstyle="->", color=GREY))
ax.plot([e1["old_style_lower_bound"], e1["exact_dim"][0]], [-0.15, -0.15], color=GREY, lw=1)
ax.text(e1["old_style_lower_bound"], -0.3, "old-style lower bound 0.3918", fontsize=9, color="#555555")
ax.set_xlim(0.385, 0.47); ax.set_ylim(-0.8, 0.9); ax.set_yticks([]); ax.set_xlabel("dimension")
ax.set_title("the 3-adic break (3/80, 3/40 on the surface)")
fig.suptitle("Fig 32. From a bracket to an exact number (DRAFT)")
plt.tight_layout(); plt.savefig("figures/fig32-exact-vs-bracket.png", dpi=130); plt.close()

# Fig 33: where the exact formula is proved
g = np.array([[r[0], r[1], r[5] if r[5] is not None else np.nan, r[6] if r[6] is not None else np.nan] for r in E["E3_region"]], dtype=float)
fig, ax = plt.subplots(figsize=(8.6, 6.6))
okm = ~np.isnan(g[:, 2])
sc = ax.scatter(g[okm, 0], g[okm, 1], c=g[okm, 3] - g[okm, 2], s=14, cmap="viridis", marker="s")
ax.scatter(g[~okm, 0], g[~okm, 1], color="#eeeeee", s=14, marker="s")
cb = plt.colorbar(sc); cb.set_label("size lost to the collision (D* minus exact size)")
ax.plot([3/80], [3/40], "*", color=RED, ms=16); ax.annotate("3-adic break\n(exact 0.42327)", (3/80, 3/40), xytext=(0.09, 0.025), color=RED, bbox=dict(fc="white", ec="none", alpha=0.85),
                                                         arrowprops=dict(arrowstyle="->", color=RED))
ax.plot([1/20], [1/10], "X", color="k", ms=11); ax.annotate("original break\n(pinned, not exact)", (1/20, 1/10), xytext=(0.10, 0.17), bbox=dict(fc="white", ec="none", alpha=0.85),
                                                         arrowprops=dict(arrowstyle="->"))
seq = E["E4_sequence"]
ax.plot([eval(s_["mu"]) for s_ in seq], [eval(s_["kappa"]) for s_ in seq], "o", color=ORANGE, ms=5, mfc="none")
ax.text(0.005, 0.45, "coloured: the exact formula is proved at every rational point here whose\n"
        "two numerators share an odd prime (a dense set); grey: hypotheses (I), (II) fail\n"
        "orange rings: exact points p/(20p+1), p/(10p+1) closing in on the original", fontsize=9, va="top")
ax.set_xlabel("μ (second map about 0);  λ = μ(2 − κ)"); ax.set_ylabel("κ (map about 2);  ν = κ/(2 − κ)")
ax.set_title("Fig 33. Where on the collision surface the exact size is proved (DRAFT)")
plt.tight_layout(); plt.savefig("figures/fig33-where-exact.png", dpi=130); plt.close()

# Fig 34: the original break's ratios and the bad windows
L0, M0 = 19/200, 1/20
fig, ax = plt.subplots(figsize=(11, 5))
pts = [(a, b, L0**a / M0**b) for a in range(1, 31) for b in range(1, 31)]
pts = [p for p in pts if 0.2 < p[2] < 3]
bN = [float(eval(v)) for v in o["bad_N"]]; bK = [float(eval(v)) for v in o["bad_K"]]
ax.axhspan(*bN, color=RED, alpha=0.15); ax.axhspan(*bK, color=RED, alpha=0.15)
ax.text(30.5, np.sqrt(bN[0] * bN[1]), "bad window\n(N against N)", color=RED, va="center", fontsize=9)
ax.text(30.5, np.sqrt(bK[0] * bK[1]), "bad window\n(K against N)", color=RED, va="center", fontsize=9)
for a, b, r in pts:
    badp = bN[0] <= r <= bN[1] or bK[0] <= r <= bK[1]
    ax.plot(a, r, "o", color=RED if badp else BLUE, ms=6 if badp else 4)
ax.axvline(12.5, color=GREEN, ls="--"); ax.text(12.3, 2.4, "runs with at most 12 copies of λ:\nno ratio hits a window",
                                                 ha="right", color=GREEN, fontsize=10)
ax.annotate("first hit: λ¹³/μ¹⁰ = 0.5257", (13, L0**13 / M0**10), xytext=(15, 0.27), arrowprops=dict(arrowstyle="->", color=RED), color=RED)
ax.set_yscale("log"); ax.set_xlim(0, 34); ax.set_xlabel("a, the power of λ left after cancelling")
ax.set_ylabel("stretch ρ = λ^a / μ^b  (log scale)")
ax.set_title("Fig 34. The original break: no shared prime, so ρ can hit a bad window once a ≥ 13 (DRAFT)")
plt.tight_layout(); plt.savefig("figures/fig34-original-windows.png", dpi=130); plt.close()
print("figures 31-34 written")
