"""Round 7 figures (README section 11): fig19 pieces and certified discs, fig20 the 3-adic step, fig21 the last step
for the original example, fig22 what is proved about the size."""
import json, random, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from fractions import Fraction as F
import plane as P
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"; PURPLE = "#7b4fa0"
COL = [BLUE, ORANGE, PURPLE]
R = json.load(open("plane.json"))

def chaos(name, n=120000, seed=1):
    ex = P.EXAMPLES[name]
    maps = [(complex(float(c[0]), float(c[1])), complex(float(l[0]), float(l[1])), j)
            for j, (c, ls) in enumerate(ex) for l in ls]
    rng = random.Random(seed); z = 0.3 + 0.3j; pts = []; cls = []
    for t in range(n + 50):
        c, l, j = maps[rng.randrange(len(maps))]; z = c + l * (z - c)
        if t >= 50: pts.append(z); cls.append(j)
    return np.array(pts), np.array(cls)

SETUP = {"ORIG": ((F(3,5), F(3,8)), F(87,100)), "MOD": ((F(29,50), F(19,50)), F(17,20))}
LAM = {"ORIG": "class 2: i/3,  class 3: 1/3", "MOD": "class 2: 3i/10,  class 3: 3/10"}

# Fig 19: the fractal by outer class, and the certified disc covers that keep the classes apart
fig, axs = plt.subplots(1, 2, figsize=(13, 6.4))
for ax, name in zip(axs, ["ORIG", "MOD"]):
    cB, Rd = SETUP[name]
    res, X, B2 = P.pingpong(name, cB, Rd, F(1, 60), F(1, 60))
    for j, cover in enumerate(X):
        for z, r in cover:
            ax.add_patch(Circle((float(z[0]), float(z[1])), float(r), fill=False, lw=0.35, color=COL[j], alpha=0.35))
    ax.add_patch(Circle((float(cB[0]), float(cB[1])), float(Rd), fill=False, lw=1.2, ls="--", color=GREY))
    pts, cls = chaos(name)
    for j in range(3):
        sel = cls == j
        ax.scatter(pts[sel].real, pts[sel].imag, s=0.15, color=COL[j])
    for (x, y), lab, j in [((0, 0), "centre 0", 0), ((1, 0), "centre 1", 1), ((0, 1), "centre i", 2)]:
        ax.plot([x], [y], "k+", ms=10); ax.text(x + 0.03, y - 0.07, lab, fontsize=9)
    ax.text(-0.32, 1.22, "dashed: the invariant disc B", color=GREY, fontsize=9)
    ax.set_aspect("equal"); ax.set_xlim(-0.35, 1.5); ax.set_ylim(-0.5, 1.3)
    gaps = res["min_gap_float"]
    ax.set_title(f"{'Original (round 2)' if name == 'ORIG' else 'Modified example (round 7)'}\n"
                 f"class 1 (blue): (1+i)/4, 3/10;  {LAM[name]}  (2 orange, 3 purple)\n"
                 f"smallest gap between covers: {min(gaps.values()):.2f} (certified: {res['certified']})", fontsize=10)
fig.suptitle("Figure 19. The plane union fractal, coloured by the class of the last map applied. Faint circles: the "
             "exact-fraction disc covers of the\nthree class pieces. No circle of one colour meets a circle of another, "
             "so the class of a word's outer run can be read off from where it sends a point.", fontsize=11.5)
plt.tight_layout(); plt.savefig("figures/fig19-plane-pieces.png", dpi=130, bbox_inches="tight"); plt.close()

# Fig 20: the 3-adic step
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
for ax, name, col in zip(axs, ["MOD", "ORIG"], [GREEN, RED]):
    d = {int(k): v for k, v in R["C3"][name].items()}; tot = sum(d.values()); xs = sorted(d)
    ax.bar(xs, [d[x] / tot for x in xs], color=col, width=0.8)
    ax.set_xlim(-20, 4); ax.set_xlabel("3-adic valuation of g′(1)")
    if name == "MOD":
        ax.set_title("Modified: every word lands on valuation 0 (proved)\n"
                     "λ₁ᵅ·λ₂⁻ᵝ has valuation −β ≠ 0, so it can't link two such points", fontsize=10.5)
        ax.annotate("", xy=(-3, 0.5), xytext=(-0.5, 0.5), arrowprops=dict(arrowstyle="->", lw=2, color=BLUE))
        ax.text(-12.5, 0.55, "multiplying by λ₁ᵅλ₂⁻ᵝ moves a\npoint left by β ≥ 1, off the bar", color=BLUE, fontsize=9)
    else:
        ax.set_title(f"Original: words spread over {len(xs)} valuations\n(so this 3-adic proof can't work there)", fontsize=10.5)
axs[0].set_ylabel("share of words")
fig.suptitle("Figure 20. The arithmetic step. Random words whose outer run is class 2 or 3, applied to the point 1, "
             "and measured 3-adically\n(how many times 3 divides the result). Classes 2 and 3 of the modified example "
             "shrink 3-adically toward their centres 1 and i.", fontsize=11.5)
plt.tight_layout(); plt.savefig("figures/fig20-three-adic.png", dpi=130, bbox_inches="tight"); plt.close()

# Fig 21: the last step for the original example: which rho = lam1^d lam2^-e can link two class-2/3 points
pts, cls = chaos("ORIG", n=60000, seed=3)
W = pts[cls > 0]; rng = np.random.default_rng(0)
a = W[rng.integers(0, len(W), 40000)]; b = W[rng.integers(0, len(W), 40000)]; rat = a / b
cB, Rd = SETUP["ORIG"]; res, X, B2 = P.pingpong("ORIG", cB, Rd, F(1, 60), F(1, 60))
WD = X[1] + X[2]
Z = np.array([complex(float(z[0]), float(z[1])) for z, r in WD]); Rr = np.array([float(r) for z, r in WD])
def inside(rho):
    D = np.abs(rho * Z[None, :] - Z[:, None]) - (abs(rho) * Rr[None, :] + Rr[:, None])
    return D.min() <= 0
fig, ax = plt.subplots(figsize=(11, 5.6))
ax.scatter(np.degrees(np.angle(rat)), np.log10(np.abs(rat)), s=0.3, color=GREY, alpha=0.5)
l1 = (1 + 1j) / 4; opened = []
for e in range(1, 40):
    for d in range(1, 50):
        rho = l1**d / 0.3**e
        x, y = np.degrees(np.angle(rho)), np.log10(abs(rho))
        if abs(y) > 0.6: continue
        if e <= 6:
            ax.plot(x, y, "o", color=GREEN, ms=7)
        elif inside(rho):
            ax.plot(x, y, "o", color=RED, ms=7); opened.append((d, e))
            if len(opened) <= 6: ax.annotate(f"({d},{e})", (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8, color=RED)
        else:
            ax.plot(x, y, ".", color=BLUE, ms=5)
ax.plot([], [], "o", color=GREEN, label="β-gap e ≤ 6: certified outside the covers (all d)")
ax.plot([], [], ".", color=BLUE, label="e ≥ 7, outside the covers")
ax.plot([], [], "o", color=RED, label="e ≥ 7, inside the covers: the only cases left open")
ax.scatter([], [], s=8, color=GREY, label="ratios w₁/w₂ of points in the class-2 and class-3 pieces")
ax.set_xlabel("angle of the ratio (degrees)"); ax.set_ylabel("log₁₀ of its size"); ax.set_ylim(-0.6, 0.6)
ax.legend(fontsize=9, loc="lower left")
fig.suptitle("Figure 21. The last step for the original example. A collision would need ρ = λ₁ᵈ·λ₂⁻ᵉ to equal a ratio of two "
             "points of the class-2/3 pieces.\nThe certified covers rule that out for every d whenever e ≤ 6, which "
             "is enough for words with at most 6 b's per class-1 run.", fontsize=11.5)
plt.tight_layout(); plt.savefig("figures/fig21-last-step-original.png", dpi=130, bbox_inches="tight"); plt.close()

# Fig 22: what is proved about the size
fig, ax = plt.subplots(figsize=(11, 4.4))
o_chain = R["C6_ORIG_chain"][0]; o_sub = R["C8_ORIG_beta_le_6"][0]; o_D = R["C5_ORIG"][1]; m_D = R["C5_MOD"]
ax.plot([o_sub], [1], "D", color=GREEN, ms=9)
ax.plot([o_chain], [1], "|", color=GREEN, ms=22, mew=2)
ax.text(o_chain, 1.2, f"{o_chain:.4f}\n(pure a-runs, round 7)", ha="center", fontsize=9)
ax.annotate(f"original: between {o_sub:.5f} and {o_D:.5f}, proved", xy=(o_D, 1.08), xytext=(1.07, 1.42),
            fontsize=10, color=GREEN, arrowprops=dict(arrowstyle="->", color=GREEN))
ax.plot([1.194], [1], "x", color=GREY, ms=9, mew=2); ax.text(1.194, 0.7, "box count 1.194\n(round 2)", ha="center", fontsize=8.5, color=GREY)
ax.plot([1.251], [1], "x", color=GREY, ms=9, mew=2); ax.text(1.251, 0.7, "if no maps commuted\n1.251", ha="center", fontsize=8.5, color=GREY)
ax.plot(m_D, [0, 0], color=GREEN, lw=12, solid_capstyle="butt"); ax.plot([m_D[0]], [0], "o", color=GREEN, ms=8)
ax.text(m_D[0] + 0.004, -0.07, f"modified: exactly D* = {m_D[0]:.5f}, proved", fontsize=10, color=GREEN)
ax.axvline(1, color="k", lw=0.8, ls=":"); ax.text(1.002, -0.42, "1 (a curve)", fontsize=8.5)
ax.set_yticks([0, 1]); ax.set_yticklabels(["modified", "original"]); ax.set_ylim(-0.5, 1.7); ax.set_xlim(0.99, 1.27)
ax.set_xlabel("Hausdorff dimension (all green ranges rest only on Hochman 2015, Corollary 1.7)")
fig.suptitle("Figure 22. The true size of the plane union fractal. The original's range is under 0.00003 wide; it closes at "
             "D* if the last open\ncases (red in Figure 21) never collide, which holds for every word up to length 11.",
             fontsize=11.5)
plt.tight_layout(); plt.savefig("figures/fig22-plane-size.png", dpi=130, bbox_inches="tight"); plt.close()
print("ok", opened[:8])
