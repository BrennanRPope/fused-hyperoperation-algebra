"""Round 8 figures (README section 12): fig23 the prime recipe of Corollary P9, fig24 which family members pass,
fig25 the dense family in a slice of parameter space, fig26 a dense-family attractor in the plane."""
import json, random, math, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from fractions import Fraction as F
import general as G
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"; PURPLE = "#7b4fa0"
COL = [BLUE, ORANGE, PURPLE, GREEN]
D = json.load(open("general.json"))

# Fig 23: the valuation pattern of a (3, 2, 1) member on the line
rng = random.Random(23)
cl, own = G.corollary_member((3, 2, 1), "Q", rng)
primes = [p for o in own for p in o]
rows, labels = [], []
for i, (c, ls) in enumerate(cl):
    for j, l in enumerate(ls):
        rows.append([G.val(l, p, "Q") for p in primes]); labels.append(f"class {i+1}, map {j+1}:  {l[0]}")
M = np.array(rows, dtype=float)
fig, ax = plt.subplots(figsize=(8.6, 5.2))
ax.imshow(M > 0, cmap="Greens", vmin=0, vmax=1.6, aspect="auto")
for r in range(M.shape[0]):
    for k in range(M.shape[1]):
        ax.text(k, r, int(M[r, k]), ha="center", va="center", fontsize=12)
ax.set_xticks(range(len(primes)), [f"p = {p}\n(class {next(i+1 for i,o in enumerate(own) if p in o)}'s prime)" for p in primes])
ax.set_yticks(range(len(labels)), labels, fontsize=10)
for b in np.cumsum([len(ls) for _, ls in cl])[:-1]: ax.axhline(b - 0.5, color="k", lw=1)
ax.set_title("Figure 23. The prime recipe (Corollary P9), one member with classes of sizes 3, 2, 1\n"
             "Each cell: how many times the prime divides the multiplier", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig23-prime-recipe.png", dpi=150); plt.close(fig)

# Fig 24: which sampled members pass; x = largest |lambda| over the Lemma A' bound, y = D*
fig, ax = plt.subplots(figsize=(8.6, 5.4))
for fam, mk, lab in [("F1", "o", "F1 plane, 3 classes"), ("F2", "s", "F2 plane, 4 classes"), ("F3", "^", "F3 line, 3 classes")]:
    for r in D["G3_" + fam]["rows"]:
        cls = [((F(c[0]), F(c[1])), [(F(l[0]), F(l[1])) for l in ls]) for c, ls in r["classes"]]
        okE, info = G.lemma_A_explicit(cls); x = info["eps"] / info["bound"]
        y = r["Dstar"][0] if "Dstar" in r else float(G.certified_Dstar(cls)[0])
        ax.scatter([x], [y], marker=mk, s=38, color=GREEN if (r["A"] and r["B"]) else RED, alpha=0.8, edgecolor="none")
for r in D["G5"]["rows"]:
    ax.scatter([r["eps"] / r["bound"]], [r["Dstar"][0]], marker="D", s=26, color=BLUE, alpha=0.8, edgecolor="none")
ax.axvline(1, color=GREY, ls="--"); ax.text(1.03, 0.42, "Lemma A′ bound:\nleft of this line,\nseparation is automatic", fontsize=9, color="#555")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GREEN, label="passes (A) and (B): no collisions, exact size"),
                   Line2D([], [], marker="o", ls="", color=RED, label="certified cover not found (test silent)"),
                   Line2D([], [], marker="D", ls="", color=BLUE, label="dense-family members, 14 class patterns"),
                   Line2D([], [], marker="o", ls="", color="k", label="circle F1, square F2, triangle F3")],
          fontsize=9, loc="upper left")
ax.set_xlabel("largest |λ| ÷ Lemma A′ bound δ/(4D + δ)"); ax.set_ylabel("formula D* (certified)")
ax.set_title("Figure 24. Every sampled member that passes the test, and where it sits", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig24-family-pass.png", dpi=150); plt.close(fig)

# Fig 25: a slice of parameter space, line, centres 0, 1, 2; class 1 = {lam, mu}, classes 2, 3 = {3/32}
def Dstar(ls_by_class):
    f = lambda s: 1 + sum(np.prod([1 - abs(l)**s for l in ls]) - 1 for ls in ls_by_class)
    return brentq(f, 1e-4, 1.5)
n = 7; B = 1 / 9; pts = []
for a in range(-int(B * 2**n), int(B * 2**n) + 1):
    for b in range(-int(B * 2**n), int(B * 2**n) + 1):
        if a == 0 or b == 0 or a % 3 == 0 or b % 3 != 0: continue
        lam, mu = a / 2**n, b / 2**n
        if abs(lam) < B and abs(mu) < B: pts.append((lam, mu, Dstar([[lam, mu], [3/32], [3/32]])))
pts = np.array(pts)
fig, ax = plt.subplots(figsize=(7.4, 6.4))
sc = ax.scatter(pts[:, 0], pts[:, 1], c=pts[:, 2], s=14, cmap="viridis")
ax.add_patch(plt.Rectangle((-B, -B), 2 * B, 2 * B, fill=False, ls="--", color=GREY))
ax.set_xlabel("λ  (class 1's first map; numerator not divisible by 3)")
ax.set_ylabel("μ  (class 1's second map; numerator divisible by 3)")
plt.colorbar(sc, label="dimension = D*, proved at every dot")
ax.set_title(f"Figure 25. The dense family in one slice (denominator 2^{n} shown; finer\n"
             "denominators fill the dashed box). Line, centres 0, 1, 2; classes 2 and 3 shrink by 3/32", fontsize=10)
ax.set_aspect("equal"); fig.tight_layout(); fig.savefig("figures/fig25-dense-slice.png", dpi=150); plt.close(fig)

# Fig 26: one plane dense-family member, coloured by outer class, with the Lemma A' discs
row = next(r for r in D["G5"]["rows"] if r["field"] == "Qi" and r["pattern"] == [2, 2, 1])
cls = [((float(F(c[0])), float(F(c[1]))), [complex(float(F(l[0])), float(F(l[1]))) for l in ls]) for c, ls in row["classes"]]
maps = [(complex(*c), l, j) for j, (c, ls) in enumerate(cls) for l in ls]
rng = random.Random(26); z = 0.2 + 0.2j; P_, C_ = [], []
for t in range(150050):
    c, l, j = maps[rng.randrange(len(maps))]; z = c + l * (z - c)
    if t >= 50: P_.append(z); C_.append(j)
P_ = np.array(P_); C_ = np.array(C_)
fig, axs = plt.subplots(1, 2, figsize=(13, 6.4))
cc = [complex(*c) for c, _ in cls]; eps = row["eps"]
Dm = max(abs(a - b) for a in cc for b in cc); Rad = Dm * (1 + eps) / (1 - eps)
ax = axs[0]
ax.add_patch(plt.Circle((cc[0].real, cc[0].imag), Rad, fill=False, ls="--", color=GREY))
for j, c in enumerate(cc):
    sel = C_ == j; ax.scatter(P_[sel].real, P_[sel].imag, s=1.5, color=COL[j])
    ax.add_patch(plt.Circle((c.real, c.imag), eps * (Rad + abs(c - cc[0])), fill=False, color=COL[j], lw=1.2))
    ax.text(c.real + 0.1, c.imag + 0.25, f"class {j+1}", color=COL[j], fontsize=10)
ax.set_xlim(-Rad - 0.1, Rad + 0.1); ax.set_ylim(-Rad - 0.1, Rad + 0.1); ax.set_aspect("equal")
ax.set_title("Whole picture: dashed = invariant disc,\ncoloured = Lemma A′'s class discs (never meet)", fontsize=10)
ax = axs[1]; sel = C_ == 0
ax.scatter(P_[sel].real, P_[sel].imag, s=0.3, color=COL[0])
for j in (1, 2):
    sel = C_ == j; ax.scatter(P_[sel].real - (cc[j] - cc[0]).real, P_[sel].imag - (cc[j] - cc[0]).imag, s=0.3, color=COL[j], alpha=0.0)
ax.set_xlim(-0.14, 0.06); ax.set_ylim(-0.13, 0.07); ax.set_aspect("equal")
ax.set_title("Zoom on class 1's piece (near centre 0): small copies\nof the whole, turned by the non-real multipliers", fontsize=10)
fig.suptitle(f"Figure 26. A plane member of the dense family (classes of sizes 2, 2, 1; proved dimension D* = {row['Dstar'][0]:.4f})", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig26-plane-member.png", dpi=150); plt.close(fig)
print("figures 23-26 written")
