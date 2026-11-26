"""fig7: the four set-operations levers, each with what the literature already has (orange) and what is ours (blue)."""
import math, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from setops_levers import E_k, draws_p, neg_mass
BLUE, ORANGE, GREY = "#1f5fbf", "#e07b00", "#8a8a8a"
plt.rcParams.update({"font.size": 13})
f, ax = plt.subplots(2, 2, figsize=(15, 11))
# A: E(k)
a = ax[0, 0]; ks = np.linspace(-0.9, 6, 400)
a.plot(ks, [E_k(k) for k in ks], color=BLUE, lw=3, label="E(k) = (1+k)^(1/k) + k/(1+k)  (general formula: not found in print)")
for k, lab in [(-0.5, "3 at k = -1/2"), (0, "e at k = 0"), (1, "5/2 at k = 1")]:
    a.plot(k, E_k(k), "o", color=BLUE, ms=10); a.annotate(lab, (k, E_k(k)), xytext=(10, 8), textcoords="offset points")
a.plot(math.e - 1, E_k(math.e - 1), "D", color=ORANGE, ms=12, label="k = e - 1: Vandervelde 2008 (known)")
a.axhline(2, color=GREY, ls=":"); a.text(4.3, 2.03, "k → ∞: 2", color=GREY)
a.set_ylim(1.8, 4.5); a.set_xlabel("dial k"); a.set_ylabel("mean number of draws")
a.set_title("1. e along the k dial: draws until x₁∘ₖx₂∘ₖ… passes 1"); a.legend(loc="upper right", fontsize=11)
# B: draws vs p
b = ax[0, 1]; ps = np.linspace(0.5, 4, 200)
b.plot(ps, [draws_p(p) for p in ps], color=ORANGE, lw=3, label="draws until Σ xᵖ passes 1 (problem: Newman & Klamkin 1959)")
for p, lab in [(1, "e (simplex)"), (2, "3.9258 (disc / sphere)"), (3, "5.161")]:
    b.plot(p, draws_p(p), "o", color=BLUE, ms=10); b.annotate(lab, (p, draws_p(p)), xytext=(10, -14), textcoords="offset points")
b.set_xlabel("ℓp exponent p"); b.set_ylabel("mean number of draws")
b.set_title("2. e along p: a Mittag-Leffler value E₁/ₚ(Γ(1+1/p))"); b.legend(loc="upper left", fontsize=11)
# C: k = -1/chi
c = ax[1, 0]; c.set_xlim(0, 10); c.set_ylim(0, 6); c.axis("off")
c.set_title("3. k = −1/χ(whole space): the dial is set by where the shapes sit")
rows = [("closed segment [0,1]", 1, "k = −1  (probability)", 4.8, True), ("box of size M", None, "k = −1/M  (overlap shrinks)", 3.3, True),
        ("open line ℝ", -1, "k = +1  (signed size)", 1.8, False)]
for name, chi, ktxt, y, closed in rows:
    c.plot([0.5, 3.5], [y, y], color=BLUE, lw=5)
    if chi == 1: c.plot([0.5, 3.5], [y, y], "o", color=BLUE, ms=11)
    if chi == -1: c.plot([0.5, 3.5], [y, y], "o", color="white", mec=BLUE, mew=3, ms=11)
    c.text(0.5, y + 0.35, name + (f",  χ = {chi:+d}" if chi else ""), fontsize=12)
    c.text(4.6, y - 0.1, ktxt, fontsize=13, color=BLUE)
c.text(0.5, 0.5, "Euler characteristic with open cells = −1: Schanuel 1991, Propp 2003 (known).\n"
       "Restricting to a set A turns the dial to −1/a: that is conditioning (known). The reading is ours.", fontsize=11, color=GREY)
# D: Bell
d = ax[1, 1]; Ss = np.linspace(1.6, 4, 25)
d.plot([2, 4], [0, 0.5], color=ORANGE, lw=3, label="(S − 2)/4, in print\n(Morris et al. 2022)")
d.plot(Ss, [neg_mass(S) for S in Ss], "o", color=BLUE, ms=7, label="our LP, 16 corner\nstrategies")
d.axvline(2, color=GREY, ls=":"); d.axvline(2 * math.sqrt(2), color=GREY, ls=":"); d.axvline(4, color=GREY, ls=":")
for x, lab in [(2, "Bell 2"), (2 * math.sqrt(2), "Tsirelson 2√2\nmass (√2−1)/2"), (4, "no-signalling 4")]:
    d.text(x - 0.04, {2: 0.33, 2 * math.sqrt(2): 0.40}.get(x, 0.20), lab, fontsize=11, color=GREY, ha="right")
d.set_xlabel("CHSH value S"); d.set_ylabel("least total negative size")
d.set_title("4. Beating Bell needs negative size (S − 2)/4"); d.legend(loc="lower right", fontsize=11)
f.suptitle("Fig 7. Four set-operations facts as levers. Blue = ours, orange = already in print", fontsize=15)
f.tight_layout(rect=[0, 0, 1, 0.96]); f.savefig("figures/fig7-set-operations-levers.png", dpi=110); print("ok")
