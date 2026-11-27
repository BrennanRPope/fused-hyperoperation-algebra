"""Round 22 pictures (DRAFT): fig71 the cone of word-shaped ratios in the coordinates 2^j (19/10)^E (1/20)^F, and why two
numbers against one give only 1 + 19 = 20; fig72 where Conjecture 21.2 stands for each pair of side sizes.
Usage: python3 figures22.py (needs explore22.json)."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

J = json.load(open("explore22.json"))

# ---------------------------------------------------------------- fig71
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
cols = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c"}
for j, off in ((3, 0.18), (2, 0.0), (1, -0.18)):
    pts = [(E, Fv) for Fv in range(j, 8) for E in range(-j, Fv - j + 1)]
    ax1.scatter([E + off for E, Fv in pts], [Fv for E, Fv in pts], s=36, color=cols[j], label=f"j = {j} (a product of {j} word ratio{'s' if j > 1 else ''})", zorder=3)
for (E, Fv, lab) in ((-1, 1, "N  (1/19)"), (0, 1, "K  (1/10)"), (-1, 2, "mN"), (0, 2, "lN"), (1, 2, "lK")):
    ax1.annotate(lab, (E - 0.18, Fv), xytext=(E - 1.6, Fv + 0.55), fontsize=8.5, color="#1f77b4",
                 arrowprops=dict(arrowstyle="-", color="#1f77b4", lw=0.7))
ax1.set_xlabel("E  (power of 19/10)"); ax1.set_ylabel("F  (power of 1/20)")
ax1.set_xlim(-4.2, 7); ax1.set_ylim(0.3, 7.6); ax1.grid(alpha=0.3); ax1.legend(fontsize=8, loc="lower right")
ax1.set_title("A. Every ratio a word can step by is 2^j (19/10)^E (1/20)^F with F >= j\nand -j <= E <= F - j (Proposition 22.2): "
              "a cone, one layer per j", fontsize=9.5)
ax2.axis("off"); ax2.set_xlim(-0.25, 10.25); ax2.set_ylim(0, 10)
def box(x, y, wdt, h, txt, col):
    ax2.add_patch(FancyBboxPatch((x, y), wdt, h, boxstyle="round,pad=0.08", fc=col, ec="#444444"))
    ax2.text(x + wdt / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.6)
box(2.6, 8.0, 4.8, 1.4, "a + a rho = b,  rho = 19^E / (2^A 5^B) in the cone:\nis 1 + rho a product of 2, 5, 19?", "#eeeeee")
box(0.0, 4.6, 2.3, 2.3, "E >= 1\n19^E + 2^A 5^B\nhas no factor\n2, 5 or 19: no", "#f8d7d7")
box(2.55, 4.6, 2.3, 2.3, "E = 0\n1 + 2^A 5^B is\n0 or 2 mod 3,\nnever 19^g: no", "#f8d7d7")
box(5.1, 4.6, 2.3, 2.3, "E < 0, A, B >= 1\n1 + 19^e 2^A 5^B\nhas no factor\n2, 5 or 19: no", "#f8d7d7")
box(7.65, 4.6, 2.3, 2.3, "rho = 1/19^e\n1 + 19^e = 2^x 5^y\nonly for e = 1\n(mod 8, LTE)", "#d9f2d9")
for x in (1.15, 3.7, 6.25, 8.8):
    ax2.annotate("", xy=(x, 7.0), xytext=(5.0, 7.95), arrowprops=dict(arrowstyle="->", color="#444444"))
box(6.4, 1.2, 3.6, 1.8, "19c + c = 20c:\nthe one balance\n(Proposition 22.1)", "#d9f2d9")
ax2.annotate("", xy=(8.2, 3.05), xytext=(8.8, 4.55), arrowprops=dict(arrowstyle="->", color="#444444"))
ax2.set_title("B. Two numbers against one, proved by hand: only 1 + 19 = 20", fontsize=9.5)
fig.suptitle("Fig 71 (DRAFT). The word-shaped ratios form a simple cone, and its smallest case of Conjecture 21.2 is proved", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig("figures/fig71-the-cone.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig72
fig, ax = plt.subplots(figsize=(9.5, 6.2))
status = {}   # (p, q) with p >= q: (colour, text)
status[(1, 1)] = ("#2ca02c", "trivial")
status[(2, 1)] = ("#2ca02c", "proved\n(22.1)")
status[(2, 2)] = ("#98df8a", "searched\n10^300")
for q in (1, 2, 3):
    status[(3, q)] = ("#ffbb78", "searched\n10^35*")
status[(4, 1)] = status[(4, 2)] = status[(4, 3)] = status[(4, 4)] = ("#ffe9cc", "searched\n10^15")
for q in range(1, 6):
    status[(5, q)] = ("#f2f2f2", "searched\n10^12")
for q in range(1, 7):
    status[(6, q)] = ("#f7f7f7", "searched\n10^10")
for (p, q), (c, t) in status.items():
    for (x, y) in {(p, q), (q, p)}:
        ax.add_patch(plt.Rectangle((x - 0.47, y - 0.47), 0.94, 0.94, color=c))
        ax.text(x, y, t, ha="center", va="center", fontsize=8)
ax.set_xlim(0.4, 6.6); ax.set_ylim(0.4, 6.6); ax.set_aspect("equal")
ax.set_xticks(range(1, 7)); ax.set_yticks(range(1, 7))
ax.set_xlabel("numbers on one side"); ax.set_ylabel("numbers on the other side")
ax.set_title("Fig 72 (DRAFT). Conjecture 21.2 by side sizes. Dark green: proved. Lighter: searched, nothing found\n"
             "(bound on the numbers; * = on the denominators d1 d2 of the two ratios). Sizes 5 and 6: the referee's search",
             fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig72-where-21-2-stands.png", dpi=130); plt.close(fig)
print("ok")
