"""DRAFT figure D3: orbit turning across the magic square. Run after checks_magic.py."""
import os, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "..", "figures")
KNOWN, OURS, GREY = "#d9822b", "#2b6cb0", "#bbbbbb"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
fig, axs = plt.subplots(1, 3, figsize=(20, 6.4), gridspec_kw={"width_ratios": [1.1, 1.15, 0.95]})
ax = axs[0]
ns = [1, 2, 3, 4]; bps = [n-1 for n in ns]; neu = [2*(n-1) for n in ns]; w = 0.38
ax.bar([n-w/2 for n in ns], bps, w, color=OURS, label="BPS bodies, opposite charges: n − 1")
ax.bar([n+w/2 for n in ns], neu, w, color="#7aa6d6", label="neutral masses: 2(n − 1)")
for n, b, c in zip(ns, bps, neu):
    ax.text(n-w/2, b+0.1, str(b), ha="center", fontweight="bold"); ax.text(n+w/2, c+0.1, str(c), ha="center", fontweight="bold")
ax.annotate("N = 8 half-BPS:\nno turn\n(Caron-Huot–Zahraee)", xy=(1-w/2, 0.05), xytext=(0.55, 2.6), color=KNOWN, fontsize=9.5, arrowprops=dict(arrowstyle="->", color=KNOWN))
ax.annotate("double copy of a\ncolour charge: 2", xy=(2+w/2+0.12, 2.3), xytext=(1.75, 4.3), fontsize=9.5, arrowprops=dict(arrowstyle="->"))
ax.annotate("plain GR: 6", xy=(4+w/2, 6.0), xytext=(3.0, 6.3), color=KNOWN, fontsize=10, arrowprops=dict(arrowstyle="->", color=KNOWN))
ax.set_xticks(ns); ax.set_xticklabels(["n = 1\na = √3", "n = 2\na = 1", "n = 3\na = 1/√3", "n = 4\na = 0"])
ax.set_xlabel("charges the body carries, n = 4/(1 + a²)   (a = how hard the dilaton pulls)")
ax.set_ylabel("turn per lap ÷ (π K²/L²)"); ax.set_ylim(0, 7.3); ax.legend(loc="upper left", fontsize=9.5)
ax.set_title("A. The turn counts the body's charges, not the algebras")

ax = axs[1]; ax.axis("off"); ax.set_xlim(-0.9, 3.1); ax.set_ylim(-0.6, 3.5)
names = ["C", "H", "O"]; Ns = {(0, 0): 2, (0, 1): 3, (0, 2): 5, (1, 1): 4, (1, 2): 6, (2, 2): 8}
notes = {(0, 0): "graviphoton only\n→ RN bodies\nBPS turn 3", (1, 1): "N = 4 half-BPS:\na = √3 or 1\n→ turn 0 or 1",
         (2, 2): "all n = 1…4\n→ turn 0 to 3\nhalf-BPS: 0"}
for i in range(3):
    for j in range(3):
        key = (min(i, j), max(i, j)); N = Ns[key]; col = OURS if key in notes else GREY
        y = 2 - i
        ax.add_patch(plt.Rectangle((j, y), 0.95, 0.95, color=col, alpha=0.85 if key in notes else 0.5))
        ax.text(j+0.475, y+0.78, f"{names[i]}⊗{names[j]}:  N = {N}", ha="center", color="white" if key in notes else "black", fontsize=9.5, fontweight="bold")
        ax.text(j+0.475, y+0.38, notes.get(key, "not worked\nout here"), ha="center", va="center", color="white" if key in notes else "black", fontsize=8.8)
for j in range(3): ax.text(j+0.475, 3.05, names[j] + "\nSYM", ha="center", fontsize=10)
for i in range(3): ax.text(-0.1, 2-i+0.475, names[i] + " SYM", ha="right", va="center", fontsize=10)
ax.text(-0.85, -0.45, "Every cell also holds the neutral double-copy masses (graviton + dilaton): turn 2 in all nine.", fontsize=9.5)
ax.set_title("B. The 4-dimensional magic square (gravity = gauge theory²)")

ax = axs[2]
pairs = [("Maxwell →\nplain GR", 1, 6), ("colour charge →\ndouble-copy mass", 1, 2), ("N = 4 BPS →\nN = 8 BPS", 0, 0)]
for k, (lab, s, d) in enumerate(pairs):
    ax.bar(k-0.18, s, 0.36, color=KNOWN if k != 1 else "#e8b07a"); ax.bar(k+0.18, d, 0.36, color=OURS)
    ax.text(k-0.18, s+0.1, str(s), ha="center", fontweight="bold"); ax.text(k+0.18, d+0.1, str(d), ha="center", fontweight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels([p[0] for p in pairs], fontsize=9.5); ax.set_ylim(0, 7.6)
ax.text(1.0, 7.1, "left: single copy (gauge theory)   right: double copy (gravity)", ha="center", fontsize=9)
ax.set_title("C. Single copy and double copy, side by side\n(single copy: a vector adds +1, a scalar −1)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figD3-orbit-turning-across-the-magic-square.png"), dpi=130); print("ok")
