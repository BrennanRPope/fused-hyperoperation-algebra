"""Figure 22: the Lee-Yang proofs refactored into one core (run after general_p_far.py)."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
D = json.load(open("general_p_far.json"))
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 11.5))
gs = fig.add_gridspec(2, 2, wspace=0.24, hspace=0.36)
# A: the refactor as a call graph
ax = fig.add_subplot(gs[0, :]); ax.set_xlim(0, 16); ax.set_ylim(0.2, 5.1); ax.axis("off")
def box(x, y, w, h, text, fc, ec, fs=10.5, color="k"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=1.4))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=color)
def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", lw=1.5, color="#444"))
box(0.2, 3.3, 4.3, 1.6, "Toolkit (cited, shared by every p)\nPólya: base a ≤ 1\nKrein: climb b → b + 1\nbreak past a = p", "#f3f3f3", GREY)
box(0.2, 0.4, 4.3, 1.9, "AlternatesBelow p\n(the one per-p input)\nG at b + 1 alternates at the\nzeros of G at b, 0 < b ≤ p − 1", "#fde8d7", ORANGE)
box(5.6, 1.6, 4.6, 2.2, "Core theorem (proved in Lean)\nholds_of_alternatesBelow:\nLee–Yang for every a ≤ p\n+ break: the rule as an iff", "#dbe7f5", BLUE)
arrow(4.5, 4.0, 5.6, 3.2); arrow(4.5, 1.3, 5.6, 2.2)
inst = [("p = 1: nothing needed", "#d7eedf", GREEN), ("p = 2: Bessel interlacing (cited)", "#d7eedf", GREEN),
        ("p = 3: computer certificate", "#d7eedf", GREEN), ("every p: open; the far zone is\nuniform, the near zone isn't yet", "#fbe3e3", RED)]
for k, (t, fc, ec) in enumerate(inst):
    y = 4.15 - 1.2*k
    box(11.3, y, 4.4, 0.85 if k < 3 else 1.0, t, fc, ec, fs=10.5)
    arrow(10.2, 2.7, 11.3, y + 0.42)
ax.text(8.0, 0.55, "p = 4 at whole a used a different tool (Laguerre's operator),\nso it isn't an instance of this core yet.",
        ha="center", fontsize=9.5, color=GREY)
ax.set_title("The refactor: one shared core, one input per roundness", fontsize=13)
# B: takeover and near-zone size vs p
ax = fig.add_subplot(gs[1, 0])
sc = D["scan"]; ps = sorted(float(p) for p in sc)
T = [max(sc[str(p)]["takeover"], 2.0) for p in ps]   # 2 = smallest s sampled
nz = [D["near_zone_zeros"][str(p)] for p in ps]
ax.plot(ps, T, "o-", color=BLUE, label="far zone takes over (exact pieces)")
ax.plot([3], [21], "s", color=PURPLE, ms=8, label="p = 3 certificate's bound (s = 21)")
for p, t, n in zip(ps, T, nz):
    ax.annotate(str(n), (p, t), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=10, color=ORANGE, fontweight="bold")
ax.text(1.1, 13, "orange numbers: zeros of G at b that\ncome before the takeover (the part\nstill boxed by hand-built checks)", fontsize=9, color=ORANGE)
ax.set_xlabel("roundness p"); ax.set_ylabel("s where all three far-zone tests hold from then on")
ax.set_ylim(-1, 30); ax.legend(fontsize=9, loc="upper left")
ax.set_title("Where the far-out argument takes over, by roundness\n(numerical, not certified)")
# C: far margin vs p
ax = fig.add_subplot(gs[1, 1])
pp = np.linspace(1.001, 8, 800); m = 1 - np.abs(np.sin(np.pi*pp/2))/pp**(pp - 1)
ax.plot(pp, m, color=GREEN, label="kink margin 1 − |sin(πp/2)| / p^(p−1)")
ax.plot(ps, [sc[str(p)]["far_margin"] for p in ps], "o", color=ORANGE, label="measured tight margin, s > 200")
for p in [2, 4, 6]: ax.axvline(p, color=GREY, lw=0.7, ls=":")
ax.text(2.05, 0.3, "dotted: even p,\nno kink, margin 1", fontsize=9, color=GREY)
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("roundness p"); ax.set_ylabel("margin far out (top step b = p − 1)"); ax.set_ylim(-0.05, 1.1)
ax.legend(fontsize=9, loc="center right")
ax.set_title("Far out, the tight test has the same margin at every p:\nthe kink share, positive for every p > 1")
fig.suptitle("Figure 22. Lee–Yang refactored: what's shared, what's per roundness (DRAFT)", fontsize=14, y=0.995)
fig.savefig("figures/fig22-lee-yang-refactor.png", dpi=110, bbox_inches="tight"); print("saved")
