"""DRAFT picture for round 16. Run after checks.py. Blue = the law, orange = Regge read-back, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 4, figsize=(21, 5.4))

a = ax[0]
styles = {"1.0": (BL, "o-", 2.5), "0.9": (GR, "s--", 1), "1.1": (GR, "^--", 1), "0.25": (GR, "v:", 1), "4.0": (GR, "D:", 1),
          "0.0": (OR, "x-.", 1), "-1.0": (OR, "+-.", 1)}
for k, cnt in R["closure"].items():
    col, st, lw = styles[k]; lab = {"1.0": "q = 1: closes (120 corners)", "0.0": "q = 0 (flat)", "-1.0": "q = −1 (open)"}.get(k, f"q = {float(k):g}")
    a.semilogy(range(len(cnt)), cnt, st, color=col, lw=lw, ms=6, label=lab)
a.axhline(120, color=BL, lw=0.8, alpha=0.4)
a.set_xlabel("generations of kicks (kick size held fixed)"); a.set_ylabel("distinct corners reached")
a.legend(fontsize=8, loc="upper left"); a.set_ylim(0.8, 8000)
a.set_title("A. With the kick fixed, the law closes space\nat one dial only, and never at q ≤ 0")

b = ax[1]
lv = np.arange(4); mk = {"0.25": ("v", 9), "1.0": ("o", 6), "4.0": ("D", 4)}
for k, vals in R["lambda_star"].items():
    m, s = mk[k]; b.plot(lv, np.array(vals) / float(k), m + "-", color=OR, ms=s, label=f"q = {float(k):g}" if k != "1.0" else "q = 1")
b.axhline(3, color=GR, ls="--", lw=1); b.text(0.1, 2.96, "Λ = 3q (the time-symmetric constraint)", color=GR, fontsize=8.5)
for L, v in zip(lv, R["lambda_star"]["1.0"]):
    b.annotate(f"{v:.3f}", (L, v), textcoords="offset points", xytext=(-6 if L < 3 else -34, 9), fontsize=8.5, color=OR)
b.set_xticks(lv); b.set_xticklabels([f"{L}\n({600 * 8**L:,} pieces)" for L in lv])
b.set_xlabel("refinement level of the law's slice"); b.set_ylabel("Λ / q that lets a slab rest on the slice")
b.legend(fontsize=8.5, title="the three dials overlap", title_fontsize=8.5); b.set_ylim(2.9, 3.5)
b.set_title("B. At a moment of rest, Regge fixes Λ\nfrom the law's slice: Λ / q → 3")

c = ax[2]
r2 = R["grown_runs"]["q=1,L=2"]; tt = np.linspace(0, 0.65, 200)
c.plot(tt, r2["throat"] * np.cosh(r2["kappa"] * (tt - r2["bounce_t"])), color=GR, lw=1.5,
       label=f"forward + inverse run of one step, k = {r2['kappa']:.4f}√q")
mk2 = {"q=0.25,L=2": ("v", 9, "q = ¼"), "q=1,L=2": ("o", 6, "q = 1"), "q=4,L=2": ("D", 4, "q = 4")}
for key, (m, s, lab) in mk2.items():
    r = R["grown_runs"][key]; q = r["q"]; A = np.array(r["a"]); t = 0.025 / math.sqrt(q) * np.arange(len(A))
    c.plot(t * math.sqrt(q), A * math.sqrt(q), m, color=OR, ms=s, mfc="none" if m != "D" else OR, label=f"glued slabs at {lab} (level 2)")
r3 = R["grown_runs"]["q=1,L=3"]
c.text(0.40, 1.003, f"level 3: k/√q = {r3['k_over_sqrtq']:.4f}\nsize at the bounce × √q = {r3['size_sqrtq']:.4f}", fontsize=8.5, color=BL)
c.set_xlabel("time × √q"); c.set_ylabel("size of the slice × √q")
c.legend(fontsize=8, loc="upper left")
c.text(0.40, 1.035, "the run balances half a step\nbefore the law's slice (Q4, Q5)", fontsize=8.5, color=GR)
c.set_title("C. Grown from the law's slice at rest, every dial\nfalls on one curve: the shift's k is √q")

d = ax[3]; d.set_xlim(0, 2); d.set_ylim(0, 2); d.axis("off")
cells = {(0, 1): ("cosh: de Sitter\nexpands, Λ = +3q", True), (1, 1): ("cos: a 4-sphere\nΛ = +3q", True),
         (0, 0): ("cos: anti-de Sitter\nΛ = 3q < 0", False), (1, 0): ("cosh: hyperbolic 4-space\nΛ = 3q < 0", False)}
for (i, j), (txt, ok) in cells.items():
    d.add_patch(plt.Rectangle((i, j), 1, 1, fc=(BL if ok else GR), alpha=0.18 if ok else 0.08, ec="k", lw=0.8))
    d.text(i + 0.5, j + 0.55, txt, ha="center", va="center", fontsize=9.5)
    d.text(i + 0.5, j + 0.15, "the law closes space" if ok else "the law builds no closed space", ha="center", fontsize=8, color=BL if ok else GR)
d.text(0.5, 2.05, "boosts: q_v = −1\n(Lorentzian time)", ha="center", va="bottom", fontsize=9)
d.text(1.5, 2.05, "rotations: q_v = +1\n(Euclidean time)", ha="center", va="bottom", fontsize=9)
d.text(-0.05, 1.5, "space\nq_s > 0", ha="right", va="center", fontsize=9); d.text(-0.05, 0.5, "space\nq_s < 0", ha="right", va="center", fontsize=9)
d.text(1, -0.12, "k² = −q_v q_s: the space dial sets the size, the velocity dial the sign", ha="center", fontsize=8.5, color=GR)
d.set_title("D. The sign of time is not in q:\nit is which end of the dial time's kicks use", pad=38)

fig.tight_layout()
fig.savefig("figures/f23_does_q_set_k.png", bbox_inches="tight")
