"""DRAFT picture for round 12. Run after checks.py. Blue = the law, orange = Regge read-back, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(16, 4.9))

a = ax[0]
for key, mk, col, lab in (("1.0", "o", OR, "q = 1: action vs 6π²"), ("0.25", "x", BL, "q = 1/4: action vs 12π² (radius 2)")):
    rows = R["sphere_from_law"][key]; tgt = 6 * math.pi**2 / math.sqrt(float(key))
    a.semilogy([r[0] for r in rows], [abs(r[2] / tgt - 1) for r in rows], mk + "-", color=col, label=lab, ms=8)
lv = [r[0] for r in R["sphere_from_law"]["1.0"]]
a.semilogy(lv, [abs(x - 1) for x in R["lambda_rho2"]], "s--", color=GR, label="with Λ: Λρ² → 1")
a.set_xticks(lv); a.set_xticklabels([f"{l}\n({r[1]:,} tets)" for l, r in zip(lv, R["sphere_from_law"]["1.0"])])
a.set_xlabel("halvings with the law's half kick"); a.set_ylabel("relative error"); a.legend(fontsize=8)
a.set_title("A. A 3-sphere sized only by the law:\nRegge's action tends to Einstein's for each q")

b = ax[1]
cols = {"1.0": "#0b3d91", "0.25": BL, "0.0": GR, "-0.25": "#f0a040", "-1.0": OR}
for key, rows in R["star"].items():
    e2 = [r[0]**2 for r in rows]; d = [r[1] for r in rows]
    b.plot(e2, d, "o-", color=cols[key], label=f"q = {float(key):+.2f}")
b.axhline(0, color=GR, lw=0.8)
b.set_xlabel("size² of the 20 law-built tetrahedra round a corner"); b.set_ylabel("deficit on an edge from the centre (rad)")
b.legend(fontsize=8, title="dial", title_fontsize=8)
b.set_title("B. Round one corner: gaps when q > 0,\noverlaps when q < 0, none at q = 0")

c = ax[2]
qs = [float(k) for k in R["star"]]; slope = [R["star"][k][-1][1] / R["star"][k][-1][0]**2 for k in R["star"]]
qq = np.linspace(-1.1, 1.1, 50); k0 = np.mean(list(R["star_constant"].values()))
c.plot(qq, k0 * qq, "-", color=GR, lw=1, label=f"deficit / size² = {k0:.4f} × q")
c.plot(qs, slope, "o", color=OR, ms=9, label="read back by Regge from the law's lengths")
c.axhline(0, color=GR, lw=0.6); c.axvline(0, color=GR, lw=0.6)
c.set_xlabel("the law's dial q"); c.set_ylabel("deficit / size²"); c.legend(fontsize=8, loc="upper left")
c.set_title("C. The curvature read back is the dial:\nstraight line through zero, sphere to saddle")
fig.tight_layout(); fig.savefig("figures/f19_pieces_sized_by_the_law.png", bbox_inches="tight")
print("f19 written")
