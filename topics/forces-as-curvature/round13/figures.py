"""DRAFT picture for round 13. Run after checks.py. Blue = the law, orange = Regge read-back, grey = reference."""
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
a.semilogy(range(len(R["growth_72"])), R["growth_72"], "o-", color=BL, ms=8, label="kicks of 72°: closes at 120 corners")
a.semilogy(range(len(R["growth_36"])), R["growth_36"], "s-", color=GR, label="kicks of 36°: never closes (×11 per round)")
a.axhline(120, color=BL, ls=":", lw=1); a.text(5, 70, "120 = the 600-cell's corners", ha="right", color=BL, fontsize=8.5)
a.set_xlabel("rounds of kicks from rest"); a.set_ylabel("corners reached"); a.legend(fontsize=8, loc="upper left")
a.set_title("A. Grow from rest with 12 kicks toward\nan icosahedron's corners")

b = ax[1]
sc = np.array(R["scan"]); cap = sc[:, 1] >= 3000
b.semilogy(sc[cap, 0], np.full(cap.sum(), 3000), "v", color=GR, ms=5, label="passed 3000 corners: never closes")
b.semilogy(sc[~cap, 0], sc[~cap, 1], "o", color=BL, ms=10, label="closes")
for d, n in sc[~cap]:
    b.text(d - 3, n * 1.6, f"{int(d)}°: {int(n)}", ha="right" if d > 170 else "center", va="bottom", fontsize=9, color=BL)
b.set_ylim(0.3, 12000); b.set_xlim(0, 190)
b.set_xlabel("kick size (degrees on the 3-sphere)"); b.set_ylabel("corners reached")
b.legend(fontsize=8, loc="lower left")
b.set_title("B. Only two kick sizes close up:\n72° and 144° (180° just flips rest)")

c = ax[2]
rows = R["regge"]; lv = [r[0] for r in rows]; tgt = 6 * math.pi**2
c.semilogy(lv, [abs(r[2] / tgt - 1) for r in rows], "o-", color=OR, label="Regge's action vs Einstein's 6π²")
c.semilogy(lv, [abs(x - 1) for x in R["lambda_rho2"]], "s-", color=BL, label="with Λ: Λρ² → 1")
c.semilogy(lv, [abs(rows[0][2] / tgt - 1) / 4**l for l in lv], ":", color=GR, label="second order (÷4 per level)")
c.set_xticks(lv); c.set_xticklabels([f"{l}\n({r[1]:,} tets)" for l, r in zip(lv, rows)])
c.set_xlabel("halvings with the law's half kick"); c.set_ylabel("relative error"); c.legend(fontsize=8)
c.set_title("C. Glue the grown corners, read back:\nEinstein's equations, curvature +1")
fig.tight_layout(); fig.savefig("figures/f20_corners_grown_by_the_law.png", bbox_inches="tight")
print("f20 written")
