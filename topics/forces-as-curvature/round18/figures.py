"""DRAFT picture for round 18. Run after checks.py. Blue = the law, orange = the mirror or the other order, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 4, figsize=(23, 6.2), gridspec_kw=dict(width_ratios=[1.35, 1, 1, 1]))

a = ax[0]
labels = ["3D law vs\nh·y·h", "plain x·y vs\nh(x⊕y)h⁻¹", "2D law vs 3D\nlaw in a plane", "2D boost vs\nplaced point", "B(u)B(v) vs\nB(u⊕v)·turn", "for contrast:\nplain x·y vs law"]
vals = [R["sandwich"]["sandwich"], R["sandwich"]["group_vs_law"], R["one_sphere"]["plane"], R["one_sphere"]["point"],
        R["one_sphere"]["composition"], R["sandwich"]["plain_vs_law"]]
cols = [BL] * 5 + [GR]
a.bar(range(6), [max(v, 1e-17) for v in vals], color=cols)
for i, v in enumerate(vals):
    a.text(i, max(v, 1e-17) * 3, f"{v:.0e}" if v < 0.1 else f"{v:.2f}", ha="center", fontsize=8.5)
a.set_yscale("log"); a.set_ylim(1e-17, 30)
a.set_xticks(range(6)); a.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
a.set_ylabel("largest mismatch over 2,000 random pairs")
a.axhline(1e-9, color=GR, ls="--", lw=0.8); a.text(-0.4, 2e-9, "rounding level of the 2D formula", color=GR, fontsize=8, ha="left")
a.set_title("A. The law's 3D sum is the 3-sphere's own product,\nand its 2D boosts are points of that same 3-sphere")

b = ax[1]
Q = np.array(R["group"]["points_wxyz"])
order = np.lexsort((np.arctan2(Q[:, 2], Q[:, 1]), -Q[:, 0])); Q = Q[order]
def qmul(a_, b_):
    a0, av = a_[..., :1], a_[..., 1:]; b0, bv = b_[..., :1], b_[..., 1:]
    return np.concatenate([a0 * b0 - (av * bv).sum(-1, keepdims=True), a0 * bv + b0 * av + np.cross(av, bv)], -1)
P = qmul(Q[:, None, :], Q[None, :, :])
im = b.imshow(P[..., 0], cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
b.set_xlabel("corner y (sorted from rest outward)"); b.set_ylabel("corner x")
cb = fig.colorbar(im, ax=b, fraction=0.046, pad=0.03); cb.set_label("where x·y lands: 1 = rest, −1 = far point", fontsize=8.5)
b.text(0.02, -0.2, f"all 14,400 products land on corners (largest miss {R['group']['closure']:.0e});\n"
       f"{R['group']['corners']} elements turn space by {R['group']['rotations']} rotations: the binary icosahedral group",
       transform=b.transAxes, fontsize=8.5, color=BL)
b.set_title("B. The 120 corners the law grew (round 13)\nform a group: its multiplication table")

c = ax[2]
L = R["transport"]
ha = np.array([r["half_area"] for r in L])
xx = np.linspace(0, 5.5, 50); c.plot(xx, xx, color=GR, ls="--", lw=1, label="half the enclosed area")
c.plot(ha, [r["law_phase"] for r in L], "o", color=BL, ms=11, mfc="none", label="phase the law's transport collects")
for r in L:
    c.annotate(f"cap {r['cap']}", (r["half_area"], r["law_phase"]), textcoords="offset points", xytext=(8, -12), fontsize=8.5, color=BL)
c.set_xlabel("half the area a loop of velocities encloses"); c.set_ylabel("phase round the loop")
c.legend(fontsize=8.5, loc="upper left")
c.set_title("C. Carrying a frame round a loop by composing\nthe law's boosts collects half the area (round 10)")

d = ax[3]
cases = [("the law's transport\nB(u)B(v)", R["transport"][1]), ("seen through\nthe parity x → −x", R["mirror"]["through_parity"][0]),
         ("same boosts,\nother order", R["mirror"]["other_order"][0])]
w = 0.36
for i, (name, r) in enumerate(cases):
    d.bar(i - w / 2, max(abs(r["along_hopf"]), 1e-17), w, color=OR, label="along round 17's 'Hopf' field (curl −2)" if i == 0 else None)
    d.bar(i + w / 2, max(abs(r["along_mirror"]), 1e-17), w, color=BL, label="along round 17's 'mirror' field (curl +2)" if i == 0 else None)
d.set_yscale("log"); d.set_ylim(1e-17, 50)
d.set_xticks(range(3)); d.set_xticklabels([n for n, _ in cases], fontsize=8.5)
d.set_ylabel("how far the frame moves along each field (cap 0.9)")
d.legend(fontsize=8.5, loc="center", title="short bars are rounding:\nthe frame never moves along that field", title_fontsize=8)
g6 = R["law_field_on_pieces"]
d.text(0.02, -0.2, "the law's own field on the pieces, curl: " + ", ".join(f"{r['curl']:+.3f}" for r in g6) + " (levels 0-2) → +2;\n"
       "the law moves along exactly one field; a mirror or the other order picks the other",
       transform=d.transAxes, fontsize=8.5, color=BL)
d.set_title("D. The law's transport follows one Hopf field\n(round 17's 'mirror'); the hand is a choice")

fig.tight_layout()
fig.savefig("figures/f25_one_three_sphere.png", bbox_inches="tight")
