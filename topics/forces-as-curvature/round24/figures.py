"""DRAFT picture for round 24. Run after checks.py. Blue = the law's mesh, orange = charges and their light, grey = smooth reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(2, 2, figsize=(16, 13)); ax = ax.ravel()

a = ax[0]
P = R["profiles"]["pair"]; B = R["profiles"]["background"]
ch = np.array(P["chi"]); ph = np.array(P["phi"]); s = (ch > 0.25) & (ch < math.pi - 0.25)
a.plot(ch[s], ph[s], ".", color=BL, ms=2.5, alpha=0.6, label="+ and − pair at the law's two poles: potential at each corner (level 2)")
c = np.linspace(0.25, math.pi - 0.25, 400); a.plot(c, 1 / np.tan(c) / (4 * math.pi), color=GR, lw=1.6, label="smooth: cot χ / (4π)")
chb = np.array(B["chi"]); phb = np.array(B["phi"]); sb = chb > 0.25
a.plot(chb[sb], phb[sb] - 0.35, ".", color=OR, ms=2.5, alpha=0.6, label="lone + charge with an even − background (shifted down 0.35)")
E_bg = lambda x: (1 - (x - math.sin(x) * math.cos(x)) / math.pi) / (4 * math.pi * math.sin(x)**2)
from scipy.integrate import quad
ref = np.array([quad(E_bg, x, math.pi / 2)[0] for x in c]); a.plot(c, ref - 0.35, color=GR, lw=1.6, ls="--", label="smooth, with background")
a.set_xlabel("angle χ from the + charge on the law's 3-sphere"); a.set_ylabel("electric potential")
a.set_xticks([0, math.pi / 2, math.pi]); a.set_xticklabels(["0 (+ charge)", "π/2", "π (antipode)"])
a.legend(fontsize=8.5, loc="upper right"); a.set_ylim(-0.6, 0.35)
a.set_title("A. Static charges on the law's mesh: Coulomb on a closed space (Q2, Q3)")

b = ax[1]
cc = np.linspace(0.001, math.pi - 0.001, 400)
b.plot(cc, np.ones_like(cc), color=BL, lw=2.4, label="pair: field × area of the sphere at χ = 1 all the way (the − charge takes it back)")
b.plot(cc, 1 - (cc - np.sin(cc) * np.cos(cc)) / np.pi, color=OR, lw=2.4, label="lone charge + background: the enclosed charge falls to 0 at the antipode")
for rw in R["pair"]:
    b.plot([0.8, 1.4, 2.0], rw["enclosed"], "o", color=BL, ms=5 + 2 * rw["level"], mfc="none")
b.text(1.62, 0.88, "on the mesh: enclosed charge = 1 to 12 digits\nfor every cap (levels 0, 1, 2)", fontsize=9, color=BL)
b.text(0.1, 0.25, "a lone charge without a background has no field at all:\nthe flux has nowhere to end on a closed space (Q1)", fontsize=9, color=GR)
b.set_ylim(0, 1.3); b.set_xticks([0, math.pi / 2, math.pi]); b.set_xticklabels(["0", "π/2", "π"])
b.set_xlabel("angle χ of the sphere around the + charge"); b.set_ylabel("charge enclosed (= field × area)")
b.legend(fontsize=8.5, loc="upper right", bbox_to_anchor=(1.0, 0.93)); b.set_title("B. Gauss's law on the law's closed space: total charge must be 0")

c_ = ax[2]
runs = R["moving"]["runs"]
for rw, al in zip(runs, (1.0, 0.75, 0.5, 0.3)):
    t = np.array(rw["t"]); E = np.array(rw["E_rad"])
    c_.plot(t, E, color=OR, lw=2, alpha=al, label=f"move takes {rw['T']}: light left behind {rw['E_rad_after']:.3g}")
    c_.axvline(rw["T"], color=OR, lw=0.6, ls=":", alpha=al)
c_.set_yscale("log"); c_.set_ylim(1e-3, 50)
c_.set_xlabel("time"); c_.set_ylabel("energy in the gradient-free part of the field (light)")
c_.text(0.98, 0.04, f"Gauss's law exact throughout (worst {max(r['gauss'] for r in runs):.0e});\nafter the stop the light's energy is kept to 1e-9", transform=c_.transAxes, ha="right", fontsize=9, color=GR)
c_.legend(fontsize=8.5, loc="upper right"); c_.set_title(f"C. A charge walking {R['moving']['path_edges']} edges of the law's mesh makes light (Q4)")

d = ax[3]
rows = R["charge_half"]["levels"]
x = np.arange(3); w = 0.35
d.bar(x - w / 2, [max(r["magnetic_source"], 1e-17) for r in rows], w, color=BL, label="magnetic source of the law's own field (flux out of any piece)")
d.bar(x + w / 2, [max(r["electric_source"], 1e-17) for r in rows], w, color=BL, alpha=0.5, label="electric source of the law's own field (divergence)")
d.set_yscale("log"); d.set_ylim(1e-18, 1)
d.set_xticks(x); d.set_xticklabels([f"mesh level {i}" for i in x])
d.axhline(R["neutral"][2]["single_charge_residual"], color=OR, lw=1.5, ls="--", label="for scale: the smallest miss when forcing one lone unit charge (level 2)")
d.text(1.0, 1e-7, "round 10's charge ½ is not a charge in the law's space:\nthe law's field has no source here. Its 2π of flux, one Dirac unit,\nis on the sphere of directions (circulation round a fibre = 2π)", ha="center", fontsize=9, color=GR)
d.set_ylabel("size of the source (log)"); d.legend(fontsize=8.5, loc="upper left")
d.set_title("D. The law's own field is sourceless in its space (Q5)")

fig.suptitle("DRAFT round 24: charge as a source of light on the law's mesh", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("figures/f31_charge_on_the_law.png", bbox_inches="tight")
