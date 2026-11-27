"""DRAFT picture for round 27. Run after checks.py. Blue = the law, red = ruled out, green = agrees, amber = open, grey = reference."""
import json, math, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BL, RD, GN, AM, GR = "#1f6fc2", "#c0392b", "#2e8b57", "#e09b00", "#777777"
R = json.load(open("results.json")); D = R["data"]; FG = R["figure"]
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig = plt.figure(figsize=(16, 13.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.08])
ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

# A. curvature: the law's two versions against data
a = ax[0]
a.axvspan(-0.012, 0, color=AM, alpha=0.13, label="closed slices of any size: the law's sign-only version (chosen, open)")
a.axvline(0, color=GR, lw=0.8)
for y, (nm, (m, s)) in zip((0.62, 0.38), (("DESI DR2 BAO + CMB", D["OK_DESI"]), ("Planck 2018 + BAO", D["OK_PlanckBAO"]))):
    a.errorbar([m], [y], xerr=[[2 * s], [2 * s]], fmt="none", ecolor=GR, elinewidth=1.2, capsize=4)
    a.errorbar([m], [y], xerr=[s], fmt="o", color="k", ms=7, elinewidth=2.6, capsize=5)
    a.text(m, y + 0.07, f"{nm}: {m} ± {s}", ha="center", fontsize=9.5)
a.annotate(f"the law read literally (q = Λ/3 sets today's curvature):\nΩ_k = −Ω_Λ = {R['G1']['Ok_naive']}, far off the left edge,\n{R['G1']['sigma_DESI']:.0f}σ from DESI, {R['G1']['sigma_PlanckBAO']:.0f}σ from Planck: ruled out",
           (-0.0118, 0.15), (-0.0105, 0.1), fontsize=9.5, color=RD, arrowprops=dict(arrowstyle="<-", color=RD, lw=2))
a.text(0.0062, 0.86, f"DESI + CMB lean open:\nP(closed) = {R['G3']['P_closed_DESI']:.3f}", fontsize=9.5, color=AM, ha="center")
a.set_xlim(-0.012, 0.012); a.set_ylim(0, 1); a.set_yticks([])
a.set_xlabel("Ω_k (spatial curvature today; negative = closed, positive = open). Thick bars 1σ, thin 2σ")
a.legend(fontsize=8.8, loc="upper left")
a.set_title("A. The one sign the law leans toward, and the data (G1, G3)")

# B. no moment of rest in our past
b = ax[1]
z = np.array(FG["g2"]["z"])
lab = {str(D["OK_DESI"][0]): ("measured, Ω_k = 0.0023", "k", 2.6), "-0.3": ("Ω_k = −0.3", GR, 1.6),
       str(round(R["G2"]["Ok_threshold"], 4)): (f"threshold Ω_k = {R['G2']['Ok_threshold']:.2f}: at rest once, at 1 + z = {R['G2']['rest_at_1pz']:.2f}", BL, 2.4),
       "-1.2": ("Ω_k = −1.2: a bounce (the law's de Sitter shape)", BL, 1.4)}
for k, v in FG["g2"]["E2"].items():
    nm, col, lw = lab[k]; ls = ":" if k == "-1.2" else "-"
    b.plot(z, v, color=col, lw=lw, ls=ls, label=nm)
b.axhline(0, color=RD, lw=1)
b.text(3.2, 0.6, "H² = 0 is the law's moment of rest,\nwhere rounds 16, 23, 25 tie q to the slices", fontsize=9.5, color=RD)
b.set_ylim(-1.5, 12); b.set_xlim(0, 6)
b.set_xlabel("redshift z (looking back in time)"); b.set_ylabel("H²/H₀² (expansion rate squared)")
b.legend(fontsize=8.8, loc="upper left")
b.set_title("B. Our universe never stood still: its past is a big bang, not the law's bounce (G2)")

# C. the ledger
c = ax[2]; c.axis("off")
rows = R["ledger"]
col = lambda v: RD if v.startswith("ruled") else (AM if v.startswith("open") else (GN if v.startswith("agrees") else GR))
W = lambda t, n: "\n".join(textwrap.wrap(t.replace("Lambda", "Λ").replace("Omega_k", "Ω_k").replace("gamma", "γ"), n))
cells = [[W(r["c"], 34), W(r["forced"], 20), W(r["verdict"], 34)] for r in rows]
tb = c.table(cellText=cells, colLabels=["candidate", "does the law force it?", "verdict"], loc="center", cellLoc="left", colWidths=[0.40, 0.24, 0.40], bbox=[0.0, 0.0, 1.0, 0.97])
tb.auto_set_font_size(False); tb.set_fontsize(8.8)
for i, r in enumerate(rows):
    tb[(i + 1, 2)].set_facecolor(col(r["verdict"])); tb[(i + 1, 2)].set_alpha(0.28)
for j in range(3): tb[(0, j)].set_facecolor("#dddddd")
c.set_title("C. Ten candidates: one open and testable; the one forced and new is ruled out (G6)", pad=4)

# D. the scales
d = ax[3]
RH = 1.0
items = [("weak force's reach (Higgs screening)", 2.5e-18 / 1.37e26, RD),
         ("Hubble radius c/H₀", 1.0, GR),
         ("edge of the observable universe", 3.3, GR),
         (f"closed universe's curvature radius, at least ({R['G4']['a0_min_hubble']:.0f}× c/H₀)", R["G4"]["a0_min_hubble"], AM),
         (f"law's lowest wave, at least ({R['G5']['lambda_min_hubble']:.0f}× c/H₀)", R["G5"]["lambda_min_hubble"], BL)]
for i, (nm, x, cl) in enumerate(items):
    d.plot([x], [i], "o", color=cl, ms=10)
    if x < 1e-10: d.text(x * 30, i, nm, va="center", fontsize=9.5, color=cl)
    else: d.text(x / 30, i, nm, va="center", ha="right", fontsize=9.5, color=cl)
d.axvspan(1e-46, 3.3, color=GN, alpha=0.06); d.text(1e-30, 4.55, "what can be observed", color=GN, fontsize=9.5)
d.set_xscale("log"); d.set_xlim(1e-46, 1e4); d.set_ylim(-0.6, 4.9); d.set_yticks([])
d.set_xlabel("length, in Hubble radii (log scale)")
d.text(1e-44, -0.4, f"the law's f = 1 twist filling such a universe: Ω ≤ {R['G4']['Omega_YM']:.0e} of the critical density (G4)", fontsize=9.5, color=BL)
d.set_title("D. What the law does fix sits beyond reach (G4, G5)")

fig.suptitle("Round 27 (DRAFT): hunting for a gravity prediction from the law", fontsize=14, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig("figures/f34_prediction_hunt.png")
print("wrote figures/f34_prediction_hunt.png")
