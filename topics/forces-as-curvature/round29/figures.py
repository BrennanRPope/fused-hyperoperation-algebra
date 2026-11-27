"""DRAFT round 29 picture: f36_dial_slope.png. Reads results.json (this folder) and ../dial-potential/results.json."""
import json, math, os, sys, textwrap
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "dial-potential"))
import quintessence as Q
R = json.load(open(os.path.join(HERE, "results.json")))
P = json.load(open(os.path.join(HERE, "..", "..", "dial-potential", "results.json")))["slope_profile"]
GREY, BLUE, GREEN, RED, ORANGE = "#777777", "#1f6fc5", "#2e8b57", "#c8372d", "#e08a1e"

fig = plt.figure(figsize=(20, 12.5))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.32, wspace=0.25)
fig.suptitle("Round 29 (DRAFT): does the law fix the slope of the dial as rolling dark energy? (audit C6)", fontsize=16)

# A: chi2 profile with band and slopes
ax = fig.add_subplot(gs[0, 0])
lam = np.array(P["lam"]); c = np.array(P["chi2"]); c0 = c[0]
ax.plot(lam, c - c0, color="k", lw=2, label="data: Δχ² vs a cosmological constant")
lo, hi = R["D4"]["band"]
ax.axvspan(lo, hi, color=GREY, alpha=0.18, label="data's 1σ band, 0.61 to 0.97 (fitted)")
ax.axhline(0, color=GREY, lw=1, ls=":")
marks = [(math.sqrt(2 / 3), GREEN, "√(2/3): chosen"),
         (math.sqrt(2), RED, "√2: slice curvature"),
         (math.sqrt(10 / 3), RED, "√(10/3): hidden sphere"),
         (7 / math.sqrt(7.5), RED, "2.56: meron on\nhidden sphere")]
for i, (v, col, lab) in enumerate(marks):
    ax.axvline(v, color=col, lw=2, ls="--")
    ax.text(v - 0.1 if abs(v - math.sqrt(2)) < 1e-9 else v + 0.035, 38, lab.replace("\n", " "), color=col, fontsize=10, va="center", rotation=90,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
ax.set_xlim(0, 2.7); ax.set_ylim(-8, 60)
ax.set_xlabel("slope λ of V ∝ exp(−λφ/M_Pl)"); ax.set_ylabel("Δχ² (Pantheon+ and DESI DR2 BAO)")
ax.set_title("A. Every slope the law forces misses the data (D4)")
ax.legend(loc="upper left", fontsize=9)

# B: w(z)
ax = fig.add_subplot(gs[0, 1])
z = Q.ZG
for lv, col, lab in [(0.0, GREY, "cosmological constant, w = −1"), (0.8177, GREEN, "best fit λ = 0.82 (≈ √(2/3))"),
                     (math.sqrt(2), RED, "λ = √2 (slice curvature)")]:
    r = next(rr for rr in R["D4"]["rows"] if abs(rr["lam"] - lv) < 1e-3)
    E, w = Q.history(r["Om"], "exp", lv)
    ax.plot(z, w, color=col, lw=2.2, label=lab, ls="-" if lv != math.sqrt(2) else "--")
ax.axhline(-1 / 3, color=ORANGE, lw=1.2, ls=":")
ax.text(2.0, -0.31, "w = −1/3: no acceleration above this", color=ORANGE, fontsize=9, ha="right", va="bottom")
ax.set_xlim(0, 2.5); ax.set_ylim(-1.05, -0.2)
ax.set_xlabel("redshift z"); ax.set_ylabel("dial's equation of state w(z)")
ax.set_title("B. The fitted dial thaws a little (D5)")
ax.legend(loc="center right", fontsize=9)
b = R["D5"]
ax.text(1.0, -0.55, f"λ = 0.82: w0 = {b['w0']:.2f}, wa = {b['wa']:.2f}", fontsize=10, color=GREEN)

# C: hidden-sphere potential
ax = fig.add_subplot(gs[0, 2])
s = np.linspace(-0.25, 1.6, 400); Ac, Bf = 1.0, 1.0
Vc, Vf = -Ac * np.exp(-5 * s), Bf * np.exp(-7 * s)
ax.plot(s, Vc, color=BLUE, lw=2, label="sphere's own curvature: −e^(−5ψ)  (D1)")
ax.plot(s, Vf, color=ORANGE, lw=2, label="law's meron on the sphere: +e^(−7ψ)  (D2)")
ax.plot(s, Vc + Vf, color="k", lw=2.5, label="together")
smin = math.log(7 / 5) / 2; vmin = -Ac * math.exp(-5 * smin) + Bf * math.exp(-7 * smin)
ax.plot([smin], [vmin], "o", color=RED, ms=9)
ax.annotate("held at fixed size,\nbut V < 0: anti-de Sitter", (smin, vmin), (0.55, -0.55), color=RED, fontsize=10,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.axhline(0, color=GREY, lw=1, ls=":")
ax.set_ylim(-1.2, 1.2); ax.set_xlim(-0.25, 1.6)
ax.set_xlabel("hidden sphere's size, ψ = ln b"); ax.set_ylabel("4D potential (arbitrary units)")
ax.set_title("C. The hidden sphere's dial has the wrong sign (D1, D2)")
ax.legend(loc="upper right", fontsize=9)

# D: ledger
ax = fig.add_subplot(gs[1, :]); ax.axis("off")
ax.set_title("D. What each reading of the dial gives", fontsize=13)
rows = [r for r in R["D4"]["rows"]]
head = ["reading of the dial", "slope λ", "forced or chosen", "Δχ² vs Λ", "w0, wa", "verdict"]
verdict = {
    "cosmological constant": "baseline",
    "best fit (dial-potential search)": "the hint DESI DR3 and Euclid will test; a fit, not a prediction",
    "conformal sqrt(2/3)": "inside the band, but nothing in the law picks it",
    "slice curvature sqrt 2 (D3)": "this is space's curvature, w = −1/3: ruled out",
    "hidden-sphere curvature sqrt(10/3) (D1), sign flipped": "real sign V < 0 cannot accelerate (D2); even flipped, ruled out",
    "meron on hidden sphere 7/sqrt(15/2) (D2)": "no history reaches today's dark energy; with curvature, anti-de Sitter",
}
cells = []
for r in rows:
    cells.append(["\n".join(textwrap.wrap(r["name"], 34)), f"{r['lam']:.3f}", "\n".join(textwrap.wrap(r["status"], 26)),
                  "—" if r["chi2"] is None else f"{r['dchi2_vs_LCDM']:+.1f}",
                  "—" if r["w0"] is None else f"{r['w0']:.2f}, {r['wa']:.2f}",
                  "\n".join(textwrap.wrap(verdict[r["name"]], 46))])
tb = ax.table(cellText=cells, colLabels=head, loc="center", cellLoc="left",
              colWidths=[0.22, 0.07, 0.17, 0.08, 0.1, 0.33])
tb.auto_set_font_size(False); tb.set_fontsize(10.5); tb.scale(1, 3.1)
colr = {"forced": "#fbe3e1", "chosen": "#e3f2e8", "fitted": "#f0f0f0", "measured": "#f0f0f0"}
for i, r in enumerate(rows, start=1):
    key = next(k for k in colr if r["status"].startswith(k))
    for j in range(len(head)):
        tb[i, j].set_facecolor(colr[key])
for j in range(len(head)):
    tb[0, j].set_text_props(weight="bold")
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "f36_dial_slope.png"), dpi=130, bbox_inches="tight")
