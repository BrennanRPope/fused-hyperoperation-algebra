import math, json
import numpy as np, matplotlib
import matplotlib.patches
from scipy.integrate import solve_ivp
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.size": 10})

def deff_of(F, r, h=1e-4):
    lr = math.log(r); a, b = abs(F(math.exp(lr + h))), abs(F(math.exp(lr - h)))
    return 1 - (math.log(a) - math.log(b)) / (2 * h)

# ---------------- h1: the map
fig, ax = plt.subplots(figsize=(12, 7.2))
rows = ["unlimited reach\n(Δ = 0, massless)", "screened\n(Δ > 0, massive carrier)", "rungs\n(Δ < 0, closed direction)", "tube\n(D_eff = 1, confined)"]
cols = ["order off (abelian:\ncarriers don't carry the charge)", "order on (non-abelian:\ncarriers carry the charge)"]
cells = {
    (0, 0): ("known", "light (electromagnetism)\n\n(gravity also has unlimited reach, but it is\nspin 2 and self-coupled: it sits outside\nthis two-column split)"),
    (0, 1): ("open", "EMPTY: a long-range non-abelian force.\nStays long-range only with many matter\nfields (SU(3): about 8+ flavours, frozen\nstrength; 17+, free); otherwise it confines."),
    (1, 0): ("open", "EMPTY among fundamental forces\n(pion exchange sits here, but it comes\nfrom the strong force). Searched as\nfifth forces, dark photon, Z'."),
    (1, 1): ("known", "weak force (W, Z)"),
    (2, 0): ("open", "EMPTY: evenly spaced heavier copies of\nlight or gravity = one hidden extra\ndirection (Kaluza–Klein). Searched in\nshort-range gravity tests."),
    (2, 1): ("open", "EMPTY: Kaluza–Klein copies of W or gluons\n(the Yang–Mills page put the\nglueball gap's ladder on this side)"),
    (3, 0): ("half", "only inside materials: magnetic flux\ntubes in superconductors.\nEMPTY among fundamental forces."),
    (3, 1): ("known", "strong force (quarks, gluons)"),
}
col = {"known": "#cfe8cf", "open": "#fde0c5", "half": "#f3f0c8"}
for (i, j), (kind, txt) in cells.items():
    y = 3 - i
    ax.add_patch(Rectangle((j, y), 1, 1, facecolor=col[kind], edgecolor="k"))
    ax.text(j + 0.5, y + 0.5, txt, ha="center", va="center", fontsize=9.5,
            fontweight="bold" if kind == "known" else "normal")
for i, rl in enumerate(rows):
    ax.text(-0.04, 3 - i + 0.5, rl, ha="right", va="center", fontsize=10)
for j, cl in enumerate(cols):
    ax.text(j + 0.5, 4.06, cl, ha="center", va="bottom", fontsize=10)
ax.set_xlim(-0.75, 2.02); ax.set_ylim(-0.05, 4.45); ax.axis("off")
ax.set_title("h1. The forces our readings allow. Green: forces we know. Orange: empty cells (all are hypotheses physicists already search for)\n"
             "The theory sets the shape of each cell's force; it sets the strength only in the rung cell, once electromagnetism is read as the closed direction (h6)", fontsize=11)
fig.tight_layout(); fig.savefig("figures/h1_force_map.png", dpi=150); plt.close(fig)

# ---------------- h2: D_eff signatures
rs = np.exp(np.linspace(math.log(0.01), math.log(3e4), 700))
Ft = lambda r: (lambda h: -(-(1/math.tanh((r+h)/2))/(r+h) + (1/math.tanh((r-h)/2))/(r-h))/(2*h))(1e-6*r)
Fy = lambda a: (lambda r: (1 + a*math.exp(-r)*(1 + r))/r**2)
mu_g = lambda gN: gN/2*(1 + math.sqrt(1 + 4/gN))
Fm = lambda r: mu_g(1/(r/30)**2)
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.plot(rs, [3]*len(rs), color="gray", lw=2, label="plain gravity or light: 3 at every distance")
ax.plot(rs, [deff_of(Ft, r) for r in rs], color="#1f77b4", lw=2.5, label="rung class (hidden extra direction, size R = 1): steps from 4 down to 3")
ax.plot(rs, [deff_of(Fy(1.0), r) for r in rs], color="#d62728", lw=2.5, label="screened class (Yukawa, α = +1, range 1): a bump")
ax.plot(rs, [deff_of(Fy(-0.3), r) for r in rs], color="#d62728", lw=1.8, ls="--", label="screened class, α = −0.3: a dip")
ax.plot(rs, [deff_of(Fm, r) for r in rs], color="#2ca02c", lw=2.5, label="flux exponent drifting to 3/2 at low field (deep MOND): falls to 2")
ax.set_xscale("log"); ax.set_ylim(1.7, 4.3)
ax.set_xlabel("distance from the source (units of the new length scale)"); ax.set_ylabel("effective dimension D_eff = 1 − d ln F / d ln r")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("h2. Each empty cell leaves its own fingerprint on the effective-dimension reading")
fig.tight_layout(); fig.savefig("figures/h2_deff_fingerprints.png", dpi=150); plt.close(fig)

# ---------------- h3: rung tower
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 2.2]})
ax = axs[0]
for n in range(0, 6):
    ax.hlines(n, -0.4, 0.4, color="k" if n == 0 else "#1f77b4", lw=3)
    ax.text(0.5, n, "rung 0: massless (the force we see)" if n == 0 else f"rung {n}: mass {n}·√(−Δ), two states (±{n})", va="center", fontsize=9)
ax.set_xlim(-0.6, 3.3); ax.set_ylim(-0.6, 5.6); ax.axis("off")
ax.set_title("A: the law's rungs, evenly spaced")
ax = axs[1]
x = np.exp(np.linspace(math.log(0.05), math.log(30), 400))
ax.plot(x, 1/np.tanh(x/2), color="#1f77b4", lw=2.5, label="all rungs summed: coth(r / 2R) (exact; equals the field on space × circle)")
ax.plot(x, 1 + 2*np.exp(-x), color="#d62728", lw=1.8, ls="--", label="lightest rung only: 1 + 2 e^(−r/R) (Yukawa with α = 2)")
ax.plot(x, 2/x, color="gray", lw=1.5, ls=":", label="4-dimensional force law 2R / r")
ax.axhline(1, color="k", lw=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.8, 50)
ax.set_xlabel("distance r / R"); ax.set_ylabel("potential relative to the plain 1/r law")
ax.legend(fontsize=8.5); ax.set_title("B: the copies add up to one extra direction at short range")
fig.suptitle("h3. The rung class: a closed direction of length 2πR gives copies with masses n/R; R = 1/√(−Δ)")
fig.tight_layout(); fig.savefig("figures/h3_rung_tower.png", dpi=150); plt.close(fig)

# ---------------- h4: roundness
th = np.linspace(0, 2*np.pi, 721)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8), subplot_kw={}, gridspec_kw={"width_ratios": [1, 1.2]})
axs[0].remove(); axp = fig.add_subplot(1, 2, 1, projection="polar")
for p, c in [(1.5, "#9467bd"), (2.0, "k"), (3.0, "#ff7f0e")]:
    X, Y = np.cos(th), np.sin(th); h = 1e-6
    V = lambda X, Y: 1/(np.abs(X)**p + np.abs(Y)**p)**(1/p)
    F = np.hypot((V(X+h, Y) - V(X-h, Y))/(2*h), (V(X, Y+h) - V(X, Y-h))/(2*h))
    axp.plot(th, F/np.mean(F), color=c, lw=2, label=f"p = {p}")
axp.set_ylim(0.8, 1.15); axp.set_yticklabels([]); axp.legend(fontsize=8, loc="lower left", bbox_to_anchor=(-0.15, -0.12))
axp.set_title("A: strength of the pull at one distance, all directions", fontsize=10)
R = json.load(open("results.json"))["E_cos4_amplitude"]
ps = sorted(float(k) for k in R)
axs[1].plot(ps, [R[str(p)] for p in ps], "o", color="#ff7f0e", ms=8)
axs[1].set_xlabel("roundness dial p (2 = round)"); axs[1].set_ylabel("size of the four-lobed (cos 4θ) part")
axs[1].set_title("B: no 1st, 2nd or 3rd harmonic; the first change is four lobes")
fig.suptitle("h4. The roundness dial on a force: away from p = 2 the pull depends on direction, with four lobes")
fig.tight_layout(); fig.savefig("figures/h4_roundness_force.png", dpi=150); plt.close(fig)

# ---------------- h5: deep MOND
r = np.exp(np.linspace(math.log(0.3), math.log(300), 400))
gN = 1/r**2
g = np.array([mu_g(v) for v in gN])
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
axs[0].plot(r, np.sqrt(gN*r), color="gray", lw=2, label="Newton: speed falls as 1/√r")
axs[0].plot(r, np.sqrt(g*r), color="#2ca02c", lw=2.5, label="flux exponent drifting to s = 3/2 at low field")
axs[0].axhline(1, color="#2ca02c", ls=":", lw=1.2); axs[0].text(30, 1.04, "flat: v⁴ = G M a₀", color="#2ca02c")
axs[0].axvline(1, color="k", ls=":", lw=1); axs[0].text(1.1, 0.25, "where the pull\nequals a₀", fontsize=9)
axs[0].set_xscale("log"); axs[0].set_xlabel("radius (units where G M = a₀ = 1)"); axs[0].set_ylabel("orbital speed")
axs[0].legend(fontsize=9); axs[0].set_title("A: orbital speed around a point mass")
sv = np.linspace(0.9, 2.2, 200)
axs[1].plot(sv, 1 + 2*(sv - 1), color="#d62728", lw=2.5)
for s_, lab, c in [(2, "s = 2: Newton / Maxwell", "#d62728"), (1.5, "s = 3/2: deep MOND,\nthe confinement border,\nunchanged by stretching space", "#2ca02c"), (1, "s = 1: tube", "k")]:
    axs[1].plot([s_], [1 + 2*(s_ - 1)], "o", color=c, ms=9); axs[1].text(s_ + 0.04, 1 + 2*(s_ - 1) - (0.15 if s_ == 1 else 0.35), lab, fontsize=9, color=c, va="top" if s_ != 1 else "bottom")
axs[1].set_ylim(0.5, 3.6); axs[1].set_xlabel("flux exponent s (field energy ∝ |D|^s)"); axs[1].set_ylabel("D_eff in 3D")
axs[1].set_title("B: where deep MOND sits on the flux-exponent dial")
fig.suptitle("h5. Not a new force but a known force changing at low field: s = 3/2 is the deep-MOND law (Bekenstein–Milgrom 1984)")
fig.tight_layout(); fig.savefig("figures/h5_deep_mond.png", dpi=150); plt.close(fig)
print("ok")

# ---------------- h6: strength from the closed direction (Part G)
fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.0))
ax = axs[0]
groups = ["two copies,\nsame charge", "two copies,\nopposite charge", "ordinary matter\n(no circle charge)"]
tens = [0.5, 0.5, 0.5]; sca = [1.5, 1.5, 0.5/3]; vec = [-2.0, 2.0, 0.0]
xs = np.arange(3); w = 0.24
ax.bar(xs - w, tens, w, color="#7f7f7f", label="gravity (tensor)")
ax.bar(xs, sca, w, color="#2ca02c", label="circle-size field (scalar)")
ax.bar(xs + w, vec, w, color="#1f77b4", label="electric (vector); negative = repels")
tot = [t + s + v for t, s, v in zip(tens, sca, vec)]
for x, t in zip(xs, tot):
    ax.text(x, 2.45, "net " + {0.0: "0", 4.0: "4"}.get(round(t, 6), "2/3"), ha="center", fontsize=10, fontweight="bold")
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(xs); ax.set_xticklabels(groups, fontsize=9)
ax.set_ylim(-2.4, 2.8); ax.set_ylabel("pull, in units where Newton between masses = 1/2")
ax.legend(fontsize=8, loc="lower right")
ax.set_title("A. One extra circle splits into three forces\n(checks G1, G2, G6)", fontsize=10.5)

ax = axs[1]
Rs = np.exp(np.linspace(math.log(1), math.log(1e3), 300))
ax.plot(Rs, 4 / Rs**2, color="#1f77b4", lw=2.5, label="our reading: α = 4 G |Δ| = 4 (l_P / R)²")
ax.axhline(1/137.036, color="#d62728", ls="--", lw=1.5, label="measured α = 1/137.036")
Rem = 2 * math.sqrt(137.036)
ax.plot([Rem], [1/137.036], "o", color="k", ms=8)
ax.annotate(f"circle radius R = {Rem:.1f} Planck lengths\n= 3.8e-34 m; first copy 5.2e17 GeV",
            (Rem, 1/137.036), xytext=(60, 0.15), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("circle radius R (Planck lengths)"); ax.set_ylabel("electric strength α")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("B. The same reading that spaces the rungs\nsets the electric strength (checks G3 to G5)", fontsize=10.5)

ax = axs[2]
ns = np.arange(1, 8)
ax.plot(ns, ns / (ns + 2), "s-", color="#2ca02c", lw=2, label="new force's strength vs gravity: n/(n+2)")
ax.plot(ns, 1 / (ns + 1), "o-", color="#9467bd", lw=2, label="light-bending γ if it reached far: 1/(n+1)")
ax.axhline(1.0, color="#d62728", lw=1.5, label="Cassini: γ = 1 to 1 part in 40,000")
ax.annotate("one circle: strength 1/3, γ = 1/2\n(21,740 standard errors off)", (1, 0.5), xytext=(2.4, 0.86),
            fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("number of hidden circles n"); ax.set_ylim(0, 1.6); ax.set_xticks(ns)
ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("C. The predicted new force must be short-range\n(the theory doesn't set the range; checks G6 to G7b)", fontsize=10.5)
fig.suptitle("h6. Strength derived in the special case: electromagnetism read as gravity's closed direction", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h6_strength_from_circle.png", dpi=150); plt.close(fig)

# ---------------- h7: what holds the circle, and the range it gives (part H)
fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))
ax = axs[0]
Rg = np.linspace(0.55, 3.0, 400)
ax.plot(Rg, 1/Rg**6, color="#1f77b4", lw=2, label="vacuum energy alone, repulsive: the circle grows")
ax.plot(Rg, -1/Rg**6, color="#d62728", lw=2, label="vacuum energy alone, attractive: the circle shrinks")
ax.plot(Rg, -6/Rg + 1/Rg**6, color="#2ca02c", lw=2.8, label="plus a negative 5D vacuum energy: held at R*")
ax.plot([1], [-5], "o", color="k", ms=7)
ax.annotate("held here; leaves a 4D vacuum\nenergy −5C/R*⁴ to cancel", (1, -5), xytext=(1.35, -2.0), fontsize=9,
            arrowprops=dict(arrowstyle="->"))
ax.set_ylim(-7, 4); ax.set_xlabel("circle radius R / R*"); ax.set_ylabel("energy of the circle (units of C/R*⁴)")
ax.legend(fontsize=8, loc="lower right")
ax.set_title("A. Its own vacuum energy can't hold the circle;\nwith a second term it can (checks H3, H4)", fontsize=10.5)

ax = axs[1]
marks = [(lP_m := 1.616255e-35, "Planck length"), (3.784e-34, "circle radius\n(part G)"),
         (3.21e-32, "new force's range\n(this hold)"), (8.4e-16, "proton"), (38.6e-6, "lab limit\n(Eöt-Wash, 38.6 µm)")]
pos = [(0.30, "center"), (-0.08, "center"), (0.08, "left"), (-0.08, "center"), (0.08, "center")]
for (v, lab_), (dy, ha_) in zip(marks, pos):
    ax.plot([math.log10(v)], [0], "o", ms=9, color="#d62728" if "new force" in lab_ else "k")
    ax.text(math.log10(v) + (0.3 if ha_ == "left" else 0), dy, lab_, ha=ha_, va="bottom" if dy > 0 else "top", fontsize=9)
ax.plot([math.log10(1.616255e-35)]*2, [0.02, 0.28], color="gray", lw=0.8)
ax.axhline(0, color="gray", lw=1)
ax.annotate("", xy=(math.log10(38.6e-6), -0.42), xytext=(math.log10(3.21e-32), -0.42),
            arrowprops=dict(arrowstyle="<->", lw=1.8))
ax.text((math.log10(38.6e-6) + math.log10(3.21e-32))/2, -0.47, "a factor of 10²⁷", ha="center", va="top", fontsize=10, fontweight="bold")
ax.set_ylim(-0.65, 0.5); ax.set_xlim(-36, -2); ax.set_yticks([])
ax.set_xlabel("distance (log10 of metres)")
ax.set_title("B. The new force's range is about 2000 Planck lengths,\nfar below anything a lab can reach (checks H5 to H7)", fontsize=10.5)

ax = axs[2]
ks = np.logspace(-30, 0.5, 300)
ax.plot(ks, (23.41**2/ks)*1.616255e-35, color="#9467bd", lw=2.5, label="range = R² / (κ l_P), circle at 23.4 l_P")
ax.axhline(38.6e-6, color="#d62728", ls="--", label="lab limit, 38.6 µm")
ax.plot([0.276], [3.21e-32], "o", color="k", ms=7)
ax.annotate("vacuum-energy hold: κ ≈ 0.3", (0.276, 3.21e-32), xytext=(1e-16, 1e-25), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.plot([2.3e-28], [38.6e-6], "s", color="#d62728", ms=7)
ax.annotate("κ ≈ 2e-28 would be needed\nto reach the lab", (2.3e-28, 38.6e-6), xytext=(1e-24, 1e-10), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("strength κ of the hold (mass = κ l_P / R²)"); ax.set_ylabel("range of the new force (m)")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("C. Any hold built only from the circle and gravity\ngives a Planck-short range (check H6b)", fontsize=10.5)
fig.suptitle("h7. What holds the circle: its own vacuum energy, with help. The new force is then real but far too short-range to see", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h7_what_holds_the_circle.png", dpi=150); plt.close(fig)

# ---------------- h8: ways to derive strength: shape and material (part S)
fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))
ax = axs[0]
names = ["circle\n(abelian, like light)", "sphere\n(non-abelian SU(2))", "needed for the weak\nand strong forces"]
vals = [4, 6, 336]
ax.bar(range(3), vals, color=["#1f77b4", "#2ca02c", "#d62728"])
for i, v in enumerate(vals):
    ax.text(i, v*1.15, f"{v}", ha="center", fontsize=11, fontweight="bold")
ax.set_yscale("log"); ax.set_ylim(1, 2000); ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("shape factor s in α = s G / size²")
ax.set_title("A. Quality: the hidden space's shape sets\nthe number in front (checks S1 to S3, S6)", fontsize=10.5)

ax = axs[1]
MZ_ = 91.1876; aem_inv_ = 127.951; s2w_ = 0.23122
a1_ = 0.6*aem_inv_*(1 - s2w_); a2_ = aem_inv_*s2w_; a3_ = 1/0.1180
mus = np.logspace(math.log10(MZ_), 19.1, 300)
for ai, bi, c_, lab_ in [(a1_, 41/10, "#7f7f7f", "hypercharge (unified scale)"), (a2_, -19/6, "#1f77b4", "weak, SU(2)"), (a3_, -7.0, "#d62728", "strong, SU(3)")]:
    ax.plot(mus, ai - bi/(2*math.pi)*np.log(mus/MZ_), color=c_, lw=2.2, label=lab_)
ax.axvline(9.71e16, color="k", ls=":", lw=1.3)
ax.text(9.71e16/1.5, 22, "weak and strong meet: 9.7e16 GeV", ha="right", fontsize=8.5)
ax.axvspan(7.3e17, 8.9e17, color="#2ca02c", alpha=0.25)
ax.text(7e17/1.4, 6, "green band: copies of a circle or sphere\nwith that strength sit at 8e17 GeV", ha="right", fontsize=8.5)
ax.set_xscale("log"); ax.set_ylim(0, 70)
ax.set_xlabel("energy (GeV)"); ax.set_ylabel("1 / strength (1/α)")
ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("B. Material: the particles each force carries\nchange its strength with energy (checks S5, S6)", fontsize=10.5)

ax = axs[2]
mus2 = np.logspace(0, 4, 300)
for nf, c_ in [(3, "#ff7f0e"), (5, "#d62728"), (6, "#8c564b"), (16, "#9467bd")]:
    b0_ = 11 - 2*nf/3
    Lam_ = MZ_*math.exp(-2*math.pi/(b0_*0.1180))
    ax.plot(mus2, 2*math.pi/(b0_*np.log(mus2/Lam_)), color=c_, lw=2.2,
            label=f"{nf} flavours: b0 = {b0_:.2f}, Λ = {Lam_*1000:.3g} MeV" if Lam_ > 1e-6 else f"{nf} flavours: b0 = {b0_:.2f}")
ax.axvline(MZ_, color="gray", ls=":", lw=1); ax.text(MZ_*1.1, 0.06, "all matched at the Z mass", fontsize=8.5)
ax.set_xscale("log"); ax.set_ylim(0.05, 0.4)
ax.set_xlabel("energy (GeV)"); ax.set_ylabel("strong strength α_s")
ax.legend(fontsize=8, loc="upper right")
ax.set_title("C. The strong force's strength comes from its material\nand one scale Λ (check S4)", fontsize=10.5)
fig.suptitle("h8. Ways to derive a force's strength: its shape (quality) sets the number in front, its particles (material) set how it runs", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h8_strength_shape_and_material.png", dpi=150); plt.close(fig)

# ---------------- h9: what carries it and what feels it (part I)
res_ = json.load(open("results.json"))
lam_ = float(res_["H_range_lP_NF2"]); Rl_ = 2*math.sqrt(137.035999084); mrad_ = float(res_["I_radion_GeV"]); M1_ = float(res_["I_copy1_GeV"])
fig, axs = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axs[0]
levels = [(1.4e4, "LHC collision energy (14 TeV)", "gray"), (4.7e13, "fastest possible inflation (H, from r < 0.036)", "#9467bd"),
          (mrad_, "the carrier: the circle's size (spin 0, no charge)", "#2ca02c")] + \
         [(n*M1_, "copies n = 1, 2, 3 (charges e, 2e, 3e)" if n == 2 else "", "#1f77b4") for n in (1, 2, 3)] + [(1.22e19, "Planck mass", "k")]
for v, lab_, c_ in levels:
    ax.hlines(math.log10(v), 0, 1, color=c_, lw=2.5)
    if lab_:
        ax.text(1.04, math.log10(v), lab_, va="center", fontsize=9, color=c_)
ax.set_xlim(0, 2.6); ax.set_ylim(3, 20); ax.set_xticks([])
ax.set_ylabel("mass or energy (log10 GeV)")
ax.text(0.5, 9, "a gap of 10¹¹:\nnothing from this\ntheory down here", ha="center", fontsize=9, color="gray")
ax.set_title("A. The ladder: one neutral carrier, then charged copies\n(checks I1, I9, I13)", fontsize=10.5)
ax = axs[1]
who = ["ordinary\nmatter", "a copy\n(charge ne)", "a light ray\n(passing)", "electric\nstrength α"]
vals_ = [1/3, 3, 0, 0]
ax.bar(range(4), vals_, color=["#7f7f7f", "#1f77b4", "#ff7f0e", "#ffffff"], edgecolor="k")
for i, v in enumerate(vals_[:3]):
    ax.text(i, v + 0.08, ["1/3", "3", "0"][i], ha="center", fontsize=11, fontweight="bold")
ax.text(3, 0.5, "changes with\nthe circle:\nα ∝ 1/R²", ha="center", fontsize=9)
ax.set_xticks(range(4)); ax.set_xticklabels(who, fontsize=9); ax.set_ylim(0, 3.6)
ax.set_ylabel("pull of the new force, relative to gravity's")
ax.set_title("B. Who feels it (checks I2, I3)", fontsize=10.5)
fig.suptitle("h9. What carries the new force and what feels it", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h9_carrier_and_couplings.png", dpi=150); plt.close(fig)

# ---------------- h10: the force and its field up close (part I)
fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))
ax = axs[0]
rr = np.exp(np.linspace(math.log(1), math.log(2e5), 600))
fl = lambda r: 1 + math.exp(-r/lam_)/3 + 4/3*(1/math.tanh(r/(2*Rl_)) - 1)
ax.plot(rr, [fl(r) for r in rr], color="#2ca02c", lw=2.5, label="gravity plus the new force (circle held)")
ax.plot(rr, [4/3/math.tanh(r/(2*Rl_)) for r in rr], color="#1f77b4", ls="--", lw=1.6, label="if the circle were free: (4/3) coth(r/2R)")
ax.axhline(1, color="gray", lw=1); ax.axhline(4/3, color="gray", lw=0.8, ls=":")
ax.axvline(Rl_, color="k", ls=":", lw=1); ax.axvline(lam_, color="#d62728", ls=":", lw=1)
ax.text(Rl_*1.1, 2.75, "circle\n23 l_P", fontsize=9); ax.text(lam_*1.1, 2.75, "range\n~2000 l_P", fontsize=9, color="#d62728")
ax.text(3e4, 1.05, "Newton", fontsize=9, color="gray"); ax.text(3e4, 1.36, "4/3 Newton", fontsize=9, color="gray")
ax.text(1.3, 2.3, "5D gravity\n(grows as R/r)", fontsize=9)
ax.set_xscale("log"); ax.set_ylim(0.9, 3.0)
ax.set_xlabel("distance between two masses (Planck lengths)"); ax.set_ylabel("V r / (G m₁ m₂): pull relative to Newton")
ax.legend(fontsize=8, loc="center right")
ax.set_title("A. The force law: 5D inside the circle, 4/3 just outside,\nNewton beyond the range (checks I4 to I4c)", fontsize=10.5)
ax = axs[1]
Ft_ = lambda r: (lambda h: -(fl(r + h)/(r + h) - fl(r - h)/(r - h))/(2*h))(1e-5*r)
ax.plot(rr, [deff_of(Ft_, r) for r in rr], color="#2ca02c", lw=2.5)
ax.axhline(3, color="gray", lw=1); ax.axhline(4, color="gray", lw=0.8, ls=":")
ax.axvline(Rl_, color="k", ls=":", lw=1); ax.axvline(lam_, color="#d62728", ls=":", lw=1)
ax.text(1.3, 4.05, "4: one extra direction", fontsize=9); ax.text(3e4, 3.03, "3", fontsize=9)
ax.text(3e3, 3.5, "small bump: the new\nforce switching off", ha="center", fontsize=9, color="#d62728")
ax.set_xscale("log"); ax.set_ylim(2.8, 4.3)
ax.set_xlabel("distance (Planck lengths)"); ax.set_ylabel("effective dimension D_eff")
ax.set_title("B. The same in the D_eff reading: a step at the circle,\na bump at the range (section 3's fingerprints together)", fontsize=10.5)
ax = axs[2]
g_ = np.linspace(-3*lam_, 3*lam_, 401); XX, YY = np.meshgrid(g_, g_); RRr = np.sqrt(XX**2 + YY**2) + lam_/50
field = np.exp(-RRr/lam_)/(RRr/lam_)
im = ax.imshow(np.log10(field), extent=[-3, 3, -3, 3], origin="lower", cmap="Greens", vmin=-3, vmax=1.7)
th_ = np.linspace(0, 2*math.pi, 200); ax.plot(np.cos(th_), np.sin(th_), color="#d62728", lw=1.5)
ax.text(0.75, 0.75, "one range", color="#d62728", fontsize=9)
ax.plot([0], [0], "ko", ms=5); ax.text(0.12, -0.25, "mass", fontsize=9)
cb = fig.colorbar(im, ax=ax, fraction=0.046); cb.set_label("log10 of how much the circle swells (relative)")
ax.set_xlabel("distance (ranges)"); ax.set_ylabel("distance (ranges)")
ax.set_title("C. The field: the circle swells near a mass,\ndR/R = (2/3) G M e^(−r/λ)/(r c²); α dips twice as much (I5)", fontsize=10.5)
fig.suptitle("h10. The new force up close, at its own range", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h10_force_up_close.png", dpi=150); plt.close(fig)

# ---------------- h11: indirect traces (part I)
fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))
ax = axs[0]
Ms = np.linspace(0, 4, 50)
ax.fill_between(Ms, 0, 2*Ms, color="#cfe8cf", label="black holes allowed: charge ≤ 2 √G × mass")
ax.plot(Ms, 2*Ms, color="k", lw=1.5, label="extreme black holes (a = √3)")
ax.plot([1, 2, 3], [2, 4, 6], "o", color="#1f77b4", ms=9, label="the copies n = 1, 2, 3")
ax.set_xlabel("mass (units of the first copy)"); ax.set_ylabel("charge (units of √G × first copy's mass)")
ax.set_xlim(0, 4); ax.set_ylim(0, 8); ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("A. Each copy sits exactly on the extreme-black-hole line,\nso like copies feel no net force (checks I7, G1)", fontsize=10.5)
ax = axs[1]
Tlev = [(1.4e4, "LHC"), (4.7e13, "inflation's H, at most"), (1.41e16, "hottest the universe could\nhave been after inflation"),
        (float(res_["I_relic_copy_Tsafe_GeV"]), "copies made only above about here\n(rough)"), (mrad_, "carrier mass"),
        (float(res_["I_settle_T_GeV"]), "carrier settles if shaken"), (M1_, "first copy")]
order_ = sorted(Tlev, key=lambda t: -t[0])
ylab = {}
for j, (v, lab_) in enumerate(order_[:5]):
    ylab[lab_] = 18.0 - j
ylab[order_[5][1]] = 12.6; ylab[order_[6][1]] = math.log10(order_[6][0])
for v, lab_ in Tlev:
    theory = ("copy" in lab_ or "carrier" in lab_ or "copies" in lab_)
    ax.plot([0], [math.log10(v)], "o", color="#d62728" if theory else "k", ms=7)
    ax.plot([0.02, 0.25], [math.log10(v), ylab[lab_]], color="gray", lw=0.7)
    ax.text(0.28, ylab[lab_], lab_.replace("\n", " "), ha="left", va="center", fontsize=8.5, color="#d62728" if theory else "k")
ax.axvline(0, color="gray", lw=1); ax.set_xlim(-0.2, 1.6); ax.set_xticks([]); ax.set_ylim(3, 18.6)
ax.set_ylabel("energy or temperature (log10 GeV)")
ax.set_title("B. Early universe: the copies are too heavy ever to have been made;\nthe carrier, if made, decays at once (I9, I11, I12)", fontsize=10.5)
ax = axs[2]
ss = np.linspace(-0.5, 6, 400)
epsf = lambda s_: 0.5*((-6*(-math.sqrt(2/3))*math.exp(-math.sqrt(2/3)*s_) + (-6*math.sqrt(2/3))*math.exp(-6*math.sqrt(2/3)*s_)) /
                       (-6*math.exp(-math.sqrt(2/3)*s_) + math.exp(-6*math.sqrt(2/3)*s_)))**2
ax.plot(ss, [epsf(s_) for s_ in ss], color="#2ca02c", lw=2.5, label="slope parameter ε of the circle's potential")
ax.axhline(1/3, color="#2ca02c", ls=":", lw=1); ax.text(4.2, 0.36, "tail: ε → 1/3", fontsize=9)
ax.axhline(0.01, color="#d62728", ls="--", label="slow roll needs ε < 0.01")
ax.set_yscale("log"); ax.set_ylim(1e-3, 1e2)
ax.set_xlabel("circle size field σ (reduced Planck units; R = R* e^(√(2/3) σ))"); ax.set_ylabel("ε")
ax.legend(fontsize=8.5, loc="upper right")
ax.set_title("C. The circle's size can't be the inflaton: its potential\nis too steep everywhere along the tail (I10)", fontsize=10.5)
fig.suptitle("h11. Indirect traces: black holes, the early universe and inflation", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h11_indirect_traces.png", dpi=150); plt.close(fig)

# ---------------- h12: the inverse (part J) ----------------
rj = json.load(open("results.json"))
fig, axs = plt.subplots(2, 2, figsize=(13.5, 10))
ax = axs[0, 0]
pairs = [("matter + matter", [1, 1/3, 0]), ("matter + copy", [1, 1, 0]), ("copy + copy\n(like charges)", [1, 3, -4])]
w = 0.25; cols_ = ["#1f77b4", "#2ca02c", "#d62728"]; labs_ = ["gravity (spin 2)", "circle's size (spin 0)", "circle's vector (spin 1)"]
for i, (nm, v) in enumerate(pairs):
    for k in range(3):
        ax.bar(i + (k - 1)*w, v[k], w, color=cols_[k], label=labs_[k] if i == 0 else None)
        if v[k] != 0:
            ax.text(i + (k - 1)*w, v[k] + (0.12 if v[k] > 0 else -0.35), f"{v[k]:+.2g}" if k != 1 or i else "+1/3", ha="center", fontsize=9)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(range(3)); ax.set_xticklabels([p[0] for p in pairs])
ax.set_ylabel("strength, in units of the pair's gravity\n(+ pulls together, − pushes apart)"); ax.set_ylim(-5, 4.2)
ax.legend(fontsize=9, loc="lower left")
ax.set_title("A. Only the spin-1 piece pushes, and only copies carry it.\nOrdinary matter never feels a push from the circle (J1, J2)", fontsize=10.5)
ax = axs[0, 1]
lam_ = float(rj["H_range_lP_NF2"]); rr = np.logspace(1.3, 5, 400)
ax.plot(rr, 1 + 3*np.exp(-rr/lam_) - 4, color="#d62728", lw=2.5, label="two like copies: 1 + 3e^(−r/λ) − 4")
ax.plot(rr, 1 + (1/3)*np.exp(-rr/lam_), color="#1f77b4", lw=2.5, label="two lumps of matter: 1 + (1/3)e^(−r/λ)")
ax.axhline(0, color="k", lw=0.8); ax.axvline(lam_, color="gray", ls=":")
ax.text(lam_*1.1, 0.6, "size carrier's range λ ≈ 2000 l_P", fontsize=9, color="gray")
ax.text(25, -2.0, "beyond λ: a push of 3 × gravity,\nwith no end to its reach\n(this is electromagnetism\nfor our circle's copies)", fontsize=9, color="#d62728", va="top")
ax.set_xscale("log"); ax.set_ylim(-3.6, 2.0)
ax.set_xlabel("distance r (Planck lengths; the tower of rungs near the circle size is left out)"); ax.set_ylabel("net force / gravity")
ax.legend(fontsize=9, loc="upper right")
ax.set_title("B. The 'inverse' is already in the circle: a short-range pull\nand an unlimited-range push, from the same circle (J3)", fontsize=10.5)
ax = axs[1, 0]
R2 = np.logspace(0, 31, 400)
ax.plot(R2, 4/R2**2, color="k", lw=2, label="coupling rule α2 = 4G/R2²")
lo, hi = float(rj["J_R2_min_lP"]), float(rj["J_R2_max_lP"])
ax.axvspan(lo, hi, color="#2ca02c", alpha=0.2, label=f"allowed with an electric circle: {lo:.1f} to {hi:.0f} l_P")
ax.axvspan(hi, 1e31, color="#d62728", alpha=0.08)
ax.text(1e11, 1e-15, "forbidden: the electric circle would be\nsmaller than the fundamental length (J5)", fontsize=9, color="#d62728")
R1_ = float(rj["G_circle_radius_lP"])
ax.plot([R1_], [1/137.036], "o", color="#1f77b4", ms=8); ax.annotate("our electric circle, α = 1/137", (R1_, 1/137.036), xytext=(1e5, 1e-2), fontsize=9, arrowprops=dict(arrowstyle="->", color="#1f77b4"), color="#1f77b4")
ax.plot([hi], [4/hi**2], "o", color="#2ca02c", ms=8); ax.annotate("weakest allowed: α2 = 1/26,400", (hi, 4/hi**2), xytext=(1e5, 1e-7), fontsize=9, arrowprops=dict(arrowstyle="->", color="#2ca02c"), color="#2ca02c")
R50_ = 50e-6/1.616255e-35; Rde_ = float(rj["J_R_dark_energy_m"])/1.616255e-35
ax.plot([R50_], [4/R50_**2], "s", color="#d62728", ms=7); ax.annotate("50 µm: my first hand estimate, 1e-61", (R50_, 4/R50_**2), xytext=(1e12, 1e-52), fontsize=9, arrowprops=dict(arrowstyle="->", color="#d62728"), color="#d62728")
ax.plot([Rde_], [4/Rde_**2], "^", color="#9467bd", ms=8); ax.annotate("9 µm: dark energy size (J7)", (Rde_, 4/Rde_**2), xytext=(1e12, 1e-42), fontsize=9, arrowprops=dict(arrowstyle="->", color="#9467bd"), color="#9467bd")
ax.axvline(38.6e-6/1.616255e-35, color="gray", ls="--"); ax.text(38.6e-6/1.616255e-35/2, 1e-28, "gravity tests:\nbelow 38.6 µm", fontsize=8.5, color="gray", ha="right")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-63, 10); ax.set_xlim(1, 1e31)
ax.set_xlabel("second circle's radius R2 (Planck lengths)"); ax.set_ylabel("strength of its push α2")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("C. A second circle: bigger means weaker, but our law caps it\nat 14 × the electric circle, so the weakest push is 1/26,400 (J4-J6)", fontsize=10.5)
ax = axs[1, 1]
Rm = np.logspace(-35, -2, 400)
c5_ = 3*1.0369277551433699/(64*math.pi**6)
ax.plot(Rm, (2*c5_)**0.25*1.973269804e-7/Rm, color="#9467bd", lw=2.5, label="a circle's vacuum energy, (energy density)^(1/4),\ntwo more fermions than bosons")
ax.axhline(float(rj["J_rhoL_quarter_meV"])*1e-3, color="k", ls="--", label="observed dark energy, 2.24 meV")
ax.axvspan(hi*1.616255e-35, 1e-2, color="#d62728", alpha=0.08)
ax.plot([Rde_*1.616255e-35], [float(rj["J_rhoL_quarter_meV"])*1e-3], "^", color="#9467bd", ms=9)
ax.annotate("crosses at 8.8 µm\n(the known 'dark dimension')", (Rde_*1.616255e-35, 2.24e-3), xytext=(1e-17, 3e-6), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.plot([R1_*1.616255e-35], [(2*c5_)**0.25*1.973269804e-7/(R1_*1.616255e-35)], "o", color="#1f77b4", ms=8)
ax.annotate("our electric circle", (R1_*1.616255e-35, (2*c5_)**0.25*1.973269804e-7/(R1_*1.616255e-35)), xytext=(1e-29, 1e26), fontsize=9, color="#1f77b4", arrowprops=dict(arrowstyle="->", color="#1f77b4"))
ax.text(1e-35, 1e1, "red: sizes our law forbids for a\nsecond circle (everything above 5e-33 m)", fontsize=9, color="#d62728")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-6, 1e28)
ax.set_xlabel("circle radius (m)"); ax.set_ylabel("energy scale (eV)")
ax.legend(fontsize=8.5, loc="upper right")
ax.set_title("D. The cosmic push: a 9 µm circle would give dark energy's size,\nbut our law can't hold one next to the electric circle (J7)", fontsize=10.5)
fig.suptitle("h12. Could the new force be the inverse, a very weak, very long-range push?", fontsize=12)
fig.tight_layout(); fig.savefig("figures/h12_the_inverse.png", dpi=150); plt.close(fig)

# ---------------- h13: the other unknown forces (part K)
rk = json.load(open("results.json"))
fig, axs = plt.subplots(2, 3, figsize=(17.5, 10.5))
ax = axs[0, 0]
r_ = np.logspace(-1.3, 1.7, 400)
ax.plot(r_, np.ones_like(r_), color="k", lw=2, label="light: unlimited reach (known)")
ax.plot(r_, (1 + r_)*np.exp(-r_), color="#ff7f0e", lw=2, label="massive new force: dies past its range")
ax.plot(r_, 1/np.tanh(r_/2) + (r_/2)/np.sinh(r_/2)**2, color="#2ca02c", lw=2, label="copies (rungs): stronger inside the circle")
ax.plot(r_, np.ones_like(r_)*1.0, color="#9467bd", lw=2, ls="--", label="long-range non-abelian, frozen strength:\nCoulomb-like at every distance")
ax.plot(r_, 0.1*r_**2, color="#8c564b", lw=2, label="tube: constant pull, never lets go")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-3, 1e3)
ax.set_xlabel("distance (in units of the force's own length: range, circle size or tube width)")
ax.set_ylabel("force × distance², relative to plain 1/r²")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("A. How each unknown force falls off with distance", fontsize=10.5)
ax = axs[0, 1]
ax.plot([3.2e-32], [1/3], "o", color="#d62728", ms=9); ax.annotate("our size carrier\n(strength 1/3, range 3e-32 m)", (3.2e-32, 1/3), xytext=(1e-30, 1e-6), fontsize=9, color="#d62728", arrowprops=dict(arrowstyle="->", color="#d62728"))
ax.fill_between([38.6e-6, 1e12], [1, 1], [1e3, 1e3], color="gray", alpha=0.25)
ax.plot([38.6e-6], [1], "s", color="k"); ax.text(38.6e-6, 3, "Eöt-Wash 2020: at strength 1,\nrange must be below 38.6 µm", fontsize=8.5, ha="center")
lm = rk["K_lambda_moon_m"]; ax.fill_between(lm, [3e-11]*2, [1e3]*2, color="#1f77b4", alpha=0.3)
ax.text(lm[0], 1e-12, "Moon's perigee:\nbelow 3e-11", fontsize=8.5, color="#1f77b4", ha="center")
ax.fill_between([1.5e11, 1e13], [rk["K_alpha_cassini_longrange"]]*2, [1e3]*2, color="#2ca02c", alpha=0.3)
ax.text(1.5e11, 2e-6, "Cassini: a scalar reaching\nbeyond the Sun, below 1e-5", fontsize=8.5, color="#2ca02c", ha="right")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-35, 1e13); ax.set_ylim(1e-14, 1e3)
ax.set_xlabel("range λ (m)"); ax.set_ylabel("strength α, relative to gravity")
ax.set_title("B. A massive new force: our law sets the range, not the strength.\nShaded: ruled out (three anchors from pass 2, not full curves) (K1)", fontsize=10.5)
ax = axs[0, 2]
for N, c_, AF in [(3, "#d62728", 16.5), (2, "#1f77b4", 11.0)]:
    lo_e = float(rk["K_BZ_edge_SU3"] if N == 3 else rk["K_BZ_edge_SU2"])
    nf = np.linspace(lo_e + 0.02, AF - 1e-3, 300)
    b0 = 11*N/3 - 2*nf/3; b1 = 34*N**2/3 - nf*(10*N/3 + (N**2 - 1)/N)
    ax.plot(nf, -4*np.pi*b0/b1, color=c_, lw=2.5, label=f"SU({N}): frozen strength α*, {lo_e:.1f} to {AF:g} flavours")
ax.plot([16], [rk["K_BZ_alpha_SU3_16"]], "o", color="#d62728"); ax.annotate("16 flavours: α* = 0.042\n(Banks–Zaks)", (16, rk["K_BZ_alpha_SU3_16"]), xytext=(12.0, 4e-3), fontsize=9, arrowprops=dict(arrowstyle="->"))
ax.text(0.4, 0.03, "left of each curve:\nthe force confines\n(tube cell, like the\nstrong force)", fontsize=8.5, color="#8c564b")
ax.axhspan(0.3, 1, color="gray", alpha=0.12); ax.text(0.3, 0.5, "above about 0.3 two loops can't be trusted", fontsize=8, color="gray")
ax.set_yscale("log"); ax.set_ylim(1e-3, 1); ax.set_xlim(0, 18)
ax.set_xlabel("number of matter flavours (the 'material')"); ax.set_ylabel("long-range strength α*")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("C. A long-range non-abelian force: in this window the material alone\nsets the strength, forgetting size and shape (two loops; K2)", fontsize=10.5)
ax = axs[1, 0]
def vtx(c, rmax=8.0):
    rhs = lambda r, s: [(1 - s[1])*s[0]/r, r*(1 - s[0]**2)/2]
    return solve_ivp(rhs, [1e-4, rmax], [c*1e-4, 1e-8/4], rtol=1e-11, atol=1e-13, dense_output=True)
lo, hi = 0.1, 2.0
for _ in range(80):
    c = 0.5*(lo + hi); s_ = vtx(c, 12.0)
    if s_.status != 0 or np.max(s_.y[0]) > 1: hi = c
    else: lo = c
s_ = vtx(0.5*(lo + hi)); rv = np.linspace(1e-3, 7, 300); fv, av = s_.sol(rv)
ax.plot(rv, fv, color="#1f77b4", lw=2.5, label="field strength |φ|/v (zero at the core)")
ax.plot(rv, av, color="#2ca02c", lw=2.5, label="enclosed flux / one flux unit")
ax.plot(rv, (1 - fv**2)/2/0.5, color="#d62728", lw=2.5, label="magnetic field (relative to centre)")
ax.text(3.0, 0.2, f"tension = π v² per flux unit\n(solved: {rk['K_vortex_tension_over_pi']:.4f} × π v²)", fontsize=9.5)
ax.set_xlabel("distance from the tube's axis (in units of 1/(e v))"); ax.set_ylim(0, 1.1)
ax.legend(fontsize=8.5, loc="center right", bbox_to_anchor=(1, 0.55))
ax.set_title("D. The abelian tube: a flux tube with a fixed tension per length.\nIts fundamental form is a cosmic string (K4)", fontsize=10.5)
ax = axs[1, 1]
lev = [(1.4e4, "LHC collision energy", "k"), (3.3e3, "Kaluza–Klein gluons ruled out below\n3.3 TeV (CMS, warped model)", "k"),
       (rk["K_string_eta_max_GeV_CMB"], "cosmic strings: breaking scale\nbelow 1.8e15 GeV (CMB)", "#8c564b"),
       (rk["K_string_eta_max_GeV_PTA"], "below 5e13 GeV (pulsar timing,\nNambu–Goto loops)", "#8c564b"),
       (rk["K_KKgluon_GeV"], "first copy of W or gluon on our circle", "#d62728"), (2*rk["K_KKgluon_GeV"], "second copy", "#d62728"), (1.22e19, "Planck energy", "k")]
for v, lab_, c_ in lev:
    ax.plot([0], [math.log10(v)], "o", color=c_, ms=7)
order = sorted(lev, key=lambda t: t[0]); ys = {}
yl = [3.0, 4.8, 13.2, 15.2, 17.2, 18.1, 19.2]
for (v, lab_, c_), yy in zip(order, yl):
    ax.plot([0.02, 0.22], [math.log10(v), yy], color="gray", lw=0.7)
    ax.text(0.25, yy, lab_, fontsize=8.5, va="center", color=c_)
ax.axvline(0, color="gray"); ax.set_xlim(-0.2, 1.7); ax.set_ylim(1.8, 20); ax.set_xticks([])
ax.set_ylabel("energy (log10 GeV)")
ax.set_title("E. Copies of W or gluons, and cosmic strings, against\nwhere experiments reach (K3, K4)", fontsize=10.5)
ax = axs[1, 2]
th_ = np.linspace(0.05, np.pi/2, 200)
ax.plot(np.degrees(th_), np.cos(th_), color="#9467bd", lw=2.5, label="mixing between our photon and the\nhidden circle's photon, ε = cos θ")
ins = ax.inset_axes([0.5, 0.45, 0.45, 0.45])
t0 = np.radians(60); ins.add_patch(matplotlib.patches.Polygon([[0, 0], [1, 0], [1 + 0.6*math.cos(t0), 0.6*math.sin(t0)], [0.6*math.cos(t0), 0.6*math.sin(t0)]], closed=True, fill=False, lw=2))
ins.annotate("", (1, 0), (0, 0), arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=2)); ins.text(0.45, -0.12, "our circle", fontsize=8, color="#1f77b4")
ins.annotate("", (0.6*math.cos(t0), 0.6*math.sin(t0)), (0, 0), arrowprops=dict(arrowstyle="->", color="#d62728", lw=2)); ins.text(-0.12, 0.45, "hidden\ncircle", fontsize=8, color="#d62728", ha="right")
ins.text(0.18, 0.08, "θ", fontsize=10); ins.set_xlim(-0.6, 1.5); ins.set_ylim(-0.25, 0.7); ins.axis("off")
ax.set_xlabel("angle θ between the two circles (degrees)"); ax.set_ylabel("mixing ε")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("F. A tilted torus: the hidden circle's photon mixes with ours,\ngiving its copies a small electric charge (K1c, K1d)", fontsize=10.5)
fig.suptitle("h13. The other unknown forces on the map, one cell at a time", fontsize=12.5)
fig.tight_layout(); fig.savefig("figures/h13_other_unknown_forces.png", dpi=150); plt.close(fig)
