import math, json
import numpy as np, matplotlib
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
    (0, 1): ("open", "EMPTY: a long-range non-abelian force.\nOnly stays long-range with many matter\nfields (SU(3): 17+ flavours); otherwise\nit confines. Searched as 'dark QCD'."),
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
