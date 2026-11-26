"""Pictures for parabolic-fingerprints/README.md. Reads identities.json, horns_*.json; height-flow numbers for e^x - 1."""
import json, os, glob, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
G, R, B, GR, O, P = "#2b8a3e", "#c92a2a", "#1c7ed6", "#868e96", "#e8590c", "#7048e8"
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name), dpi=130); plt.close(fig); print("wrote", name)
H = {os.path.basename(f)[6:-5]: json.load(open(f)) for f in glob.glob(os.path.join(HERE, "horns_*.json"))}

# ---- fig 1: the fingerprint map
pts = [  # (resit, size, label, colour)
    (1/3, 64.1531845185, "rung step e^x − 1\n= critical Poisson branching\n= tetration at the kiss base e^(1/e)", R),
    (1/3, H["poissonlike"]["size_float"], "law (1/3, 1/2, 0, 1/6)\n(same σ² and residue as Poisson)", O),
    (1.0, H["binary"]["size_float"], "binary splitting\n= diamond lattice at Q = 0\n= Mandelbrot cusp = logistic r = 1", B),
    (H["logistic3_0"]["resit_float"], H["logistic3_0"]["size_float"], "type-I intermittency\n(logistic³ at r = 1 + √8)", P),
    (1/9, H["negbin2"]["size_float"], "negative binomial r = 2", G),
    (2/3, H["binom3"]["size_float"], "binomial(3, 1/3)", G),
]
fig, ax = plt.subplots(figsize=(10.5, 5.6))
offs = [(0.06, -0.15), (-0.36, 0.35), (-0.5, -1.05), (0.06, -0.55), (0.04, 0.35), (0.04, -0.65)]
for (rs, sz, lab, c), (dx, dy) in zip(pts, offs):
    ax.plot(rs, math.log10(sz), "o", color=c, ms=10)
    ax.annotate(lab, (rs, math.log10(sz)), xytext=(rs + dx, math.log10(sz) + dy), fontsize=8, color=c)
ax.plot(0, -1.6, "v", color=GR, ms=12)
ax.annotate("Möbius x/(1 − x): horn size 0\n= geometric branching\n= Hu–Rudnick fixed point of intermittency\n(the only class with a continuous flow)",
            (0, -1.6), xytext=(0.05, -1.45), fontsize=8, color=GR)
ax.set_xlabel("iterative residue (formal fingerprint; scheme-free two-loop number)")
ax.set_ylabel("log10 horn-map size √|A₁A₋₁| (analytic fingerprint)")
ax.set_xlim(-0.1, 1.45); ax.set_ylim(-2, 9)
ax.set_title("Fig 1. Every system here creeps into its fixed point like 1/n, but the two fingerprints sort them apart.\n"
             "Only critical Poisson branching and the kiss-base tetration share the rung step's dot.", fontsize=10)
save(fig, "fig1-fingerprint-map.png")

# ---- fig 2: residue across the binomial / Poisson / negative binomial family
fig, ax = plt.subplots(figsize=(9, 4.6))
sig = np.linspace(0.45, 2.2, 300)
ax.plot(sig, (2 - sig)/(3*sig), color="k", lw=1)
named = [(0.5, "binomial(2, ½)\n= binary splitting", B), (2/3, "binomial(3, ⅓)", G), (1.0, "Poisson(1)\n= rung step", R),
         (1.5, "neg. binomial r = 2", G), (4/3, "neg. binomial r = 3", G), (2.0, "geometric\n= Möbius", GR)]
for sg, lab, c in named:
    ax.plot(sg, (2 - sg)/(3*sg), "o", color=c, ms=8)
    ax.annotate(lab, (sg, (2 - sg)/(3*sg)), xytext=(sg + 0.04, (2 - sg)/(3*sg) + 0.06), fontsize=8, color=c)
ax.axhline(0, color=GR, lw=0.5)
ax.axvline(1, color=R, lw=0.5, ls=":")
ax.set_ylim(-0.09, 1.2)
ax.text(0.47, -0.065, "under-dispersed (binomial)", fontsize=8, color=GR)
ax.text(1.35, -0.065, "over-dispersed (negative binomial)", fontsize=8, color=GR)
ax.set_xlabel("offspring variance σ² (mean fixed at 1)")
ax.set_ylabel("iterative residue")
ax.set_title("Fig 2. Along the binomial–Poisson–negative-binomial family the residue is (2 − σ²)/(3σ²);\n"
             "the rung step's 1/3 is exactly the Poisson point where the two halves meet", fontsize=10)
save(fig, "fig2-branching-residue.png")

# ---- fig 3: where the residue shows up: the survival chance of a family line
mp.mp.dps = 30
laws = [("geometric (residue 0)", lambda s: 1/(2 - s), 2, 0, GR),
        ("Poisson = rung step (residue 1/3)", lambda s: mp.exp(s - 1), 1, mp.mpf(1)/3, R),
        ("binary splitting (residue 1)", lambda s: (s*s + 1)/2, 1, 1, B)]
fig, ax = plt.subplots(figsize=(9, 4.6))
fits = {}
for lab, f, s2, rho, c in laws:
    q, ns, ys, x = [], [], [], mp.mpf(0)
    n_max = 20000
    for n in range(1, n_max + 1):
        x = f(x)
        if n >= 50 and (n % 50 == 0):
            qn = 1 - x
            ns.append(n); ys.append(float(n*n*(qn - mp.mpf(2)/(s2*n))))
    ln = np.log(ns); ys = np.array(ys)
    sel = np.array(ns) >= 2000
    nn = np.array(ns, dtype=float)[sel]
    Amat = np.vstack([ln[sel], np.ones(sel.sum()), ln[sel]/nn, 1/nn, ln[sel]**2/nn]).T   # next corrections ~ (log n)^j / n
    slope = np.linalg.lstsq(Amat, ys[sel], rcond=None)[0][0]
    fits[lab] = (slope, float(-2*rho/s2))
    ax.plot(ln, ys, color=c, label=f"{lab}: slope {slope:.4f} (predicted {float(-2*rho/s2):.4f})")
ax.set_xlabel("log n (generations)")
ax.set_ylabel("n² × (survival chance − 2/(σ² n))")
ax.legend(fontsize=8)
ax.set_title("Fig 3. All three lines survive like 2/(σ² n); the residue is the slope of the next correction, −2·residue/σ²", fontsize=10)
save(fig, "fig3-survival-correction.png")
json.dump({k: {"slope": v[0], "predicted": v[1]} for k, v in fits.items()}, open(os.path.join(HERE, "survival_fits.json"), "w"), indent=1)

# ---- fig 4: doubling (the intermittency renormalization) washes both fingerprints out
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
ks = np.arange(0, 9)
for rs, lab, c in [(1/3, "rung step", R), (1.0, "binary / diamond Q = 0", B), (51/49, "intermittency", P)]:
    axs[0].plot(ks, rs/2.0**ks, "o-", color=c, label=lab)
axs[0].set_yscale("log"); axs[0].set_xlabel("doublings k (map f^(2^k), rescaled)"); axs[0].set_ylabel("iterative residue")
axs[0].legend(fontsize=8)
axs[0].set_title("(a) each doubling halves the residue:\nan irrelevant direction with eigenvalue 1/2", fontsize=10)
A1 = 0.0890584364122
for Y, c in [(1.0, R), (2.0, "#f08c00")]:
    axs[1].plot(ks, [math.log10(A1/2**k) - 2*math.pi*(2**k)*Y/math.log(10) for k in ks], "o-", color=c,
                label=f"rung step, height Im W = {Y:g} above the axis")
axs[1].set_xlabel("doublings k"); axs[1].set_ylabel("log10 size of the first horn wave")
axs[1].legend(fontsize=8)
axs[1].set_title("(b) the horn wave survives only at frequency 2^k,\nso it dies faster than any power (from h₂ = h(2W − c)/2 + c)", fontsize=10)
fig.suptitle("Fig 4. Under doubling every map here flows to the Möbius fixed point (Hu–Rudnick): the fingerprints are the non-universal part", fontsize=10)
save(fig, "fig4-doubling.png")
