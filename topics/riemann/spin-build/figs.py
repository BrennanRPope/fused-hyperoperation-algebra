import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("spinbuild.json")); M = np.load("mixture.npz"); P = np.load("primon.npz")
blue, red, grey, green, purple, orange = "#2b6cb0", "#c53030", "#718096", "#2f855a", "#6b46c1", "#dd6b20"
Z = np.array([2 * g for g in P["gam"]])

# ------------------------------------------------------------------ Figure A: simple spins and a ferromagnetic link
fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
a = ax[0]
a.plot(Z[Z < 200], 0 * Z[Z < 200], "o", color="k", ms=5, label="zeros of ζ's spin (z = 2γ)")
x1 = Z[0]
for m, col, lab in [(1, "k", "x₁ = 28.27"), (2, orange, "2x₁ (needed by a uniform spin)"), (3, red, "3x₁ (needed by an Ising spin)")]:
    a.axvline(m * x1, color=col, ls="--" if m > 1 else ":", lw=1.4, label=lab)
r = d["three_x_rule"][0]
a.annotate(f"misses by {r['miss_2x']:.2f}", (2 * x1, 0.15), (2 * x1 + 4, 0.55), color=orange, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=orange))
a.annotate(f"misses by {r['miss_3x']:.2f}", (3 * x1, -0.15), (3 * x1 + 4, -0.55), color=red, fontsize=9,
           arrowprops=dict(arrowstyle="->", color=red))
a.set_ylim(-1, 1); a.set_yticks([]); a.set_xlim(0, 200); a.set_xlabel("z")
a.set_title("A1  Independent simple spins can't build it\nfirst zero x₁ would force 2x₁ or 3x₁ to be a zero")
a.legend(fontsize=8, loc="upper right")

a = ax[1]
aa, bb = 1.0, 0.37
Js = np.linspace(-0.8, 0.8, 81); cnt = []
zs = np.linspace(0, 40, 20001)
for J in Js:
    f = math.exp(J) * np.cos((aa + bb) * zs) + math.exp(-J) * np.cos((aa - bb) * zs)
    cnt.append(int(np.sum(f[:-1] * f[1:] < 0)))
a.plot(Js, cnt, color=blue, lw=2)
a.axhline((aa + bb) * 40 / math.pi, color=grey, ls="--", lw=1); a.text(-0.78, (aa + bb) * 40 / math.pi + 0.4, "all zeros real", color=grey, fontsize=9)
a.axvline(0, color="k", lw=0.8)
a.text(0.05, 12, "ferromagnetic\nJ ≥ 0", color=green, fontsize=10); a.text(-0.75, 12, "antiferromagnetic\nJ < 0", color=red, fontsize=10)
a.set_xlabel("coupling J between two Ising spins (amplitudes 1 and 0.37)"); a.set_ylabel("real zeros in 0 < z < 40")
a.set_title("A2  One ferromagnetic link: all zeros real for J ≥ 0\nand they no longer follow the 3x rule")

a = ax[2]
u = M["u"]
a.plot(u, M["Phi"], color="k", lw=2.4, label="ζ's spin Φ(u) (positive, symmetric)")
for n, col in zip([1, 2, 3, 4, 5], [red, orange, green, blue, purple]):
    a.plot(u, M[f"term{n}"], color=col, lw=1.1, label=f"n = {n} term: n^(−1/2) φ₁(u + ½ ln {n})")
a.axhline(0, color=grey, lw=0.6)
a.set_xlim(-0.8, 0.8)
a.set_xlabel("spin value u"); a.set_ylabel("weight")
a.set_title("A3  ζ's spin is a SUM over n of shifted copies of one bump\nthe bump dips below zero; only the total is a spin")
a.legend(fontsize=7.5, loc="upper left")
plt.tight_layout(); plt.savefig("figures/figA_simple_spins.png", dpi=130); plt.close()

# ------------------------------------------------------------------ Figure B: the symmetric primon gas
fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
a = ax[0]
sg = np.linspace(0.5, 2.0, 300)
a.fill_between(sg, -(sg - 1).clip(0), (sg - 1).clip(0), color=green, alpha=0.25, label="the spin measure exists here (|Re h| < σ − 1)")
a.plot(sg, sg - 0.5, color=blue, lw=2, label="zeros of ζ's factors if RH: |Re h| = σ − ½")
a.plot(sg, -(sg - 0.5), color=blue, lw=2)
beta = 0.7
a.plot(sg, np.abs(sg - beta), color=red, ls="--", lw=1.6, label=f"a zero off the line (Re ρ = {beta}) would sit here")
a.plot(sg, -np.abs(sg - beta), color=red, ls="--", lw=1.6)
a.axvspan(0.5, 1.0, color=grey, alpha=0.12)
a.axvline(1.0, color="k", lw=1); a.text(1.02, -1.0, "σ = 1: phase transition\n(pole of ζ)", fontsize=9)
a.text(0.55, -1.25, "no measure here\n(½ < σ ≤ 1)", fontsize=9, color=grey)
a.plot([0.5], [0], "k*", ms=13); a.text(0.53, -0.42, "RH: blue lines\nmeet at σ = ½", fontsize=9)
a.set_xlabel("inverse temperature σ"); a.set_ylabel("Re h (field)")
a.set_xlim(0.5, 2.0); a.set_ylim(-1.5, 1.5)
a.set_title("B1  Symmetric primon gas: RH ⇔ Lee–Yang at σ = ½\nbut the gas only exists for σ > 1")
a.legend(fontsize=7.5, loc="upper left")

a = ax[1]
t = P["t"]
a.plot(t, P["zeta2"], color="k", lw=2, label="log |ζ(½ + it)|²  (has zeros)")
for Pv, col in [(10, orange), (100, green), (10000, blue)]:
    a.plot(t, P[f"logZ_{Pv}"], color=col, lw=1.1, label=f"primes ≤ {Pv}: log Z_P(it), never −∞")
for g in P["gam"]:
    if 10 < g < 50: a.axvline(g, color=red, lw=0.6, ls=":")
a.set_ylim(-14, 6)
a.set_xlabel("imaginary field h = i t"); a.set_ylabel("log of the partition function")
a.set_title("B2  Every finite gas is zero-free on the line\nzeros (red) appear only in the limit, which diverges at σ = ½")
a.legend(fontsize=7.5, loc="lower left")

a = ax[2]
rows = ["independent\nIsing / uniform", "ζ's spin Φ\n(Newman)", "symmetric primon\ngas, σ = ½", "Knauf's\nspin chain"]
cols = ["symmetric\nspins", "ferromagnetic\nlinks", "field linear\nin the total", "measure\nexists"]
grid = [["✓", "✓", "✓", "✓"], ["✓", "?", "✓", "✓"], ["✓", "✓", "✓", "✗"], ["✓", "✓", "✗", "✓"]]
colr = {"✓": green, "✗": red, "?": orange}
for i, row in enumerate(grid):
    for j, v in enumerate(row):
        a.add_patch(plt.Rectangle((j, -i), 1, 1, color=colr[v], alpha=0.3))
        a.text(j + 0.5, -i + 0.5, v, ha="center", va="center", fontsize=20, color=colr[v])
for i, rlab in enumerate(rows):
    a.text(-0.1, -i + 0.5, rlab, ha="right", va="center", fontsize=9)
for j, clab in enumerate(cols):
    a.text(j + 0.5, 1.15, clab, ha="center", va="bottom", fontsize=9)
a.text(2, -3.45, "row 1 fails the 3x rule instead;  ? = this is RH", ha="center", fontsize=9, color=grey)
a.set_xlim(-1.6, 4.1); a.set_ylim(-3.6, 1.9); a.axis("off")
a.set_title("B3  The Lee–Yang theorem needs all four.\nEach known ζ model misses a different one")
plt.tight_layout(); plt.savefig("figures/figB_primon_gas.png", dpi=130); plt.close()
print("ok")
