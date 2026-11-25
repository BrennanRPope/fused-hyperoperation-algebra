"""Pictures for topics/why-i/kappa-poincare (Nov 25, 2026). Run: python3 figures.py"""
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
circ = lambda a, b, k: a + b + k*a*b

# k1. the identification: massless momenta in kappa-Poincare compose by o_k with k = -1/kappa
K = 1.0; e = np.linspace(0, 4.5, 300); pm = lambda p0: K*(1 - np.exp(-p0/K))
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax[0].plot(e, e, color=GREY, ls=":", label="special relativity: |p| = E (no ceiling)")
ax[0].plot(e, pm(e), color=BLUE, lw=2.4, label="kappa-Poincare: |p| = kappa (1 - e^{-E/kappa})")
ax[0].axhline(K, color=RED, lw=1.2, ls="--"); ax[0].text(3.0, K + 0.06, "wall of the circle law at -1/k = kappa", color=RED, fontsize=9)
e1, e2 = 0.8, 1.3
for x, lab, col in [(e1, "a", ORANGE), (e2, "b", GREEN), (e1 + e2, "a o_k b", PURPLE)]:
    ax[0].plot([x, x], [0, pm(x)], color=col, lw=1.5); ax[0].plot(x, pm(x), "o", color=col); ax[0].text(x + 0.05, pm(x) - 0.12, lab, color=col)
ax[0].set_ylim(0, 1.6); ax[0].set_xlabel("energy E (units of the Planck energy kappa)"); ax[0].set_ylabel("momentum |p| of a massless particle")
ax[0].set_title("energies add; momenta compose by a + b - ab/kappa"); ax[0].legend(fontsize=8.5, loc="upper left")
a = np.linspace(-0.5, 0.95, 200)
for b, col in [(0.3, ORANGE), (0.6, GREEN), (0.9, PURPLE)]:
    ax[1].plot(a, circ(a, b, -1/K), color=col, lw=2, label=f"b = {b} kappa")
ax[1].plot(a, a + 0.6, color=GREY, ls=":", label="plain addition, b = 0.6 kappa")
ax[1].axhline(K, color=RED, ls="--", lw=1.2)
ax[1].set_xlabel("momentum a (units of kappa)"); ax[1].set_ylabel("composed momentum a o_k b, k = -1/kappa")
ax[1].set_title("the circle law with k = -1/kappa: nothing passes kappa"); ax[1].legend(fontsize=8.5)
save(fig, "k1-massless-momenta-are-a-circle-law")

# k2. placements on the repo's (k, d) plane
fig, ax = plt.subplots(figsize=(8.6, 6.4))
kk = np.linspace(-2.2, 2.2, 400)
ax.fill_between(kk, kk**2/4, 2.6, color=GREEN, alpha=.10); ax.fill_between(kk, -2.0, kk**2/4, color=ORANGE, alpha=.10)
ax.plot(kk, kk**2/4, color=PURPLE, lw=2.4)
ax.text(-2.1, 2.3, "kappa > 0: no wall (tangent addition)", color=GREEN, fontsize=9)
ax.text(-2.1, -1.85, "kappa < 0: two walls (velocity addition)", color=ORANGE, fontsize=9)
ax.text(1.25, 0.62, "k^2 = 4d:\nwalls merge\n(kappa = 0)", color=PURPLE, fontsize=9)
pts = [(0, -1, "velocity addition\n(k = 0, d = -1/c^2)", ORANGE), (1, 0, "repo relativity coordinate a = 2u/(1-u)\n(k = 1 at c = 1, d = 0)", ORANGE),
       (-1, 0, "kappa-Poincare massless momenta\n(k = -1/kappa, d = 0)", BLUE), (-1.2, 0.36, "Magueijo-Smolin energy law\n(k = -2l, d = l^2)", PURPLE),
       (0, 0, "plain addition: bicrossproduct\nenergies (p0 + q0)", GREY)]
offs = [(0.08, -0.3), (-0.75, -0.38), (-1.15, -0.38), (-0.95, 0.25), (0.08, 0.18)]
for (x, y, t, col), (dx, dy) in zip(pts, offs):
    ax.plot(x, y, "o", ms=9, color=col, mec="k"); ax.text(x + dx, y + dy, t, fontsize=8.3, color=col)
ax.axhline(0, color="k", lw=.5); ax.axvline(0, color="k", lw=.5)
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.0, 2.6); ax.set_xlabel("k"); ax.set_ylabel("d")
ax.set_title("where known laws sit on the repo's plane (a+b+kab)/(1-dab)\n(placements depend on the energy coordinate chosen; only the matrix law is physical)")
save(fig, "k2-placements-on-the-plane")

# k3. the literature line, and where the circle law sits
items = [(1947, "Snyder: momentum space is de Sitter;\ncomposition is non-associative", RED),
         (1991, "Lukierski, Ruegg, Nowicki, Tolstoy (1991); Lukierski, Nowicki, Ruegg (1992):\nkappa-Poincare (Planck dial kappa)", BLUE),
         (1994, "Majid-Ruegg: bicrossproduct basis;\np + e^{-p0/kappa} q  (= our affine chain)", BLUE),
         (2000, "Amelino-Camelia: DSR, a second\nobserver-independent scale", GREEN),
         (2002, "Magueijo-Smolin (DSR2); Judes-Visser:\nE(+)E' = U^-1(U(E)+U(E')) (= a point on our plane)", PURPLE),
         (2003, "Kowalski-Glikman (2002), Kowalski-Glikman-Nowak (2003): DSR bases are\ncoordinates on de Sitter momentum space", GREEN),
         (2011, "Amelino-Camelia, Freidel, Kowalski-Glikman, Smolin:\nrelative locality; torsion = order gap, curvature = non-associativity", GREEN)]
fig, ax = plt.subplots(figsize=(12, 5.4)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.plot([0.07, 0.07], [0.04, 0.93], color="k", lw=1)
for i, (yr, t, col) in enumerate(items):
    y = 0.9 - i*0.135
    ax.plot(0.07, y, "o", color=col, ms=8); ax.text(0.0, y, str(yr), va="center", fontsize=10, weight="bold", color=col)
    ax.text(0.1, y, t.replace("\n", " "), va="center", fontsize=9.5, color=col)
ax.set_title("kappa-Poincare and DSR in the literature (blue and purple: where our circle law already appears)")
save(fig, "k3-literature-line")

# k4. the geometry of each composition law
rows = [("scalar circle law a o_k b", "yes", "yes", "none", "none", "trivially: it's one number", GREY),
        ("affine / matrix chain\n= kappa-Poincare", "yes", "no", "1/kappa\n(= order gap)", "none", "yes: exactly (K2)", BLUE),
        ("Magueijo-Smolin energy law", "yes", "yes", "none", "none", "a point on the repo plane,\non k^2 = 4d (K3)", PURPLE),
        ("Snyder", "no", "no", "none at\nleading order", "yes", "no: chains are associative;\nonly the full rotation group (K6)", RED)]
fig, ax = plt.subplots(figsize=(13.5, 4.2)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
hx = [0.0, 0.24, 0.35, 0.46, 0.6, 0.7]
for x, h in zip(hx, ["composition law", "associative?", "commutative?", "torsion", "curvature", "can a chain of weighted circles give it?"]):
    ax.text(x, 0.92, h, weight="bold", fontsize=9.5)
for i, r in enumerate(rows):
    y = 0.74 - i*0.2
    for x, t in zip(hx, r[:6]): ax.text(x, y, t, fontsize=9, va="center", color=r[6] if x in (0.0, 0.7) else "k")
ax.set_title("momentum-space geometry read off each composition law (relative-locality dictionary)")
save(fig, "k4-geometry-of-each-law")

# k5. scorecard
rows = [("Already done?", "yes: the circle law, with k = -1/kappa, is how kappa-Poincare composes massless momenta (1991-94)", BLUE),
        ("The full version", "the affine stretch-and-shift chain is the whole kappa-Poincare momentum group AN(3)", BLUE),
        ("DSR2", "Magueijo-Smolin's energy law is a point on the repo's plane, on the merged-wall line", PURPLE),
        ("Out of reach", "Snyder's non-associative law: chains are always associative", RED),
        ("What's ours", "only placements: the k = -1/kappa reading, the wall = maximum momentum, DSR2 on k^2 = 4d", ORANGE),
        ("Caveat", "a one-number law fits every DSR model after a change of coordinate; the physics is in the matrix law", GREY)]
fig, ax = plt.subplots(figsize=(12.5, 3.7)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for i, (q, a_, col) in enumerate(rows):
    y = 0.9 - i*0.16
    ax.add_patch(plt.Rectangle((0.0, y - 0.055), 0.012, 0.11, color=col))
    ax.text(0.025, y, q, va="center", fontsize=10, weight="bold"); ax.text(0.17, y, a_, va="center", fontsize=10, color=col)
ax.set_title("the circle law against kappa-Poincare and doubly special relativity")
save(fig, "k5-scorecard")
print("figures written")
