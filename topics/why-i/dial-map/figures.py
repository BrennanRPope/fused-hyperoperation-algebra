"""Pictures for topics/why-i/dial-map (Nov 25, 2026). Run: python3 figures.py"""
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
def box(ax, x, y, text, col, w=1.9, h=0.62, fs=8.5):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.04", fc="white", ec=col, lw=1.8))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color="k")

# m1. the Bronstein (cGh) cube: what is proved along each edge
P = {"": (0, 0), "c": (4, 0), "G": (1.6, 1.4), "h": (0, 3.2), "cG": (5.6, 1.4), "ch": (4, 3.2), "Gh": (1.6, 4.6), "cGh": (5.6, 4.6)}
lab = {"": "Newtonian\nmechanics", "c": "special\nrelativity", "G": "Newtonian\ngravity", "h": "quantum\nmechanics",
       "cG": "general\nrelativity", "ch": "quantum field\ntheory", "Gh": "non-relativistic\nquantum gravity", "cGh": "theory of\neverything?"}
edges = [("", "c", "proved"), ("", "h", "proved"), ("", "G", "limit"), ("c", "cG", "limit"), ("G", "cG", "limit"),
         ("c", "ch", "proved"), ("h", "ch", "proved"), ("h", "Gh", "open"), ("G", "Gh", "open"), ("cG", "cGh", "open"),
         ("ch", "cGh", "open"), ("Gh", "cGh", "open")]
sty = {"proved": dict(color=GREEN, lw=3), "limit": dict(color=BLUE, lw=2, ls="-."), "open": dict(color=GREY, lw=1.5, ls=":")}
fig, ax = plt.subplots(figsize=(10, 6.6)); ax.axis("off")
for a, b, s in edges:
    (x1, y1), (x2, y2) = P[a], P[b]
    if "G" in b and "G" not in a:      # the G edges: highlight
        ax.plot([x1, x2], [y1, y2], color=RED, lw=9, alpha=0.18, solid_capstyle="round")
    ax.plot([x1, x2], [y1, y2], **sty[s])
for kk, (x, y) in P.items():
    box(ax, x, y, lab[kk], GREEN if kk in ("", "c", "h", "ch") else (BLUE if kk in ("G", "cG") else GREY), w=1.55, h=0.7)
ax.annotate("", xy=(5.2, -0.6), xytext=(-0.3, -0.6), arrowprops=dict(arrowstyle="->", color="k")); ax.text(2.4, -0.85, "1/c on", ha="center")
ax.annotate("", xy=(-0.95, 3.6), xytext=(-0.95, -0.2), arrowprops=dict(arrowstyle="->", color="k")); ax.text(-1.15, 1.7, "hbar on", rotation=90, va="center")
ax.annotate("", xy=(3.5, 1.15), xytext=(2.3, 0.3), arrowprops=dict(arrowstyle="->", color="k")); ax.text(2.55, 0.85, "G on", rotation=35)
ax.plot([], [], **sty["proved"], label="downhill limit proved, and uphill is a stable deformation (Inonu-Wigner; Vey, Bayen et al.; Mendes)")
ax.plot([], [], **sty["limit"], label="downhill limit proved, uphill no stability theorem (Ehlers 1981 for c -> infinity in gravity)")
ax.plot([], [], **sty["open"], label="open: quantum gravity (Bronstein 1936)")
ax.plot([], [], color=RED, lw=9, alpha=0.18, label="the G edges (red glow): no edge going up is a forced, stable step; the uniqueness results fix the FORM of gravity, not G")
ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.02), fontsize=8.3, frameon=False)
ax.set_xlim(-1.6, 6.8); ax.set_ylim(-1.0, 5.3)
ax.set_title("The Bronstein-Zelmanov-Okun cube: a classification by which constants are switched on\n(G edges highlighted)")
save(fig, "m1-bronstein-cube")

# m2. our dial lattice, layer by layer, with the laws each step switches on
fig, ax = plt.subplots(figsize=(13, 7.6)); ax.axis("off"); ax.set_xlim(-0.6, 13.5); ax.set_ylim(-0.6, 8.3)
nodes = {
 "add": (1.2, 0.5, "all dials off: plain +\nPoisson e^x, Shannon, Galilean v-addition,\nsuperposition", GREY),
 "k": (1.2, 2.5, "k on: one wall at -1/k\nat most m copies / multisets, inclusion-exclusion,\nTsallis, primon gases (k = -1, +1)", ORANGE),
 "ki": (1.2, 4.5, "k imaginary: unit phases\ncheckerboard amplitudes (k = -i),\ngauge phases", ORANGE),
 "kd2": (4.6, 2.5, "k, d with kappa < 0: two walls\nspecial-relativity velocity addition\n(d = -1/c^2)", ORANGE),
 "kd0": (4.6, 4.5, "k, d with kappa > 0: no wall\nrotation / tangent addition", ORANGE),
 "torus": (4.6, 6.5, "beyond rational laws: torus\n(elliptic formal group; pinch -> walls)", GREY),
 "k2": (8.2, 2.5, "matrices + kappa2 = -1/c^2\nLorentz boosts, Thomas-Wigner rotation,\nMaxwell invariance, Doppler lift", GREEN),
 "k1": (8.2, 4.5, "matrices + kappa1 (curvature)\nde Sitter / anti-de Sitter, Newton-Hooke;\ntwist = curvature x area", GREEN),
 "hbar": (11.6, 2.5, "order dial hbar on\nHeisenberg phase, Weyl / Moyal;\ni and hbar put in by hand", GREEN),
 "zoom": (11.6, 0.5, "zoom dial s on (new input: size)\nMoran dimension, complex dimensions", PURPLE),
 "G": (11.6, 6.5, "G: no forcing dial\nbinding = k = -G/(r c^2) only to 1st order;\nbootstrap sums to another ruler", RED),
}
for kk, (x, y, t, col) in nodes.items(): box(ax, x, y, t, col, w=3.15, h=1.05, fs=7.8)
arrows = [("add", "k", "coordinate change"), ("k", "ki", "k -> i k"), ("k", "kd2", "add d"), ("ki", "kd0", "add d"), ("kd0", "torus", ""),
          ("kd2", "k2", "compose non-commuting\nmotions: stable step"), ("kd0", "k1", "stable step"), ("k2", "k1", ""),
          ("add", "zoom", "new input"), ("k2", "hbar", "functions + order:\nstable step")]
for a, b, t in arrows:
    x1, y1 = nodes[a][:2]; x2, y2 = nodes[b][:2]
    if y1 == y2:   # horizontal: span the gap between boxes, label underneath
        sx = 1 if x2 > x1 else -1
        ax.annotate("", xy=(x2 - sx*1.62, y2), xytext=(x1 + sx*1.62, y1), arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
        if t: ax.text((x1 + x2)/2, y1 - 0.62, t, fontsize=7.5, color=GREY, ha="center", va="top")
    elif x1 == x2:
        sy = 1 if y2 > y1 else -1
        ax.annotate("", xy=(x2, y2 - sy*0.55), xytext=(x1, y1 + sy*0.55), arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
        if t: ax.text(x1 + 0.1, (y1 + y2)/2, t, fontsize=7.5, color=GREY, ha="left", va="center")
    else:
        ax.annotate("", xy=(x2 - 1.62, y2 - 0.3), xytext=(x1 + 1.62, y1 + 0.3), arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
        if t: ax.text((x1 + x2)/2, (y1 + y2)/2, t, fontsize=7.5, color=GREY, ha="center")
for x, t, col in [(1.2, "scalar dial k\n(not a deformation)", ORANGE), (4.6, "scalar dials k, d\n(still one kappa)", ORANGE),
                  (8.2, "Cayley-Klein dials\n(matrices; stable deformations)", GREEN), (11.6, "other dials", GREY)]:
    ax.text(x, 7.75, t, ha="center", fontsize=9.5, color=col, weight="bold")
ax.set_title("The framework's dials, layer by layer: what each one switches on (orange: change of ruler; green: published stable deformation)")
save(fig, "m2-dial-lattice")

# m3. line-up: cube axes, Cayley-Klein dials, our dials
rows = [("1/c", "kappa2 = -1/c^2 (signature)", "d at k = 0, or k = 2/c, for one boost; full match only with matrices", GREEN, "match; stable step"),
        ("G", "none (kappa1 is curvature / Lambda, not G)", "k = -G/(r c^2) for binding, first order only", RED, "no forcing match"),
        ("(Planck scale, DSR)", "kappa1 of momentum space", "affine stretch-and-shift = kappa-Poincare (known)", BLUE, "match; hypothesis, not forced"),
        ("hbar", "none (Cayley-Klein is classical kinematics)", "order dial (Moyal), put in by hand", GREEN, "match; stable step"),
        ("(Lambda, not on the cube)", "kappa1 (curvature)", "k for one time-shift motion; full match with matrices", BLUE, "match; stable (de Sitter rigid)"),
        ("none", "none", "k alone, d alone: rulers; zoom: a new size input", ORANGE, "not deformations")]
fig, ax = plt.subplots(figsize=(13, 4.5)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for x, h in [(0.01, "Bronstein cube axis"), (0.17, "Cayley-Klein dial"), (0.45, "our dial"), (0.83, "verdict")]:
    ax.text(x, 0.93, h, weight="bold", fontsize=10)
for i, (a, b, c, col, v) in enumerate(rows):
    y = 0.80 - i*0.135
    ax.text(0.01, y, a, fontsize=9.5); ax.text(0.17, y, b, fontsize=9.5); ax.text(0.45, y, c, fontsize=9.5)
    ax.text(0.83, y, v, fontsize=9.5, color=col, weight="bold")
ax.set_title("Lining up the dials")
save(fig, "m3-line-up")

# m4. stability test: dimension of H^2 (room to deform)
names = ["Galilei", "Bargmann\n(with mass)", "Poincare", "de Sitter", "anti-de Sitter"]
h11 = [5, 6, 2, 0, 0]; h31 = [4, None, 1, 0, 0]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
x = np.arange(5)
ax[0].bar(x - .18, h11, .34, color=[ORANGE, ORANGE, BLUE, GREEN, GREEN], label="1+1")
ax[0].bar([0 + .18, 2 + .18, 3 + .18, 4 + .18], [4, 1, 0, 0], .34, color=[ORANGE, BLUE, GREEN, GREEN], alpha=.55, hatch="//", label="3+1")
for xi, v in zip(x, h11): ax[0].text(xi - .18, v + .1, str(v), ha="center", fontsize=8)
for xi, v in zip([0, 2, 3, 4], [4, 1, 0, 0]): ax[0].text(xi + .18, v + .1, str(v), ha="center", fontsize=8)
ax[0].set_xticks(x); ax[0].set_xticklabels(names, fontsize=8.5); ax[0].set_ylabel("dim H^2(g, g): room to deform")
ax[0].legend(fontsize=8); ax[0].set_title("turn on 1/c^2, then curvature: the room to deform runs out\n(0 = rigid = stable)")
ax[1].bar([0, 1, 2], [9, 5, 0], color=[ORANGE, BLUE, GREY])
for xi, v in zip([0, 1, 2], [9, 5, 0]): ax[1].text(xi, v + .15, str(v), ha="center", fontsize=8)
ax[1].set_xticks([0, 1, 2]); ax[1].set_xticklabels(["x, p, 1 commuting", "Heisenberg\n[x, p] = hbar", "one-parameter law\n(our k or d)"], fontsize=8.5)
ax[1].set_title("quantum side: the 3-dim piece is not the whole story (Moyal is);\nour scalar k has nothing to deform")
save(fig, "m4-stability")

# m5. scorecard
rows = [("What the cube is", "a classification by which of c, G, hbar a theory keeps (Bronstein, Zelmanov, Okun)"),
        ("What is proved", "limits downhill; uphill, relativity and quantum mechanics are stable deformations"),
        ("Not proved", "anything about G uphill; the quantum-gravity corner is open"),
        ("Gravity forcing", "bootstrap and Lovelock fix GR's form from spin-2 + universal coupling; not G, not via algebra"),
        ("Our matches", "kappa2 = -1/c^2 (relativity) and the order dial hbar (Moyal): both stable steps"),
        ("Our non-matches", "k and d alone are changes of ruler; zoom is a new input; no G dial"),
        ("Brennan's layers", "right shape: each switched-on dial adds a layer of laws, but only green steps are deformations")]
fig, ax = plt.subplots(figsize=(12.5, 4.0)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for i, (q, a_) in enumerate(rows):
    y = 0.92 - i*0.135; col = [GREY, GREEN, RED, RED, GREEN, ORANGE, BLUE][i]
    ax.add_patch(plt.Rectangle((0.0, y - 0.055), 0.012, 0.11, color=col))
    ax.text(0.025, y, q, va="center", fontsize=10, weight="bold"); ax.text(0.2, y, a_, va="center", fontsize=10, color=col)
ax.set_title("the dial map against the cube")
save(fig, "m5-scorecard")

# m6. the gravity forcing results: what each assumes, what it forces, and whether a chain can say it
rows = [("Gupta 1954, Kraichnan 1955,\nFeynman lectures, Deser 1970", "a massless spin-2 field that couples\nto all energy, its own included",
         "the full nonlinear Einstein equations\n(Deser: in one step)", "field theory", "no: needs fields on spacetime\n(scalar toy only, G3)", RED),
        ("Weinberg 1964-65\n(soft gravitons)", "Lorentz-invariant S-matrix with a\nmassless spin-2 particle",
         "one universal coupling = equivalence\nprinciple", "S-matrix", "no", RED),
        ("Boulware-Deser 1975", "same as the bootstrap, any order", "GR is the unique consistent completion", "field theory", "no", RED),
        ("Lovelock 1971", "4D, metric only, 2nd-order equations,\nsymmetric and divergence-free",
         "Einstein tensor + Lambda g\n(G and Lambda left free)", "differential geometry", "no", RED),
        ("Bacry-Levy-Leblond 1968;\nLevy-Nahas 1967", "a kinematic Lie algebra",
         "Poincare deforms only to de Sitter / AdS:\nLambda is the stable dial", "algebra (stability)", "yes: matrix chains carry kappa1\n= Lambda/3 (checked, S3, S8)", GREEN),
        ("MacDowell-Mansouri 1977", "a de Sitter gauge connection, broken\nto Lorentz", "GR + Lambda + topological term;\nthe combination G x Lambda set by the gauge coupling", "gauge theory", "partly: the dS algebra yes,\nthe field action no", ORANGE),
        ("kappa-Poincare (Lukierski et al.\n1991); DSR", "Hopf-algebra deformation with a\nPlanck energy kappa", "curved (de Sitter) momentum space", "hypothesis (Fermi:\nE_QG > 1.2 E_Planck)", "yes: it IS the affine\nstretch-and-shift law (G4)", BLUE)]
fig, ax = plt.subplots(figsize=(15, 7.4)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
hx = [0.0, 0.19, 0.42, 0.63, 0.76]
for x, h in zip(hx, ["result", "assumes", "forces", "kind of forcing", "can a chain of weighted circles say it?"]):
    ax.text(x, 0.96, h, weight="bold", fontsize=10)
for i, (a, b, c, d, e, col) in enumerate(rows):
    y = 0.865 - i*0.128
    ax.add_patch(plt.Rectangle((0.745, y - 0.05), 0.006, 0.1, color=col))
    for x, t in zip(hx, [a, b, c, d]): ax.text(x, y, t, fontsize=8.6, va="center")
    ax.text(0.76, y, e, fontsize=8.6, va="center", color=col, weight="bold")
ax.set_title("What forces gravity's form, and whether the framework can express it (red: no; green / blue: yes, known structure)")
save(fig, "m6-gravity-forcing")

# m7. two concrete gravity matches, side by side
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
alpha = 0.5; m = np.linspace(0, 1.6, 200)
Linv = lambda y: (-1 + np.sqrt(1 + 4*alpha*y))/(2*alpha)
L = lambda x: x + alpha*x**2
adm = Linv(L(m) + L(m)); circ = 2*m - 2*alpha*m*m; newt = 2*m
ax[0].plot(m, newt, color=GREY, ls=":", label="no gravity: m + m")
ax[0].plot(m, adm, color=BLUE, lw=2.2, label="ADM shell law (exact self-energy)")
ax[0].plot(m, circ, color=ORANGE, lw=1.8, ls="--", label="weighted circle, k = -G/R")
ax[0].set_ylim(0, 3.3); ax[0].set_xlabel("each mass m (units with G/2R = 0.5)"); ax[0].set_ylabel("combined mass")
ax[0].set_title("binding energy: the circle law matches to 2nd order,\nthen turns over; the exact law never does")
ax[0].legend(fontsize=8.5)
kap = 1.0; p0 = np.linspace(-2.5, 2.5, 200)
ax[1].plot(p0, np.exp(-p0/kap), color=BLUE, lw=2.2, label="kappa-Poincare stretch e^{-p0/kappa}")
ax[1].plot(p0, 1 + 0*p0, color=GREY, ls=":", label="Poincare (kappa -> infinity): no stretch")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel("energy p0 (units of the Planck dial kappa)"); ax[1].set_ylabel("lift 1 + k a")
ax[1].set_title("the Planck dial: kappa-Poincare's momentum law is the\naffine stretch-and-shift law with 1 + k a = e^{-p0/kappa}")
ax[1].legend(fontsize=8.5)
save(fig, "m7-two-gravity-matches")
print("figures written")
