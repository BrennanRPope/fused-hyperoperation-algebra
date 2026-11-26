"""Pictures for ns-same-shape (DRAFT). Run checks.py first (reads results.json)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
R = json.load(open("results.json")); out = "figures/"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def box(ax, x, y, w, h, text, fc, fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=fc, ec="k", lw=1))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, wrap=True)
def arrow(ax, a, b, txt="", c="k", ls="-", rad=0.0, fs=8):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14, color=c, ls=ls, lw=1.5,
                                 connectionstyle=f"arc3,rad={rad}"))
    if txt: ax.text((a[0]+b[0])/2, (a[1]+b[1])/2+0.04, txt, ha="center", fontsize=fs, color=c)

# ---- Fig 1: the shape, named ----
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
ax = axs[0]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.set_title("(a) Two sectors, one cross term")
box(ax, 0.03, 0.55, 0.38, 0.3, "right-handed waves\nkeep their own\nΣ k|u₊|²  (positive)", "#d8ecff")
box(ax, 0.59, 0.55, 0.38, 0.3, "left-handed waves\nkeep their own\nΣ k|u₋|²  (positive)", "#ffe1d6")
arrow(ax, (0.42, 0.66), (0.58, 0.66), "", "crimson"); arrow(ax, (0.58, 0.74), (0.42, 0.74), "", "crimson")
ax.text(0.5, 0.92, "cross term ε (right meets left)", ha="center", color="crimson", fontsize=9)
box(ax, 0.15, 0.08, 0.7, 0.3, "ε = 0: both extra quantities survive → smooth\n(Biferale–Titi 2013, proved)\nε > 0: only energy survives → open\n(our shell toy: blows up for every ε tried)", "#f4f4f4", 8.5)
ax = axs[1]; ax.set_title("(b) Energy is too weak at small scales")
l = np.logspace(0, 3, 50)
ax.loglog(l, l**0.0, "g-", lw=2, label="2D: energy scale-free (critical)")
ax.loglog(l, l**-1.0, "r-", lw=2, label="3D: energy shrinks as λ⁻¹ (supercritical)")
ax.set_xlabel("zoom factor λ (smaller scales →)"); ax.set_ylabel("energy left to control the zoomed flow")
ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=.3)
ax = axs[2]; ax.set_title("(c) Self-feeding growth: Y′ = Y² vs the rung line")
t = np.linspace(0, 0.495, 300)
ax.semilogy(t, 1/(0.5-t), "r-", lw=2, label="Y′ = Y² (stretching): wall at t = 1/2")
ax.semilogy(t, np.exp(np.log(2)*np.exp(t)), "b-", lw=2, label="Y′ = Y log Y (2D, rung 3): lives forever")
ax.semilogy(t, 2*np.exp(t), "g-", lw=2, label="Y′ = Y (rung 2)")
ax.set_xlabel("time (all start at Y = 2)"); ax.set_ylabel("vorticity size Y"); ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.suptitle("The shape of the Navier–Stokes danger: a cross term between two sectors breaks a second invariant, leaving only a supercritical energy, and the loss shows up as self-feeding growth", fontsize=10.5)
fig.tight_layout(); fig.savefig(out+"fig1_the_shape.png"); plt.close(fig)

# ---- Fig 2: scale entanglement ----
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("(a) A field as a chain of scale 'qubits'")
for i in range(8):
    ax.add_patch(Circle((0.08+0.12*i, 0.6), 0.04, fc="#cfe3ff", ec="k"))
    ax.text(0.08+0.12*i, 0.6, f"b{i+1}", ha="center", va="center", fontsize=8)
    if i < 7: ax.plot([0.12+0.12*i, 0.16+0.12*i], [0.6, 0.6], "k-")
ax.plot([0.5, 0.5], [0.45, 0.75], "r--", lw=2); ax.text(0.5, 0.79, "cut", color="r", ha="center")
ax.text(0.08, 0.48, "coarsest\nbit", ha="center", fontsize=8); ax.text(0.92, 0.48, "finest\nbit", ha="center", fontsize=8)
ax.text(0.5, 0.22, "Bit j of the position says which half at scale 2⁻ʲ.\nEntanglement across the cut = how much the fine\nscales depend on the coarse ones (Schmidt spectrum).\nThis is the 'interscale correlation' of Gourianov et al. 2022.",
        ha="center", fontsize=8.5)
ax = axs[1]; ax.set_title("(b) Entanglement across each scale cut (2²⁰ points)")
cols = ["#2a9d8f", "#264653", "#e76f51", "#f4a261", "#8d99ae", "#6a4c93"]
for (name, S), c in zip(R["entropy_by_cut_n20"].items(), cols):
    ax.plot(range(1, len(S)+1), S, "-o", ms=3, color=c, label=name)
ax.set_xlabel("cut after bit j (coarse → fine)"); ax.set_ylabel("entanglement entropy (bits)"); ax.grid(alpha=.3); ax.legend(fontsize=7)
ax = axs[2]; ax.set_title("(c) Rank needed (error 10⁻⁶) as scales are added")
g = R["max_rank_vs_n"]; nn = g["n"]
for (name, c) in zip([k for k in g if k != "n"], cols):
    ax.semilogy(nn, g[name], "-o", color=c, label=name)
ax.set_xlabel("number of scales n (grid 2ⁿ)"); ax.set_ylabel("largest bond dimension"); ax.grid(alpha=.3); ax.legend(fontsize=7)
fig.suptitle("Singular but self-similar profiles are as low-entanglement as smooth ones; what costs rank is disorder across scales, and viscosity caps it", fontsize=10.5)
fig.tight_layout(); fig.savefig(out+"fig2_scale_entanglement.png"); plt.close(fig)

# ---- Fig 3: multiplicative cascade ----
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
ax = axs[0]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("(a) Mandelbrot cascade = ∘₁ down a tree")
def tree(x, y, w, d):
    if d == 0: return
    for s in (-1, 1):
        ax.plot([x, x+s*w/2], [y, y-0.22], "k-", lw=1)
        ax.text(x+s*w/4-0.03*s, y-0.1, "W", fontsize=7, color="purple")
        tree(x+s*w/2, y-0.22, w/2, d-1)
tree(0.5, 0.95, 0.5, 3)
ax.text(0.5, 0.06, "mass of a leaf = Π W along its path,  W = 1 + x\nso each step is mult-add: (1+a)(1+b) − 1 = a ∘₁ b", ha="center", fontsize=8.5)
ax = axs[1]; ax.set_title("(b) Typical total mass vs depth")
for key, c in zip(R["cascade"], ["#2a9d8f", "#457b9d", "k", "#e76f51", "#9d0208"]):
    lab = f"σ = {float(key):.2f}" + ("  (critical line)" if abs(float(key)-R["sigma_critical"]) < 1e-3 else "")
    ax.semilogy(R["cascade"][key]["median"], "-o", ms=3, color=c, label=lab)
ax.set_xlabel("depth (number of scales)"); ax.set_ylabel("median total mass"); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.text(0.03, 0.62, "Kahane–Peyrière (1976):\nsurvives iff σ²/2 < log 2", transform=ax.transAxes, fontsize=8)
ax = axs[2]; ax.set_title("(c) Where the mass sits after 12 scales")
gg = np.random.default_rng(3)
for s, c, off in [(0.6, "#2a9d8f", 0), (2.0, "#9d0208", 1)]:
    m = np.ones(1)
    for d in range(12):
        m = np.repeat(m, 2)*np.exp(s*gg.standard_normal(2*len(m)) - s**2/2)/2
    ax.plot(np.linspace(0, 1, len(m)), m/m.max() + 1.15*off, color=c, lw=0.6, label=f"σ = {s}")
ax.set_yticks([]); ax.set_xlabel("position"); ax.legend(fontsize=8)
fig.suptitle("A cascade made of ∘₁ steps has an exact critical line; past it, the typical mass collapses onto a dust of points (intermittency, not blow-up)", fontsize=10.5)
fig.tight_layout(); fig.savefig(out+"fig3_multiplicative_cascade.png"); plt.close(fig)

# ---- Fig 4: chirality family ----
L = 32; tt = np.linspace(0, 2*np.pi, L, endpoint=False); X, Y = np.meshgrid(tt, tt, indexing="ij"); Z = 0.7
def abc(k, A=1.0, B=0.7, C=0.4):
    return np.array([A*np.sin(k*Z)+C*np.cos(k*Y), B*np.sin(k*X)+A*np.cos(k*Z), C*np.sin(k*Y)+B*np.cos(k*X)])
u1, u2 = abc(1), abc(2)
same = np.linalg.norm(np.cross(u1, 1*u1, axis=0), axis=0)
cross = np.linalg.norm(np.cross(u1, 2*u2, axis=0)+np.cross(u2, 1*u1, axis=0), axis=0)
fig = plt.figure(figsize=(15, 4.8)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.7])
for i, (F, ttl) in enumerate([(same, "(a) one handed wave with itself:\nnonlinearity |u×ω| = 0 everywhere"),
                              (cross, "(b) two waves, different twist:\n|cross term| = |λ₂−λ₁|·|u₁×u₂|")]):
    ax = fig.add_subplot(gs[i]); im = ax.imshow(F.T, origin="lower", extent=[0, 2*np.pi, 0, 2*np.pi], vmin=0, vmax=cross.max(), cmap="magma")
    ax.set_title(ttl, fontsize=9.5); ax.set_xlabel("x"); ax.set_ylabel("y"); plt.colorbar(im, ax=ax, fraction=0.046)
ax = fig.add_subplot(gs[2]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("(c) The single-chirality family (solid = known link, dashed = our reading)")
box(ax, 0.02, 0.62, 0.28, 0.24, "homochiral 3D Navier–Stokes\n2D-like inverse cascade (2012)\nproved regular (2013)", "#d8ecff", 8)
box(ax, 0.36, 0.62, 0.28, 0.24, "2D Euler\nregular (Wolibner, Yudovich)\ninfinitely many invariants ∫F(ω)", "#e8f5e9", 8)
box(ax, 0.70, 0.62, 0.28, 0.24, "self-dual Yang–Mills / gravity\nintegrable; all-same-helicity\ntree amplitudes vanish", "#fff3cd", 8)
box(ax, 0.36, 0.12, 0.28, 0.24, "area-preserving\ndiffeomorphisms\n(the shared algebra)", "#f4f4f4", 8)
arrow(ax, (0.31, 0.74), (0.355, 0.74), "", "k")
arrow(ax, (0.50, 0.61), (0.50, 0.37), "Arnold 1966", "k")
arrow(ax, (0.84, 0.61), (0.62, 0.30), "Monteiro–O'Connell 2011", "k", rad=-0.2)
arrow(ax, (0.16, 0.61), (0.38, 0.30), "ours?", "crimson", "--", rad=0.2)
ax.text(0.5, 0.03, "Open question this suggests: does homochiral flow carry an infinite family of invariants like 2D Euler?", ha="center", fontsize=8.5, color="crimson")
fig.tight_layout(); fig.savefig(out+"fig4_chirality_family.png"); plt.close(fig)

# ---- Fig 5: spacetime as turbulence ----
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
ax = axs[0]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("(a) Fluid/gravity: the dictionary (known)")
rows = [("Einstein's equations near a horizon", "Navier–Stokes (incompressible limit)"), ("horizon area", "entropy of the fluid"),
        ("ripples on the horizon", "eddies"), ("horizon in steady turbulence", "fractal-like, dimension d + 4/3"),
        ("black string thinning (5D)", "fluid thread pinching (Rayleigh–Plateau)"), ("cosmic censorship", "regularity (no naked blow-up)")]
ax.text(0.25, 0.94, "gravity", ha="center", weight="bold"); ax.text(0.75, 0.94, "fluid", ha="center", weight="bold")
for i, (g_, f_) in enumerate(rows):
    y = 0.82-0.14*i; ax.text(0.25, y, g_, ha="center", fontsize=8); ax.text(0.75, y, f_, ha="center", fontsize=8)
    ax.annotate("", xy=(0.52, y+0.005), xytext=(0.48, y+0.005), arrowprops=dict(arrowstyle="<->"))
ax = axs[1]; ax.set_title("(b) Black-string cascade (Lehner–Pretorius 2010, sketch)")
yy = 0
for gen in range(5):
    r = 0.5**gen; n = 2**gen; xs = (np.arange(n)+0.5)/n
    for xc in xs:
        ax.add_patch(Circle((xc, -1.2*gen), 0.18*r, fc="k"))
    ax.plot([0, 1], [-1.2*gen]*2, "k-", lw=3*0.4**gen)
    ax.text(1.04, -1.2*gen, f"gen {gen}: string ×{0.4**gen:.2f}", fontsize=7, va="center")
ax.set_xlim(-0.05, 1.45); ax.set_ylim(-5.5, 0.6); ax.axis("off")
ax.text(0.5, -5.4, "ever-smaller beads on ever-thinner strings; the string\nreaches zero width in finite time (numerical, not proved)", ha="center", fontsize=8)
ax = axs[2]; ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("(c) Does the cascade end in a singularity?")
cells = [("AdS₄ (Einstein + matter)", "cascade → black hole\nproved for massless Vlasov\n(Moschidis 2022)", "#ffd6d6"),
         ("AdS₃ (our union-rule space)", "cascade, but energy is below\nthe black-hole threshold → stays\nregular (numerics, 2013)", "#e8f5e9"),
         ("2D Navier–Stokes", "cascade (inverse) → regular\n(proved)", "#e8f5e9"),
         ("3D Navier–Stokes", "cascade (forward) → ???\n(the Millennium problem)", "#fff3cd")]
for i, (a_, b_, c_) in enumerate(cells):
    y = 0.76-0.24*i; box(ax, 0.0, y, 0.38, 0.18, a_, "#f4f4f4", 8); box(ax, 0.42, y, 0.56, 0.18, b_, c_, 8)
fig.suptitle("Spacetime can be a turbulent flow: near a horizon Einstein becomes Navier–Stokes, and gravity has its own cascades, one with a proved collapse and one protected by an energy threshold", fontsize=10.5)
fig.tight_layout(); fig.savefig(out+"fig5_spacetime_as_turbulence.png"); plt.close(fig)

# ---- Fig 6: the map ----
cands = ["Tao's averaged NS (2016)", "Single-chirality family\n(self-dual YM, 2D Euler)", "Scale entanglement\n(tensor networks)",
         "Multiplicative cascades\n(Mandelbrot, ∘₁)", "Supercritical defocusing NLS,\ncompressible Euler (MRRS)", "Boussinesq / axisymmetric\nEuler (Chen–Hou)",
         "NLS energy cascade\non a torus (CKSTT)", "Spacetime turbulence\n(fluid/gravity, AdS)"]
feats = ["energy\nsupercritical", "two sectors +\ncross term", "cascade to\nsmall scales", "blow-up\nproved there?"]
M = np.array([[1, 0.5, 1, 1], [0, 1, 0, 0], [0, 0.5, 1, 0], [0, 0, 1, 0], [1, 0, 0.5, 1], [1, 0.5, 1, 1], [0, 0, 1, 0], [1, 0.5, 1, 1]])
verd = ["proves the fine\nstructure matters", "LEAD: question", "LEAD: microscope,\nnot a proof", "picture only", "known; shows shape\ncan break", "known; same family", "known; slow version", "LEAD: tools\n+ AdS₃ contrast"]
vc = ["#cfe8cf", "#ffe08a", "#ffe08a", "#e0e0e0", "#cfe8cf", "#cfe8cf", "#e0e0e0", "#ffe08a"]
fig, ax = plt.subplots(figsize=(12, 6.2)); ax.set_xlim(-0.5, 5.6); ax.set_ylim(-0.5, len(cands)-0.5); ax.invert_yaxis()
for i in range(len(cands)):
    ax.text(-0.45, i, cands[i], va="center", ha="right", fontsize=8.5)
    for j in range(4):
        v = M[i, j]; ax.add_patch(plt.Rectangle((j-0.4, i-0.4), 0.8, 0.8, fc={1: "#2a9d8f", 0.5: "#a8dadc", 0: "#f1f1f1"}[v], ec="w"))
    ax.add_patch(FancyBboxPatch((4.0, i-0.38), 1.5, 0.76, boxstyle="round,pad=0.02", fc=vc[i], ec="k", lw=0.6))
    ax.text(4.75, i, verd[i], ha="center", va="center", fontsize=7.5)
for j, f in enumerate(feats): ax.text(j, -0.75, f, ha="center", va="bottom", fontsize=8.5)
ax.text(4.75, -0.75, "verdict", ha="center", va="bottom", fontsize=8.5, weight="bold")
ax.axis("off"); ax.set_title("Same-shape candidates against the three features (dark = has it, light = partly)", pad=40)
fig.tight_layout(); fig.savefig(out+"fig6_map.png", bbox_inches="tight"); plt.close(fig)
print("figures written")
