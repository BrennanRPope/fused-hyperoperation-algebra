"""Pictures for closed-leads-fresh-look (DRAFT). Every curve is from a formula or the checks; nothing fitted.
Run: python3 figures.py  -> figures/fig1..fig6 (PNG)."""
import numpy as np, matplotlib, math
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, Circle
from scipy import special
import sort_data as S
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
GREEN, ORANGE, BLUE, GREY, RED, PURPLE = "#2a9d55", "#e07b1f", "#2f6fb5", "#9a9a9a", "#c8382f", "#7a4fb5"
def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---------------- fig1: the sort board
verdict_col = {"F2": GREEN, "F3": GREEN, "F4": GREEN, "F1": BLUE, "F5": ORANGE, "answered": BLUE, "flag": PURPLE, "": GREY}
fig, ax = plt.subplots(figsize=(11, 6.2))
cols = {"K": 0, "O": 1, "U": 2}
rng = np.random.default_rng(3)
for li in range(6):
    for pv in "KOU":
        rows = [r for r in S.L if r["layer"] == li and r["prov"] == pv]
        for j, r in enumerate(rows):
            x = cols[pv] + 0.08 + (j % 9)*0.095
            y = 5 - li + 0.28 - (j//9)*0.22
            c = verdict_col[r["here"]]
            ax.scatter(x, y, s=95 if r["here"] else 55, color=c, edgecolor="k" if r["here"].startswith("F") else "none", zorder=3)
            if r["here"].startswith("F"):
                ax.annotate(r["here"], (x, y), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8, weight="bold")
ax.set_xlim(0, 3); ax.set_ylim(-0.6, 6)
ax.set_xticks([0.5, 1.5, 2.5]); ax.set_xticklabels(["only other people's methods\n(known math restated,\nor a standard outside test)",
                                                     "only our oldest tool\n(a ∘ₖ placement or fit)",
                                                     "our newer tools\nalready applied"])
ax.set_yticks([5.15 - i for i in range(6)]); ax.set_yticklabels(S.LAYERS)
for xx in (1, 2): ax.axvline(xx, color="#ddd", lw=1)
for yy in range(6): ax.axhline(yy - 0.5 + 0.5, color="#eee", lw=0.6)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="", color=GREEN, mec="k", ms=9, label="reworked here: opens a lead"),
     Line2D([], [], marker="o", ls="", color=ORANGE, mec="k", ms=9, label="reworked here: changes the reading only"),
     Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="settled by a later page of ours (F1 checked here)"),
     Line2D([], [], marker="o", ls="", color=PURPLE, ms=8, label="untried with our newer tools: flagged"),
     Line2D([], [], marker="o", ls="", color=GREY, ms=7, label="stays closed")]
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False, fontsize=9)
ax.set_title("Fig 1. 100 rows of closed leads, most fundamental at the top, sorted by whose tools were ever used on them")
save(fig, "fig1-sort-board")

# ---------------- fig2: F1 scale
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(wspace=0.45))
ax = axs[0]
K, D = np.meshgrid(np.linspace(-2, 2, 400), np.linspace(-1.5, 1.5, 300))
R = 1/3 - 4*D/(3*K**2 + 1e-12)
im = ax.contourf(K, D, np.clip(R, -2, 1), levels=30, cmap="coolwarm")
for cc in (-0.5, 0.25, 0.5, 1.0):
    kk = np.linspace(-2, 2, 200); ax.plot(kk, cc*kk**2, color="k", lw=0.8, ls="--")
ax.plot(np.linspace(-2, 2, 2), [0, 0], color=GREEN, lw=3)
ax.annotate("d = 0: only the sign of k\nsurvives rescaling", (1.2, 0), xytext=(0.1, -1.35), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.annotate("rescaling slides a setting\nalong a parabola d = c·k²", (-1.3, 0.25*1.69), xytext=(-1.95, 1.05), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_xlabel("homothety k"); ax.set_ylabel("curvature dial d")
fig.colorbar(im, ax=ax, label="iterative residue 1/3 − 4d/(3k²)")
ax.set_title("(a) The residue is constant along each parabola")
ax = axs[1]
cs = np.linspace(-1.5, 1.5, 300)
ax.plot(cs, 1/3 - 4*cs/3, color=BLUE, lw=2)
ax.scatter([0], [1/3], color=GREEN, zorder=3, s=60); ax.annotate("cube, e^x − 1: 1/3", (0, 1/3), xytext=(0.2, 0.9), arrowprops=dict(arrowstyle="->"))
ax.scatter([0.25], [0], color=ORANGE, zorder=3, s=60); ax.annotate("parabola class (Δ = 0): 0", (0.25, 0), xytext=(0.45, -0.9), arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("the scale-free number d/k² (a reading)"); ax.set_ylabel("iterative residue")
ax.set_title("(b) A whole real number survives once d is in")
fig.suptitle("Fig 2 (F1). k alone keeps only its sign; k with d keeps the number d/k², and the height-flow residue reads it", y=1.02)
save(fig, "fig2-scale-free-number")

# ---------------- fig3: F2 order / history
fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6))
ax = axs[0]
base = 1000
ax.add_patch(Rectangle((0, 0), 1, base, color="#cfe0f3")); ax.text(0.5, base/2, "1000", ha="center")
ax.add_patch(Rectangle((0, base), 1, 100, color=BLUE, alpha=0.7)); ax.text(0.5, base+50, "+10%: 100", ha="center", color="w", fontsize=8)
ax.add_patch(Rectangle((0, base+100), 1, 500, color="#f3d9bf")); ax.text(0.5, base+350, "+500", ha="center")
ax.text(0.5, 1680, "raise, then grant\n= 1600", ha="center")
ax.add_patch(Rectangle((1.4, 0), 1, base, color="#cfe0f3")); ax.text(1.9, base/2, "1000", ha="center")
ax.add_patch(Rectangle((1.4, base), 1, 500, color="#f3d9bf")); ax.text(1.9, base+250, "+500", ha="center")
ax.add_patch(Rectangle((1.4, base+500), 1, 100, color=BLUE, alpha=0.7)); ax.text(1.9, base+550, "10% of 1000", ha="center", color="w", fontsize=8)
ax.add_patch(Rectangle((1.4, base+600), 1, 50, color=ORANGE)); ax.text(2.45, base+625, "← 10% of 500 = 50:\n   the corner", fontsize=8, va="center")
ax.text(1.9, 1720, "grant, then raise\n= 1650", ha="center")
ax.set_xlim(-0.1, 3.3); ax.set_ylim(0, 1850); ax.set_xticks([]); ax.set_ylabel("budget")
ax.set_title("(a) Two reforms, two orders")
ax = axs[1]
import itertools
rng = np.random.default_rng(20261126)
kd = 0.4
same = [(1 + kd*s, s) for s in rng.uniform(-0.5, 1, 5)]
diff = [(1 + a_, c_) for a_, c_ in zip(rng.uniform(-0.3, 0.6, 5), rng.uniform(-1, 1, 5))]
def comp(seq):
    A, C = 1.0, 0.0
    for al, c in seq: A, C = al*A, al*C + c
    return C
od = [comp(p) for p in itertools.permutations(diff)]
os_ = [comp(p) for p in itertools.permutations(same)]
ax.hist(od, bins=30, color=GREEN, alpha=0.8, label="different centres: 120 outcomes")
ax.axvline(os_[0], color=BLUE, lw=3, label="shared centre: all 120 orders agree")
ax.set_xlabel("final outcome (shift) after five reforms"); ax.set_ylabel("orders")
ax.legend(fontsize=8, frameon=False); ax.set_title("(b) History appears exactly when centres differ")
ax = axs[2]
# upper half plane: points (C, A); geodesics; commutator loop
th = np.linspace(0, np.pi, 200)
for c0, r0 in ((0, 1), (0, 2), (1.5, 1.2), (-1.5, 1.0)):
    ax.plot(c0 + r0*np.cos(th), r0*np.sin(th), color="#ccc", lw=0.8)
for xv in (-2, -1, 0, 1, 2): ax.plot([xv, xv], [0.05, 3], color="#eee", lw=0.6)
f = (1.5, 0.4); g = (0.8, -0.6)  # (stretch A, shift C)
def act(m, p): A, C = p; a, c = m; return (a*A, a*C + c)
def inv(m): a, c = m; return (1/a, -c/a)
P0 = (1.0, 0.0); pts = [P0]
for m in (inv(g), inv(f), g, f):
    pts.append(act(m, pts[-1]))
xs = [p[1] for p in pts]; ys = [p[0] for p in pts]
ax.plot(xs, ys, "-o", color=GREEN, ms=4)
ax.annotate("start: do nothing", (0, 1), xytext=(0.6, 0.25), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.annotate("back after the loop:\nshifted by the area", (xs[-1], ys[-1]), xytext=(-2.4, 2.4), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.set_xlim(-2.6, 2.6); ax.set_ylim(0, 3)
ax.set_xlabel("shift C"); ax.set_ylabel("stretch A (> 0)")
ax.set_title("(c) Reforms as points of the hyperbolic plane")
fig.suptitle("Fig 3 (F2). Composition gives history: the cost of order is an oriented area, and the space of reforms is curved", y=1.02)
save(fig, "fig3-history-from-composition")

# ---------------- fig4: F3 elliptic
fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6))
ax = axs[0]
u = np.linspace(-4, 4, 400)
for m, c in zip((0, 0.5, 0.9, 0.99, 1.0), plt.cm.viridis(np.linspace(0, 0.9, 5))):
    ax.plot(u, special.ellipj(u, m)[0], color=c, lw=2, label=f"m = {m}")
ax.text(-3.9, 1.12, "m = 0: sin (round side, d > 0)", color=plt.cm.viridis(0), fontsize=8)
ax.text(-3.9, -1.25, "m = 1: tanh (relativity side, d < 0)", color=plt.cm.viridis(0.9), fontsize=8)
ax.set_ylim(-1.4, 1.4); ax.set_xlabel("u"); ax.set_ylabel("sn(u | m)")
ax.legend(fontsize=7, frameon=False, loc="lower right"); ax.set_title("(a) Jacobi's modulus runs from sin to tanh")
ax = axs[1]
th = np.linspace(0, 2*np.pi, 800)
for de, c, lab in ((0.0, BLUE, "Edwards parameter 0: the circle (our conic)"), (-1.0, GREEN, "−1: Euler–Gauss lemniscatic curve"),
                   (-8.0, PURPLE, "−8"), (0.6, ORANGE, "+0.6")):
    # x^2 + y^2 = 1 + de x^2 y^2 in polar form: r^2 = 1 + de r^4 c^2 s^2  -> solve
    cs2 = (np.cos(th)*np.sin(th))**2
    with np.errstate(invalid="ignore", divide="ignore"):
        r2 = np.where(np.abs(de*cs2) < 1e-12, 1.0, (1 - np.sqrt(1 - 4*de*cs2))/(2*de*cs2 + 1e-300))
    r = np.sqrt(r2)
    ax.plot(r*np.cos(th), r*np.sin(th), color=c, lw=2, label=lab)
ax.set_aspect("equal"); ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1))
ax.set_title("(b) x² + y² = 1 + (Edwards parameter)·x²y²")
ax = axs[2]
ms = np.linspace(0.001, 0.999, 500)
j = 16*(1 + 14*ms + ms**2)**3/(ms*(1 - ms)**4)
ax.semilogy(ms, j, color=GREEN, lw=2)
ax.axvspan(0, 0.02, color=BLUE, alpha=0.2); ax.axvspan(0.98, 1, color=BLUE, alpha=0.2)
ax.text(0.03, 3e6, "our plane dial:\nj = ∞ at both ends\n(nodal cubics)", fontsize=8)
ax.set_xlabel("Jacobi's modulus m (= Edwards parameter at a = 1)"); ax.set_ylabel("j-invariant (log scale)")
ax.set_title("(c) Elliptic curves live between the two ends")
fig.suptitle("Fig 4 (F3). Elliptic functions are not a fractional rung: they are one parameter off our plane dial, joining its round and hyperbolic sides", y=1.02)
save(fig, "fig4-elliptic-off-the-plane")

# ---------------- fig5: F4 tunnelling
def delta_M(g, kw):
    be = g/(2*kw); return np.array([[1 - 1j*be, -1j*be], [1j*be, 1 + 1j*be]])
def prop(kw, Lx): return np.array([[np.exp(1j*kw*Lx), 0], [0, np.exp(-1j*kw*Lx)]])
def odds(M): T = 1/abs(M[1, 1])**2; return (1 - T)/T
fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6))
ax = axs[0]
# hyperbolic triangle in the Poincare disc: two sides D1, D2 with turn phi
def disc_point(dist, ang): r = np.tanh(dist/2); return r*np.cos(ang), r*np.sin(ang)
ax.add_patch(Circle((0, 0), 1, fill=False, color="#bbb"))
D1, D2, phi = 1.4, 1.1, 2.1
P1 = disc_point(D1, 0)
# second side: from P1 at turn phi, approximate with Mobius translate
z1 = complex(*P1)
w = np.tanh(D2/2)*np.exp(1j*(np.pi - phi))
z2 = (w + z1)/(1 + np.conj(z1)*w)
def geod(a, b, n=60):
    # geodesic in disc via Mobius to origin
    t = np.linspace(0, 1, n); m = (b - a)/(1 - np.conj(a)*b)
    r = np.tanh(np.arctanh(abs(m))*t)*np.exp(1j*np.angle(m))
    return (r + a)/(1 + np.conj(a)*r)
for (a, b), c, lab in (((0, z1), BLUE, "barrier 1 (odds t₁): length D₁"), ((z1, z2), ORANGE, "barrier 2 (odds t₂): length D₂"), ((0, z2), GREEN, "both: length D")):
    g = geod(complex(a), complex(b)); ax.plot(g.real, g.imag, color=c, lw=2.5, label=lab)
ax.annotate("turn = loop phase φ", (z1.real, z1.imag), xytext=(0.1, -0.6), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.text(-0.95, 0.85, "1 + 2t = cosh D,  t = sinh²(D/2)", fontsize=8)
ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.set_title("(a) Two barriers make a hyperbolic triangle")
ax = axs[1]
kw = 1.0; M1, M2 = delta_M(1.3, kw), delta_M(0.8, kw); u_, v_ = odds(M1), odds(M2)
Ls = np.linspace(0, np.pi, 600)
rho = np.array([odds(M2 @ prop(kw, Lx) @ M1) for Lx in Ls])
ax.plot(Ls, rho, color="k", lw=1.5, label="two barriers: odds t = R/T")
ax.axhline(u_ + v_ + 2*u_*v_, color=GREEN, lw=2, label="average = t₁ ∘₂ t₂ (mean resistance)")
ax.axhline(u_ + v_ + u_*v_, color=BLUE, lw=2, ls="--", label="log-average = t₁ ∘₁ t₂ (typical)")
ax.axhline(u_ + v_, color=ORANGE, lw=2, ls=":", label="from mean transmission = t₁ + t₂")
ax.set_xlabel("spacing × wavenumber (one full turn of φ)"); ax.set_ylabel("odds of reflection t")
ax.legend(fontsize=7, frameon=False); ax.set_title("(b) One pair, three averages, three settings of k")
ax = axs[2]
g = 1.0; u1 = odds(delta_M(g, kw)); Ns = np.arange(1, 11)
ax.semilogy(Ns, ((1 + 2*u1)**Ns - 1)/2, color=GREEN, lw=2, label="mean: ∘₂ power ((1+2t₁)ᴺ − 1)/2")
ax.semilogy(Ns, (1 + u1)**Ns - 1, color=BLUE, lw=2, ls="--", label="typical: ∘₁ power (1+t₁)ᴺ − 1")
rng = np.random.default_rng(5); mc_m, mc_t = [], []
Mg = delta_M(g, kw)
for N in Ns:
    vals = []
    for _ in range(6000):
        M = Mg.copy()
        for _j in range(1, N): M = Mg @ prop(kw, rng.uniform(0, np.pi)) @ M
        vals.append(odds(M))
    vals = np.array(vals); mc_m.append(vals.mean()); mc_t.append(np.exp(np.log1p(vals).mean()) - 1)
ax.semilogy(Ns, mc_m, "o", color=GREEN, ms=5, label="simulated barriers, mean")
ax.semilogy(Ns, mc_t, "s", color=BLUE, ms=5, label="simulated barriers, typical")
ax.set_xlabel("number of barriers N"); ax.set_ylabel("odds of reflection"); ax.legend(fontsize=7, frameon=False)
ax.set_title("(c) Many barriers: two growth rates (localization)")
fig.suptitle("Fig 5 (F4). Tunnelling through separated barriers: exact laws on the odds dial, k = 2 for the mean and k = 1 for the typical value", y=1.02)
save(fig, "fig5-tunnelling-hyperbolic")

# ---------------- fig6: F5 acids
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
acids = {"citric": (3.13, 4.76, 6.40), "phosphoric": (2.14, 7.20, 12.34), "arsenic": (2.21, 6.93, 11.52)}
for nm, (p1, p2, p3) in acids.items():
    ax.scatter(p2 - p1, p3 - p2, s=60, color=BLUE, zorder=3); ax.annotate(nm, (p2 - p1, p3 - p2), xytext=(5, -10), textcoords="offset points")
xx = np.linspace(0, 6, 10)
ax.plot(xx, xx, color=GREEN, lw=2, label="pair-only weights, uniform: equal gaps")
ax.axvline(math.log10(3), color=RED, ls="--", lw=2, label="rank-only fold k^(m−1): lower gap = log 3")
ax.set_xlabel("pKa₂ − pKa₁"); ax.set_ylabel("pKa₃ − pKa₂"); ax.legend(fontsize=8, frameon=False)
ax.set_title("(a) The acids pick the pair-only corner weights")
ax = axs[1]
ys = np.array([0.0, 0.7, 1.6])
ax.plot(ys, [0, 0, 0], "o", color=BLUE, ms=9)
for i, j, h in ((0, 1, 0.25), (1, 2, 0.25), (0, 2, 0.55)):
    ax.annotate("", (ys[j], 0.02), xytext=(ys[i], 0.02), arrowprops=dict(arrowstyle="-", connectionstyle=f"arc3,rad=-{h}", color=GREY))
ax.text(0.8, 0.5, "on a line, distances add:\nartanh√w₁₃ = artanh√w₁₂ + artanh√w₂₃", ha="center", fontsize=8)
tri = np.array([[2.6, 0], [3.6, 0], [3.1, 0.866]])
ax.add_patch(Polygon(tri, fill=False, color=GREEN, lw=2)); ax.plot(tri[:, 0], tri[:, 1], "o", color=GREEN, ms=9)
ax.text(3.1, -0.35, "three equal pair factors need\nthree equal distances: a plane", ha="center", fontsize=8)
ax.set_xlim(-0.3, 4); ax.set_ylim(-0.6, 1.1); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) Cauchy pair factors from points log p")
fig.suptitle("Fig 6 (F5). Triprotic acids: the old page called them 'against the dial'; they sit in our pair-only corner-weight family", y=1.02)
save(fig, "fig6-acids-corner-weights")
print("figures written")
