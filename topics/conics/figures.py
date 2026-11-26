import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12.5, "figure.dpi": 130, "savefig.bbox": "tight"})
BLUE, ORANGE, GREY, MAG, GREEN = "#1f5fbf", "#e07b00", "#8a8f99", "#e0115f", "#2e8b57"
EC, PC, HC = "#1f5fbf", "#e0115f", "#e07b00"   # ellipse, parabola, hyperbola colours
def ctype(k, d):
    D = k*k - 4*d
    return "parabola" if abs(D) < 1e-9 else ("ellipse" if D < 0 else "hyperbola")
COL = {"ellipse": EC, "parabola": PC, "hyperbola": HC}
def draw_conic(ax, k, d, lim=2.2, lvl_rhs="y", color=None, lw=2, n=500):
    x = np.linspace(-lim, lim, n); y = np.linspace(-lim, lim, n); X, Y = np.meshgrid(x, y)
    Z = X**2 + k*X*Y + d*Y**2 - (Y if lvl_rhs == "y" else lvl_rhs)
    ax.contour(X, Y, Z, [0], colors=[color or COL[ctype(k, d)]], linewidths=lw)
def P(s, k, d):
    q = 1 + k*s + d*s*s; return np.array([s/q, s*s/q])
F = lambda a, b, k, d: (a + b + k*a*b)/(1 - d*a*b)

# ---------- Fig 1: the dial plane is a map of conics
fig, ax = plt.subplots(figsize=(10.5, 7.2))
kk = np.linspace(-3.4, 3.4, 400)
ax.fill_between(kk, kk**2/4, 4.2, color=EC, alpha=.10)
ax.fill_between(kk, -2.6, kk**2/4, color=HC, alpha=.10)
ax.plot(kk, kk**2/4, color=PC, lw=3)
for kap in (-1.5, -0.75, 0.75, 1.5):
    ax.plot(kk, kap + kk**2/4, color=GREY, lw=.8, ls="--")
ax.axhline(0, color="k", lw=.6); ax.axvline(0, color="k", lw=.6)
ax.plot(kk, 0*kk, color=GREY, lw=2.5, alpha=.5)
ax.text(-3.3, 3.85, "no wall: ELLIPSE (circle at k=0, d=1)", color=EC, fontsize=12, weight="bold")
ax.annotate("one double wall:\nPARABOLA on the curve d = k²/4", (-2.9, 2.1025), xytext=(-3.3, 2.55), color=PC, fontsize=11.5, weight="bold", arrowprops=dict(arrowstyle="->", color=PC))
ax.text(-3.3, -2.45, "two walls: HYPERBOLA", color=HC, fontsize=12, weight="bold")
ax.text(3.35, -0.08, "d = 0: the ∘ₖ line (grey band), touching the\nparabola at its tip: hyperbolas except k = 0", color="#555", fontsize=9.5, ha="right", va="top")
ax.text(-3.35, 1.6+ 2.89-0.1, "", fontsize=1)
ax.text(3.35, 4.05, "dashed: lines of equal geometry κ = d − k²/4\n(also parabolas)", color="#555", fontsize=9, ha="right", va="top")
pts = [(0, 1, "circle\n(0, 1)", 2.0, 1), (1.3, 2.7, "ellipse\n(1.3, 2.7)", 1.0, 1), (0, 0, "addition\n(0, 0)", 2.0, 1), (2, 1, "parabola\n(2, 1)", 2.0, 1),
       (0, -1, "relativity\n(0, −1/c²), c = 1", 2.0, 1), (-1, 0, "probability\nk = −1", 2.0, -1), (-2.4, 0.9, "hyperbola", 2.0, 1), (2.6, -1.6, "hyperbola", 2.0, -1)]
for (k0, d0, lab, lm, side) in pts:
    ax.plot(k0, d0, "o", color="k", ms=5, zorder=5)
    ins = ax.inset_axes([k0 - .42, d0 + .12, .84, .84], transform=ax.transData)
    draw_conic(ins, k0, d0, lim=lm, lw=1.8)
    ins.plot(0, 0, "o", color="k", ms=2.5)
    ins.set_xticks([]); ins.set_yticks([]); ins.set_aspect("equal")
    for sp_ in ins.spines.values(): sp_.set_color("#bbb")
    ax.text(k0 + side*.47, d0 + .5, lab, fontsize=8.5, va="center", ha="left" if side > 0 else "right")
ax.set_xlim(-3.4, 3.4); ax.set_ylim(-2.6, 4.2)
ax.set_xlabel("k  (the first dial)"); ax.set_ylabel("d  (the plane dial)")
ax.set_title("Figure 1. Every setting (k, d) of the plane law is one conic, x² + kxy + dy² = y. Its type is the dial's region.", loc="left")
fig.savefig("figures/fig1-dial-plane-is-a-map-of-conics.png"); plt.close(fig)

# ---------- Fig 2: the chord rule on each conic is our law
fig, axs = plt.subplots(1, 3, figsize=(14, 5))
cases = [(0, 1, 0.55, 0.9, "ellipse  (k, d) = (0, 1)\nchord rule = adding angles"),
         (0, 0, 0.6, 0.8, "parabola  (k, d) = (0, 0)\nchord rule = plain addition"),
         (0, -1, 0.3, 0.6, "hyperbola  (k, d) = (0, −1)\nchord rule = adding velocities")]
for ax, (k0, d0, s, t, ttl) in zip(axs, cases):
    lim = {1: 1.3, 0: 1.8, -1: 2.2}[d0]; draw_conic(ax, k0, d0, lim=lim + .5)
    O = np.zeros(2); Pp, Qq = P(s, k0, d0), P(t, k0, d0); r = F(s, t, k0, d0); Rr = P(r, k0, d0)
    ax.plot(*zip(Pp, Qq), color="k", lw=1.6)
    v = (Qq - Pp)/np.linalg.norm(Qq - Pp); seg = np.array([O - .3*v, Rr + .2*v])
    ax.plot(seg[:, 0], seg[:, 1], color=GREEN, lw=1.6, ls="--")
    ax.axhline(0, color=GREY, lw=.8)
    for pt, lab, c in [(O, "O  (identity, a = 0)", "k"), (Pp, f"P  (a = {s})", "k"), (Qq, f"Q  (a = {t})", "k"),
                       (Rr, f"P∘Q  (a = {r:.3f})", GREEN)]:
        ax.plot(*pt, "o", color=c, ms=6, zorder=6); ax.annotate(lab, pt, xytext=((6, -14) if lab.startswith("O") else (-118, 6) if c == GREEN else (8, -6) if lab.startswith("P ") else (8, 2)), textcoords="offset points", fontsize=9.5, color=c)
    if d0 < 0:
        xx = np.linspace(-lim, lim, 2)
        for w in (1, -1):   # asymptote directions y = w x through centre (0,-1/2)
            ax.plot(xx, w*xx - 0.5, color=HC, lw=.9, ls=":")
        ax.text(-2.1, -1.25, "asymptotes =\nthe walls (a = ±1)", color=HC, fontsize=9.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim + .6, lim + .6); ax.set_aspect("equal"); ax.set_title(ttl)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Figure 2. Draw the chord PQ, then the parallel line through O: it lands on P∘Q. In the slope coordinate a this is exactly "
             "(a + b + kab)/(1 − dab).", x=0.01, ha="left", y=1.07, fontsize=12.5)
fig.savefig("figures/fig2-chord-rule-is-the-law.png"); plt.close(fig)

# ---------- Fig 3: the levels of the plane number's size: all central conics, degenerate ones included
fig, axs = plt.subplots(3, 3, figsize=(11, 11))
rowlab = [("κ > 0 (no wall)", 1), ("κ = 0 (double wall)", 0), ("κ < 0 (two walls)", -1)]
names = {(1, 1): "ellipse (circle)", (1, 0): "a single point", (1, -1): "nothing real\n(imaginary ellipse)",
         (0, 1): "two parallel lines", (0, 0): "a double line", (0, -1): "nothing real",
         (-1, 1): "hyperbola", (-1, 0): "two crossing lines\n= the walls (light cone)", (-1, -1): "the other hyperbola"}
x = np.linspace(-2, 2, 600); X, Y = np.meshgrid(x, x)
for i, (rl, kap) in enumerate(rowlab):
    for j, lev in enumerate((1, 0, -1)):
        ax = axs[i, j]; Z = X**2 + kap*Y**2
        c = [EC, PC, HC][i]
        if kap == 1 and lev == 0: ax.plot(0, 0, "o", color=c, ms=9)
        elif kap == 0 and lev == 0: ax.axvline(0, color=c, lw=5, alpha=.8)
        elif kap == -1 and lev == 0: ax.plot(x, x, color=c, lw=2.2); ax.plot(x, -x, color=c, lw=2.2)
        elif lev == -1 and kap >= 0: ax.text(0, 0, "∅", fontsize=34, ha="center", va="center", color=GREY)
        else: ax.contour(X, Y, Z - lev, [0], colors=[c], linewidths=2.2)
        if kap == -1 and lev != 0:
            ax.plot(x, x, color=c, lw=.7, ls=":"); ax.plot(x, -x, color=c, lw=.7, ls=":")
        ax.set_title(names[(kap, lev)], fontsize=11); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.axhline(0, color="#ddd", lw=.6, zorder=0); ax.axvline(0, color="#ddd", lw=.6, zorder=0)
        if j == 0: ax.set_ylabel(rl, fontsize=12, weight="bold", color=c)
        if i == 0: ax.annotate(["size N = 1", "size N = 0", "size N = −1"][j], (0.5, 1.25), xycoords="axes fraction", ha="center", fontsize=12.5, weight="bold")
fig.suptitle("Figure 3. The size N = X² + κY² of a plane number (after completing the square). Its three levels give every central conic:\n"
             "size 1 is the 'unit circle', size 0 is where numbers can't be divided by (the walls).", x=0.01, ha="left", y=1.005, fontsize=12.5)
fig.savefig("figures/fig3-levels-give-every-central-conic.png"); plt.close(fig)

# ---------- Fig 4: the ruler is Cayley's distance with the walls as the absolute
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
def ruler(a, k, d):
    kap = d - k*k/4; T = a/(1 + k*a/2)
    if kap > 1e-12: return np.arctan(np.sqrt(kap)*T)/np.sqrt(kap)
    if kap < -1e-12: return np.arctanh(np.sqrt(-kap)*T)/np.sqrt(-kap)
    return T
for ax, (k0, d0, ttl, c) in zip(axs, [(0, 1, "no real wall: angle\n(Cayley's elliptic measure)", EC),
                                      (2, 1, "one double wall at a = −1\n(parabolic measure)", PC),
                                      (-1, 0, "two walls at a = 1 and ∞: log of a cross-ratio\n(hyperbolic measure; this is ∘₋₁, probability)", HC)]):
    a = np.linspace(-4, 4, 4001)
    with np.errstate(all="ignore"):
        L = ruler(a, k0, d0)
    if d0 == 0:
        with np.errstate(all="ignore"): L = np.log(1 + k0*a)/k0
    L[np.abs(np.gradient(L)) > 1] = np.nan
    ax.plot(a, L, color=c, lw=2.4)
    D = k0*k0 - 4*d0
    if D >= 0:
        ws = np.roots([d0, k0, 1]) if d0 != 0 else [-1/k0]
        for w in np.real(ws): ax.axvline(w, color=c, ls=":", lw=1.4); ax.text(w + .08, 2.6, f"wall\na = {w:g}", color=c, fontsize=9)
    ax.axhline(0, color=GREY, lw=.6); ax.axvline(0, color=GREY, lw=.6)
    ax.set_ylim(-3.2, 3.2); ax.set_xlim(-4, 4); ax.set_title(ttl); ax.set_xlabel("a (the quantity composed)")
axs[0].set_ylabel("ruler: the distance from 0 to a")
fig.suptitle("Figure 4. The law's ruler is Cayley's projective distance (1859), with the walls playing the absolute conic. "
             "The ruler runs to infinity exactly at the walls.", x=0.01, ha="left", y=1.07, fontsize=12.5)
fig.savefig("figures/fig4-ruler-is-cayley-distance.png"); plt.close(fig)

# ---------- Fig 5: Cayley–Klein / Hilbert: ellipse, parabola and both branches of a hyperbola are one disc
fig, axs = plt.subplots(1, 3, figsize=(14.5, 5))
th = np.linspace(0, 2*np.pi, 400)
# hyperbolic circles of radius rho about centres, drawn in the Klein disc
def klein_circle(c, rho, n=300):
    # point at hyperbolic distance rho from c in direction phi (Klein model via hyperboloid)
    cz = np.array(c); g = 1/np.sqrt(1 - cz @ cz)
    H0 = np.array([g, g*cz[0], g*cz[1]])    # hyperboloid point
    # orthonormal tangent basis at H0 (Minkowski)
    e1 = np.array([0, 1, 0.]); e2 = np.array([0, 0, 1.])
    def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def proj(v): return v + mdot(v, H0)*H0
    t1 = proj(e1); t1 /= np.sqrt(mdot(t1, t1)); t2 = proj(e2); t2 -= mdot(t2, t1)*t1; t2 /= np.sqrt(mdot(t2, t2))
    out = []
    for phi in np.linspace(0, 2*np.pi, n):
        v = np.cos(phi)*t1 + np.sin(phi)*t2; Hp = np.cosh(rho)*H0 + np.sinh(rho)*v
        out.append(Hp[1:]/Hp[0])
    return np.array(out)
cents = [(0, 0), (0.5, 0.2), (-0.3, -0.5), (0.1, 0.7), (-0.65, 0.3)]
ax = axs[0]; ax.plot(np.cos(th), np.sin(th), color=EC, lw=2.5)
for c in cents: K = klein_circle(c, 0.5); ax.plot(K[:, 0], K[:, 1], color=GREY, lw=1.2)
for (p1, p2) in [((-1, 0), (0.6, 0.8)), ((-0.8, -0.6), (0.9, -0.43)), ((0, -1), (0.2, 0.98))]:
    ax.plot(*zip(p1, p2), color="k", lw=.9)
ax.set_title("circle (any ellipse): Hilbert's distance inside\n= the hyperbolic plane (Klein model)"); ax.set_aspect("equal"); ax.axis("off")
ax = axs[1]
inv = lambda z: np.array([ (z[:, 0]*(1 + (1 + z[:, 1])/(1 - z[:, 1])))/2 if False else z[:, 0]/(1 - z[:, 1]), (1 + z[:, 1])/(1 - z[:, 1]) ]).T
xx = np.linspace(-2.2, 2.2, 400); ax.plot(xx, xx**2, color=PC, lw=2.5)
for c in [(0, -0.4), (0.4, -0.2), (-0.35, -0.1), (0, 0.25), (0.45, 0.35), (-0.5, -0.6)]:
    K = inv(klein_circle(c, 0.35)); ax.plot(K[:, 0], K[:, 1], color=GREY, lw=1.2)
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-0.3, 4.8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("parabola: the same disc, touching infinity once.\nSame circles (grey), same geometry")
ax = axs[2]; yy = np.linspace(1, 4.2, 300)
for sgn in (1, -1): ax.plot(np.sqrt(yy**2 - 1), sgn*yy, color=HC, lw=2.5); ax.plot(-np.sqrt(yy**2 - 1), sgn*yy, color=HC, lw=2.5)
inv2 = lambda z: np.array([z[:, 0]/z[:, 1], 1/z[:, 1]]).T
for c in [(0, 0.6), (0.3, -0.62), (-0.4, 0.5), (0.35, 0.45), (-0.2, -0.7)]:
    K = klein_circle(c, 0.3)
    if np.all(K[:, 1] > 0.02) or np.all(K[:, 1] < -0.02):
        Kh = inv2(K); ax.plot(Kh[:, 0], Kh[:, 1], color=GREY, lw=1.2)
ax.plot([-4, 4], [4, -4], color=HC, lw=.7, ls=":"); ax.plot([-4, 4], [-4, 4], color=HC, lw=.7, ls=":")
ax.text(0, 0, "the line at infinity\ncuts the disc here", ha="center", va="center", fontsize=9, color="#555")
ax.set_xlim(-3.3, 3.3); ax.set_ylim(-4.2, 4.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("hyperbola: the insides of BOTH branches are the\ntwo halves of one disc; one branch alone is not")
fig.suptitle("Figure 5. Projectively, ellipse, parabola and hyperbola are one conic and give one geometry. They differ only in how they meet "
             "the line at infinity: 0, 1 or 2 points, which is our wall count.", x=0.01, ha="left", y=1.07, fontsize=12.5)
fig.savefig("figures/fig5-one-conic-three-looks.png"); plt.close(fig)

# ---------- Fig 6: the dial is the conic's affine curvature, not its eccentricity; Kepler
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
ax = axs[0]
for lam, ls in [(1, "-"), (2, "--"), (3, ":")]:
    draw_conic(ax, 0.5*lam, 1.0*lam*lam, lim=1.5, color=EC, lw=2)
ax.plot([], [], color=EC, ls="-", label="(k, d) = (0.5, 1)"); ax.plot([], [], color=EC, ls="--", label="(1, 4)"); ax.plot([], [], color=EC, ls=":", label="(1.5, 9)")
ax.legend(fontsize=9, loc="upper left"); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlim(-.9, .9); ax.set_ylim(-.15, 1.2)
ax.set_title("one law (a → λa), three different ellipses:\neccentricity is a choice of ruler")
ax = axs[1]
dd = np.linspace(-1.2, 1.8, 600); kap = dd
def ecc(k, d):
    ev = np.linalg.eigvalsh([[1, k/2], [k/2, d]]); lp, lm = ev.max(), ev.min()
    if abs(lm) < 1e-12: return 1.0
    return np.sqrt(1 - lm/lp) if lm > 0 else np.sqrt(1 + abs(lm)/lp)
e = np.array([ecc(0, x) for x in dd])
ax.plot(dd, 2**(4/3)*kap, color=MAG, lw=2.6, label="affine curvature = 2^(4/3) · κ (exact)")
ax.plot(dd, e, color=GREY, lw=2, label="eccentricity e (k = 0 slice)")
ax.axvline(0, color=PC, lw=1, ls=":"); ax.axhline(0, color="k", lw=.5); ax.axhline(1, color=GREY, lw=.5, ls=":")
ax.text(0.04, -2.6, "parabola:\nκ = 0, e = 1", color=PC, fontsize=9)
ax.fill_between(dd, -3.2, 3.2, where=dd > 0, color=EC, alpha=.07); ax.fill_between(dd, -3.2, 3.2, where=dd < 0, color=HC, alpha=.07)
ax.set_ylim(-3.2, 3.2); ax.set_xlabel("d  (with k = 0, so κ = d)"); ax.legend(fontsize=8.8, loc="upper left")
ax.set_title("the dial is the conic's affine curvature,\nsign and size, straight through the parabola")
ax = axs[2]
mu, Lm = 1.0, 1.0
phi = np.linspace(-np.pi, np.pi, 2000)
for E, c, lab in [(-0.3, EC, "E < 0: ellipse"), (0.0, PC, "E = 0: parabola"), (0.4, HC, "E > 0: hyperbola")]:
    ee = np.sqrt(1 + 2*E*Lm**2/mu**2); r = (Lm**2/mu)/(1 + ee*np.cos(phi))
    r[(1 + ee*np.cos(phi)) <= 0.03] = np.nan
    ax.plot(r*np.cos(phi), r*np.sin(phi), color=c, lw=2.2, label=f"{lab} (e = {ee:.2f})")
ax.plot(0, 0, "*", color="k", ms=12); ax.set_xlim(-4.5, 1.5); ax.set_ylim(-3.2, 3.2); ax.set_aspect("equal")
ax.legend(fontsize=8.8, loc="lower left"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Kepler, same angular momentum: energy (the plane dial)\npicks the type; e² − 1 = 2EL²/μ²")
fig.suptitle("Figure 6. What the dial measures on a conic: its type and its affine curvature. Ordinary eccentricity needs one more number.",
             x=0.01, ha="left", y=1.07, fontsize=12.5)
fig.savefig("figures/fig6-affine-curvature-not-eccentricity.png"); plt.close(fig)

# ---------- Fig 7: the other dials each touch the conics at one point
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
ax = axs[0]
e1 = np.array([1, -1, 0])/np.sqrt(2); e2 = np.array([1, 1, -2])/np.sqrt(6)
for T, c in [(0.05, GREY), (0.5, "#6b8fd6"), (5, BLUE), (100, MAG)]:
    pts = []
    for a_ in np.linspace(0, 2*np.pi, 400):
        dv = np.cos(a_)*e1 + np.sin(a_)*e2; lo, hi = 0, 60.
        for _ in range(80):
            m = (lo + hi)/2; f = T*np.log(np.sum(np.exp(m*dv/T))) - T*np.log(3) - 1
            lo, hi = (m, hi) if f < 0 else (lo, m)
        pts.append(lo*np.array([np.cos(a_), np.sin(a_)]))
    pts = np.array(pts); pts /= np.abs(pts).max()
    ax.plot(pts[:, 0], pts[:, 1], color=c, lw=2, label=f"T = {T}")
ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=9, loc="lower left")
ax.set_title("roundness (temperature): cold = triangle (3 lines),\nhot = circle. The conic is only the hot limit")
ax = axs[1]
th = np.linspace(0, 2*np.pi, 800)
for p, c in [(1, GREY), (1.5, "#9aa"), (2, MAG), (3, "#6b8fd6"), (12, BLUE)]:
    xs = np.sign(np.cos(th))*np.abs(np.cos(th))**(2/p); ys = np.sign(np.sin(th))*np.abs(np.sin(th))**(2/p)
    ax.plot(xs, ys, color=c, lw=2.6 if p == 2 else 1.6, label=f"p = {p}" + ("  (the only conic)" if p == 2 else ""))
ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=9, loc="center", bbox_to_anchor=(0.5, 0.5))
ax.set_title("ℓ^p dial (cross-polytope → sphere → cube):\nconic only at p = 2; ends are pairs of parallel lines")
ax = axs[2]; ax.axis("off")
ax.add_patch(Rectangle((0.02, 0.32), 0.96, 0.34, color="#f1f3f8", transform=ax.transAxes))
ax.annotate("", (0.88, 0.49), (0.18, 0.49), xycoords="axes fraction", arrowprops=dict(arrowstyle="-|>", lw=2, color=GREY))
ins = ax.inset_axes([0.02, 0.36, 0.26, 0.26]); ins.plot(np.cos(th), 0.7*np.sin(th), color=EC, lw=2); ins.set_aspect("equal"); ins.axis("off")
ins2 = ax.inset_axes([0.74, 0.36, 0.26, 0.26]); ins2.fill([0, 1, .5, 0], [0, 0, .87, 0], fill=False, ec=GREY, lw=2); ins2.set_aspect("equal"); ins2.axis("off")
ax.text(0.15, 0.25, "round end:\nevery conic\n(Pick form = 0)", ha="center", fontsize=10, color=EC, transform=ax.transAxes)
ax.text(0.87, 0.25, "cornered end:\ntriangle", ha="center", fontsize=10, color=GREY, transform=ax.transAxes)
ax.text(0.5, 0.75, "shape dial (Benzécri space, up to projective maps)", ha="center", fontsize=11, transform=ax.transAxes)
ax.text(0.5, 0.06, "ellipse, parabola and hyperbola are all the same point here", ha="center", fontsize=10, transform=ax.transAxes, style="italic")
ax.set_title("shape dial: the conics are one end of it")
fig.suptitle("Figure 7. The roundness dials and the shape dial each meet the conic family in just one place: the circle.", x=0.01, ha="left", y=1.07, fontsize=12.5)
fig.savefig("figures/fig7-other-dials-touch-once.png"); plt.close(fig)
print("ok")

# ---------- Fig 0: sorting board, where every conic lives
rows = [
 ("circle", "e", "(0, 1) in the plain ruler; every\nno-wall dial after completing the square", "only with imaginary k", "p = 2; the hot limit\nof temperature", "round end (Pick = 0)"),
 ("ellipse", "e", "inside the parabola d > k²/4", "only with imaginary k", "—", "round end (same point)"),
 ("parabola", "p", "on the curve d = k²/4; (0, 0) is\naddition. The curve itself is a parabola", "k = 0 (plain addition)", "—", "round end (same point)"),
 ("hyperbola", "h", "outside, d < k²/4: relativity\n(0, −1/c²), probability (−1, 0)", "every real k ≠ 0;\nwall −1/k = an asymptote", "—", "round end (both branches)"),
 ("2 parallel lines", "pl", "size-1 level of a flat dial", "k = 0: x² = 1", "ends of ℓ^p are\npairs of these", "—"),
 ("2 crossing lines", "cl", "size-0 level of a two-wall dial:\nthe walls (light cone)", "x = 0 and x + ky = 0", "—", "—"),
 ("double line", "dl", "size-0 level of a flat dial\n(dual numbers' zero divisors)", "k = 0", "—", "—"),
 ("single point", "pt", "size-0 level of a no-wall dial:\nonly 0 can't be divided by", "—", "—", "—"),
 ("nothing real", "no", "size −1 of a no-wall dial;\nCayley–Klein's elliptic absolute", "—", "—", "—"),
]
fig = plt.figure(figsize=(15, 11.2)); 
cols = ["conic", "plane dial (k, d)", "k alone (d = 0)", "roundness dials", "shape dial"]
xs = [0.0, 0.17, 0.48, 0.67, 0.84]; ws = [0.17, 0.31, 0.19, 0.17, 0.16]
H = 1/(len(rows) + 1.7)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
for j, c in enumerate(cols): ax.text(xs[j] + 0.01, 1 - 0.6*H, c, fontsize=13, weight="bold", va="center")
for i, (name, g, pd, kk_, rd, sd) in enumerate(rows):
    y0 = 1 - (i + 1.2)*H - H
    ax.add_patch(Rectangle((0, y0), 1, H, color="#f4f6fb" if i % 2 == 0 else "white", zorder=0))
    col = {"e": EC, "p": PC, "h": HC}.get(g, "#555")
    gi = fig.add_axes([0.005, y0 + 0.1*H, 0.07, 0.8*H]); gi.set_aspect("equal"); gi.axis("off"); gi.set_xlim(-1.3, 1.3); gi.set_ylim(-1.3, 1.3)
    t_ = np.linspace(0, 2*np.pi, 200); u_ = np.linspace(-1.2, 1.2, 100)
    if g == "e": gi.plot(np.cos(t_), (0.95 if name == "circle" else 0.55)*np.sin(t_), color=col, lw=2)
    elif g == "p": gi.plot(u_, u_**2 - 1, color=col, lw=2)
    elif g == "h":
        v_ = np.linspace(-1, 1, 100); gi.plot(np.cosh(v_)*0.6, np.sinh(v_), color=col, lw=2); gi.plot(-np.cosh(v_)*0.6, np.sinh(v_), color=col, lw=2)
    elif g == "pl": gi.plot([-.5, -.5], [-1.2, 1.2], color=col, lw=2); gi.plot([.5, .5], [-1.2, 1.2], color=col, lw=2)
    elif g == "cl": gi.plot([-1.1, 1.1], [-1.1, 1.1], color=col, lw=2); gi.plot([-1.1, 1.1], [1.1, -1.1], color=col, lw=2)
    elif g == "dl": gi.plot([0, 0], [-1.2, 1.2], color=col, lw=5)
    elif g == "pt": gi.plot(0, 0, "o", color=col, ms=7)
    else: gi.text(0, 0, "∅", fontsize=20, ha="center", va="center", color=GREY)
    ax.text(0.08, y0 + H/2, name, fontsize=12, weight="bold", va="center", color=col)
    for j, txt in enumerate([pd, kk_, rd, sd], start=1):
        ax.text(xs[j] + 0.01, y0 + H/2, txt, fontsize=10.6, va="center", color="#222" if txt != "—" else "#bbb")
ax.text(0.01, 0.012, "Figure 0. Where every conic lives. Colours: blue = no wall (ellipse type), magenta = one double wall (parabola type), "
        "orange = two walls (hyperbola type).", fontsize=11.5)
fig.savefig("figures/fig0-sorting-board.png"); plt.close(fig)
print("board ok")
