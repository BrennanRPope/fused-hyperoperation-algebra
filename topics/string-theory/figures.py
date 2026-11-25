"""Figures for topics/string-theory (Nov 25, 2026). Run: python3 figures.py -> figures/*.png"""
import numpy as np, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Polygon, Rectangle
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREEN, RED, AMBER, GRAY, PURPLE = "#2b6cb0", "#2f855a", "#c53030", "#b7791f", "#718096", "#6b46c1"
def save(fig, name): fig.savefig(os.path.join(OUT, name), dpi=140, bbox_inches="tight"); plt.close(fig)
def box(ax, x, y, w, h, text, color, fs=10):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=color, ec="none", alpha=0.18))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc="none", ec=color, lw=1.5))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, wrap=True)
def arrow(ax, x1, y1, x2, y2, color=GRAY):
    ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="->", color=color, lw=1.5))

# s1: what string theory rests on
fig, ax = plt.subplots(figsize=(12, 6.2)); ax.set_xlim(0, 12); ax.set_ylim(0, 6.6); ax.axis("off")
ax.text(0.2, 6.35, "What string theory rests on: inputs (left), consistency demands (middle), what comes out (right)", fontsize=13, weight="bold")
ins = [("A string sweeps a surface;\naction = its area\n(Nambu-Goto, Polyakov)", BLUE), ("Two dials: string length\nl_s = sqrt(alpha')\nand coupling g_s = e^phi", BLUE),
       ("Quantum mechanics +\nspecial relativity\n(assumed, not derived)", BLUE), ("Supersymmetry on the\nworldsheet (for fermions\nand no tachyon)", AMBER)]
for i, (t, c) in enumerate(ins): box(ax, 0.2, 4.7 - 1.45*i, 3.2, 1.15, t, c)
mid = [("Worldsheet conformal\ninvariance survives quantizing\n-> D = 26 (or 10)", GREEN), ("No ghosts, no anomalies\n(Goddard-Thorn 1972;\nGreen-Schwarz 1984)", GREEN),
       ("Modular invariance of\nloop surfaces\n(cuts off the UV)", GREEN), ("Beta functions = 0\n(conformal invariance\non a curved background)", GREEN)]
for i, (t, c) in enumerate(mid): box(ax, 4.4, 4.7 - 1.45*i, 3.2, 1.15, t, c)
outs = [("A massless spin-2 state\n(the graviton) is always\nin the spectrum", PURPLE), ("Einstein's equations\n+ corrections in alpha'\n(Callan-Friedan-Martinec-Perry 1985)", PURPLE),
        ("G = (numbers) g_s^2 l_s^8 / V6\n(derived, not put in)", PURPLE), ("Gauge fields, chiral fermions;\nbut a landscape of vacua\n(conjectured ~10^500)", RED)]
for i, (t, c) in enumerate(outs): box(ax, 8.6, 4.7 - 1.45*i, 3.2, 1.15, t, c)
for i in range(4):
    arrow(ax, 3.45, 5.27 - 1.45*i, 4.35, 5.27 - 1.45*i); arrow(ax, 7.65, 5.27 - 1.45*i, 8.55, 5.27 - 1.45*i)
ax.text(0.2, 0.05, "Blue: assumed inputs.  Amber: extra assumption.  Green: consistency demands (proved in perturbation theory).  Purple: derived outputs.  Red: open or conjectured.", fontsize=9, color=GRAY)
save(fig, "s1-what-string-theory-rests-on.png")

# s2: critical dimension
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
n = np.arange(1, 30); ax = axs[0]
ax.plot(n, np.cumsum(n), "o-", color=GRAY, ms=3, label="plain partial sums 1+2+...+n (diverge)")
ax.axhline(-1/12, color=RED, lw=2, label="zeta-regularized value: -1/12")
ax.set_yscale("symlog", linthresh=0.05); ax.set_xlabel("number of vibration modes kept, n"); ax.set_ylabel("sum of mode numbers")
ax.set_title("Zero-point energy of the string's modes"); ax.legend(fontsize=9, loc="center right")
ax = axs[1]; Ds = np.arange(4, 30); ax.plot(Ds, (Ds - 2)/24, color=BLUE, label="bosonic: a = (D-2)/24")
ax.plot(Ds, (Ds - 2)*(1/24 + 1/48), color=GREEN, label="superstring NS: a = (D-2)/16")
ax.axhline(1, color=BLUE, ls=":"); ax.axhline(0.5, color=GREEN, ls=":")
ax.plot([26], [1], "o", color=BLUE, ms=9); ax.plot([10], [0.5], "o", color=GREEN, ms=9)
ax.annotate("D = 26: a = 1, the level-1\nvector is massless (needed\nfor Lorentz invariance)", (26, 1), (15, 1.25), arrowprops=dict(arrowstyle="->", color=BLUE), fontsize=9)
ax.annotate("D = 10: a = 1/2", (10, 0.5), (11.5, 0.2), arrowprops=dict(arrowstyle="->", color=GREEN), fontsize=9)
ax.set_xlabel("spacetime dimension D"); ax.set_ylabel("normal-ordering constant a"); ax.set_title("Why 26 and 10 (checks A1-A3)"); ax.legend(fontsize=9, loc="upper left")
save(fig, "s2-critical-dimension.png")

# s3: graviton forced
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.2]})
ax = axs[0]; Mx = np.zeros((24, 24))
for i in range(24):
    for j in range(24): Mx[i, j] = 0 if i == j else (1 if i < j else 2)
ax.imshow(Mx, cmap=matplotlib.colors.ListedColormap(["#d6bcfa", "#9ae6b4", "#fbd38d"]))
ax.set_title("Closed string, lowest massless level:\nleft mover (row) x right mover (column) = 24 x 24"); ax.set_xticks([]); ax.set_yticks([])
ax.text(17, 5, "symmetric\npart", ha="center", fontsize=9); ax.text(6, 18, "antisymmetric\npart", ha="center", fontsize=9)
ax = axs[1]; parts = [("graviton\n(symmetric, traceless)", 299, PURPLE), ("B-field\n(antisymmetric)", 276, GREEN), ("dilaton\n(trace)", 1, AMBER)]
ax.barh([0, 1, 2], [p[1] for p in parts], color=[p[2] for p in parts]); ax.set_yticks([0, 1, 2]); ax.set_yticklabels([p[0] for p in parts])
for i, p in enumerate(parts): ax.text(p[1] + 5, i, str(p[1]), va="center")
ax.invert_yaxis(); ax.set_xlabel("number of states (total 576)")
ax.set_title("A massless spin-2 particle is unavoidable (check A4);\nthe dilaton's value is the coupling g_s"); 
save(fig, "s3-graviton-is-forced.png")

# s4: G derived
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axs[0]; gs = np.logspace(-3, 0, 200); ls = np.logspace(-1, 1, 200); GS, LS = np.meshgrid(gs, ls)
Gv = (2*np.pi)**7*GS**2*LS**8/(16*np.pi)
cs = ax.contour(GS, LS, np.log10(Gv), levels=np.arange(-4, 12, 2), cmap="viridis"); ax.clabel(cs, fmt=lambda v: f"G=1e{int(v)}", fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("string coupling g_s"); ax.set_ylabel("string length l_s (any unit)")
ax.set_title("Newton's constant in 10D = (2pi)^7 g_s^2 l_s^8 / 16pi\nlines of equal G (check A5)")
ax = axs[1]; N = np.arange(1, 60); ax.loglog(N, np.pi/(2*N**2), "o", color=PURPLE, ms=4)
ax.set_xlabel("number of colours N in the gauge theory on the boundary"); ax.set_ylabel("G5 / L^3  (Newton's constant in AdS units)")
ax.set_title("AdS/CFT: G = pi / (2 N^2) in AdS units (check A7)\nthe best-known 'G dial'")
ax.text(3, 1e-3, "each handle on a worldsheet\ncosts g_s^2 ~ 1/N^2,\nthe same power as G", fontsize=9, color=PURPLE)
save(fig, "s4-G-is-derived.png")

# b1: T-duality spectrum against log R
fig, ax = plt.subplots(figsize=(11, 5)); u = np.linspace(-2, 2, 400); R = np.exp(u)
for nn, w, col, lab in [(1, 0, BLUE, "momentum mode n=1 (mass 1/R)"), (0, 1, GREEN, "winding mode w=1 (mass R/alpha')"), (2, 0, BLUE, None), (0, 2, GREEN, None)]:
    ax.plot(u, np.sqrt((nn/R)**2 + (w*R)**2), color=col, lw=2 if nn + w == 1 else 1, ls="-" if nn + w == 1 else "--", label=lab)
m = np.sqrt(np.maximum((1/R)**2 + R**2 - 2, 0)); ax.plot(u, m, color=RED, lw=2, label="n = w = 1 with N = 1, Ntilde = 0: massless only at R = sqrt(alpha')")
ax.axvline(0, color=GRAY, ls=":"); ax.text(0.05, 6.3, "self-dual radius\nR = sqrt(alpha')\n= mult-add identity a = 0", fontsize=9)
ax.annotate("", (1.5, 5.8), (-1.5, 5.8), arrowprops=dict(arrowstyle="<->", color=GRAY)); ax.text(-0.9, 6.0, "T-duality: u -> -u (mirror)", fontsize=9)
ax.set_ylim(0, 7.5); ax.set_xlabel("u = log(R / sqrt(alpha'))  = log(1 + a), the mult-add flat coordinate at k = 1")
ax.set_ylabel("mass (units of 1/sqrt(alpha'))"); ax.set_title("T-duality: the spectrum is the same at R and alpha'/R (checks T1, T2, T5, T6)"); ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
save(fig, "b1-t-duality-mirror.png")

# b2: radius change is a boost
fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
for ax, Rv, title in [(axs[0], 0.7, "R = 0.7 sqrt(alpha')"), (axs[1], 1.6, "R = 1.6 sqrt(alpha')")]:
    xs = np.linspace(-4, 4, 400)
    for c in (-4, -2, 2, 4):
        if c > 0: ax.plot(np.sqrt(xs**2 + c), xs, color=GRAY, lw=0.6); ax.plot(-np.sqrt(xs**2 + c), xs, color=GRAY, lw=0.6)
        else: ax.plot(xs, np.sqrt(xs**2 - c), color=GRAY, lw=0.6); ax.plot(xs, -np.sqrt(xs**2 - c), color=GRAY, lw=0.6)
    for nn in range(-3, 4):
        for w in range(-3, 4):
            pl, pr = (nn/Rv + w*Rv)/np.sqrt(2), (nn/Rv - w*Rv)/np.sqrt(2)
            col = PURPLE if nn*w > 0 else (GREEN if nn*w < 0 else BLUE)
            ax.plot(pl, pr, "o", color=col, ms=5)
            if (nn, w) in [(1, 1), (1, 0), (0, 1)]: ax.annotate(f"(n,w)=({nn},{w})", (pl, pr), (pl + 0.15, pr + 0.25), fontsize=8)
    ax.plot([-4, 4], [-4, 4], color=AMBER, lw=0.8, ls="--"); ax.plot([-4, 4], [4, -4], color=AMBER, lw=0.8, ls="--")
    ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal"); ax.set_xlabel("p_L"); ax.set_ylabel("p_R"); ax.set_title(title)
fig.suptitle("Changing the radius slides every lattice point along its hyperbola p_L^2 - p_R^2 = 2nw: a Lorentz boost with rapidity log R (check T4).\n"
             "The radius is the hyperbolic Cayley-Klein dial; T-duality is the mirror across the light line (n <-> w).", fontsize=10)
save(fig, "b2-radius-is-a-boost.png")

# b3: half-plane: Borel, S-duality, fundamental domain, dial lines
fig, ax = plt.subplots(figsize=(10.5, 6.2)); th = np.linspace(0, np.pi, 300)
x = np.linspace(-0.5, 0.5, 200); top = 3.0
ax.fill_between(x, np.sqrt(1 - x**2), top, color=AMBER, alpha=0.25, label="fundamental domain F of SL(2,Z):\none copy of every torus shape tau")
for kk, col in [(1.0, BLUE), (0.5, BLUE), (-0.5, GREEN), (-1.0, GREEN)]:
    ss = np.linspace(-1/kk, 3, 200) if kk > 0 else np.linspace(-3, -1/kk, 200)
    ax.plot(ss, 1 + kk*ss, color=col, lw=1.2); ax.plot([-1/kk], [0], "o", color=col, ms=5)
    ax.text(-1/kk, -0.18, f"wall -1/k\n(k={kk:g})", ha="center", fontsize=8, color=col)
ax.axhline(1, color=GRAY, lw=1, ls="--"); ax.text(1.9, 1.05, "k = 0: shift C0 (axion)", fontsize=8, color=GRAY)
ax.axvline(0, color=PURPLE, lw=1.5); ax.text(0.04, 2.75, "pure multiplication:\nchange g_s at C0 = 0", fontsize=8, color=PURPLE)
ax.plot(np.cos(th), np.sin(th), color=RED, lw=1); ax.plot([0], [1], "*", color=RED, ms=14)
ax.annotate("tau = i: dial's 'do nothing' point\n= S-duality's self-dual coupling g_s = 1", (0, 1), (0.9, 1.9), arrowprops=dict(arrowstyle="->", color=RED), fontsize=9)
ax.annotate("", (0, 0.45), (0, 2.2), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.5)); ax.text(-1.4, 0.45, "S: g_s <-> 1/g_s\n(= group inverse\non this line only)", fontsize=8, color=RED)
ax.set_xlim(-3, 3); ax.set_ylim(-0.4, top); ax.set_xlabel("Re tau = s (shift; for IIB the axion C0; for a torus B)"); ax.set_ylabel("Im tau = lambda (stretch; for IIB 1/g_s)")
ax.set_title("One half-plane, three readings: the order dial's group, the IIB coupling tau, the torus modulus (checks S1-S5)"); ax.legend(loc="upper left", fontsize=8)
save(fig, "b3-half-plane-dualities.png")

# b4: modular invariance cuts the UV
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]; x = np.linspace(-0.5, 0.5, 200)
ax.fill_between(x, 0, 2.5, color=GRAY, alpha=0.15, label="strip: a field-theory loop would integrate here")
ax.fill_between(x, np.sqrt(1 - x**2), 2.5, color=AMBER, alpha=0.5, label="F: string theory integrates only here")
ax.fill_between(x, 0, 0.25, color=RED, alpha=0.3, label="UV region (tau2 -> 0): never reached")
for g in [lambda t: t + 1, lambda t: t - 1]:
    pass
ax.set_xlim(-0.6, 0.6); ax.set_ylim(0, 2.5); ax.set_xlabel("tau1"); ax.set_ylabel("tau2 (loop 'proper time')"); ax.legend(fontsize=8, loc="upper right")
ax.set_title("Modular invariance removes the ultraviolet\n(area of F = pi/3; strip area diverges; check M2)")
ax = axs[1]
from numpy import exp, pi
def eta(tau, terms=200):
    q = np.exp(2j*pi*tau); n = np.arange(1, terms + 1); return np.exp(2j*pi*tau/24)*np.prod(1 - q[..., None]**n, axis=-1)
X, Y = np.meshgrid(np.linspace(-1.5, 1.5, 300), np.linspace(0.25, 2.2, 200)); T = X + 1j*Y
Z = np.log10(Y**12*np.abs(eta(T))**48)
im = ax.pcolormesh(X, Y, Z, cmap="magma", shading="auto"); plt.colorbar(im, ax=ax, label="log10 of tau2^12 |eta|^48")
for s0 in (-1, 0, 1): ax.plot(s0 + np.cos(th), np.sin(th), color="white", lw=0.7)
for s0 in (-1.5, -0.5, 0.5, 1.5): ax.axvline(s0, color="white", lw=0.7, ls=":")
ax.set_xlabel("tau1"); ax.set_ylabel("tau2"); ax.set_title("The torus weight repeats under tau -> tau + 1\nand tau -> -1/tau (check M1)")
save(fig, "b4-modular-invariance-cuts-the-uv.png")

# b5: Hagedorn
Nmax = 1500; p24 = [1] + [0]*Nmax
for nn in range(1, Nmax + 1):
    for _ in range(24):
        for j in range(nn, Nmax + 1): p24[j] += p24[j - nn]
import math
Ns = np.arange(2, Nmax + 1); lp = np.array([math.log(p24[N]) for N in Ns])
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axs[0]; ax.plot(np.sqrt(Ns - 1), lp, color=BLUE, label="log(number of states) at mass M = sqrt(N-1)")
ax.plot(np.sqrt(Ns - 1), 4*np.pi*np.sqrt(Ns) - 6.75*np.log(Ns) - 1.6, "--", color=RED, label="4 pi sqrt(alpha') M - (27/4) log N + const")
ax.set_xlabel("mass M (units of 1/sqrt(alpha'))"); ax.set_ylabel("log of state count"); ax.legend(fontsize=9); ax.set_title("States grow exponentially with mass (check M3)")
ax = axs[1]; M = np.sqrt(Ns - 1)
for T, col in [(0.9/(4*np.pi), GREEN), (1.0/(4*np.pi), AMBER), (1.1/(4*np.pi), RED)]:
    ax.plot(M, lp - M/T, color=col, label=f"T = {T*4*np.pi:.1f} T_H")
ax.set_xlabel("mass M"); ax.set_ylabel("log(states x Boltzmann weight)"); ax.set_title("Hagedorn wall T_H = 1/(4 pi sqrt(alpha')):\nbelow it the sum converges, above it diverges (M4)")
ax.legend(fontsize=9)
save(fig, "b5-hagedorn-wall.png")

# b6: minimum length and the (k,d) plane
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]; dp = np.logspace(-1.3, 1.3, 400); ax.loglog(dp, 1/(2*dp), color=GRAY, ls="--", label="Heisenberg: hbar/(2 dp)")
ax.loglog(dp, dp/2, color=GREEN, ls="--", label="string stretching: alpha' dp/(2 hbar)"); ax.loglog(dp, 1/(2*dp) + dp/2, color=BLUE, lw=2, label="string GUP (sum)")
ax.plot([1], [1], "o", color=RED); ax.annotate("smallest dx = sqrt(alpha') = l_s\nat dp = hbar/sqrt(alpha')", (1, 1), (1.6, 0.2), arrowprops=dict(arrowstyle="->", color=RED), fontsize=9)
ax.set_xlabel("momentum spread dp (units hbar/l_s)"); ax.set_ylabel("position spread dx (units l_s)"); ax.legend(fontsize=8); ax.set_title("Minimum length; the curve is mirror-symmetric\nunder dp -> hbar^2/(alpha' dp) (checks L1, L2)")
ax = axs[1]; kk = np.linspace(-2.2, 2.2, 300)
for kap, col in [(-0.5, BLUE), (0, GRAY), (0.5, GREEN)]: ax.plot(kk, kap + kk**2/4, color=col, lw=1, label=f"kappa = {kap:g}")
ax.fill_between(kk, kk**2/4, 2.2, color=GREEN, alpha=0.1); ax.fill_between(kk, -1.5, kk**2/4, color=BLUE, alpha=0.08)
ax.plot([0], [1], "s", color=GREEN, ms=9); ax.annotate("Kempf-Mangano-Mann minimum length\n(0, beta): elliptic, no wall,\nflat coordinate bounded", (0, 1), (0.3, 1.6), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.plot([-1], [0], "o", color=BLUE, ms=9); ax.annotate("kappa-Poincare maximum momentum\n(-1/kappa_P, 0): hyperbolic, wall", (-1, 0), (-2.1, -1.2), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.plot([0], [-1], "^", color=PURPLE, ms=9); ax.annotate("velocity addition (0, -1/c^2)", (0, -1), (0.4, -0.8), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.5, 2.2); ax.set_xlabel("k"); ax.set_ylabel("d"); ax.legend(fontsize=8, loc="lower right")
ax.set_title("On the repo's (k, d) plane, minimum length and\nmaximum momentum sit on opposite sides (L3-L6, ours)")
save(fig, "b6-minimum-length-on-the-plane.png")

# b7: Veneziano softness
from scipy.special import loggamma
ss = np.arange(4, 300, 2) + 0.5
la = np.array([(loggamma(complex(-1 - s)) + loggamma(complex(-1 + s/2)) - loggamma(complex(-2 - s/2))).real for s in ss])
fig, ax = plt.subplots(figsize=(10, 4.6))
ax.plot(ss, la, color=BLUE, lw=2, label="string (Veneziano) amplitude at fixed angle, log|A|")
ax.plot(ss, -np.log(2)*ss + la[-1] + np.log(2)*ss[-1], "--", color=RED, label="slope -ln 2 per unit alpha' s: exponential fall")
ax.plot(ss, -2*np.log(ss) + la[0] + 2*np.log(ss[0]), ":", color=GRAY, lw=2, label="a field theory: power-law fall (e.g. 1/s^2)")
ax.set_xlabel("collision energy squared, alpha' s"); ax.set_ylabel("log |amplitude|"); ax.legend(fontsize=9)
ax.set_title("Strings are soft at high energy: hard scattering dies exponentially, the root of 'no distance below l_s' (check L7)")
save(fig, "b7-strings-are-soft.png")

# b8: genus dial
fig, ax = plt.subplots(figsize=(12, 4.6)); ax.set_xlim(0, 12); ax.set_ylim(0, 4.6); ax.axis("off")
def surf(cx, cy, g):
    ax.add_patch(Ellipse((cx, cy), 2.4 if g else 1.6, 1.3 if g else 1.6, fc="#bee3f8", ec=BLUE, lw=1.5))
    for h in range(g):
        hx = cx + (h - (g - 1)/2)*0.85; ax.add_patch(Ellipse((hx, cy), 0.5, 0.22, fc="white", ec=BLUE, lw=1.2))
for i, (g, w) in enumerate([(0, "g_s^-2"), (1, "g_s^0"), (2, "g_s^+2")]):
    surf(1.6 + 3.3*i, 3.1, g); ax.text(1.6 + 3.3*i, 1.95, f"genus {g}: chi = {2 - 2*g}\nweight {w}", ha="center", fontsize=10)
ax.text(10.3, 3.5, "each handle\ncosts g_s^2", fontsize=11, color=PURPLE, ha="center"); ax.text(10.3, 2.7, "G ~ g_s^2 l_s^8", fontsize=11, color=PURPLE, ha="center")
ax.text(0.2, 1.0, "Matrix version ('t Hooft 1974, Harer-Zagier 1986): (1/N) E tr X^4 = 2 + 1/N^2,  (1/N) E tr X^6 = 5 + 10/N^2, ...  "
        "planar diagrams + 1/N^2 per handle (checks N1, N2).", fontsize=10)
ax.text(0.2, 0.45, "The weight comes from the worldsheet's own Einstein term: (phi/4pi) x integral of curvature = phi x chi (Gauss-Bonnet, N4). "
        "So spacetime's G is built from the worldsheet's 'Newton constant'.", fontsize=10)
ax.set_title("The genus dial: g_s counts joinings and splittings; 1/N^2 does the same for matrices", fontsize=13)
save(fig, "b8-genus-dial.png")

# b9: borrow scorecard
rows = [("T-duality R <-> alpha'/R", "shape matches: mult-add (k=1) inverse, u -> -u", "borrowable framing", GREEN),
        ("Radius modulus", "rapidity of an O(1,1) boost = flat coordinate", "borrowable framing", GREEN),
        ("S-duality, SL(2,Z)", "dial group = Borel part; S itself is outside", "needs full SL(2) + lattice", AMBER),
        ("Modular invariance", "integrate over H/SL(2,Z): UV cut off", "needs the lattice (same gap as path-integral page)", AMBER),
        ("Minimum length / GUP", "elliptic point (0, beta) vs kappa-P wall", "ours, coordinate-dependent", AMBER),
        ("Hagedorn wall", "same mechanism as the primon gas wall", "analogy, known (Julia 1990)", AMBER),
        ("alpha', g_s as dials", "alpha': contraction-like; g_s: asymptotic series", "partly analogy", AMBER),
        ("G from g_s^2 l_s^8 (or 1/N^2)", "needs a length dial + a handle-counting dial", "lead: 1/N of matrix chains", PURPLE),
        ("Worldsheet CFT, D = 26/10, SUSY,\nanomalies, beta = Einstein, branes", "field theory on surfaces", "chains can't express", RED)]
fig, ax = plt.subplots(figsize=(13, 6)); ax.set_xlim(0, 13); ax.set_ylim(0, len(rows) + 1); ax.axis("off")
for j, h in enumerate(["String theory idea", "How it meets our framework", "Verdict"]): ax.text([0.1, 4.3, 9.2][j], len(rows) + 0.4, h, weight="bold")
for i, (a, b, c, col) in enumerate(rows):
    y = len(rows) - i - 0.4; ax.add_patch(Rectangle((0, y - 0.42), 13, 0.84, fc=col, alpha=0.12, ec="none"))
    ax.text(0.1, y, a, va="center", fontsize=9); ax.text(4.3, y, b, va="center", fontsize=9); ax.text(9.2, y, c, va="center", fontsize=9, color=col, weight="bold")
ax.set_title("What we could borrow from string theory, judged on its own merits", fontsize=13)
save(fig, "b9-borrow-scorecard.png")
print("figures written:", sorted(os.listdir(OUT)))
