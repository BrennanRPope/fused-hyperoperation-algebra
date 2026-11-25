"""Figures for why-i/e-i-pi/asymptotic-safety. Run: python3 figures.py (numpy only)."""
import numpy as np, os, sys
from math import gamma as Gamma, pi, erf
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from svgkit import *
_lines = Fig.lines
Fig.lines = lambda self, x, y, s, **kw: _lines(self, x, y, "\n".join(s) if isinstance(s, list) else s, **kw)
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True); P = lambda n: os.path.join(OUT, n)

def rk4(f, v, dt, n):
    v = np.array(v, float); out = [v.copy()]
    for _ in range(n):
        k1 = f(v); k2 = f(v + dt/2*k1); k3 = f(v + dt/2*k2); k4 = f(v + dt*k3)
        v = v + dt/6*(k1 + 2*k2 + 2*k3 + k4); out.append(v.copy())
    return np.array(out)

# ---------- a1: the unit lattice ----------
f = Fig(1180, 600); f.title("Units form a lattice; setting ħ = c = 1 collapses it to one line, and G sits two steps below a pure number", size=17)
ax = Axes(f, 80, 110, 1020, 60, (-3.5, 3.5), (0, 1))
f.line(ax.X(-3.5), ax.Y(0.5), ax.X(3.5), ax.Y(0.5), color=DARK, width=2)
for v in range(-3, 4):
    f.line(ax.X(v), ax.Y(0.5) - 6, ax.X(v), ax.Y(0.5) + 6, color=DARK); f.text(ax.X(v), ax.Y(0.5) - 14, f"{v:+d}" if v else "0", size=13, anchor="middle")
f.text(ax.X(-3.5), ax.Y(0.5) - 34, "mass dimension (ħ = c = 1)", size=13, color=GRAY)
pts = [(-2, "G", RED, 0), (2, "Λ", BLUE, 0), (1, "probe scale k", GREEN, 0), (0, "ħ, c", GRAY, 0), (0, "g = G k²", PURPLE, 1), (0, "G Λ", ORANGE, 2)]
for v, lab, col, row in pts:
    f.circle(ax.X(v), ax.Y(0.5), 7, color=col, fill=col)
    f.text(ax.X(v), ax.Y(0.5) + 32 + 22*row, lab, size=14, anchor="middle", color=col, weight="bold")
f.line(ax.X(-2), ax.Y(0.5) + 92, ax.X(-0.08), ax.Y(0.5) + 92, color=PURPLE, width=1.5, arrow=True)
f.text(ax.X(-1), ax.Y(0.5) + 110, "multiply by k² (two steps up)", size=12, anchor="middle", color=PURPLE)
# index panel
f.text(80, 345, "Do three constants reach every unit with whole powers? The lattice index says (check U1–U2):", size=14, weight="bold")
rows = [("c, ħ, m (a mass)", "index 1", "every unit is a whole-power product: no roots (Compton units)", GREEN),
        ("c, ħ, G", "index 2", "a length needs a square root: ℓ_P = √(ħG/c³)", RED),
        ("c, ħ, Λ", "index 2", "again a square root: L = √(3/Λ)", RED)]
for i, (a, b2, cc, col) in enumerate(rows):
    y = 380 + i*48
    f.rect(80, y, 1020, 40, color=col, fill=col, opacity=0.08, rx=6)
    f.text(100, y + 26, a, size=14, weight="bold"); f.text(330, y + 26, b2, size=14, color=col, weight="bold"); f.text(450, y + 26, cc, size=13)
f.text(80, 545, "Reading (ours): gravity's scales sit at 'odd' places in the unit lattice, which is why G = c³ℓ²/ħ comes as a square (g-dials S7).", size=13, color=GRAY)
f.text(80, 568, "Known: units as a free abelian group (Buckingham π, the Lean formalisation of dimensional analysis); the √ in Planck units.", size=12, color=GRAY)
f.save(P("a1-unit-lattice.svg"))

# ---------- a2: different units at different levels ----------
gs, g0 = 0.707, 1e-6
t = np.linspace(-2, 14, 400)
g = gs/(1 + (gs/g0 - 1)*np.exp(-2*t))
Gk = g/(g0*np.exp(2*t))
ds = 4/(1 + g/gs)
x = t/np.log(10)
f = Fig(1180, 760); f.title("Brennan's idea in the literature: G has different units at different levels (toy flow, shape only)", size=17)
ax = Axes(f, 110, 80, 700, 280, (x[0], x[-1]), (-7, 0.8))
ax.frame(xticks=[], yticks=[(0, "1"), (-2, "10⁻²"), (-4, "10⁻⁴"), (-6, "10⁻⁶")], ylabel="G_k / G₀")
ax.plot(x, np.log10(Gk), color=RED, width=3)
xc = np.log10(np.sqrt(gs/g0))
ax.vline(xc, color=GRAY, dash="4,4")
ax.text(x[0] + 0.1, 0.25, "classical: G constant, units length² (in ħ = c = 1)", size=12, color=RED)
ax.text(xc + 0.25, -4.6, "fixed point: G_k = g*/k², falls as the", size=12, color=RED)
ax.text(xc + 0.25, -5.2, "probe scale grows; g = G k² is a pure number", size=12, color=RED)
ax2 = Axes(f, 110, 420, 700, 240, (x[0], x[-1]), (1.7, 4.3))
ax2.frame(xticks=[(x[0], "low energy"), (xc, "Planck scale"), (x[-1], "high energy")], yticks=[2, 3, 4], ylabel="spectral dimension", xlabel="log of the probe scale k")
ax2.plot(x, ds, color=BLUE, width=3); ax2.vline(xc, color=GRAY, dash="4,4")
ax2.text(x[0] + 0.1, 4.08, "4: ordinary spacetime", size=12, color=BLUE); ax2.text(x[-1] - 0.1, 2.12, "2: where G is a pure number", size=12, color=BLUE, anchor="end")
f.lines(850, 110, ["Known (Lauscher–Reuter 2005):", "at the fixed point the graviton's", "anomalous dimension is η_N = −2.", "That cancels G's two units and", "makes the propagator 1/p⁴, so a", "random walk sees dimension 2.", "", "Two spacetimes, one theory:", "4D with dimensionful G at large", "distances, 2D with dimensionless G", "at short ones. The Planck length", "sets where one turns into the other.", "", "Reading (ours): the 2 units G loses", "are exactly two steps of the unit", "lattice (picture a1)."], size=13)
f.text(110, 735, "Toy β(g) = 2g − 2g²/g*, g* = 0.707 (the Einstein–Hilbert value, check A6); d_s = 2d/(2 − η_N) (check R3). Checks A1, A7, D1, D6.", size=12, color=GRAY)
f.save(P("a2-units-change-with-scale.svg"))

# ---------- a3: the flow is mult-add ----------
f = Fig(1180, 560); f.title("In the variable 1/g, running G is our weighted circle operation: equal RG steps, a wall at the fixed point", size=17)
ax = Axes(f, 80, 110, 1020, 80, (0, 12), (0, 1))
f.text(80, 95, "Gravity (G has units): u = 1/g, step u → u ∘ₖ s with k = −g*", size=14, weight="bold", color=PURPLE)
f.line(ax.X(0), ax.Y(0.4), ax.X(12), ax.Y(0.4), color=DARK, width=2)
uw = 1/0.707
f.line(ax.X(uw), ax.Y(0.0), ax.X(uw), ax.Y(1.0), color=RED, width=3)
f.text(ax.X(uw), ax.Y(0.0) + 20, "wall −1/k = 1/g*: the UV fixed point", size=13, color=RED, anchor="middle")
u = 11.5
for i in range(7):
    un = uw + (u - uw)*np.exp(-0.35*2)
    f.line(ax.X(u), ax.Y(0.55), ax.X(un) + 3, ax.Y(0.55), color=PURPLE, width=2, arrow=True)
    f.circle(ax.X(u), ax.Y(0.4), 4, color=PURPLE, fill=PURPLE); u = un
f.text(ax.X(11.8), ax.Y(0.75), "← toward high energy (each arrow = the same RG time)", size=12, anchor="end", color=GRAY)
f.text(ax.X(12), ax.Y(0.4) + 34, "u → ∞: g → 0, the Gaussian point (free gravitons)", size=12, anchor="end", color=GRAY)
ax2 = Axes(f, 80, 310, 1020, 80, (-6, 6), (0, 1))
f.text(80, 295, "A dimensionless coupling (QCD-like, check R2): 1/α steps by plain addition, k = 0, no wall", size=14, weight="bold", color=GREEN)
f.line(ax2.X(-6), ax2.Y(0.4), ax2.X(6), ax2.Y(0.4), color=DARK, width=2)
for i in range(7):
    a0 = -5 + 1.5*i
    f.line(ax2.X(a0), ax2.Y(0.55), ax2.X(a0 + 1.5) - 3, ax2.Y(0.55), color=GREEN, width=2, arrow=True)
    f.circle(ax2.X(a0), ax2.Y(0.4), 4, color=GREEN, fill=GREEN)
f.lines(80, 440, ["Ours (checks A2–A5): with β(g) = (d − 2) g − w g², u = 1/g obeys du/dt = −(d − 2) u + w, so each RG step is u ∘ₖ s(t) = u + s + k u s with k = −g*.",
                  "Steps compose by the same rule, the wall is the UV fixed point, and RG time is the flat coordinate log|1 + ku|/k.",
                  "The rate d − 2 is G's canonical units: as d → 2 the dial k = −g* → 0 and gravity's flow becomes the plain-addition flow below.",
                  "Limit: exact for this one-loop shape; the full Einstein–Hilbert flow is two-dimensional (picture a5)."], size=13)
f.save(P("a3-flow-is-mult-add.svg"))

# ---------- a4: the d - 2 dial ----------
f = Fig(1180, 520); f.title("Spacetime dimension minus two is the dial strength: the wall merges with the Gaussian point as d → 2", size=17)
ax = Axes(f, 100, 80, 480, 340, (0, 2.2), (0, 0.8))
ax.frame(xticks=[0, 0.5, 1, 1.5, 2], yticks=[0, 0.2, 0.4, 0.6, 0.8], xlabel="d − 2 (G's mass dimension, sign flipped)", ylabel="fixed point g* = |k|")
e = np.linspace(0, 2.1, 50); ax.plot(e, e*0.707/2, color=RED, width=3)
ax.dot(2, 0.707, r=6, color=RED); ax.text(1.95, 0.74, "d = 4", size=12, anchor="end", color=RED)
ax.dot(0, 0, r=6, color=GREEN); ax.text(0.05, 0.04, "d = 2: G is a pure number, k = 0", size=12, color=GREEN)
ax2 = Axes(f, 680, 80, 440, 340, (-2, 14), (0, 1.1))
ax2.frame(xticks=[], yticks=[0, 0.5, 1], xlabel="RG time (log k)", ylabel="g / g*")
for ee, col, yl in [(2, RED, 0.62), (1, ORANGE, 0.45), (0.3, PURPLE, 0.28)]:
    tt = np.linspace(-2, 14, 300); gg = 1/(1 + 999*np.exp(-ee*(tt + 2)))
    ax2.plot(tt, gg, color=col, width=2.5); ax2.text(13.8, yl if ee != 2 else 1.03, f"d − 2 = {ee:g}", size=12, anchor="end", color=col)
ax2.hline(1, color=GRAY, dash="4,4")
f.text(100, 470, "Known: Weinberg's 2 + ε fixed point has g* = O(ε) (1979). Reading (ours): ε = d − 2 is the strength of the mult-add dial on G's flow, and the", size=13, color=GRAY)
f.text(100, 492, "UV fixed point's η_N = 2 − d is the flow undoing exactly that much dimension: at the fixed point gravity behaves as if d were 2. Check A5.", size=13, color=GRAY)
f.save(P("a4-dimension-is-the-dial.svg"))

# ---------- a5: Einstein-Hilbert flow plane ----------
d = 4; pre = (4*pi)**(1 - d/2)
Phi = lambda p, n, w: 1/Gamma(n + 1)/(1 + w)**p
Phit = lambda p, n, w: 1/Gamma(n + 2)/(1 + w)**p
B1 = lambda l: pre/3*(d*(d+1)*Phi(1, d/2-1, -2*l) - 6*d*(d-1)*Phi(2, d/2, -2*l) - 4*d*Phi(1, d/2-1, 0) - 24*Phi(2, d/2, 0))
B2 = lambda l: -pre/6*(d*(d+1)*Phit(1, d/2-1, -2*l) - 6*d*(d-1)*Phit(2, d/2, -2*l))
eta = lambda g, l: g*B1(l)/(1 - g*B2(l))
def beta(v):
    g, l = v; e = eta(g, l)
    bl = -(2 - e)*l + 0.5*g*pre*(2*d*(d+1)*Phi(1, d/2, -2*l) - 8*d*Phi(1, d/2, 0) - d*(d+1)*e*Phit(1, d/2, -2*l))
    return np.array([(d - 2 + e)*g, bl])
def jac(v, h=1e-7):
    J = np.zeros((2, 2))
    for jj in range(2):
        dv = np.zeros(2); dv[jj] = h; J[:, jj] = (beta(v + dv) - beta(v - dv))/(2*h)
    return J
v = np.array([0.7, 0.2])
for _ in range(60): v = v - np.linalg.solve(jac(v), beta(v))
w, Vv = np.linalg.eig(jac(v))
f = Fig(1180, 640); f.title("The real (Einstein–Hilbert) flow: every trajectory spirals out of one fixed point, then runs to classical gravity", size=17)
ax = Axes(f, 100, 70, 640, 480, (-0.3, 0.5), (0, 1.0))
ax.frame(xticks=[-0.2, 0, 0.2, 0.4], yticks=[0, 0.25, 0.5, 0.75, 1.0], xlabel="λ = Λ/k² (cosmological constant in probe units)", ylabel="g = G k²")
cols = [BLUE, BLUE, PURPLE, ORANGE, ORANGE, BLUE]
for ph, col in zip([0.0, 0.8, 1.6, 2.2, 2.6, 3.6], cols):
    v0 = v + 1e-4*np.real(Vv[:, 0]*np.exp(1j*ph))
    T = rk4(beta, v0, -0.005, 4000)
    ok = (T[:, 1] < 0.48) & (T[:, 1] > -0.3) & (T[:, 0] < 1.0) & (T[:, 0] > 0)
    jump = np.r_[np.hypot(np.diff(T[:, 0]), np.diff(T[:, 1])) > 0.02, True]
    ok = ok & ~np.cumsum(jump[::-1])[::-1].astype(bool) | ok & (np.cumsum(jump) == 0)
    bad = np.where(~ok)[0]; stop = bad[0] if len(bad) else len(T)
    ax.plot(T[:stop, 1], T[:stop, 0], color=col, width=1.8)
ax.dot(v[1], v[0], r=7, color=RED); f.line(ax.X(v[1]) + 6, ax.Y(v[0]) - 6, ax.X(0.33), ax.Y(0.86), color=RED, width=1); ax.text(0.335, 0.865, "fixed point (0.707, 0.193)", size=12, color=RED)
ax.dot(0, 0, r=7, color=GREEN); ax.text(0.01, 0.03, "Gaussian point", size=12, color=GREEN)
ax.vline(0.5, color=GRAY, dash="3,3"); ax.text(0.49, 0.95, "λ = ½", size=11, color=GRAY, anchor="end")
f.lines(780, 100, ["Reproduced (check A6): the optimised-cutoff", "fixed point and exponents θ = 1.475 ± 3.043 i.", "",
                   "Complex θ means the flow into the fixed point", "both shrinks and turns: z → e^(−θt) z (A8).", "That is e and i together (reading ours).", "",
                   "Leaving the fixed point toward low energy,", "each trajectory ends where g·λ = GΛ is", "constant: its label is the g-dials number", "GħΛ/c³ (A9; Reuter–Weyer 2004).", "",
                   "Blue: Λ < 0. Orange: Λ > 0 (these hit the", "λ = ½ wall unless tuned close to the", "separatrix, as our universe's must be).", "Purple: near the separatrix."], size=13)
f.text(100, 610, "Flow equations: Reuter–Saueressig form with the Litim cutoff, d = 4; truncated, and fixed-point values depend on the scheme.", size=12, color=GRAY)
f.save(P("a5-einstein-hilbert-flow.svg"))

# ---------- a6: dimension flows are straight lines in one coordinate ----------
f = Fig(1180, 660); f.title("Three published dimension flows are exactly straight in one coordinate: two walls, scale as rapidity", size=17)
axL = Axes(f, 90, 80, 470, 400, (-3, 3), (1.7, 4.3))
axL.frame(xticks=[(-3, "short"), (0, ""), (3, "long")], yticks=[2, 3, 4], xlabel="log scale (each curve centred)", ylabel="dimension")
s = np.linspace(-3, 3, 200)
ys_as = 2 + 2/(1 + 2*np.exp(-2*s))            # toy AS: d_s between 2 and 4
D0, a_ = 4.02 - 119/54, 4.02
ys_cdt = D0 + (a_ - D0)*(np.exp(s*np.log(10))/(1 + np.exp(s*np.log(10))))
alp = 0.5; y_ = np.exp(s*np.log(10)*(1 - alp))*alp; dH = (alp + y_)/(1 + y_)
ys_cal = 2 + 2*(dH - alp)/(1 - alp)
ys_erf = np.array([3 + erf(z) for z in s])
for ys, col, lab in [(ys_as, PURPLE, "asymptotic safety (toy flow)"), (ys_cdt, BLUE, "CDT fit 4.02 − 119/(54 + σ)"), (ys_cal, GREEN, "Calcagni multiscale (rescaled to 2–4)"), (ys_erf, GRAY, "an erf-shaped curve (not from physics)")]:
    axL.plot(s, ys, color=col, width=2.5)
axR = Axes(f, 660, 80, 470, 400, (-3, 3), (-8, 8))
axR.frame(xticks=[(-3, "short"), (3, "long")], yticks=[-8, -4, 0, 4, 8], xlabel="log scale", ylabel="log[(d − d_UV)/(d_IR − d)]")
for ys, col, lo, hi in [(ys_as, PURPLE, 2, 4), (ys_cdt, BLUE, D0, a_), (ys_cal, GREEN, 2, 4), (ys_erf, GRAY, 2, 4)]:
    r = np.log(np.clip((ys - lo), 1e-12, None)/np.clip((hi - ys), 1e-12, None)); mk = np.abs(r) < 7.9
    axR.plot(s[mk], r[mk], color=col, width=2.5, dash="6,4" if col == GRAY else None)
axR.text(-2.9, 6.8, "straight: AS toy, CDT, Calcagni", size=12, color=DARK); axR.text(2.9, -7.2, "grey: bends (D5)", size=12, color=GRAY, anchor="end")
yl = 565
for i, (col, lab) in enumerate([(PURPLE, "asymptotic safety toy (D1)"), (BLUE, "CDT fit, Ambjørn–Jurkiewicz–Loll 2005 (D2)"), (GREEN, "Calcagni's multiscale measure (D3)"), (GRAY, "erf-shaped two-plateau curve (D5)")]):
    xx = 90 + (i % 2)*540; yy = yl + (i//2)*24
    f.line(xx, yy - 4, xx + 30, yy - 4, color=col, width=3); f.text(xx + 40, yy, lab, size=13, color=col)
f.text(90, 632, "Ours: each is (A + B sᵖ)/(1 + C sᵖ); its 'position between the two walls' composes by relativistic velocity addition as the scale moves (D4).", size=12, color=GRAY)
f.save(P("a6-dimension-flow-is-rapidity.svg"))

# ---------- a7: deformed covariance ----------
f = Fig(1180, 600); f.title("Deformed general covariance (loop quantum gravity): a lift whose wall is a signature change, with G setting where", size=17)
ax = Axes(f, 100, 80, 620, 380, (0, 1.05), (-1.15, 1.15))
ax.frame(xticks=[(0, "0"), (0.5, "ρ_c / 2"), (1, "ρ_c")], yticks=[-1, 0, 1], xlabel="energy density ρ (ρ_c ≈ 0.41 Planck density)", ylabel="")
r = np.linspace(0, 1, 200)
ax.hline(0, color=DARK, width=1)
ax.plot(r[r <= 0.5], 1 - 2*r[r <= 0.5], color=BLUE, width=3); ax.plot(r[r >= 0.5], 1 - 2*r[r >= 0.5], color=RED, width=3)
ax.plot(r, 4*r*(1 - r), color=GREEN, width=2, dash="6,4")
ax.dot(0.5, 0, r=6, color=DARK); ax.dot(1, 0, r=6, color=GREEN)
ax.text(0.03, -0.3, "β = 1 − 2ρ/ρ_c > 0: Lorentzian", size=13, color=BLUE)
ax.text(0.03, -0.75, "β < 0: Euclidean (four space directions, no time)", size=13, color=RED)
ax.text(0.52, 0.08, "wall: β = 0, the Carroll point", size=12, color=DARK)
ax.text(0.56, 1.06, "H² ∝ ρ(1 − ρ/ρ_c), bounce at ρ_c", size=12, color=GREEN)
cx, cy, R = 900, 230, 110
f.circle(cx, cy, R, color=GRAY); f.line(cx - R - 20, cy, cx + R + 20, cy, color="#bbb"); f.line(cx, cy - R - 20, cx, cy + R + 20, color="#bbb")
for ang, col in [(0.0, BLUE), (np.pi/2, DARK), (2.4, RED)]:
    px, py = cx + R*np.cos(ang), cy - R*np.sin(ang); f.line(cx, cy, px, py, color=col, width=2); f.circle(px, py, 5, color=col, fill=col)
f.text(cx + R + 12, cy - 6, "δK = 0", size=12, color=BLUE); f.text(cx + 8, cy - R - 12, "quarter turn: β = 0", size=12, color=DARK)
f.text(cx - R*0.74 - 8, cy - R*0.68 - 10, "Euclidean", size=12, color=RED, anchor="end")
f.text(cx, cy + R + 40, "β = cos(2δK): the real part of a turn", size=13, anchor="middle"); f.text(cx, cy + R + 60, "the signature flips at a quarter turn", size=13, anchor="middle")
f.lines(100, 548, ["Known: Bojowald–Paily 2012; Bojowald–Mielczarek 2015 (β, signature change); Ashtekar–Pawlowski–Singh (ρ_c, bounce). Physical status debated.",
                   "Reading (ours): β is the lift 1 + kρ with k = −2/ρ_c ∝ G²ħ. The slot it deforms is the signature (c's direction); G only places the wall. Checks B1–B4."], size=12, color=GRAY)
f.save(P("a7-deformed-covariance.svg"))

# ---------- a8: verdict board ----------
f = Fig(1180, 640); f.title("Is G a deformation of a spacetime made of anomalous dimensions? Claim by claim", size=17)
rows = [("Units form a lattice; natural units collapse it", "known", GREEN),
        ("A quantity's units can change with scale (anomalous dimension)", "known", GREEN),
        ("At the fixed point η_N = −2: G is dimensionless and spacetime looks 2D", "known, if the fixed point exists", GREEN),
        ("G (with Λ) labels the trajectory that flows from that spacetime to 4D GR", "known RG reading", GREEN),
        ("The label is one pure number, GħΛ/c³, the same as in g-dials", "known in substance; link ours", ORANGE),
        ("That flow is a deformation in the Inönü–Wigner / rigidity sense", "a stretch: no algebra, no H² test", RED),
        ("In 1/g the one-loop flow is exactly mult-add; the wall is the fixed point", "ours, exact for that shape", PURPLE),
        ("Dimension flows (AS, CDT, Calcagni) are exact two-wall rapidity curves", "ours, exact for the published forms", PURPLE),
        ("Deformed covariance bends a G-dependent structure", "known; but its dial is the signature", ORANGE),
        ("Either supplies the area law's 1/4", "no", RED)]
for i, (claim, st, col) in enumerate(rows):
    y = 70 + i*54
    if i % 2 == 0: f.rect(30, y, 1120, 54, color="none", fill="#f5f6f8")
    f.text(45, y + 33, claim, size=14)
    f.rect(800, y + 12, 340, 30, color=col, fill=col, opacity=0.18, rx=6); f.text(970, y + 32, st, size=13, anchor="middle", color=DARK)
f.text(45, 625, "Green: established in the literature. Orange: established pieces with a caveat or our link. Purple: ours. Red: not supported.", size=12, color=GRAY)
f.save(P("a8-verdict-board.svg"))
print("ok")
