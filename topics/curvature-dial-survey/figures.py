"""Pictures for the curvature-dial survey. Run: python3 figures.py -> figures/*.png
Colour key: blue = our reading, orange = already known, grey = endpoints / background."""
import json, math, cmath, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

BLUE, ORANGE, GREY, RED = "#1f5fbf", "#e07b00", "#8a8a8a", "#c0392b"
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
                     "axes.spines.top": False, "axes.spines.right": False})
os.makedirs("figures", exist_ok=True)
def save(fig, name):
    fig.tight_layout(); fig.savefig(f"figures/{name}.png", dpi=130); plt.close(fig); print("figures/" + name + ".png")
def law(k, d, a, b): return (a + b + k*a*b)/(1 - d*a*b)

# ---------- hero ----------
fig, (A, B) = plt.subplots(1, 2, figsize=(14, 6.2))
dn = np.linspace(-2, -0.02, 300)
A.plot(dn, 1/np.sqrt(-dn), color=ORANGE, lw=3); A.plot(dn, -1/np.sqrt(-dn), color=ORANGE, lw=3)
A.text(-1.95, 2.4, "two walls ±1/√(−d)\n(light speed at d = −1/c²)", color=ORANGE, fontsize=13)
dp = np.linspace(0.06, 2, 300)
A.plot(dp, np.pi/(2*np.sqrt(dp)), color=BLUE, lw=3)
A.text(0.25, 4.4, "no wall: a turn,\nquarter turn π/(2√d)", color=BLUE, fontsize=13)
A.axvline(0, color=GREY, lw=1.5, ls="--")
for xx, lab in [(-1, "Einstein\nd = −1/c²"), (0, "Galileo\nd = 0"), (1, "tangent\nd = +1")]:
    A.plot([xx], [0], "o", color="k", ms=11); A.annotate(lab, (xx, 0), (xx - .35, -1.9), fontsize=13)
A.set_xlim(-2.05, 2.05); A.set_ylim(-3, 6); A.set_xlabel("curvature dial d   (k = 0)"); A.set_ylabel("where the law runs out")
A.set_title("A. Turning d at k = 0: three stops")
kk = np.linspace(-3, 3, 300)
B.plot(kk, kk**2/4, color=GREY, lw=2.5); B.text(-2.95, 1.55, "parabola Δ = 0", color=GREY, fontsize=12)
B.axvline(0, color=BLUE, lw=3, alpha=.6); B.axvline(1, color=BLUE, lw=2, ls="--", alpha=.6)
B.plot(kk, kk - 1, color=ORANGE, lw=2); B.text(-2.9, -3.6, "Hirzebruch χ_y line\nd = k − 1", color=ORANGE, fontsize=12)
pts = [(0, -1, "Einstein = signature (0, −1)", (-2.6, -.75)), (0, 0, "Galileo (0, 0)", (-1.75, .15)), (0, 1, "tangent (0, 1)", (-1.75, 1.15)),
       (1, 0, "circle operation = Todd (1, 0)", (1.1, -.55)), (2, 1, "Euler (2, 1)", (2.1, .6))]
for px, py, lab, pos in pts:
    B.plot([px], [py], "o", color="k", ms=9); B.annotate(lab, (px, py), pos, fontsize=11)
for N in [3, 4, 5, 6, 8]:
    B.plot([1], [1/(4*math.cos(math.pi/N)**2)], "s", color=BLUE, ms=7)
B.annotate("finite clocks on k = 1\n(section N)", (1, .5), (1.4, 2.3), fontsize=11, color=BLUE, arrowprops=dict(arrowstyle="-", color=BLUE))
B.text(-.9, 2.6, "turning (Δ < 0)", fontsize=12); B.text(-2.9, -1.4, "two walls\n(Δ > 0)", fontsize=12)
B.set_xlim(-3, 3); B.set_ylim(-4.2, 3.2); B.set_xlabel("homothety k"); B.set_ylabel("curvature dial d")
B.set_title("B. The survey lines k = 0 and k = 1")
save(fig, "fig0-hero-curvature-dial")

# ---------- sorting board ----------
L = json.load(open("laws.json"))
order = ["slides", "plausible", "partly", "disguise", "self", "no"]
cols = {"slides": BLUE, "plausible": BLUE, "partly": BLUE, "disguise": GREY, "self": "k", "no": GREY}
fig, ax = plt.subplots(figsize=(15, 11)); ax.axis("off")
xpos = {"slides": 0, "plausible": .25, "partly": .25, "disguise": .25, "self": .25, "no": .5}
ystart = {"slides": .95, "no": .95}
_y = .95
for v in ["plausible", "partly", "disguise", "self"]:
    ystart[v] = _y; _y -= .045 + .046*sum(r["verdict"] == v for r in L["rows"]) + .05
for v in order:
    rows = [r for r in L["rows"] if r["verdict"] == v]
    x0, y0 = xpos[v], ystart[v]
    ax.text(x0, y0, f"{v} ({len(rows)})", fontsize=15, weight="bold", color=cols[v], transform=ax.transAxes)
    for i, r in enumerate(rows):
        xi = x0 + (0.25 if (v == "no" and i >= 19) else 0); yi = y0 - .045 - .046*(i % 19 if v == "no" else i)
        name = f"{r['n']}. {r['law']}"; name = name if len(name) < 34 else name[:32] + "…"
        tag = f" [{r['candidate']}]" if r["candidate"] else ""
        ax.text(xi, yi, name + tag, fontsize=10.5, color=cols[v], transform=ax.transAxes)
fig.suptitle("All 65 law groups, sorted by what turning d does (k held fixed)", fontsize=17, y=.995)
save(fig, "fig-sorting-board")

# ---------- C1 ----------
fig, (A, B) = plt.subplots(1, 2, figsize=(14, 5.6))
aa = np.linspace(-3, 3, 2000)
for dd, col, lab in [(-1, ORANGE, "d = −1 (Einstein, c = 1)"), (-.25, ORANGE, "d = −1/4"), (0, GREY, "d = 0 (Galileo)"),
                     (.25, BLUE, "d = 1/4"), (1, BLUE, "d = 1 (tangent)")]:
    yv = law(0, dd, aa, .5); yv[np.abs(1 - dd*aa*.5) < 0.02] = np.nan
    A.plot(aa, yv, color=col, lw=2.5 if dd in (-1, 0, 1) else 1.4, label=lab)
A.axhline(1, color=ORANGE, ls=":", lw=1); A.axhline(-1, color=ORANGE, ls=":", lw=1)
A.set_ylim(-4, 4); A.set_xlabel("a"); A.set_ylabel("a ⊕ 0.5"); A.legend(loc="lower right", fontsize=10.5)
A.set_title("A. One law, three kinds (k = 0)")
dd = np.linspace(-1.5, .249, 400); disc = np.sqrt(1 - 4*dd)
r1 = np.where(dd != 0, (-1 + disc)/(2*np.where(dd == 0, 1, dd)), -1); r2 = (-1 - disc)/(2*dd)
B.plot(dd, -1/(1 + disc)*2, color=ORANGE, lw=3, label="near wall (−1/k at d = 0)")
farw = (-1 - disc)/(2*dd); farw[np.abs(dd) < 0.02] = np.nan
B.plot(dd, farw, color=BLUE, lw=3, label="far wall (from infinity at d = 0)")
B.axvline(.25, color=GREY, ls="--"); B.text(.27, -9, "parabola\nd = 1/4:\nwalls merge", color=GREY, fontsize=12)
B.axvline(0, color=GREY, lw=.8); B.set_ylim(-12, 12); B.set_xlim(-1.5, .6)
B.set_xlabel("d   (k = 1)"); B.set_ylabel("wall position"); B.legend(fontsize=11, loc="upper left")
B.set_title("B. At k = 1 the kind changes at d = k²/4")
save(fig, "fig1-galileo-einstein-tangent")

# ---------- C2 ----------
def U(n, k, d):
    u0, u1 = 0*d, 1 + 0*d
    for _ in range(n - 1): u0, u1 = u1, k*u1 - d*u0
    return u1
fig, ax = plt.subplots(figsize=(11, 6))
dd = np.linspace(-.6, 1.3, 800)
for n, col in zip(range(2, 8), plt.cm.Blues(np.linspace(.45, 1, 6))):
    ax.plot(dd, U(n, 1, dd), color=col, lw=2.4, label=f"U_{n}: weight of a {n}-element corner")
for N in [3, 4, 5, 6]:
    z = 1/(4*math.cos(math.pi/N)**2); ax.plot([z], [0], "o", color=RED, ms=9); ax.annotate(f"U_{N} = 0\nd = {z:.3f}", (z, 0), (z - .02 + (.12 if N == 4 else 0) - (.25 if N == 6 else 0), {3: -1.1, 4: -1.7, 5: -1.1, 6: -1.7}[N]), fontsize=11, color=RED)
ax.axvline(.25, color=GREY, ls="--"); ax.text(.0, 2.6, "d = 0: every weight is k^(n−1) = 1\n(the cube's corners)", fontsize=12, color=GREY)
ax.axhline(0, color=GREY, lw=.8); ax.set_ylim(-2.2, 3.2); ax.set_xlabel("curvature dial d   (k = 1)"); ax.set_ylabel("corner weight U_n(k, d)")
ax.legend(fontsize=10, loc="upper right"); ax.set_title("C2. Corner weights become Lucas numbers")
save(fig, "fig2-lucas-corner-weights")

# ---------- C3 ----------
fig, ax = plt.subplots(figsize=(11, 5.6))
u0, v0 = .5 + .2j, -.1 + .6j; uv = np.conj(u0)*v0
dd = np.linspace(-1.2, 1.2, 400)
om = 2*np.arctan2(dd*uv.imag, 1 - dd*uv.real)
ax.plot(dd, np.degrees(om), color=BLUE, lw=3)
ax.axhline(0, color=GREY, lw=.8); ax.axvline(0, color=GREY, lw=.8)
ax.text(-1.15, 25, "d < 0: Thomas–Wigner,\ngeodetic, Lense–Thirring", color=ORANGE, fontsize=13)
ax.text(.35, -45, "d > 0: Foucault, Berry,\nLarmor, quaternions", color=ORANGE, fontsize=13)
ax.annotate("d = 0: Galilean boosts,\nno twist at all", (0, 0), (-.9, -50), fontsize=13, color=GREY, arrowprops=dict(arrowstyle="->", color=GREY))
ax.set_xlabel("curvature dial d"); ax.set_ylabel("twist ω (degrees)")
ax.set_title("C3. tan(ω/2) = d(u × v)/(1 − d u·v): the twist turns on with d")
save(fig, "fig3-holonomy-turns-on")

# ---------- C4 ----------
fig, (A, B) = plt.subplots(1, 2, figsize=(14, 6))
th = np.linspace(0, 2*np.pi, 400)
for dd, col in [(1, BLUE), (.3, BLUE), (-.3, ORANGE), (-1, ORANGE)]:
    if dd > 0: A.plot(np.cos(th), np.sin(th)/np.sqrt(dd), color=col, lw=2.4 if dd == 1 else 1.4)
    else:
        s = np.linspace(-1.6, 1.6, 200)
        for sg in (1, -1): A.plot(sg*np.cosh(s), np.sinh(s)/np.sqrt(-dd), color=col, lw=2.4 if dd == -1 else 1.4)
A.axvline(1, color=GREY, lw=2.5); A.axvline(-1, color=GREY, lw=2.5)
A.text(-2.95, -3.85, "grey, d = 0: two parallel lines\n(dual numbers, Galilean shear)", backgroundcolor="white", color=GREY, fontsize=11)
A.text(-.45, .15, "d > 0\ncircle", color=BLUE, fontsize=12); A.text(2.05, 1.0, "d < 0\nhyperbola", color=ORANGE, fontsize=12)
A.set_xlim(-3, 3); A.set_ylim(-4, 4); A.set_aspect("equal"); A.set_xlabel("x"); A.set_ylabel("y")
A.set_title("A. Unit 'circle' x² + d y² = 1")
for E, col, lab in [(-.2, BLUE, "bound E < 0 (d > 0 side)"), (0, GREY, "parabola E = 0 (flat)"), (.3, ORANGE, "escape E > 0 (d < 0 side)")]:
    v0_ = math.sqrt(2*(E + 1))
    s = solve_ivp(lambda tt, z: [z[2], z[3], -z[0]/np.hypot(z[0], z[1])**3, -z[1]/np.hypot(z[0], z[1])**3],
                  [0, 60], [1, 0, 0, v0_], rtol=1e-10, atol=1e-12, max_step=.05)
    B.plot(s.y[2], s.y[3], color=col, lw=2.5, label=lab)
B.plot(math.sqrt(.6)*np.cos(th), math.sqrt(.6)*np.sin(th), color=ORANGE, ls=":", lw=2)
B.text(-.75, -.82, "wall: |v| = v∞", color=ORANGE, fontsize=12)
B.set_aspect("equal"); B.set_xlabel("v_x"); B.set_ylabel("v_y"); B.legend(fontsize=10, loc="upper left", bbox_to_anchor=(1.0, 1.0))
B.set_title("B. Kepler hodographs: curvature −2E")
save(fig, "fig4-number-planes-kepler")

# ---------- C5 ----------
p, a1 = 13, 1
fig, ax = plt.subplots(figsize=(11, 5.4))
for a2 in range(p):
    n = sum(1 for X in range(p) for Y in range(p) if (Y*Y + a1*X*Y - X**3 - a2*X*X) % p == 0)
    D = (a1*a1 + 4*a2) % p
    col = GREY if D == 0 else (ORANGE if pow(D, (p - 1)//2, p) == 1 else BLUE)
    ax.bar(a2, n, color=col)
ax.axhline(p, color="k", lw=.8)
ax.text(.2, 14.6, "blue: non-split node, p + 1 (turning)     orange: split node, p − 1 (two walls)     grey: cusp, p",
        fontsize=11.5)
ax.set_ylim(10, 15.5); ax.set_xlabel("−d = a₂ (mod 13), with k = −a₁ = −1"); ax.set_ylabel("smooth points over F₁₃")
ax.set_title("C5. Pinched cubics y² + xy = x³ + a₂x²: the d dial picks the kind")
save(fig, "fig5-singular-cubics")

# ---------- C6 ----------
fig, ax = plt.subplots(figsize=(11, 5.6))
dd = np.linspace(-2, .2499, 500); disc = np.sqrt(1 - 4*dd)
with np.errstate(divide="ignore", invalid="ignore"):
    near = 2/(1 + disc); far = (1 + disc)/(2*dd)
ax.plot(dd, near, color=ORANGE, lw=3, label="anchor that is 'certainty' (1 at d = 0)")
ax.plot(dd[dd < -0.01], far[dd < -0.01], color=BLUE, lw=3, label="second anchor (at infinity when d = 0)")
ax.plot(dd[dd > 0.01], far[dd > 0.01], color=BLUE, lw=3)
ax.axvline(.25, color=GREY, ls="--"); ax.text(.27, -6, "one anchor\n(d = 1/4)\nthen none", color=GREY, fontsize=12)
ax.plot([0], [1], "o", color="k", ms=9); ax.annotate("union rule x + y − xy", (0, 1), (-1.6, 4), fontsize=12, arrowprops=dict(arrowstyle="->"))
ax.set_ylim(-8, 8); ax.set_xlim(-2, .7); ax.set_xlabel("curvature dial d   (k = −1)"); ax.set_ylabel("anchor (fixed point of every move)")
ax.legend(fontsize=11, loc="lower left"); ax.set_title("C6. Probability needs two anchors: it survives only while Δ > 0")
save(fig, "fig6-anchors")

# ---------- C7 ----------
fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
for ax, dd, title in zip(axs, [-1, -.05, .5], ["d = −1: a light cone", "d → 0⁻: the cone opens flat", "d > 0: no cone (Laplace)"]):
    tt = np.linspace(-1, 1, 10)
    if dd < 0:
        c = 1/math.sqrt(-dd); xx = np.linspace(-3, 3, 10)
        ax.fill_between(xx, np.abs(xx)/c, 3, color=ORANGE, alpha=.25); ax.fill_between(xx, -3, -np.abs(xx)/c, color=ORANGE, alpha=.25)
        ax.plot(xx, xx/c, color=ORANGE, lw=2); ax.plot(xx, -xx/c, color=ORANGE, lw=2)
    else:
        th = np.linspace(0, 2*np.pi, 200)
        for r in [.8, 1.6, 2.4]: ax.plot(r*np.cos(th), r*np.sin(th)/math.sqrt(dd), color=BLUE, lw=1.6)
    ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_xlabel("x"); ax.set_ylabel("t"); ax.set_title(title, fontsize=13)
fig.suptitle("C7. The operator ∂ₓ² + d ∂ₜ²: waves need d < 0", fontsize=15)
save(fig, "fig7-light-cone-needs-walls")

# ---------- C8 ----------
fig, ax = plt.subplots(figsize=(11, 5.6))
for seg in [np.linspace(-3, .245, 300), np.linspace(.255, 3, 300)]:
    ax.plot(seg, 1/(1 - 4*seg), color=BLUE, lw=3)
ax.axvline(.25, color=GREY, ls="--"); ax.axhline(0, color=GREY, lw=.8)
for dd, s, lab in [(0, 1, "s = 1: the circle operation"), (1/8, 2, "s = 2"), (.5, -1, "s = −1: a rewind\n(same residue as log(1 + x))")]:
    ax.plot([dd], [s], "o", color="k", ms=8); ax.annotate(lab, (dd, s), (dd + .25, s + .6), fontsize=12)
ax.text(-2.9, 2.5, "rung height s = k²/Δ\nresidue Δ/(3k²) = 1/(3s)", color=BLUE, fontsize=13)
ax.set_ylim(-4, 5); ax.set_xlabel("curvature dial d   (k = 1)"); ax.set_ylabel("height s")
ax.set_title("C8. At fixed k, d reads as a height (formal match only)")
save(fig, "fig8-height-reading")

# ---------- C9 ----------
def yval(kk, dd):
    if dd == 0: return 0j
    yy = [1 - (kk + sg*cmath.sqrt(kk*kk - 4*dd))/(2*dd)*kk for sg in (1, -1)]
    return min(yy, key=lambda z: abs(z))
fig, (A, B) = plt.subplots(1, 2, figsize=(14, 6))
dd = np.r_[np.linspace(-6, -.01, 300), np.linspace(.001, .25, 200)]
A.plot(dd, [yval(1, x).real for x in dd], color=BLUE, lw=3)
for x0, y0, lab in [(0, 0, "Todd (y = 0)"), (.25, -1, "Euler (y = −1)")]:
    A.plot([x0], [y0], "o", color=ORANGE, ms=10); A.annotate(lab, (x0, y0), (x0 - 2.6, y0 + .25), fontsize=13, color=ORANGE)
A.text(-5.8, .45, "→ signature (y = 1) as d → −∞", color=ORANGE, fontsize=12)
A.set_xlabel("d   (k = 1, two-wall side)"); A.set_ylabel("χ_y parameter y"); A.set_title("A. Real y up to the parabola")
th = np.linspace(0, 2*np.pi, 300); B.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
dds = np.linspace(.2501, 8, 300); yy = np.array([yval(1, x) for x in dds])
B.plot(yy.real, np.abs(yy.imag), color=BLUE, lw=4, label="d from 1/4 up")
B.text(-.75, -1.35, "(y and its conjugate are\nthe same genus; upper half shown)", fontsize=11, color=GREY)
for N in [3, 4, 5, 6, 8]:
    z = -np.exp(-2j*np.pi/N); B.plot(z.real, z.imag, "s", color=ORANGE, ms=9)
    B.annotate(f"−ζ_{N}", (z.real, z.imag), (z.real*1.18 - .1, z.imag*1.18), fontsize=12, color=ORANGE)
B.plot([-1], [0], "o", color="k", ms=8); B.annotate("Euler, d = 1/4", (-1, 0), (-.8, -.25), fontsize=11)
B.set_aspect("equal"); B.set_xlim(-1.5, 1.5); B.set_ylim(-1.5, 1.5); B.set_title("B. Turning side: |y| = 1")
B.set_xlabel("Re y"); B.set_ylabel("Im y")
save(fig, "fig9-genus-along-d")

# ---------- C10 ----------
fig, ax = plt.subplots(figsize=(11, 5.4))
dp = np.linspace(.05, 3, 300); ax.plot(dp, np.pi/np.sqrt(dp), color=BLUE, lw=3, label="ruler length π/√d (one full turn of the slope)")
ax.axvspan(-1.5, 0, color=ORANGE, alpha=.12); ax.text(-1.4, 10, "d ≤ 0: the ruler never closes;\nexp(θu) never returns to −1", color=ORANGE, fontsize=13)
ax.plot([1], [np.pi], "o", color="k", ms=9); ax.annotate("d = 1: π, and e^{iπ} = −1", (1, np.pi), (1.3, 6), fontsize=13, arrowprops=dict(arrowstyle="->"))
ax.set_xlim(-1.5, 3); ax.set_ylim(0, 14); ax.set_xlabel("curvature dial d   (k = 0)"); ax.set_ylabel("length of the ruler ∫ dx/(1 + d x²)")
ax.legend(fontsize=11); ax.set_title("C10. π and Euler's identity live only on the turning side")
save(fig, "fig10-pi-and-i")

# ---------- C11 ----------
fig, ax = plt.subplots(figsize=(11, 5.6))
dd = np.linspace(-1, 1.2, 600); disc = np.sqrt((1 - 4*dd).astype(complex))
r1 = (1 + disc)/2; r2 = (1 - disc)/2
ax.plot(dd, r1.real, color=BLUE, lw=3, label="sample point α (real part)"); ax.plot(dd, r2.real, color=BLUE, lw=3, ls="--", label="sample point β (real part)")
ax.plot(dd, r1.imag, color=ORANGE, lw=2.4, label="imaginary part (complex step)"); ax.plot(dd, r2.imag, color=ORANGE, lw=2.4)
ax.axvline(.25, color=GREY, ls="--")
ax.text(-.95, 1.3, "Δ > 0: two real samples\n= finite difference", fontsize=12); ax.text(.27, 1.3, "Δ = 0: one point + slope\n= jet / dual numbers", fontsize=12, color=GREY)
ax.text(.55, -.9, "Δ < 0: complex pair\n= complex step", fontsize=12, color=ORANGE)
ax.set_xlabel("curvature dial d   (k = 1; d = 0 is the jet-fields page's 'two points k apart')"); ax.set_ylabel("where f is sampled (× b)")
ax.legend(fontsize=10.5, loc="lower left"); ax.set_title("C11. f(a + bu) in R[u]/(u² − ku + d)")
save(fig, "fig11-jets-finite-difference-complex-step")

# ---------- C12 ----------
fig, ax = plt.subplots(figsize=(11, 5.6))
xx = np.linspace(-6, 6, 800)
for dd, col, lab in [(-1, ORANGE, "d = −1: Yukawa, reach 1"), (-.2, ORANGE, "d = −0.2: Yukawa, reach 2.2"), (0, GREY, "d = 0: Coulomb (1D: −|x|/2)"),
                     (.2, BLUE, "d = 0.2: oscillating"), (1, BLUE, "d = 1: oscillating")]:
    if dd < 0: m = math.sqrt(-dd); g = np.exp(-m*np.abs(xx))/(2*m); g -= g.max()
    elif dd == 0: g = -np.abs(xx)/2
    else: g = -np.sin(math.sqrt(dd)*np.abs(xx))/(2*math.sqrt(dd))
    ax.plot(xx, g, color=col, lw=2.6 if dd in (-1, 0, 1) else 1.4, label=lab)
ax.set_xlabel("distance x"); ax.set_ylabel("field G(x) − G(0)"); ax.legend(fontsize=10.5, loc="lower center")
ax.set_title("C12. −w'' − d w = δ: carrier mass² = −d at k = 0")
save(fig, "fig12-force-reach")

# ---------- C13 ----------
fig, ax = plt.subplots(figsize=(11, 5.4))
KK = np.linspace(-1.3, 1.3, 10); ax.plot(KK, -2*KK, color=BLUE, lw=3, label="averaged law: k = −2d")
for K, lab, col in [(-1, "barrier, SU(1,1)\nB₋₁(u): k = 2", ORANGE), (0, "free flight\nk = 0, add", GREY), (1, "beam splitter, SU(2)\nB₊₁(u): k = −2", ORANGE)]:
    ax.plot([K], [-2*K], "o", color=col, ms=12); ax.annotate(lab, (K, -2*K), (K + .08, -2*K + .35), fontsize=12)
ax.axhline(0, color=GREY, lw=.8); ax.axvline(0, color=GREY, lw=.8)
ax.set_xlabel("curvature dial d of the amplitudes"); ax.set_ylabel("homothety k of the averaged probabilities")
ax.legend(fontsize=11, loc="upper right"); ax.set_title("C13. One twist B_d(u): barrier, free, beam splitter")
save(fig, "fig13-tunnelling")

# ---------- C14 ----------
fig, ax = plt.subplots(figsize=(10, 5.8))
s0, P0 = 2.0, .8; D0 = s0**2 - 4*P0
ts = np.linspace(-.6, .3, 300)
for tt in ts:
    disc = cmath.sqrt(s0**2 - 4*(P0 - 2*tt)); z1, z2 = (s0 + disc)/2, (s0 - disc)/2
    col = plt.cm.coolwarm((tt - ts[0])/(ts[-1] - ts[0]))
    ax.plot([z1.real, z2.real], [z1.imag, z2.imag], "o", color=col, ms=3)
ax.annotate(f"collision at t = −Δ/8 = {-D0/8:.2f}", (1, 0), (1.1, .5), fontsize=12, arrowprops=dict(arrowstyle="->"))
ax.text(.05, -.75, "blue: rewind (t < 0), the pair turns complex\nred: forward, the pair spreads on the line", fontsize=12)
ax.set_xlabel("Re z"); ax.set_ylabel("Im z"); ax.set_aspect("equal")
ax.set_title("C14. Heat flow on z² − sz + p: sum kept, product slides")
save(fig, "fig14-riemann-pair")

# ---------- C15 ----------
fig, ax = plt.subplots(figsize=(10, 5.8))
kap = np.linspace(-1.6, 1.6, 300)
ax.plot(kap, kap**2, color=ORANGE, lw=3, label="the d road: conic laws, ε = κ² (κ = d)")
ax.plot([-1, 1], [1, 1], color=ORANGE, lw=3, ls="--", label="the modulus road: ε = 1, κ = 1 − 2m")
ax.plot([-1, 1], [1, 1], "o", color=BLUE, ms=12)
ax.annotate("m = 1: tanh\n(separatrix)", (-1, 1), (-1.55, 1.6), fontsize=12, color=BLUE)
ax.annotate("m = 0: tan\n(small swings)", (1, 1), (.95, 1.6), fontsize=12, color=BLUE)
ax.annotate("d = 0: Galileo", (0, 0), (.2, .35), fontsize=12, arrowprops=dict(arrowstyle="->"))
ax.annotate("m = 1/2: square torus\n(90° swing)", (0, 1), (-.55, 1.25), fontsize=12)
ax.set_xlabel("κ (coefficient of 2t²)"); ax.set_ylabel("ε (coefficient of t⁴)"); ax.set_ylim(-.2, 2.2)
ax.legend(fontsize=10.5, loc="lower right"); ax.set_title("C15. log' = 1/√(1 + 2κt² + εt⁴): two roads, two meetings")
save(fig, "fig15-elliptic-modulus-road")

# ---------- N1 ----------
fig, (A, B) = plt.subplots(1, 2, figsize=(14, 6))
th = np.linspace(0, 2*np.pi, 300); A.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
for N, col in zip([3, 4, 5, 6], plt.cm.Blues(np.linspace(.5, 1, 4))):
    ang = 2*np.pi*np.arange(N + 1)/N
    A.plot(np.cos(ang)*(1 - .05*N/6), np.sin(ang)*(1 - .05*N/6), "-o", color=col, lw=2, ms=6,
           label=f"N = {N}: d/k² = {1/(4*math.cos(math.pi/N)**2):.3f}")
A.plot([1], [0], "*", color="k", ms=16); A.annotate("identity 0", (1, 0), (1.05, .15), fontsize=12)
A.set_aspect("equal"); A.axis("off"); A.legend(fontsize=11, loc="lower left", bbox_to_anchor=(-.1, -.12))
A.set_title("A. Repeating 'add infinity' closes after N steps")
for N in range(3, 25):
    for j in range(1, N//2 + 1):
        if math.gcd(j, N) != 1: continue
        z = 1/(4*math.cos(math.pi*j/N)**2)
        if z < 3: B.plot([z], [N], "o", color=BLUE if j == 1 else ORANGE, ms=6)
B.axvline(.25, color=GREY, ls="--"); B.text(.32, 2.8, "parabola d = k²/4", color=GREY, fontsize=12)
B.set_xlabel("d / k²   (turning side)"); B.set_ylabel("clock size N"); B.set_xlim(0, 2.1)
B.text(.9, 5, "blue: j = 1 (d/k² = 1/B_N, Beraha)\norange: other j", fontsize=12)
B.set_title("B. Where the clocks sit")
save(fig, "fig16-finite-clocks")
