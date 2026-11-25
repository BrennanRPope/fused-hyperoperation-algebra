"""Figures for topics/why-i/e-i-pi/chain-resolves. Run: python3 figures.py (numpy only)."""
import numpy as np, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from svgkit import *
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True); P = lambda n: os.path.join(OUT, n)
def Scl(x, y, t): return ((x*x + y*y)*np.cos(t) - 2*x*y)/(2*np.sin(t))
def xcl(x0, x1, T, t): return (x0*np.sin(T - t) + x1*np.sin(t))/np.sin(T)
x0, x1 = 0.0, 1.0
def Smid(xm, T): return Scl(x0, xm, T/2) + Scl(xm, x1, T/2)

# r1 / r2: middle-point integrand, short and long time
for name, T, title in [("r1-short-time-both-dials-resolve.svg", 0.6*np.pi, "Short trip (T = 0.6 half-periods): both dials pick out the classical path"),
                       ("r2-long-time-only-imaginary-resolves.svg", 1.5*np.pi, "Long trip (T = 1.5 half-periods): the real dial blows up; the imaginary dial still picks the classical path")]:
    f = Fig(1180, 530); f.title(title, size=17)
    xs = xcl(x0, x1, T, T/2); hb = 0.02
    g = np.linspace(xs - 2.5, xs + 2.5, 3000); dS = Smid(g, T) - Smid(xs, T)
    ax = Axes(f, 70, 70, 480, 360, (xs - 2.5, xs + 2.5), (0, 1.2) if T < np.pi else (0, 6))
    ax.frame(xticks=[(xs, "classical"), (xs - 2, "−2"), (xs + 2, "+2")], yticks=[0, 1] if T < np.pi else [0, 1, 2, 4, 6],
             xlabel="where the path is at the halfway time (offset from classical)", ylabel="weight e^(−S/ħ) (relative)")
    wr = np.minimum(np.exp(-dS/hb), 5.95)
    ax.plot(g, wr, color=BLUE, width=2.5)
    ax.vline(xs, color=GREEN, width=1.5, dash="4,3")
    if T < np.pi:
        ax.text(xs + 0.25, 0.9, "real dial: a sharp peak", size=13, color=BLUE); ax.text(xs + 0.25, 0.8, "at the classical path (Laplace)", size=13, color=BLUE)
    else:
        ax.text(xs - 2.35, 5.5, "real dial: weight GROWS away from", size=13, color=BLUE); ax.text(xs - 2.35, 5.15, "the classical path → the sum diverges", size=13, color=BLUE)
        ax.text(xs + 0.15, 0.4, "a dip, not a peak", size=12, color=BLUE); ax.text(xs + 0.3, 4.6, "(cut off at 6; it keeps growing)", size=11, color=GRAY)
    bx = Axes(f, 650, 70, 480, 360, (xs - 2.5, xs + 2.5), (-1.3, 1.3))
    bx.frame(xticks=[(xs, "classical"), (xs - 2, "−2"), (xs + 2, "+2")], yticks=[-1, 0, 1], xlabel="where the path is at the halfway time", ylabel="Re e^(iS/ħ)")
    hb2 = 0.02
    bx.plot(g, np.cos(dS/hb2), color=ORANGE, width=1.2)
    bx.vline(xs, color=GREEN, width=1.5, dash="4,3")
    bx.text(xs - 2.4, 1.15, "imaginary dial: fast wiggles cancel everywhere", size=12, color=ORANGE)
    bx.text(xs - 2.4, -1.22, "except the slow stretch at the classical path (stationary phase)", size=12, color=ORANGE)
    f.text(70, 505, "Harmonic oscillator from x = 0 to x = 1, ħ = 0.02. Classical point marked green. Checks R1–R5 run the same test at T = 0.6 and 1.4 half-periods.", size=13, color=GRAY)
    f.save(P(name))

# r3: Morse index and allowed dials
f = Fig(1180, 520); f.title("Why only the imaginary dial works for every trip length", size=18)
def hessian(T, N=200):
    e = T/N; n = N - 1
    return (np.diag(np.full(n, 2/e)) - np.diag(np.full(n - 1, 1/e), 1) - np.diag(np.full(n - 1, 1/e), -1)) - np.eye(n)*e
Ts = np.linspace(0.2, 3.4, 120)*np.pi
lows = np.array([np.sort(np.linalg.eigvalsh(hessian(t)))[:4] for t in Ts])*(200/Ts[:, None])   # = (n pi/T)^2 - 1 in the continuum
ax = Axes(f, 80, 70, 520, 360, (0.2, 3.4), (-1.2, 6))
ax.frame(xticks=[0.5, 1, 1.5, 2, 2.5, 3], yticks=[-1, 0, 2, 4, 6], xlabel="trip length T (in half-periods, π/ω)", ylabel="curvature (nπ/T)² − 1 of the action")
ax.hline(0, color=DARK, width=1)
for j, col in enumerate([RED, PURPLE, BLUE, GREEN]): ax.plot(Ts/np.pi, lows[:, j], color=col, width=2.2)
for n in [1, 2, 3]: ax.vline(n, color="#bbb", width=1, dash="4,3")
ax.text(1.75, -0.45, "below 0: the classical path is a MAXIMUM", size=12, color=RED); ax.text(1.75, -0.8, "in that direction (one more per half-period)", size=12, color=RED)
ax.text(1.05, 5.4, "conjugate points (caustics) at T = π, 2π, 3π", size=12, color=GRAY)
# dial plane
for (xo, tt, short) in [(700, "T < π (no caustic yet)", True), (950, "T > π", False)]:
    bx = Axes(f, xo, 110, 200, 200, (-1.3, 1.3), (-1.3, 1.3))
    bx.frame(xticks=[], yticks=[], box=True)
    if short: f.rect(bx.X(0), bx.y, bx.w/2, bx.h, color="none", fill=GREEN, opacity=0.2)
    bx.hline(0, color="#aaa", width=1); bx.vline(0, color=GREEN if not short else "#aaa", width=5 if not short else 1)
    if not short: bx.dot(0, 1, 5, GREEN); bx.dot(0, -1, 5, GREEN)
    f.text(xo + 100, 100, tt, size=13, anchor="middle", weight="bold")
    f.text(xo + 100, 335, "allowed dials κ (weight e^(−κS/ħ))", size=12, anchor="middle")
    f.text(xo + 100, 352, "whole right half-plane" if short else "only the imaginary axis", size=12, anchor="middle", color=GREEN)
    bx.text(1.0, 0.08, "1", size=11); bx.text(0.08, 1.05, "i", size=11); bx.text(0.08, -1.2, "−i", size=11)
f.text(700, 400, "Left half-plane is always out: wiggly paths have huge kinetic action.", size=12)
f.text(700, 418, "Past the first caustic the right half is out too. Only Re κ = 0 survives.", size=12)
f.text(700, 436, "At each caustic the imaginary dial picks up a quarter turn e^(−iπ/2) (Maslov).", size=12)
f.text(80, 490, "Morse index theorem (Morse 1934); 'action is not least' (Gray and Taylor 2007). Checks R6–R9. The dial-plane reading is ours.", size=12, color=GRAY)
f.save(P("r3-caustics-leave-only-the-imaginary-dial.svg"))

# r4: checkerboard reversibility
def step(psi, k, eps, mu):
    if np.isclose(k, -1): c_, t_ = 1 - eps*mu, eps*mu
    else: c_, t_ = np.cos(eps*mu), 1j*np.sin(eps*mu)
    turned = np.stack([c_*psi[:, 0] + t_*psi[:, 1], c_*psi[:, 1] + t_*psi[:, 0]], axis=1)
    return np.stack([np.roll(turned[:, 0], 1), np.roll(turned[:, 1], -1)], axis=1)
n, eps, mu, N = 401, 0.05, 1.0, 300
def run(k, start):
    psi = np.zeros((n, 2), complex); psi[start, 0] = 1; out = [psi]
    for _ in range(N): psi = step(psi, k, eps, mu); out.append(psi)
    return np.array(out)
U1, U2 = run(-1j, 200), run(-1j, 202); K1, K2 = run(-1, 200), run(-1, 202)
dU = np.linalg.norm((U1 - U2).reshape(N + 1, -1), axis=1); dK = np.abs(K1 - K2).reshape(N + 1, -1).sum(axis=1)
f = Fig(1180, 520); f.title("Two starts two sites apart, run through the checkerboard chain: does the chain remember which was which?", size=16)
ax = Axes(f, 80, 70, 560, 360, (0, N), (0, 1.15))
ax.frame(xticks=[0, 100, 200, 300], yticks=[0, 0.5, 1], xlabel="steps", ylabel="how distinguishable (1 = fully)")
ax.plot(np.arange(N + 1), dU/dU[0], color=ORANGE, width=3); ax.plot(np.arange(N + 1), dK/dK[0], color=BLUE, width=3)
ax.text(140, 1.05, "k = −i (Feynman): exactly 1 at every step", size=13, color=ORANGE)
ax.text(60, 0.25, "k = −1 (Kac): total probability kept,", size=13, color=BLUE); ax.text(60, 0.18, "but the two histories blur together", size=13, color=BLUE)
p1 = np.abs(K1[-1]).sum(axis=1); p2 = np.abs(K2[-1]).sum(axis=1)
s = np.arange(n); m1 = (p1 > 0) & (s >= 150) & (s <= 250); m2 = (p2 > 0) & (s >= 150) & (s <= 250)   # each start lives on one sublattice
ym = 1.25*max(p1[m1].max(), p2[m2].max())
bx = Axes(f, 720, 70, 400, 360, (150, 250), (0, ym))
bx.frame(xticks=[(160, "−40"), (200, "0"), (240, "+40")], yticks=[0, round(ym/2, 2)], xlabel="site, after 300 steps", ylabel="probability")
bx.plot(s[m1], p1[m1], color=BLUE, width=2.5); bx.plot(s[m2], p2[m2], color=PURPLE, width=2, dash="5,4")
bx.text(153, 0.93*ym, "Kac, after 300 steps: start 200 (solid)", size=12, color=BLUE)
bx.text(153, 0.86*ym, "and start 202 (dashed) nearly coincide", size=12, color=PURPLE)
f.text(80, 490, "Turn rate 1, step 0.05, as on the checkerboard page. Checks R15–R18. Known physics (unitary vs stochastic); the side-by-side test is ours.", size=12, color=GRAY)
f.save(P("r4-checkerboard-remembers-only-at-imaginary-dial.svg"))

# r5: two routes
f = Fig(1180, 520); f.title("Two classical routes to the same place: the real dial keeps one, the imaginary dial keeps both", size=17)
xg = np.linspace(-2.5, 2.5, 4000); S1 = (xg + 1)**2/2; S2 = (xg - 1)**2/2
ax = Axes(f, 80, 70, 480, 340, (-2.5, 2.5), (0, 3))
ax.frame(xticks=[-2, -1, 0, 1, 2], yticks=[0, 1, 2, 3], xlabel="where you arrive", ylabel="effective action")
ax.plot(xg, np.minimum(S1, 3), color="#bbb", width=1.5, dash="4,3"); ax.plot(xg, np.minimum(S2, 3), color="#bbb", width=1.5, dash="4,3")
ax.text(1.55, 1.75, "route 2", size=11, color=GRAY); ax.text(-2.0, 1.75, "route 1", size=11, color=GRAY)
hb = 0.05; Sr = -hb*np.log(np.exp(-S1/hb) + np.exp(-S2/hb))
ax.plot(xg, Sr, color=BLUE, width=3)
ax.text(-2.4, 2.75, "real dial: −ħ log(sum) → the smaller action", size=12, color=BLUE)
ax.text(-2.4, 2.55, "a kink at 0; the other route is forgotten", size=12, color=BLUE)
bx = Axes(f, 650, 70, 480, 340, (-2.5, 2.5), (0, 4.4))
bx.frame(xticks=[-2, -1, 0, 1, 2], yticks=[0, 1, 2, 3, 4], xlabel="where you arrive", ylabel="|e^(iS₁/ħ) + e^(iS₂/ħ)|²")
hb = 0.08; I = np.abs(np.exp(1j*S1/hb) + np.exp(1j*S2/hb))**2
bx.plot(xg, I, color=ORANGE, width=1.3)
bx.text(-2.4, 4.2, "imaginary dial: fringes, both routes kept", size=12, color=ORANGE)
f.text(80, 495, "Real dial = Hopf–Cole / Hopf–Lax: viscosity solutions of Hamilton–Jacobi (Crandall–Lions 1983). Imaginary dial = interference (Feynman). Checks R13–R14.", size=12, color=GRAY)
f.save(P("r5-two-routes.svg"))

# r6 scorecard
f = Fig(1180, 420); f.title("Scorecard: what resolving the chain needs, by dial", size=18)
rows = [("Resolve onto a classical path, short trips", "yes", "yes"),
        ("Resolve onto the classical path at every trip length (stationary, not least)", "no", "yes"),
        ("Give Newton's mechanics, not the upside-down potential", "no", "yes"),
        ("Compose without losing anything (reversible steps)", "no", "yes"),
        ("Keep both routes where paths cross (caustics)", "no", "yes")]
f.text(60, 80, "Requirement", size=14, weight="bold"); f.text(820, 80, "real dial", size=14, weight="bold", anchor="middle", color=BLUE); f.text(1000, 80, "imaginary dial", size=14, weight="bold", anchor="middle", color=ORANGE)
for i, (r, a, b) in enumerate(rows):
    y = 110 + i*52
    if i % 2 == 0: f.rect(50, y - 4, 1080, 50, color="none", fill="#f5f6f8")
    f.text(60, y + 27, r, size=14)
    for xx, v in [(820, a), (1000, b)]:
        col = GREEN if v == "yes" else RED
        f.rect(xx - 45, y + 8, 90, 28, color=col, fill=col, rx=6, opacity=0.2); f.text(xx, y + 27, v, size=14, anchor="middle")
f.text(60, 395, "Row 1 is 'allows'. Rows 2 and 4 each force Re κ = 0 on their own. Row 3 says the real dial resolves onto the wrong classical mechanics.", size=13, color=GRAY)
f.save(P("r6-scorecard.svg"))
print("ok")
