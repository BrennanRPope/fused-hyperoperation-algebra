"""Figures for asymptotic-safety/beyond-one-loop. Run: python3 figures.py (numpy + PIL via ../../svgkit.py)."""
import numpy as np, os, sys
from math import pi
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from svgkit import *
from flows import *
_lines = Fig.lines
Fig.lines = lambda self, x, y, s, **kw: _lines(self, x, y, "\n".join(s) if isinstance(s, list) else s, **kw)
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True); P = lambda n: os.path.join(OUT, n)

gs0 = root(lambda x: 2 + x*B1_0/(1 - x*B2_0), 0.5, 5.0); th0 = 23/11
akf = lambda G: ak_two_loop(G, 2.0); Gs2 = root(lambda G: akf(G)/G, 1e-9, 200.0)
wq, b2q = 1.0, 0.18; brf = lambda x: br_full(x, wq, wq + b2q, b2q); gsb = 1/(wq + b2q)
flows1 = [("one loop / Bonanno–Reuter closed form", lambda x: one_loop(x, 1.0), 1.0, BLUE),
          ("EH, Litim cutoff, λ = 0  (θ = 23/11)", litim0, gs0, GREEN),
          ("Bonanno–Reuter's own β  (ω′/ω = 1.18)", brf, gsb, ORANGE),
          ("two loops in 2 + ε, taken at ε = 2", akf, Gs2, RED)]

# ---------- b1: flow lines in the flat coordinate ----------
f = Fig(1180, 620); f.title("In the circle operation's flat coordinate, only the one-loop shape is a straight line", size=17)
ax = Axes(f, 90, 80, 470, 400, (0, 6), (-9, 8))
ax.frame(xticks=range(0, 7), yticks=range(-8, 9, 2), xlabel="RG time t = log k (toward the UV →)", ylabel="flat coordinate  log |1 − g*/g|", grid=True)
for name, fb, gs, col in flows1:
    T = 6.0; n = 1200; tr = rk4(lambda v: np.array([fb(v[0])]), [1e-3*gs], T/n, n); g = tr[:, 0]
    ax.plot(np.arange(n + 1)*T/n, np.log(np.abs(1 - gs/g)), color=col, width=2.5)
ax.text(0.15, -7.4, "one loop: slope exactly −(d − 2) = −2 all the way", size=12, color=BLUE)
ax.text(0.15, -8.3, "the others start at slope −2 and bend to −θ at the wall", size=12, color=GRAY)
ax2 = Axes(f, 660, 80, 470, 400, (0, 1), (1.8, 3.0))
ax2.frame(xticks=[0, 0.25, 0.5, 0.75, 1], yticks=[2.0, 2.2, 2.4, 2.6, 2.8, 3.0], xlabel="g / g*   (Gaussian point → UV wall)", ylabel="local rate  −d/dt log|1 − g*/g|", grid=True)
for name, fb, gs, col in flows1:
    x = np.linspace(1e-4, 0.999, 400)*gs; rate = fb(x)*(gs/x**2)/np.abs(1 - gs/x)*np.sign(gs - x)
    ax2.plot(x/gs, rate, color=col, width=2.5)
for i, (name, fb, gs, col) in enumerate(flows1):
    f.rect(660, 515 + i*0 , 0, 0); f.line(90 + (i % 2)*540, 528 + (i//2)*24, 120 + (i % 2)*540, 528 + (i//2)*24, color=col, width=3)
    f.text(128 + (i % 2)*540, 533 + (i//2)*24, name, size=13)
f.lines(90, 588, ["A circle-op flow (a o_k b = a + b + kab in u = 1/g) is exactly a flat line on the right. Flat at 2 needs θ = d − 2 at the wall. Checks P1–P4, T2.",
                  "Known: the four β functions (Bonanno–Reuter 2000; Litim 2004 / Reuter–Saueressig; Aida–Kitazawa 1997). The flat-coordinate test is ours."], size=12, color=GRAY)
f.save(P("b1-flat-coordinate.svg"))

# ---------- b2: one number decides ----------
f = Fig(1180, 600); f.title("One coordinate-free number decides it: the UV exponent θ divided by G's units d − 2", size=17)
rows = [("One-loop shape (Weinberg 1979)", 1.0, BLUE, "exact circle op"),
        ("Bonanno–Reuter closed form G0/(1 + ωG0k²)", 1.0, BLUE, "exact circle op (it is the one-loop solution)"),
        ("Large N matter, g direction (Donà–Eichhorn–Percacci)", 1.0, BLUE, "exact; and Λ/G is a second one at rate d"),
        ("EH, Litim cutoff, λ = 0", 23/22, GREEN, "23/22: off by 4.5%, set by B₂"),
        ("Bonanno–Reuter's own β (their cutoff)", 1.18, ORANGE, "ω′/ω = 1.18: off by 18%"),
        ("Two loops in 2 + ε, at ε = 2 (Aida–Kitazawa)", 2.8264/2, RED, "1.41: off at order ε"),
        ("Full EH (g, λ), pure gravity", None, PURPLE, "complex 1.475 ± 3.043i: no circle-op form at all")]
ax = Axes(f, 470, 70, 560, 420, (0.9, 1.5), (0, len(rows)))
for i, (lab, val, col, note) in enumerate(rows):
    y = len(rows) - i - 0.5
    f.text(455, ax.Y(y) + 5, lab, size=13, anchor="end")
    if val is not None:
        f.rect(ax.X(0.9), ax.Y(y) - 13, ax.X(val) - ax.X(0.9), 26, color=col, fill=col, opacity=0.35, rx=4)
        f.text(ax.X(val) + 8, ax.Y(y) + 5, note, size=12, color=DARK)
    else:
        f.text(ax.X(1.01), ax.Y(y) + 5, note, size=12, color=col, weight="bold")
f.line(ax.X(1.0), ax.y - 5, ax.X(1.0), ax.y + ax.h, color=DARK, width=2, dash="6,4")
f.text(ax.X(1.0), ax.y - 10, "= 1: circle op holds with both ends as true walls", size=12, anchor="middle")
for t in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]:
    f.text(ax.X(t), ax.y + ax.h + 20, f"{t:g}", size=12, anchor="middle", color="#444")
f.text(ax.X(1.2), ax.y + ax.h + 42, "θ / (d − 2)", size=13, anchor="middle")
f.lines(60, 560, ["Why this number: a circle-op (Möbius) flow always has exponents +r and −r at its two walls, whatever the dials k and d (check T1). Smooth changes",
                  "of coordinate cannot change exponents. So θ_UV = d − 2 is necessary and, on the interval, sufficient (1-D linearisation). Checks P1–P6, T1–T6, L5."], size=12, color=GRAY)
f.save(P("b2-one-number.svg"))

# ---------- b3: Lambda as the second dial ----------
f = Fig(1180, 640); f.title("Does Λ join as a second dial? At large N yes, as a second circle op; with graviton loops, no", size=17)
D = 4; aa, bb = 1.0, 0.5
btN = lambda v: np.array([(D - 2)*v[0] - aa*v[0]**2, -2*v[1] - aa*v[0]*v[1] + bb*v[0]])
gsN, lsN = 2/aa, (bb/D)*(2/aa)
ax = Axes(f, 80, 90, 460, 390, (0, 2.6), (-0.6, 1.2))
ax.frame(xticks=[0, 0.5, 1, 1.5, 2, 2.5], yticks=[-0.5, 0, 0.5, 1.0], xlabel="g (Newton coupling)", ylabel="λ (cosmological constant)", grid=True)
for l0 in np.linspace(-0.5, 1.1, 9):
    for g0 in (0.02, 2.55):
        tr = rk4(btN, [g0, l0], 0.004, 2500); ax.plot(tr[:, 0], tr[:, 1], color=BLUE, width=1.6, opacity=0.8)
ax.dot(0, 0, 6, color=DARK); ax.dot(gsN, lsN, 7, color=RED)
ax.text(0.05, -0.08, "Gaussian", size=12); ax.text(gsN + 0.05, lsN + 0.06, "UV point: θ = 2 and 4", size=12, color=RED)
f.text(310, 72, "Large N (matter loops only)", size=15, anchor="middle", weight="bold", color=BLUE)
f.lines(80, 536, ["In u = 1/g and v = λ/g: two separate circle-op flows,", "rates d − 2 = 2 (G's units) and d = 4 (vacuum energy's units)."], size=12)
beh = eh_beta(); v0 = newton(beh, [0.7, 0.2])
ax = Axes(f, 650, 90, 460, 390, (0, 1.2), (-0.3, 0.45))
ax.frame(xticks=[0, 0.3, 0.6, 0.9, 1.2], yticks=[-0.2, 0, 0.2, 0.4], xlabel="g", ylabel="λ", grid=True)
J = jac(beh, v0); w, V = np.linalg.eig(J)
for ph in np.linspace(0, 2*pi, 7)[:-1]:
    z0 = v0 + np.array([0.25*np.cos(ph), 0.12*np.sin(ph)])
    tr = rk4(beh, z0, 0.002, 2500); inside = (tr[:, 0] > 0) & (tr[:, 0] < 1.2) & (tr[:, 1] > -0.3) & (tr[:, 1] < 0.45)
    ax.plot(tr[inside, 0], tr[inside, 1], color=PURPLE, width=1.6, opacity=0.85)
for ph in (0.0, 1.0, 2.0, 3.0, 4.2, 5.2):
    tr = rk4(beh, v0 + 1e-4*np.real(V[:, 0]*np.exp(1j*ph)), -0.01, 1400)
    out = ~((tr[:, 0] > 0) & (tr[:, 0] < 1.2) & (tr[:, 1] > -0.3) & (tr[:, 1] < 0.45)); jump = np.r_[False, np.linalg.norm(np.diff(tr, axis=0), axis=1) > 0.01]; bad = out | jump | ~np.all(np.isfinite(tr), axis=1); stop = np.argmax(bad) if bad.any() else len(tr)
    ax.plot(tr[:stop, 0], tr[:stop, 1], color=GRAY, width=1.2, opacity=0.7)
ax.dot(0, 0, 6, color=DARK); ax.dot(v0[0], v0[1], 7, color=RED)
ax.text(v0[0] + 0.03, v0[1] + 0.035, "UV point: θ = 1.475 ± 3.043i", size=12, color=RED)
f.text(880, 72, "Pure gravity, EH truncation (Litim)", size=15, anchor="middle", weight="bold", color=PURPLE)
f.lines(650, 536, ["Real exponents (2, −2) at the Gaussian point, a complex pair at the UV point:",
                   "no projective (two-coupled-circle-op) flow can do both (checks L3, L5)."], size=12)
f.lines(80, 600, ["Known: the large-N exponents {4, 2} and λ* = 3/8 (Donà–Eichhorn–Percacci 2013, Table II); the EH spiral (Reuter 1998; Litim 2004).",
                  "Ours: the decoupled circle ops in (1/g, λ/g) with rates d − 2 and d (L1, L4); the no-go for the pure-gravity flow (L5). Purple: flows into the UV point; gray: away from it."], size=12, color=GRAY)
f.save(P("b3-lambda-second-dial.svg"))

# ---------- b4: exponent path in N_V ----------
f = Fig(1180, 620); f.title("Adding matter walks the spiral's exponents onto the circle-op values {2, 4}; graviton loops are the 1/N break", size=17)
NV = np.concatenate([np.arange(0, 60, 0.5), np.logspace(np.log10(60), 6, 200)]); vcur = v0.copy(); TH = []
for N in NV:
    bt = eh_beta(-4*N, 2*N); vcur = newton(bt, vcur, 30); TH.append(np.sort_complex(-np.linalg.eigvals(jac(bt, vcur))))
TH = np.array(TH)
ax = Axes(f, 110, 70, 640, 460, (1, 4.5), (-3.5, 3.5))
ax.frame(xticks=[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5], yticks=[-3, -2, -1, 0, 1, 2, 3], xlabel="Re θ", ylabel="Im θ", grid=True)
ax.hline(0, color="#bbb")
cx = np.any(np.abs(TH.imag) > 1e-9, axis=1); ic = np.where(cx)[0][-1] + 1
up = TH[:ic].real[:, 0] + 1j*np.abs(TH[:ic].imag[:, 0])
ax.plot(up.real, up.imag, color=PURPLE, width=2.5); ax.plot(up.real, -up.imag, color=PURPLE, width=2.5)
lo, hi = np.sort(TH[ic:].real, axis=1).T
ax.plot(np.r_[up.real[-1], lo], np.r_[0, 0*lo], color=PURPLE, width=2.5); ax.plot(np.r_[up.real[-1], hi], np.r_[0, 0*hi], color=PURPLE, width=2.5)
for N, lab in [(0, "N_V = 0 (pure gravity)"), (10, "10"), (100, "100"), (300, "300")]:
    i = np.argmin(np.abs(NV - N)); z = TH[i, 1]; ax.dot(z.real, z.imag, 5, color=PURPLE); ax.dot(z.real, -z.imag, 5, color=PURPLE)
    ax.text(z.real + 0.06, abs(z.imag) + 0.1, lab, size=12, color=PURPLE)
for N, lab in [(1e4, "10⁴"), (1e6, "10⁶")]:
    i = np.argmin(np.abs(NV - N)); [ax.dot(z.real, 0, 4, color=PURPLE) for z in TH[i]]
ax.dot(2, 0, 8, color=BLUE, fill="none"); ax.dot(4, 0, 8, color=BLUE, fill="none")
ax.text(2.0, -0.35, "2 = d − 2", size=12, color=BLUE, anchor="middle"); ax.text(4.0, -0.35, "4 = d", size=12, color=BLUE, anchor="middle")
ireal = np.argmax(np.all(np.abs(TH.imag) < 1e-9, axis=1)); ax.text(TH[ireal, 0].real, 0.3, f"real from N_V ≈ {NV[ireal]:.0f}", size=12, color=DARK, anchor="middle")
cpr = [(1.475, 3.043), (1.689, 2.486), (2.205, 3.214), (2.425, 1.270), (1.752, 2.069)]
for a, b in cpr: ax.dot(a, b, 4, color=GRAY); ax.dot(a, -b, 4, color=GRAY)
ax.text(2.45, 1.45, "pure gravity, five cutoffs", size=11, color=GRAY); ax.text(2.45, 1.2, "(CPR 2008 Table II): all complex", size=11, color=GRAY)
f.lines(790, 110, ["How to read it", "", "Purple: the UV exponents of the EH flow", "with N_V gauge fields added.", "",
                   "Blue rings: the values a pair of circle-op", "flows needs, θ = {d − 2, d} = {2, 4}.", "",
                   "Pure gravity sits far off the axis.", "Matter loops pull the pair down to the", "axis and onto the rings; the leftover", "shrinks like 1100/N_V here.", "",
                   "So the circle-op structure is exactly", "the large-N limit, and graviton loops", "are what break it."], size=13)
f.lines(110, 578, ["Illustrative hybrid: Litim graviton terms (parent page) + Donà–Eichhorn–Percacci matter coefficients (their eqs. 35–36). Endpoint {4, 2} and λ* = 3/8 agree with their Table II.",
                   "Gray dots: published pure-gravity EH exponents for five cutoffs (Codello–Percacci–Rahmede 2008). Checks L8, L9."], size=12, color=GRAY)
f.save(P("b4-exponent-path.svg"))

# ---------- b5: verdict board ----------
f = Fig(1180, 660); f.title("Does G's running stay the circle operation beyond one loop? Claim by claim", size=17)
rows = [("Bonanno–Reuter's G(k) = G0/(1 + ωG0k²) is exactly the circle op, k = −g*", "exact; it is the one-loop solution", BLUE),
        ("The test: θ_UV/(d − 2) = 1, coordinate-free", "ours (standard Möbius + linearisation)", PURPLE),
        ("EH, Litim, λ = 0: θ = 23/11, ratio 23/22", "broken by 4.5%, set by B₂ (ours)", ORANGE),
        ("Two loops in 2 + ε: θ = ε + (3/5)ε² + …", "broken at order ε (ours from AK's β)", RED),
        ("A new flat coordinate makes it a circle op again", "only by turning the wall into a corner", ORANGE),
        ("The plane dial d can absorb the correction", "no: d keeps the exponents ±r", RED),
        ("Large N: 1/g and Λ/G are two decoupled circle ops, rates d − 2 and d", "ours, exact (exponents known: DEP)", GREEN),
        ("Pure gravity (g, λ): Λ as a second circle-op dial", "impossible: complex pair (ours)", RED),
        ("A new testable prediction", "none: it is Sakharov-type running, rewritten", RED),
        ("Worth keeping as a lens", "yes: a clean 1/N-and-B₂ diagnostic", GREEN)]
for i, (claim, st, col) in enumerate(rows):
    y = 62 + i*54
    if i % 2 == 0: f.rect(30, y, 1120, 54, color="none", fill="#f5f6f8")
    f.text(45, y + 33, claim, size=14)
    f.rect(790, y + 12, 350, 30, color=col, fill=col, opacity=0.18, rx=6); f.text(965, y + 32, st, size=12, anchor="middle", color=DARK)
f.text(45, 638, "Blue/green: holds exactly. Purple: our test. Orange: holds only in part. Red: does not hold.", size=12, color=GRAY)
f.save(P("b5-verdict-board.svg"))
print("ok")
