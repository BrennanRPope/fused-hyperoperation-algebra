"""Figures for topics/why-i/e-i-pi. Run: python3 figures.py  (numpy only; writes SVG to figures/)."""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from svgkit import *
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"); os.makedirs(OUT, exist_ok=True)
P = lambda n: os.path.join(OUT, n)

# ---------------------------------------------------------------- f1 one dial, three constants
f = Fig(1100, 560); f.title("One dial, three constants: where the lift 1 + kx goes as x runs from 0")
ax = Axes(f, 70, 70, 430, 430, (-1.8, 3.0), (-1.8, 3.0))
ax.frame(xticks=[-1, 0, 1, 2, (np.e, "e")], yticks=[-1, 0, 1, 2], xlabel="real part of the lift", ylabel="imaginary part of the lift")
ax.hline(0, color="#ccc", width=1); ax.vline(0, color="#ccc", width=1)
x = np.linspace(0, 2*np.pi, 400)
ax.plot(np.cos(x), np.sin(x), color=ORANGE, width=2.5)
xr = np.linspace(0, 1.05, 50); ax.plot(np.exp(xr), 0*xr, color=BLUE, width=4)
xs = np.linspace(0, 2*np.pi, 400); k = 0.06 + 1j; L = np.exp(k*xs)
ax.plot(L.real, L.imag, color=PURPLE, width=1.8, dash="5,4")
ax.dot(1, 0, 6, DARK); ax.text(1.05, -0.2, "start (x = 0)", size=12)
ax.dot(np.e, 0, 6, BLUE); ax.text(np.e - 0.05, 0.15, "e", size=20, color=BLUE, weight="bold", anchor="middle")
ax.text(1.05, 0.42, "real dial k = 1:", size=12, color=BLUE); ax.text(1.05, 0.25, "at x = 1 the lift is e", size=12, color=BLUE)
ax.dot(0, 1, 6, ORANGE); ax.text(0.1, 1.12, "i  (x = π/2)", size=15, color=ORANGE, weight="bold")
ax.dot(-1, 0, 6, ORANGE); ax.text(-1.55, 0.18, "−1  (x = π)", size=15, color=ORANGE, weight="bold")
ax.text(-1.7, -1.45, "imaginary dial k = i: stays on the circle,", size=12, color=ORANGE)
ax.text(-1.7, -1.63, "halfway round at x = π, home at x = 2π", size=12, color=ORANGE)
ax.text(0.9, 2.55, "mixed dial k = 0.06 + i: a spiral", size=12, color=PURPLE)
ax.text(0.9, 2.35, "(grows and turns: the size is not kept)", size=12, color=PURPLE)
# right panel: the k plane
bx = Axes(f, 610, 70, 430, 430, (-2.2, 2.2), (-2.2, 2.2))
bx.frame(xticks=[-2, -1, 0, 1, 2], yticks=[(-1, "−i"), (1, "i"), (-2, "−2i"), (2, "2i")], xlabel="real part of the dial k", ylabel="imaginary part of k")
bx.hline(0, color=BLUE, width=4); bx.vline(0, color=ORANGE, width=4)
bx.text(1.1, -0.25, "growth / boosts", size=13, color=BLUE); bx.text(1.1, -0.45, "(e lives here)", size=12, color=BLUE)
bx.text(0.1, 1.85, "turns: size kept", size=13, color=ORANGE); bx.text(0.1, 1.62, "(i, and period 2π/|k|)", size=12, color=ORANGE)
for kk, lab in [(1, "k = 1"), (-1, "k = −1  (Kac, heat)"), (1j, "k = i"), (-1j, "k = −i  (Feynman)")]:
    bx.dot(kk.real if isinstance(kk, complex) else kk, kk.imag if isinstance(kk, complex) else 0, 6, DARK)
bx.text(1.0, 0.3, "k = 1", size=11); bx.text(-2.1, 0.3, "k = −1 (Kac, heat)", size=11)
bx.text(0.1, 1.08, "k = i", size=11); bx.text(0.1, -0.95, "k = −i (Feynman)", size=11)
bx.text(-2.05, -1.6, "everything off the two axes: spirals", size=12, color=PURPLE)
bx.text(-2.05, -1.8, "sign of Im k = which way round", size=12, color=GRAY); bx.text(-2.05, -2.0, "(the orientation, i.e. the twist)", size=12, color=GRAY)
f.save(P("f1-one-dial-three-constants.svg"))

# ---------------------------------------------------------------- f2 Euler by chaining
f = Fig(1100, 520); f.title("Chaining small steps: e on the real line, π and −1 round the wall (the circle operation, wall at −1)")
ax = Axes(f, 70, 60, 600, 400, (-3.0, 2.0), (-0.5, 2.0))
ax.frame(xticks=[-3, -2, -1, 0, 1, (np.e - 1, "e−1")], yticks=[0, 1, 2], xlabel="real part", ylabel="imaginary part")
ax.hline(0, color="#ccc", width=1)
t = np.linspace(0, np.pi, 200); ax.plot(-1 + np.cos(t), np.sin(t), color="#bbb", width=7)
ax.dot(-1, 0, 6, RED); ax.text(-1, -0.25, "wall −1", size=13, color=RED, anchor="middle")
cols = {12: PURPLE, 40: BLUE, 400: ORANGE}
for n, col in cols.items():
    pts = [(1 + 1j*np.pi/n)**m - 1 for m in range(n + 1)]
    ax.plot([p.real for p in pts], [p.imag for p in pts], color=col, width=2)
    if n <= 40:
        for p in pts: ax.dot(p.real, p.imag, 2.5, col)
    ax.dot(pts[-1].real, pts[-1].imag, 5, col)
for n, yy in [(1, 0.32), (2, 0.32), (4, 0.32), (16, 0.32), (400, 0.32)]:
    v = (1 + 1/n)**n - 1; ax.dot(v, 0, 5, BLUE if n < 400 else ORANGE)
    if n < 400: ax.text(v - (0.03 if n == 16 else 0), -0.25, str(n), size=11, anchor="middle", color=GRAY)
ax.dot(-2, 0, 7, DARK); ax.text(-2.0, -0.25, "−2", size=13, anchor="middle", weight="bold")
ax.text(-2.0, -0.42, "(lift −1)", size=12, anchor="middle")
ax.text(0.1, 0.45, "n steps of 1/n on the real line", size=12, color=BLUE)
ax.text(0.1, 0.28, "(n = 1, 2, 4, 16, 400) → e − 1", size=12, color=BLUE)
ax.text(-2.9, 1.85, "n steps of iπ/n turn round the wall", size=12)
ax.text(-2.9, 1.68, "and land ever closer to −2", size=12)
ax.text(-2.9, 1.40, "n = 12", size=13, color=PURPLE, weight="bold"); ax.text(-2.9, 1.22, "n = 40", size=13, color=BLUE, weight="bold"); ax.text(-2.9, 1.04, "n = 400", size=13, color=ORANGE, weight="bold")
ax.text(-0.3, 1.25, "grey: the exact half-turn", size=12, color=GRAY)
f.lines(720, 90, "What the picture says\n\n• The same rule, a ∘ b = a + b + ab, makes\n  both constants. Only the direction of the\n  steps differs.\n\n• Real steps compound: 1/n, n times → e − 1.\n\n• Imaginary steps turn about the wall −1:\n  iπ/n, n times → −2, so the lift\n  (1 + iπ/n)ⁿ → −1.  That is e^(iπ) = −1.\n\n• π is how far you go round the wall to\n  reach the far side; 2π brings you home.\n\nAll of this is Euler (1748). The wall reading\nis ours, from the complex-analysis and pi pages.", size=14)
f.save(P("f2-euler-by-chaining.svg"))

# ---------------------------------------------------------------- f3 i across physics
rows = [
 ("Waves and phasors", "a steady quarter turn of a (cos, sin) pair", "removable", "k imaginary: the lift turns"),
 ("AC circuits", "90° lag of current in coils and capacitors", "removable", "same turn, one per frequency"),
 ("Optics, Fourier", "a shift becomes a phase; size-keeping modes", "removable (Hartley)", "keep-size road"),
 ("Rotations, spinors", "an oriented plane (bivector) squares to −1", "i is a real geometric object", "the twist: oriented area"),
 ("Electromagnetism", "E ↔ cB quarter turn; boost = turn by iφ", "removable", "k real ↔ k imaginary flip"),
 ("Relativity", "ict, Wick rotation: the sign of the metric", "removable (ict dropped)", "Cayley–Klein κ sign"),
 ("Classical mechanics", "phase-space turn; symplectic area + length", "removable", "twist + length → J"),
 ("Causality", "future = upper half-plane (Kramers–Kronig)", "needed in the proof", "the wall's frequency split"),
 ("Heat, temperature", "temperature = period in imaginary time", "removable, but deep", "period 2π of the turn"),
 ("Quantum mechanics", "amplitudes compose; i in the Schrödinger eq.", "essential (tested 2022)", "keep-size + local tomography"),
]
f = Fig(1180, 640); f.title("Where i shows up in physics, what it encodes, and whether you can take it out")
X0, Y0, rh = 30, 70, 52
heads = [("Area", X0), ("What i encodes", X0 + 200), ("Can it be removed?", X0 + 590), ("Dial reading (ours where new)", X0 + 860)]
for h, xx in heads: f.text(xx, Y0, h, size=14, weight="bold")
for i, (a, b, c, d) in enumerate(rows):
    y = Y0 + 22 + i*rh
    if i % 2 == 0: f.rect(X0 - 10, y - 4, 1130, rh, color="none", fill="#f5f6f8")
    f.text(X0, y + 28, a, size=14, weight="bold")
    f.text(X0 + 200, y + 28, b, size=13)
    col = GREEN if c.startswith("removable") and "deep" not in c else (RED if c.startswith("essential") else AMBER)
    f.rect(X0 + 590, y + 10, 250, 28, color=col, fill=col, rx=6, opacity=0.18)
    f.text(X0 + 600, y + 29, c, size=13, color=DARK)
    f.text(X0 + 860, y + 28, d, size=13, color=PURPLE)
f.text(X0, 630, "Green: a real formulation exists and says the same.  Amber: i can be renamed but not removed (it becomes a geometric object or a proof tool).  Red: experiment rules out the real version.", size=12, color=GRAY)
f.save(P("f3-i-across-physics.svg"))

# ---------------------------------------------------------------- f4 boost = imaginary turn
f = Fig(1100, 520); f.title("A turn by angle θ traces a circle; a boost is a turn by the imaginary angle iφ and traces a hyperbola")
ax = Axes(f, 70, 70, 400, 380, (-2.2, 2.2), (-2.2, 2.2))
ax.frame(xticks=[-2, -1, 0, 1, 2], yticks=[-2, -1, 0, 1, 2], xlabel="E_y", ylabel="cB_z")
th = np.linspace(0, 2*np.pi, 300); ax.plot(np.cos(th), np.sin(th), color=ORANGE, width=3)
for t in [0, 0.5, 1.0]: ax.dot(np.cos(t), np.sin(t), 4, ORANGE)
ax.text(-0.95, 1.85, "duality turn F → e^(−iθ)F", size=13, color=ORANGE)
ax.text(-0.95, 1.62, "E_y² + c²B_z² is kept", size=12, color=ORANGE)
ph = np.linspace(-1.4, 1.4, 300)
ax.plot(np.cosh(ph), -np.sinh(ph), color=BLUE, width=3); ax.plot(-np.cosh(ph), np.sinh(ph), color=BLUE, width=3, dash="5,4")
for p in [0, 0.5, 1.0]: ax.dot(np.cosh(p), -np.sinh(p), 4, BLUE)
ax.text(0.2, -1.9, "boost of rapidity φ", size=13, color=BLUE); ax.text(0.2, -2.1, "E_y² − c²B_z² is kept", size=12, color=BLUE)
bx = Axes(f, 600, 70, 400, 380, (-2.2, 2.2), (-2.2, 2.2))
bx.frame(xticks=[-2, -1, 0, 1, 2], yticks=[-2, -1, 0, 1, 2], xlabel="x", ylabel="ct")
bx.plot([-2.2, 2.2], [-2.2, 2.2], color="#ccc", width=1, dash="3,3"); bx.plot([-2.2, 2.2], [2.2, -2.2], color="#ccc", width=1, dash="3,3")
bx.plot(np.cosh(ph), np.sinh(ph), color=BLUE, width=3)
bx.plot(-np.cosh(ph), -np.sinh(ph), color=PURPLE, width=3)
bx.dot(1, 0, 5, BLUE); bx.dot(-1, 0, 5, PURPLE)
bx.arrow(0.95, 0.15, -0.85, 0.15, color=RED)
bx.text(-0.8, 0.45, "rapidity iπ: (t, x) → (−t, −x)", size=13, color=RED)
bx.text(-0.8, 0.25, "the opposite wedge (Euler for boosts)", size=12, color=RED)
bx.text(0.2, 1.9, "boosts slide points", size=12, color=BLUE); bx.text(0.2, 1.7, "along the right wedge", size=12, color=BLUE)
bx.text(-1.55, -1.75, "x₄ = ict turns this into a circle;", size=12, color=GRAY)
bx.text(-1.55, -1.95, "dropped in 1973 (MTW, Box 2.1)", size=12, color=GRAY)
f.save(P("f4-boost-is-an-imaginary-turn.svg"))

# ---------------------------------------------------------------- f5 Kramers-Kronig
w00, gam = 1.0, 0.3
chi = lambda w: 1/(w00**2 - w**2 - 1j*gam*w)
f = Fig(1100, 500); f.title("Causality is analyticity: the real part of a causal response is fixed by its imaginary part")
ax = Axes(f, 80, 60, 620, 380, (0, 3), (-2, 3.6))
ax.frame(xticks=[0, 0.5, 1, 1.5, 2, 2.5, 3], yticks=[-2, -1, 0, 1, 2, 3], xlabel="frequency ω (resonance at 1)", ylabel="response χ(ω)")
ax.hline(0, color="#ccc", width=1)
w = np.linspace(0, 3, 600)
ax.plot(w, chi(w).real, color=BLUE, width=2.5); ax.plot(w, chi(w).imag, color=ORANGE, width=2.5)
hh = 0.004; grid = np.arange(-200, 200, hh) + hh/2
def kk(w0): g = grid; return np.sum(chi(g).imag/(g - w0))*hh/np.pi
for wt in np.arange(0.1, 3.0, 0.2):
    wt = np.round(wt/hh)*hh; ax.dot(wt, kk(wt), 4, DARK)
ax.text(1.25, 3.2, "Im χ (absorption)", size=13, color=ORANGE); ax.text(1.25, -1.5, "Re χ (dispersion)", size=13, color=BLUE)
ax.text(1.9, 1.0, "• = Re χ rebuilt from Im χ alone", size=13); ax.text(1.9, 0.75, "  (Kramers–Kronig)", size=13)
bx = Axes(f, 780, 90, 280, 280, (-1.6, 1.6), (-1.6, 1.6))
bx.frame(xticks=[-1, 0, 1], yticks=[-1, 0, 1], xlabel="Re ω", ylabel="Im ω")
bx.f.rect(bx.x, bx.y, bx.w, bx.h/2, color="none", fill=GREEN, opacity=0.12)
bx.hline(0, color=DARK, width=1)
poles = np.roots([-1, -1j*gam, w00**2])
for p in poles: bx.f.text(bx.X(p.real), bx.Y(p.imag) + 5, "×", size=20, anchor="middle", color=RED, weight="bold")
bx.text(-1.5, 1.35, "upper half: no poles", size=12, color=GREEN); bx.text(-1.5, 1.15, "= nothing before the cause", size=12, color=GREEN)
bx.text(-1.5, -1.4, "poles at Im ω = −γ/2", size=12, color=RED)
f.text(780, 440, "Flip the sign of i and the 'future' half-plane flips:", size=12, color=GRAY)
f.text(780, 458, "the same recipe then returns −Re χ (check C17).", size=12, color=GRAY)
f.save(P("f5-causality-kramers-kronig.svg"))

# ---------------------------------------------------------------- f6 local tomography as the circle operation
f = Fig(1180, 520); f.title("Local tomography is the circle operation of parameter counts: K_AB − 1 = a + b + ab, with a = K_A − 1")
def grid(x0, y0, a, b, cell, hidden=0, missing=0, label="", sub=""):
    for i in range(b + 1):
        for j in range(a + 1):
            col = "#cccccc" if i == 0 and j == 0 else (BLUE if i == 0 else (ORANGE if j == 0 else PURPLE))
            idx_from_end = (b - i)*a + (a - j)          # counts correlation cells from the bottom-right
            if i > 0 and j > 0 and idx_from_end < missing:
                f.rect(x0 + j*cell + 1, y0 + i*cell + 1, cell - 4, cell - 4, color=RED, width=2)
            else:
                f.rect(x0 + j*cell, y0 + i*cell, cell - 2, cell - 2, color="white", fill=col, opacity=0.85)
    f.text(x0, y0 - 30, label, size=16, weight="bold"); f.text(x0, y0 - 12, sub, size=12, color=GRAY)
    return x0 + (a + 1)*cell, y0 + (b + 1)*cell
xe, ye = grid(50, 120, 3, 3, 46, label="Complex (qubits)", sub="3 + 3 + 3·3 = 15 = K_AB − 1  ✓ exact")
xe, ye = grid(430, 120, 2, 2, 46, label="Real (rebits)", sub="2 + 2 + 2·2 = 8, but K_AB − 1 = 9")
f.rect(430 + 3*46 + 14, 120 + 2*46, 44, 44, color=RED, fill=RED, opacity=0.85)
f.text(430 + 3*46 + 64, 120 + 2*46 + 20, "J⊗J: a joint quarter-turn", size=12, color=RED)
f.text(430 + 3*46 + 64, 120 + 2*46 + 36, "no local test can see", size=12, color=RED)
x0 = 800; cell = 26
xe, ye = grid(x0, 120, 5, 5, cell, missing=8, label="Quaternionic", sub="5 + 5 + 5·5 = 35 needed, only 27 exist")
f.text(x0, 120 + 6*cell + 24, "red outlines: 8 needed cells that don't exist,", size=12, color=RED)
f.text(x0, 120 + 6*cell + 42, "so there is no standard way to combine", size=12, color=RED)
yl = 420
for col, s, xx in [("#cccccc", "normalization (the 1)", 50), (BLUE, "A's own parameters (a)", 260), (ORANGE, "B's own parameters (b)", 480), (PURPLE, "correlations (a·b)", 700)]:
    f.rect(xx, yl, 18, 18, color="white", fill=col); f.text(xx + 26, yl + 14, s, size=13)
f.text(50, 470, "Hardy (2001), Wootters (1990): K_AB = K_A·K_B. Read as (1 + a)(1 + b) = 1 + (a + b + ab), it is the weighted circle operation at k = 1 on parameter counts.", size=13)
f.text(50, 492, "The grid is the known Fano/Bloch split of a two-qubit state (15 = 3 + 3 + 9). Only the reading as the circle operation is ours.", size=13, color=GRAY)
f.save(P("f6-local-tomography-is-mult-add.svg"))

# ---------------------------------------------------------------- f7 axiom scorecard
items = [
 ("Hardy 2001", [("1 Probabilities as frequencies", "-"), ("2 Simplicity (K depends on N, least K)", "-"), ("3 Subspaces behave like smaller systems", "-"),
                 ("4 Composite systems: N_AB = N_A N_B, K_AB = K_A K_B", "Y"), ("5 Continuity: continuous reversible maps between pure states", "~")]),
 ("Masanes–Müller 2011", [("1 Finiteness", "-"), ("2 Local tomography", "Y"), ("3 Equivalence of subspaces", "-"), ("4 Symmetry (reversible maps reach every pure state)", "~"), ("5 All measurements allowed", "-")]),
 ("Chiribella–D'Ariano–Perinotti 2011", [("Causality", "-"), ("Perfect distinguishability", "-"), ("Ideal compression", "-"), ("Local distinguishability (= local tomography)", "Y"), ("Pure conditioning", "-"), ("Purification", "-")]),
]
f = Fig(1180, 620); f.title("The published axioms that derive complex quantum theory, checked against the framework")
y = 70
for name, ax_ in items:
    f.text(40, y, name, size=15, weight="bold"); y += 24
    for s, st in ax_:
        col, lab = {"Y": (GREEN, "expresses it exactly"), "~": (AMBER, "supplies the turning it needs"), "-": (GRAY, "not addressed")}[st]
        f.rect(60, y - 14, 220, 20, color=col, fill=col, rx=5, opacity=0.2); f.text(70, y + 1, lab, size=12)
        f.text(300, y + 1, s, size=13); y += 24
    y += 10
f.text(40, y + 6, "No axiom is violated: the framework says nothing about states, measurements or probabilities, so most axioms are simply outside it.", size=13)
f.text(40, y + 26, "Green: K_AB − 1 = (K_A − 1) ∘ (K_B − 1) (check B1). Amber: continuous reversible maps need a turning (imaginary) dial (checks A4, B7).", size=13, color=GRAY)
f.save(P("f7-axiom-scorecard.svg"))

# ---------------------------------------------------------------- f8 temperature as a period
f = Fig(1180, 540); f.title("Hawking's temperature: the one formula where e, i, π and ħ, c, G all meet")
cx, cy, R = 200, 260, 150
f.circle(cx, cy, R, color=GREEN, fill=GREEN, opacity=0.12, width=2)
for a in np.linspace(0, 2*np.pi, 13)[:-1]: f.line(cx, cy, cx + R*np.cos(a), cy + R*np.sin(a), color="#9bc", width=1)
for r in [50, 100]: f.circle(cx, cy, r, color="#9bc", width=1)
f.text(cx, cy + R + 30, "period 2π: smooth (no tip)", size=14, anchor="middle", color=GREEN, weight="bold")
f.text(cx, 80, "Euclidean horizon: time made imaginary", size=14, anchor="middle")
f.text(cx, 98, "becomes an angle round the horizon", size=14, anchor="middle")
# cone
x2, y2 = 520, 300
f.poly([(x2 - 110, y2 + 90), (x2, y2 - 120), (x2 + 110, y2 + 90)], color=RED, width=2, fill=RED, opacity=0.12, close=True)
f.ellipse(x2, y2 + 90, 110, 24, color=RED, width=2)
f.line(x2 - 110, y2 + 90, x2, y2 - 120, color=RED, width=2); f.line(x2 + 110, y2 + 90, x2, y2 - 120, color=RED, width=2)
f.text(x2 + 12, y2 - 118, "tip", size=13, color=RED)
f.text(cx + 8, cy - 8, "horizon", size=12, color=DARK)
f.line(cx + R + 8, cy, cx + R + 8, cy - 40, color=GREEN, width=2, arrow=True); f.text(cx + R + 14, cy - 18, "imaginary time τ", size=12, color=GREEN)
f.text(x2, y2 + 150, "any other period: a cone tip", size=14, anchor="middle", color=RED, weight="bold")
f.text(x2, y2 + 170, "(the horizon would be singular)", size=12, anchor="middle", color=RED)
X = 720
f.text(X, 120, "T  =  ħ c³ / (8π G M k_B)", size=26, weight="bold")
f.lines(X, 170, "• i : time t → −iτ (Wick); a thermal state is a\n   period in imaginary time (KMS, check C19)\n\n• π : the 2π is one full turn, forced by\n   smoothness at the horizon (check C21)\n\n• e : the thermal weight e^(−E/k_B T) is\n   evolution e^(−iHt/ħ) at imaginary time\n\n• ħ, c, G, k_B all appear; G only through the\n   horizon's surface gravity κ = c⁴/(4GM)\n\nKnown: Hawking 1975, Gibbons–Hawking 1977,\nUnruh 1976, Bisognano–Wichmann 1975.", size=14)
f.save(P("f8-temperature-is-a-turn.svg"))

# ---------------------------------------------------------------- f9 QG dial ledger
f = Fig(1180, 515); f.title("Which dials a quantum-gravity picture needs, read through e, i and π")
boxes = [
 (40, 80, BLUE, "c  →  signature dial", "Cayley–Klein κ₂ = −1/c². Wick rotation flips\nthis dial from real (boosts) to imaginary\n(turns): k → ik. Already in the framework."),
 (620, 80, ORANGE, "ħ  →  order dial", "Sets the size of the imaginary dial (k = −i/ħ)\nand the imaginary-time period βħ.\nDirection forced (keep size); size not."),
 (40, 270, GREEN, "Λ  →  curvature dial", "Cayley–Klein κ₁: de Sitter / anti-de Sitter.\nIn de Sitter the horizon also has a 2π period\n(Gibbons–Hawking temperature). In the framework."),
 (620, 270, RED, "G  →  no dial", "Enters Hawking's T only through the horizon:\nκ = c⁴/(4GM). No composition law forces it\n(dial-map page). This is the open piece."),
]
for x0, y0, col, head, body in boxes:
    f.rect(x0, y0, 520, 150, color=col, fill=col, rx=12, opacity=0.10); f.rect(x0, y0, 520, 150, color=col, rx=12, width=2)
    f.text(x0 + 20, y0 + 34, head, size=20, weight="bold", color=col)
    f.lines(x0 + 20, y0 + 70, body, size=14)
f.text(40, 470, "e, i, π tell us the needed kinds of dial: one that can be turned imaginary (c, via Wick), one that sets its size (ħ), and a compact period 2π that any horizon must honour.", size=13)
f.text(40, 492, "They do not supply G. The common proof (Euler) links the dials' directions, not their sizes.", size=13, color=GRAY)
f.save(P("f9-quantum-gravity-dial-ledger.svg"))
print("wrote", sorted(os.listdir(OUT)))
