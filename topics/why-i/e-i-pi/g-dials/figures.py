"""Figures for topics/why-i/e-i-pi/g-dials. Run: python3 figures.py (numpy only)."""
import numpy as np, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from svgkit import *
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True); P = lambda n: os.path.join(OUT, n)
G, hbar, c, kB, Lam, mp, Msun = 6.67430e-11, 1.054571817e-34, 2.99792458e8, 1.380649e-23, 1.1056e-52, 1.67262192e-27, 1.98847e30

# g1: G needs a partner -- pure numbers that contain G
f = Fig(1180, 560); f.title("G has units, so it can't be a dial by itself: every pure number with G in it pairs G with something else", size=17)
items = [(np.log10(hbar*G*Lam/c**3), "G ħ Λ / c³", "Planck length² × cosmological curvature (MacDowell–Mansouri, Freidel–Starodubtsev)", RED),
         (np.log10(G*mp**2/(hbar*c)), "G m_p² / (ħ c)", "gravity's 'fine-structure constant' for a proton", PURPLE),
         (np.log10(13e3/1.2209e19), "E / E_P", "LHC collision energy over the Planck energy (κ-Poincaré dial)", BLUE),
         (np.log10(G*5.972e24/(6.371e6*c**2)), "G M / (R c²)", "Earth's compactness", GREEN),
         (np.log10(G*Msun/(6.957e8*c**2)), "G M / (R c²)", "the Sun's compactness", GREEN),
         (np.log10(0.2), "G M / (R c²)", "a neutron star", GREEN),
         (np.log10(0.5), "G M / (R c²) = ½", "a black hole: the horizon (the wall)", DARK)]
ax = Axes(f, 60, 90, 1060, 60, (-125, 2), (0, 1))
f.line(ax.X(-125), ax.Y(0.5), ax.X(2), ax.Y(0.5), color=DARK, width=2)
for v in [-120, -100, -80, -60, -40, -20, 0]:
    f.line(ax.X(v), ax.Y(0.5) - 6, ax.X(v), ax.Y(0.5) + 6, color=DARK); f.text(ax.X(v), ax.Y(0.5) + 26, f"10^{v}" if v else "1", size=12, anchor="middle")
for i, (v, lab, desc, col) in enumerate(items):
    y = 210 + i*44
    f.circle(ax.X(v), ax.Y(0.5), 6, color=col, fill=col)
    f.line(ax.X(v), ax.Y(0.5) + 8, ax.X(v), y - 6, color=col, width=1, dash="3,3")
    xt = ax.X(v) + 10 if v < -60 else ax.X(v) - 10
    num = "" if v > -0.5 else f"  ≈ 10^{v:.0f}"
    f.text(xt, y + 4, f"{lab}{num}: {desc}", size=13, color=col, anchor="start" if v < -60 else "end")
f.text(60, 535, "Dimension count (check D1): with ħ, c and Λ there is exactly one pure number containing G, namely GħΛ/c³ ≈ 2.9×10⁻¹²².", size=13, color=GRAY)
f.save(P("g1-g-needs-a-partner.svg"))

# g2: scorecard
f = Fig(1180, 690); f.title("Candidate dials for G, judged as dials: a continuous knob with a stable endpoint, like c and ħ", size=17)
cols = [("Candidate", 40), ("what it is a dial on", 300), ("stable end", 590), ("form forced", 720), ("size fixed", 850), ("fits our dial", 980)]
for t, x in cols: f.text(x, 75, t, size=13, weight="bold")
rows = [("Two-scale stable algebra (Mendes)", "whether positions commute (length ℓ)", "yes", "yes", "no", "yes"),
        ("G ħ Λ/c³ (MacDowell–Mansouri, F–S)", "distance from topological BF theory", "yes", "yes", "no", "partly"),
        ("Spin-2 self-coupling √G (Deser, B–H)", "how gravitons couple to themselves", "yes", "yes", "no", "no"),
        ("κ-Poincaré / relative locality", "how momenta combine (E_P)", "yes", "no", "bound", "yes"),
        ("Barbero–Immirzi γ", "the Hodge-star (i-like) direction", "yes", "yes", "count", "yes"),
        ("Jacobson / Sakharov", "entropy per area η = c³/4Għ", "yes", "yes", "no", "partly"),
        ("Asymptotic safety g = Gk²", "energy scale (running G)", "two", "flow", "no", "weak"),
        ("Minimal length (GUP β)", "the commutator [x, p]", "yes", "no", "no", "partly"),
        ("Newton–Cartan / Carroll", "c (G rides along)", "yes", "—", "—", "c's dial"),
        ("Teleparallel couplings", "which torsion combination", "yes", "no", "no", "no")]
tone = {"yes": GREEN, "no": RED, "partly": ORANGE, "weak": ORANGE, "bound": ORANGE, "count": ORANGE, "two": GREEN, "flow": ORANGE, "—": GRAY, "c's dial": GRAY}
for i, r in enumerate(rows):
    y = 95 + i*52
    if i % 2 == 0: f.rect(30, y, 1120, 52, color="none", fill="#f5f6f8")
    f.text(40, y + 31, r[0], size=13, weight="bold" if i < 3 else "normal"); f.text(300, y + 31, r[1], size=12)
    for (t, x), v in zip(cols[2:], r[2:]):
        f.rect(x - 5, y + 13, 100, 26, color=tone[v], fill=tone[v], rx=6, opacity=0.2); f.text(x + 45, y + 31, v, size=12, anchor="middle")
f.text(40, 640, "'bound': E_P > 1.2 E_Planck for the linear effect (Fermi GRB 090510). 'count': γ = 0.2375 or 0.274 from black-hole state counting.", size=12, color=GRAY)
f.text(40, 660, "'two': a Gaussian point and an ultraviolet fixed point. The first three rows meet in one number (next picture). Checks in checks.py.", size=12, color=GRAY)
f.save(P("g2-scorecard.svg"))

# g3: the stability square
f = Fig(1180, 600); f.title("The same stability test that forces c and ħ forces two more scales: a length ℓ and a curvature radius L", size=17)
pts = {(0, 0): (230, 470), (1, 0): (230, 160), (0, 1): (640, 470), (1, 1): (640, 160)}
lab = {(0, 0): ("Poincaré + Heisenberg", "can still deform: dim H² = 7 (3+1: 3)", RED),
       (1, 0): ("+ minimal length ℓ", "positions stop commuting; dim H² = 1", ORANGE),
       (0, 1): ("+ curvature radius L", "curved momenta / de Sitter; dim H² = 1", ORANGE),
       (1, 1): ("both: rigid", "dim H² = 0 (1+1 and 3+1)", GREEN)}
for (x1, y1), (x2, y2) in [((230, 432), (230, 198)), ((640, 432), (640, 198)), ((380, 470), (490, 470)), ((380, 160), (490, 160))]:
    f.line(x1, y1, x2, y2, color="#999", width=2.5, arrow=True)
for k, (x, y) in pts.items():
    t1, t2, col = lab[k]; f.rect(x - 150, y - 38, 300, 76, color=col, fill=col, rx=10, opacity=0.15)
    f.text(x, y - 6, t1, size=15, anchor="middle", weight="bold", color=col); f.text(x, y + 18, t2, size=12, anchor="middle")
f.text(215, 330, "turn on ℓ", size=12, color=GRAY, rotate=-90); f.text(435, 455, "turn on 1/L", size=12, color=GRAY, anchor="middle")
f.lines(830, 150, "With ħ fixed by [x, p] = iħ, the rigid algebra\ncarries ONE pure number: ℓ² / L².\n\nIf ℓ is the Planck length and L the de Sitter\nradius, ℓ²/L² = G ħ Λ / (3c³) ≈ 10⁻¹²²,\nthe same number MacDowell–Mansouri and\nFreidel–Starodubtsev expand gravity in.\n\nThen G = c³ ℓ² / ħ: G is ℓ's square in\nother units, not a dial of its own.", size=14)
f.text(830, 420, "Known: the algebra and its rigidity", size=13, color=GREEN); f.text(830, 440, "(Mendes 1994; Chryssomalakos–Okon 2004).", size=13, color=GREEN)
f.text(830, 470, "Ours: the meeting with G ħ Λ/c³.", size=13, color=ORANGE); f.text(830, 490, "Open: is ℓ really the Planck length?", size=13, color=RED)
f.text(60, 575, "dim H²(g, g) computed in checks S2 and S4 (Chevalley–Eilenberg, numpy). H² = 0 means no small change of the brackets is possible: the algebra is stable.", size=12, color=GRAY)
f.save(P("g3-stability-forces-two-scales.svg"))

# g4: Immirzi dial planes
f = Fig(1180, 560); f.title("The Barbero–Immirzi parameter is a dial on the Hodge star, which plays i: its self-dual value sits on the wall", size=17)
for xo, ttl, wall, col in [(80, "Lorentzian (★² = −1)", "i", RED), (620, "Euclidean (★² = +1)", "1", BLUE)]:
    ax = Axes(f, xo, 90, 380, 380, (-2.4, 2.4), (-2.4, 2.4))
    ax.frame(xticks=[], yticks=[], box=True); ax.hline(0, color="#aaa", width=1); ax.vline(0, color="#aaa", width=1)
    f.text(xo + 190, 80, ttl, size=15, anchor="middle", weight="bold")
    ax.plot([-2.4, 2.4], [0, 0], color=GREEN, width=4)
    ax.text(1.0, 0.12, "real 1/γ: Barbero (always invertible)", size=11, color=GREEN) if wall == "i" else ax.text(-2.3, 0.15, "real 1/γ", size=11, color=GREEN)
    ax.dot(0, 0, 6, DARK); ax.text(0.08 if wall == "i" else -0.75, -0.3 if wall == "i" else -1.3, "1/γ = 0: Palatini", size=11)
    if wall == "i":
        for s in [1, -1]: ax.dot(0, s, 8, col); ax.text(0.12, s + 0.1, ("1/γ = " + ("−i" if s > 0 else "i") + ": self-dual (Ashtekar)"), size=11, color=col)
        ax.arrow(1.6, -1.45, 2.3, -1.45, color=PURPLE); ax.text(0.3, -1.75, "black-hole count γ = 0.2375 is at 1/γ = 4.2 →", size=11, color=PURPLE)
    else:
        for s in [1, -1]: ax.dot(s, 0, 8, col)
        ax.text(0.6, 0.95, "wall 1/γ = ±1", size=11, color=col); ax.text(0.6, 0.7, "UV fixed point γ² = 1", size=11, color=col)
        ax.text(0.6, 0.45, "(Benedetti–Speziale 2011)", size=11, color=col)
        ax.arrow(2.1, -0.6, 1.1, -0.1, color=col); ax.arrow(0.2, -0.6, 0.9, -0.1, color=col); ax.text(0.2, -0.9, "renormalization flow", size=11, color=col)
    f.text(xo + 190, 495, f"det(1 + ★/γ) = (1 {'+' if wall == 'i' else '−'} 1/γ²)³ is zero only at 1/γ = ±{wall}", size=13, anchor="middle")
f.text(520, 300, "Wick", size=13, anchor="middle", color=GRAY); f.line(480, 280, 560, 280, color=GRAY, arrow=True)
f.text(60, 540, "Known: γ = ±i (Lorentzian) / ±1 (Euclidean) give self-dual variables; checks I1–I4. Ours: reading those values as the wall of the lift 1 + ★/γ.", size=12, color=GRAY)
f.save(P("g4-immirzi-wall.svg"))

# g5: state counting fixes gamma
f = Fig(1180, 540); f.title("Black-hole state counting picks γ: where each counting rule crosses 1", size=18)
js = np.arange(1, 4001)/2; rj = np.sqrt(js*(js + 1)); gg = np.linspace(0.05, 0.45, 400)
F1 = np.array([np.sum(2*np.exp(-2*np.pi*g*rj)) for g in gg]); F2 = np.array([np.sum((2*js + 1)*np.exp(-2*np.pi*g*rj)) for g in gg]); F3 = 2*np.exp(-2*np.pi*gg*np.sqrt(3)/2)
ax = Axes(f, 90, 70, 640, 380, (0.05, 0.45), (0, 3))
ax.frame(xticks=[0.1, 0.2, 0.3, 0.4], yticks=[0, 1, 2, 3], xlabel="Barbero–Immirzi parameter γ", ylabel="sum over puncture spins")
ax.hline(1, color=DARK, width=1.2, dash="5,4")
for Fy, col in [(F1, BLUE), (F2, ORANGE), (F3, GREEN)]:
    m = Fy <= 3; ax.plot(gg[m], Fy[m], color=col, width=2.5)
for g0, col in [(0.2375, BLUE), (0.2741, ORANGE), (0.1274, GREEN)]: ax.vline(g0, color=col, width=1, dash="3,3"); ax.dot(g0, 1, 5, col)
f.lines(770, 120, "spin ½ only:  γ = ln2/(π√3) = 0.127\n(Ashtekar–Baez–Corichi–Krasnov 1998)", size=13, color=GREEN)
f.lines(770, 200, "all spins, weight 2:  γ = 0.2375\n(Domagala–Lewandowski; Meissner 2004)", size=13, color=BLUE)
f.lines(770, 280, "all spins, weight 2j + 1:  γ = 0.274\n(Ghosh–Mitra)", size=13, color=ORANGE)
f.lines(770, 360, "Each makes the entropy A/4 in Planck units.\nThe value depends on the counting rule:\nγ is fixed by counting, not by a symmetry.", size=13)
f.text(90, 520, "Checks I5–I6 solve each equation; the roots reproduce the published values.", size=12, color=GRAY)
f.save(P("g5-counting-fixes-gamma.svg"))

# g6: three dials with stable endpoints
f = Fig(1180, 520); f.title("Three genuine dials containing G: each starts from a stable theory, the form of the change is forced, its size is free", size=17)
rows = [("free gravitons on flat space", "√G (self-coupling)", "general relativity (+Λ), unique", "Gupta 1954, Deser 1970; Boulanger–Damour–Gualtieri–Henneaux 2001"),
        ("topological BF theory (no local gravitons)", "α = G ħ Λ / c³", "general relativity + Λ", "MacDowell–Mansouri 1977; Freidel–Starodubtsev 2005"),
        ("Poincaré + Heisenberg (can still deform)", "ℓ and L", "the rigid two-scale algebra", "Mendes 1994; Chryssomalakos–Okon 2004")]
for i, (a, d, b, src) in enumerate(rows):
    y = 90 + i*130
    f.rect(40, y, 330, 70, color=BLUE, fill=BLUE, rx=10, opacity=0.12); f.text(205, y + 41, a, size=13, anchor="middle")
    f.line(380, y + 35, 760, y + 35, color=DARK, width=2.5, arrow=True); f.text(570, y + 25, "turn " + d, size=14, anchor="middle", weight="bold", color=ORANGE)
    f.text(570, y + 58, "form forced, size free", size=12, anchor="middle", color=GRAY)
    f.rect(770, y, 370, 70, color=GREEN, fill=GREEN, rx=10, opacity=0.12); f.text(955, y + 41, b, size=13, anchor="middle")
    f.text(40, y + 92, src, size=11, color=GRAY)
f.text(40, 490, "This is the same standing as ħ: deformation quantization's form is unique up to equivalence (Moyal), its size is not. None of the three fixes the value of G.", size=13)
f.save(P("g6-three-dials-with-stable-ends.svg"))

# g7: horizon as wall
f = Fig(1180, 520); f.title("Where G sits in our dial picture: Schwarzschild's g_tt = 1 + kx with k = −2G/c², and the horizon is the wall", size=17)
xs = np.linspace(0, 1.5, 400); gtt = 1 - xs
ax = Axes(f, 90, 70, 600, 370, (0, 1.5), (-0.6, 1.1))
ax.frame(xticks=[(0, "0"), (0.5, "½"), (1, "1: horizon"), (1.5, "")], yticks=[-0.5, 0, 0.5, 1], xlabel="2GM/(rc²)  (mass per radius, in the dial's units)", ylabel="g_tt (the lift)")
ax.hline(0, color=DARK, width=1); ax.plot(xs[xs <= 1], gtt[xs <= 1], color=BLUE, width=3); ax.plot(xs[xs >= 1], gtt[xs >= 1], color=RED, width=3)
ax.vline(1, color=DARK, width=1.5, dash="5,4"); ax.dot(1, 0, 6, DARK)
ax.text(0.3, 0.9, "outside: time is time", size=13, color=BLUE); ax.text(1.03, -0.35, "inside: time and", size=13, color=RED); ax.text(1.03, -0.45, "radius trade roles", size=13, color=RED)
ax.text(0.03, 0.15, "slope at the wall → surface gravity c⁴/4GM", size=12); ax.text(0.03, 0.06, "→ with the 2π period, Hawking's T", size=12)
f.lines(730, 120, "Reading (ours): G/c² is the dial that turns mass\ninto length, and the horizon is where the lift\n1 + kx reaches its wall −1/k.\n\nPast the wall the lift goes negative, the same\nflip as mult-add past its wall (why-i, part 2).\n\nLimit: masses do not compose by mult-add\n(dial-map G1–G2), and redshift factors that do\ncompose that way have dial 1/c², not G.\nSo this places G; it does not dial it.", size=13)
f.text(90, 495, "Checks W1–W3.", size=12, color=GRAY)
f.save(P("g7-horizon-is-a-wall.svg"))

# g8: asymptotic safety toy flow
f = Fig(1180, 500); f.title("Asymptotic safety: G runs between two fixed points (toy one-loop shape)", size=18)
gs = 0.7; t = np.linspace(-6, 6, 600); gt = gs/(1 + np.exp(-2*t))
ax = Axes(f, 90, 70, 560, 340, (-6, 6), (0, 0.85))
ax.frame(xticks=[(-5, "low energy"), (0, ""), (5, "high energy")], yticks=[0, 0.35, 0.7], xlabel="log of the energy scale k", ylabel="g = G k² (dimensionless)")
ax.plot(t, gt, color=PURPLE, width=3); ax.hline(gs, color=GREEN, width=1.2, dash="5,4"); ax.hline(0, color=BLUE, width=1.2)
ax.text(-0.2, 0.75, "ultraviolet fixed point g* (scheme-dependent)", size=12, color=GREEN); ax.text(-5.8, 0.09, "Gaussian point: ordinary gravity", size=12, color=BLUE)
f.lines(690, 120, "A genuine dial, but the knob is the energy\nscale, not a contraction like c → ∞.\n\nWeak fit (ours): in the variable 1/g the flow\nmaps are stretch-and-shift maps x → λx + s,\nthe order dial's affine group, with fixed point\n1/g*. Any flow of this one-loop shape does that.", size=13)
f.text(90, 470, "Toy β(g) = 2g − 2g²/g* with g* = 0.7, shape only. Weinberg 1979; Reuter 1998. Checks A1–A2.", size=12, color=GRAY)
f.save(P("g8-running-g.svg"))
print("ok")
