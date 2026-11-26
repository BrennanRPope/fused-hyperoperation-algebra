"""Figures for newton-g/README.md. Run: python3 figures.py (numpy + PIL via svgkit)."""
import numpy as np
from svgkit import *
from gshapes import *

# ---------- n1: what matched, read as induced share ----------
f = Fig(1000, 460); f.title("What matched at one loop: the share of gravity's stiffness 1/G that is induced")
a = Axes(f, 80, 70, 400, 300, (-3, 3), (0, 1))
t = np.linspace(-3, 3, 400); R = R_logistic(t)
a.plot(np.r_[t, t[::-1]], np.r_[R, np.zeros_like(R)], color=BLUE, fill="#d6e4fb", width=0.5, close=True)
a.plot(np.r_[t, t[::-1]], np.r_[R, np.ones_like(R)], color=GRAY, fill="#eeeeee", width=0.5, close=True)
a.plot(t, R, color=BLUE, width=2.5)
a.frame(xticks=[-3, -2, -1, 0, 1, 2, 3], yticks=[0, 0.5, 1], xlabel="RG time t = ln k (energy grows to the right)", ylabel="share  R = g/g*")
a.text(-2.8, 0.82, "bare part 1/G0", size=13, color="#555"); a.text(0.6, 0.28, "induced part w k^2", size=13, color=BLUE)
a.text(-2.9, 0.06, "IR: all bare", size=12, color=DARK); a.text(0.9, 0.12, "UV end: all induced", size=12, color=DARK)
f.lines(80, 418, "One loop = Sakharov: 1/G = 1/G0 + w k^2.  R/(1-R) = G0 w k^2,\nso the odds grow by e^2 per e-fold (a union / odds rule, ours).", size=12)
b = Axes(f, 580, 70, 360, 300, (0, 1), (1.9, 2.15))
Rv = np.linspace(0.001, 0.999, 300)
b.plot(Rv, 2 + 0*Rv, color=GRAY, width=2.5, dash="6,4")
b.plot(Rv, 1/(0.5*(1 - Rv) + Rv/theta_litim), color=RED, width=2.5)
b.frame(xticks=[0, 0.5, 1], yticks=[1.9, 2.0, 23/11, 2.15], fmt="{:.3f}", xlabel="share R = g/g*", ylabel="growth rate of log-odds")
b.text(0.3, 1.975, "one loop: always 2 (G's units)", size=12, color="#555")
b.text(0.35, 2.10, "Litim, beyond one loop: 2 -> 23/11", size=12, color=RED)
f.lines(580, 418, "Beyond one loop the graviton's own loops make the rate drift:\nthe odds rule breaks by 1/22 at the fixed-point end.", size=12)
f.save("figures/n1-what-matched.svg")

# ---------- n2: branching ladder ----------
f = Fig(1000, 400); f.title("One scheme-free number, rho, places every 'double fixed point' flow on one ladder")
a = Axes(f, 70, 150, 860, 60, (-0.05, 1.3), (0, 1))
a.f.line(a.X(-0.05), a.Y(0.5), a.X(1.3), a.Y(0.5), color="#999", width=2)
for v in [0, 0.25, 0.5, 0.75, 1, 1.25]:
    f.line(a.X(v), a.Y(0.5) - 6, a.X(v), a.Y(0.5) + 6, color="#999"); f.text(a.X(v), a.Y(0.5) + 24, f"{v:g}", size=12, anchor="middle", color="#555")
pts = [(0, "circle op = geometric offspring", "one loop exactly (Mobius)", GREEN, -1),
       (1/3, "E(x) = e^x - 1 = Poisson(1) branching", "our height flow", BLUE, 1),
       (0.6, "2D gravity, c = 0 (t = ln k)", "Aida-Kitazawa two loops", RED, -1),
       (1.0, "binary splitting", "0 or 2 children", PURPLE, 1),
       (1.2, "2D gravity if t = ln k^2", "same physics, other clock", AMBER, -1)]
for v, l1, l2, c, side in pts:
    a.dot(v, 0.5, r=7, color=c)
    y = a.Y(0.5) + (-60 if side < 0 else 62)
    f.line(a.X(v), a.Y(0.5) + (-8 if side < 0 else 8), a.X(v), y + (14 if side < 0 else -16), color=c, width=1)
    f.text(a.X(v), y, l1, size=12, anchor="middle", color=c, weight="bold"); f.text(a.X(v), y + 15, l2, size=11, anchor="middle", color="#555")
f.text(500, 312, "rho = -b1/b0^2 of the flow; for critical branching rho = 1 - 2 f'''(1) / (3 sigma^4)  (ours to state; parts known)", size=12, anchor="middle")
f.text(500, 334, "Exponent at the other fixed point: theta = eps + rho eps^2, with eps the rate at the free end (d - 2 for gravity)", size=12, anchor="middle")
f.text(500, 360, "Same shape everywhere; the SIZE of rho only compares once you fix the clock (one height = ? e-folds of energy)", size=12, anchor="middle", color=RED)
f.save("figures/n2-branching-ladder.svg")

# ---------- n3: numerology trap ----------
f = Fig(1000, 440); f.title("Can E's rho = 1/3 be made to equal gravity's? Only by choosing a free number")
a = Axes(f, 90, 60, 560, 300, (-40, 20), (0, 1.2))
c = np.linspace(-40, 12.4, 400); a.plot(c, 15/(25 - c), color=RED, width=2.5)
a.hline(1/3, color=BLUE, dash="6,4"); a.hline(2/3, color=BLUE, dash="2,4")
a.dot(-20, 1/3, color=BLUE); a.dot(2.5, 2/3, color=BLUE); a.dot(0, 0.6, color=RED)
a.frame(xticks=[(v, str(v)) for v in [-40, -30, -20, -10, 0, 10, 20]], yticks=[0, 1/3, 0.6, 2/3, 1], fmt="{:.2f}", xlabel="matter central charge c", ylabel="rho")
a.text(-39, 0.36, "E, one height = 1 e-fold (r = e)", size=12, color=BLUE); a.text(-39, 0.70, "E, one height = 2 e-folds (r = e^2)", size=12, color=BLUE)
a.text(-19, 0.25, "needs c = -20", size=12, color=BLUE); a.text(3.5, 0.74, "needs c = 2.5", size=12, color=BLUE)
a.text(1, 0.53, "pure gravity: 3/5", size=12, color=RED); a.text(-12, 0.95, "gravity: rho = 15/(25 - c)", size=13, color=RED)
f.lines(690, 90, "Free numbers available:\n- the clock r (energy ratio per height)\n- matter content c\n- for branching: the offspring law\n\nOne target number (3/5).\nEach knob alone hits it.\nNone lands on a natural value.\n\nVerdict: no test, so\nno evidence either way.", size=13)
f.save("figures/n3-numerology.svg")

# ---------- n4: shape ----------
p = theta_litim/2; tt = np.linspace(-4, 4, 801)
RL = align(R_litim, tt); RR = align(lambda t: R_richards(t, p), tt); RG = align(R_logistic, tt)
f = Fig(1000, 440); f.title("Same shape, different mechanics: G's flow beyond one loop vs two textbook curves")
a = Axes(f, 80, 60, 420, 300, (-4, 4), (0, 1))
a.plot(tt, RG, color=GRAY, width=2, dash="6,4"); a.plot(tt, RR, color=ORANGE, width=2, dash="3,3"); a.plot(tt, RL, color=BLUE, width=2.5)
a.frame(xticks=[-4, -2, 0, 2, 4], yticks=[0, 0.5, 1], xlabel="RG time t (aligned at R = 1/2)", ylabel="R = g/g*")
a.text(-3.8, 0.9, "blue: Litim cutoff, exact", size=12, color=BLUE); a.text(-3.8, 0.83, "orange: Richards, p = 23/22", size=12, color=ORANGE); a.text(-3.8, 0.76, "gray: logistic (one loop)", size=12, color=GRAY)
b = Axes(f, 590, 60, 360, 300, (-4, 4), (-0.8, 0.8))
b.hline(0, color="#bbb"); b.plot(tt, 100*(RG - RL), color=GRAY, width=2, dash="6,4"); b.plot(tt, 100*(RR - RL), color=ORANGE, width=2.5)
b.frame(xticks=[-4, -2, 0, 2, 4], yticks=[-0.8, -0.4, 0, 0.4, 0.8], xlabel="RG time t", ylabel="difference from Litim, %")
f.lines(80, 410, "Richards matches both end rates (2 and 23/11) and stays within 0.2%; logistic within 0.7%.\nBut Richards needs g^(45/22) in its beta: not allowed for perturbative gravity. Shape yes, mechanism no.", size=12)
f.save("figures/n4-shape.svg")

# ---------- n5: fluid strand ----------
f = Fig(1000, 500); f.title("Spacetime as a fluid: G sets the units of every transport number, so no ratio can fix it")
a = Axes(f, 70, 70, 330, 300, (0, 4), (0, 10))
cols = [BLUE, GREEN, PURPLE, ORANGE, AMBER]; tops = 0
for i, (h, lab) in enumerate([(1.0, "1/G0 (bare)"), (2.5, "field species 1"), (3.0, "field species 2"), (1.5, "field species 3"), (1.2, "...")]):
    f.rect(a.X(1), a.Y(tops + h), a.X(2) - a.X(1), a.Y(tops) - a.Y(tops + h), color="white", fill=cols[i], opacity=0.85)
    f.text(a.X(2.15), a.Y(tops + h/2) + 4, lab, size=12, color=cols[i]); tops += h
a.frame(xticks=[], yticks=[], ylabel="stiffness 1/G  (each ~ kappa^2)", box=False)
f.lines(40, 440, "Sakharov: species add in 1/G (a supertrace over fields,\nfermions with a minus sign; Visser eq. 21). Each piece\nscales with the cutoff kappa^2, the 'molecule' size: free.", size=12)
b = Axes(f, 560, 70, 380, 300, (-2, 2), (-4, 2.5))
lg = np.linspace(-2, 2, 50)
b.plot(lg, -np.log10(16*np.pi) - lg, color=BLUE, width=2.5); b.plot(lg, -np.log10(4) - lg, color=GREEN, width=2.5)
b.plot(lg, np.log10(1/(4*np.pi)) + 0*lg, color=RED, width=2.5, dash="6,4")
b.frame(xticks=[-2, -1, 0, 1, 2], yticks=[-4, -3, -2, -1, 0, 1, 2], xlabel="log10 G (Planck units)", ylabel="log10 of quantity")
b.text(-1.9, 1.9, "entropy density s = 1/(4G)", size=12, color=GREEN); b.text(-1.9, -0.3, "viscosity eta = 1/(16 pi G)", size=12, color=BLUE)
b.text(0.1, -0.85, "eta/s = 1/(4 pi): flat", size=12, color=RED)
f.lines(560, 440, "Membrane paradigm / Jacobson: G is the inverse 'density of\nhorizon bits'. Every dimensionless ratio (1/4, 1/4 pi) is G-free.", size=12)
f.save("figures/n5-fluid.svg")

# ---------- n6: verdict ----------
rows = [("One-loop match = circle op", "exact identity; true for any quadratic flow", GREEN, "recovered"),
        ("...is Sakharov induced gravity", "g/g* = induced share of 1/G; odds x e^2 per e-fold", GREEN, "ours (reading)"),
        ("Breaks past one loop", "4.5% (Litim), order eps (2+eps), spiral (full EH)", AMBER, "recovered"),
        ("Branching ladder rho", "circle op 0, E 1/3, 2D gravity 3/5, binary 1", GREEN, "ours to state"),
        ("Does E's 1/3 match gravity?", "only with a free clock or c = -20: numerology", RED, "no test"),
        ("Richards / power-mean dial", "right shape (0.2%), wrong mechanism (non-analytic)", AMBER, "ours"),
        ("Fluid strand fixes G's value?", "no: eta/s and S/A are G-free; G = molecule size", RED, "known"),
        ("New testable prediction", "none found", RED, "")]
f = Fig(1000, 80 + 46*len(rows)); f.title("Newton's G lead, reopened: verdict board")
for i, (a1, a2, c, st) in enumerate(rows):
    y = 60 + 46*i
    f.rect(30, y, 940, 38, color="#ddd", fill="#fafafa", rx=6); f.circle(52, y + 19, 8, color=c, fill=c)
    f.text(72, y + 24, a1, size=14, weight="bold"); f.text(420, y + 24, a2, size=13, color="#333"); f.text(960, y + 24, st, size=12, anchor="end", color="#777", italic=True)
f.save("figures/n6-verdict.svg")
print("figures written")
