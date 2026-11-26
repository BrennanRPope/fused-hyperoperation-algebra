"""Figures for section 7 of newton-g/README.md (second look, Nov 26, 2026). Run: python3 figures2.py"""
import json, math
import numpy as np
from svgkit import *

R2 = json.load(open("results2.json"))
rho = lambda q: 1/3 - 4*q/3          # q = d/k^2

# ---------- n7: rho across the plane law's curvature dial ----------
f = Fig(1000, 500); f.title("The plane law's rung step: its second-loop number rho = 1/3 - 4d/(3k^2) runs over every value")
a = Axes(f, 90, 60, 560, 340, (-0.6, 0.6), (-0.5, 1.2))
f.rect(a.X(-0.6), a.y, a.X(0.25) - a.X(-0.6), a.h, color="none", fill="#fdebd6", opacity=0.6)
f.rect(a.X(0.25), a.y, a.X(0.6) - a.X(0.25), a.h, color="none", fill="#dbe9fb", opacity=0.6)
q = np.linspace(-0.6, 0.6, 200); a.plot(q, rho(q), color=DARK, width=2.5)
a.hline(0, color="#999", dash="3,3"); a.vline(0.25, color=PURPLE, width=2)
for qq, rr, c, lab, dx, dy in [(0.25, 0, PURPLE, "parabola: Moebius, rho = 0", -0.27, -0.17),
                                (0, 1/3, BLUE, "d = 0: e^x - 1, rho = 1/3", 0.02, 0.05),
                                (-0.2, 0.6, RED, "gravity's 3/5 needs d = -k^2/5", 0.02, 0.04),
                                (-0.5, 1.0, PURPLE, "rho = 1 (binary splitting's value)", 0.02, 0.06)]:
    a.dot(qq, rr, r=6, color=c); a.text(qq + dx, rr + dy, lab, size=12, color=c, weight="bold")
a.frame(xticks=[-0.6, -0.4, -0.2, 0, 0.25, 0.4, 0.6], yticks=[-0.5, -1/3, 0, 1/3, 0.6, 1, 1.2], fmt="{:.2f}",
        xlabel="curvature dial d / k^2  (only the ratio matters)", ylabel="rho")
a.text(-0.58, -0.42, "hyperbola side (two walls): rho > 0", size=12, color=AMBER)
a.text(0.27, 1.1, "ellipse side", size=12, color=BLUE); a.text(0.27, 1.03, "(spheres): rho < 0", size=12, color=BLUE)
f.lines(690, 80, "What the picture says\n\n- rho is a reading of the base dials:\n  rho = 1/3 - 4d/(3k^2)  (= discriminant/(3k^2))\n- its sign is the conic's type\n- sliding the shape dial k with d = 0\n  keeps rho = 1/3 (a rescaling)\n- 3/5 is reachable, but only by\n  choosing d = -k^2/5. No named point\n  of the plane sits there.\n- every sphere-like (ellipse) point\n  has the wrong sign for gravity\n\nChecks P1-P8.", size=13)
f.save("figures/n7-curvature-dial-rho.svg")

# ---------- n8: the clock is the base dial r ----------
f = Fig(1000, 460); f.title("The clock: E's rho in gravity's units is (ln r)/3, where r is the energy ratio per height step")
a = Axes(f, 90, 60, 560, 300, (1, 10), (0, 0.85))
r = np.linspace(1.01, 10, 300); a.plot(r, np.log(r)/3, color=BLUE, width=2.5)
a.hline(0.6, color=RED, dash="6,4"); a.hline(0, color=PURPLE, width=2)
for name, rv in [("2", 2), ("e", math.e), ("6", 6), ("e^2", math.e**2), ("e^(9/5)", math.exp(1.8))]:
    v = R2["rho_E_by_clock"][name]; a.dot(rv, v, r=5, color=BLUE if name != "e^(9/5)" else RED)
    dx, dy = {"6": (-1.35, 0.04), "e^(9/5)": (0.15, -0.06)}.get(name, (0.12, -0.05))
    a.text(rv + dx, v + dy, f"r = {name}: {v:.3f}", size=11, color=DARK)
a.text(1.2, 0.63, "gravity, Aida-Kitazawa, t = ln k: 3/5", size=12, color=RED)
a.text(1.2, 0.03, "Moebius step: rho = 0 at every clock", size=12, color=PURPLE)
a.frame(xticks=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], yticks=[0, 0.2, 1/3, 0.6, 0.8], fmt="{:.2f}",
        xlabel="similarity ratio r (energy ratio per height step; a base dial)", ylabel="rho per e-fold")
f.lines(690, 80, "One height step = ln r e-folds.\nMatch needs r = e^(9/5) = 6.05:\none height = 9/5 e-folds.\n\nr = 6 (an integer blocking\nratio) lands 0.45% short.\nThat is one free number hitting\none target: numerology.\n\nThe design rule: r is an input.\nNo reading (discriminant,\ncurvature, holonomy) can set it.\n\nChecks Q1-Q4.", size=13)
f.save("figures/n8-clock-r.svg")

# ---------- n9: design-rule audit board ----------
f = Fig(1000, 470); f.title("Design-rule audit of the first pass: base dials in, readings out")
rows = [("What the first pass wrote", "Base dial (input)", "Reading (output)", "Fix"),
        ("k = -g*", "w, dimension n: k = -(n-2)/w", "g* = 12pi/23, the wall", "reworded (A1)"),
        ("R = g/g* 'probability'", "odds t = G0 w k^2", "share R = t/(1+t)", "reworded (A3)"),
        ("Richards p = 23/22", "l^p exponent p", "Litim exponent theta = 23/11", "p set from a reading: a fit (A2)"),
        ("rho compared to 3/5", "k, d, clock r, matter c", "rho = 1/3 - 4d/(3k^2)", "rho is an output (A4)"),
        ("'one height = 5/9 e-fold'", "clock r", "rho per e-fold = (ln r)/3", "inverted; correct: 9/5 (Q2)")]
x0, y0, cw, rh = 30, 60, [230, 250, 250, 210], 58
for i, row in enumerate(rows):
    x = x0
    for j, cell in enumerate(row):
        fill = "#eeeeee" if i == 0 else ("#dbe9fb" if j == 1 else "#fdebd6" if j == 2 else "white")
        f.rect(x, y0 + i*rh, cw[j], rh, color="#999", fill=fill)
        f.text(x + 8, y0 + i*rh + 34, cell, size=13 if i else 14, weight="bold" if i == 0 else "normal",
               color=RED if (j == 3 and i in (3, 5)) else DARK)
        x += cw[j]
f.lines(30, 425, "Blue: what you choose. Orange: what you read off. Two first-pass slips: a fit was presented beside a reading,\nand the clock factor was written upside down. Neither changes the verdict.", size=13)
f.save("figures/n9-design-audit.svg")

# ---------- n10: fingerprint plane ----------
f = Fig(1000, 470); f.title("Two fingerprints: a clock change slides the residue, but cannot move E onto the flow line")
a = Axes(f, 90, 60, 560, 320, (-0.1, 1.15), (-0.5, 8))
pts = [("Moebius (circle op, one loop)", 0, 0, PURPLE), ("negative binomial r = 2", 1/9, 0.768, GRAY), ("E = Poisson(1)", 1/3, 64.153, BLUE),
       ("binomial(3, 1/3)", 2/3, 37646, GRAY), ("binary splitting", 1, 2.235e7, GRAY)]
a.hline(0, color=PURPLE, width=3)
for name, rr, h, c in pts:
    y = math.log10(1 + h); a.dot(rr, y, r=6, color=c)
    dx, dy = {"Moebius (circle op, one loop)": (-0.08, -0.45), "binary splitting": (-0.32, 0.1), "negative binomial r = 2": (0.03, 0.3)}.get(name, (0.02, 0.25))
    a.text(rr + dx, y + dy, name, size=12, color=c, weight="bold")
a.dot(0.6, 0, r=7, color=RED); a.text(0.50, -0.42, "gravity's 3/5 if its flow is analytic", size=12, color=RED)
yE = math.log10(65.153)
a.plot([0.08, 1.1], [yE, yE], color=BLUE, width=1.5, dash="5,4")
a.arrow(0.36, yE + 0.6, 0.62, yE + 0.6, color=BLUE)
a.text(0.40, yE + 0.85, "change of clock: residue x ln r", size=12, color=BLUE)
a.text(0.5, yE - 0.45, "horn stays non-zero (not to scale)", size=11, color=BLUE)
a.frame(xticks=[0, 1/3, 0.6, 1], yticks=[(0, "0"), (2, "100"), (4, "10^4"), (6, "10^6"), (8, "10^8")], fmt="{:.2f}",
        xlabel="iterative residue (= rho)", ylabel="horn-map size (log scale)")
a.text(0.75, 0.25, "purple line: maps inside a continuous flow", size=12, color=PURPLE)
f.lines(690, 80, "Horn sizes from ../parabolic-fingerprints.\n\nA continuous flow has horn size 0.\nE has 64.153 and no flow\n(../fractional-steps).\nIts iterates share one Fatou\ncoordinate, so they stay off\nthe line (standard, Ecalle).\n\nIf gravity's beta were an analytic\nflow, no clock could make E match\nit, even when the residues agree.\nCaveat: perturbative beta series\ndiverge, so 'analytic' is not proven.", size=13)
f.save("figures/n10-fingerprint-plane.svg")
print("ok")
