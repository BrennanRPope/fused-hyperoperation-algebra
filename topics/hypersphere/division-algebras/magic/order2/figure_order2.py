"""DRAFT figure D6: the pair rule at the next order. Run after checks_order2.py."""
import os, sys, json; sys.path.insert(0, '..')
import numpy as np, sympy as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from poly import poly_F
from precess import u, L, E
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "..", "..", "figures")
KNOWN, OURS, GREY, GREEN, PURPLE = "#d9822b", "#2b6cb0", "#bbbbbb", "#2f855a", "#805ad5"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
one = sp.Integer(1)
ws, ss = [one, 0, 0, 0], [-1, 1, 1, 1]          # visitor: one opposite charge
bgs = [("one charge", [one, 0, 0, 0], "k"), ("two charges", [one, one, 0, 0], GREEN),
       ("three charges", [one, one, one, 0], OURS), ("four charges", [one]*4, PURPLE)]
def orbit(ks, laps=10):
    P = sp.expand(poly_F(ks, ws, ss) + L**2*u**2); cs = [P.coeff(u, m) for m in range(5)]
    K = cs[1].subs(E, 1)/2; Lv = 3*K; Ev = sp.sqrt(1 - sp.Rational(3, 5)*K**2/Lv**2)
    c = [float(sp.N(x.subs(E, Ev))) for x in cs]; Lf = float(Lv); c[2] -= Lf**2
    F = np.polynomial.Polynomial(c); dF = F.deriv()
    r = F.roots(); r = np.sort(r[(abs(r.imag) < 1e-12) & (r.real > 0)].real); u2 = r[1]
    sol = solve_ivp(lambda p, y: [y[1], dF(y[0])/(2*Lf**2)], [0, 2*np.pi*laps], [u2, 0.0], rtol=1e-11, atol=1e-13, max_step=0.01)
    ph, uu = sol.t, sol.y[0]; return np.cos(ph)/uu, np.sin(ph)/uu
fig = plt.figure(figsize=(22, 12.5)); gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], hspace=0.38, wspace=0.22)
notes = ["closes exactly,\nat every order", "an exact precessing ellipse:\nturn known in closed form", "turns more:\ntriples join at second order", "turns most:\nall four join at third order"]
for i, (name, ks, col) in enumerate(bgs):
    ax = fig.add_subplot(gs[0, i]); x, y = orbit(ks); s = np.max(np.hypot(x, y))
    ax.plot(x/s, y/s, color=col, lw=1.2); ax.plot(0, 0, "k*", ms=12)
    ax.set_aspect("equal"); ax.set_xlim(-1.12, 1.12); ax.set_ylim(-1.12, 1.12); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{'ABCD'[i]}. Black hole with {name}", fontsize=12.5); ax.set_xlabel(notes[i], fontsize=11)
fig.text(0.5, 0.965, "Ten laps of the same kind of visitor in strong gravity (exact orbits, not approximations)", ha="center", fontsize=14)
# E. second-order turn split into pairs compounded and triples
ax = fig.add_subplot(gs[1, 0:2])
rows = {r["bg"]: r for r in json.load(open(os.path.join(HERE, "results_order2.json")))["data"]["order2_cases"] if r["probe"] == "anti1"}
labs, comp, trip = [], [], []
for key, lab in [("2", "two charges"), ("3", "three charges"), ("4", "four charges")]:
    r = rows[key]; T = sp.Rational(r["T"]); K = sp.Rational(r["K"]); p3 = sp.Rational(r["triples"])
    labs.append(lab); comp.append(float(sp.Rational(3, 4)*T**2)); trip.append(float(3*p3/K**3))
xs = np.arange(3)
ax.bar(xs, comp, 0.55, color=GREEN, label="pairs, compounded:  ¾ × (first-order turn)²")
ax.bar(xs, trip, 0.55, bottom=comp, color=OURS, label="new at this order: triples of charges")
for xi, a, b in zip(xs, comp, trip):
    ax.text(xi, a + b + 0.25, f"{sp.nsimplify(a)} + {sp.nsimplify(b)} = {sp.nsimplify(a + b)}", ha="center", fontsize=11)
ax.set_xticks(xs, labs); ax.set_ylabel("second-order turn coefficient"); ax.set_ylim(0, 13.5); ax.legend(loc="upper left", fontsize=10.5)
ax.set_title("E. At second order the pairs come back squared, and triples of charges add new turning")
# F. which charges act at which order
ax = fig.add_subplot(gs[1, 2:4]); ax.axis("off")
ax.set_title("F. The rule, order by order (a small visitor around a four-charge black hole)")
tab = [["order", "what sets the turn", "status"],
       ["1st", "pairs of charges", "ours (section 9); one charge's zero: Caron-Huot & Zahraee"],
       ["2nd", "pairs (squared) + triples", "ours: exact formula, 30 cases checked"],
       ["3rd", "+ the product of all four", "ours: read off the exact equation, checked numerically"],
       ["all", "one charge: never turns\ntwo charges: closed-form rosette", "ours: one exact orbit equation for every case"]]
t = ax.table(cellText=tab[1:], colLabels=tab[0], loc="center", cellLoc="left", colWidths=[0.1, 0.38, 0.52])
t.auto_set_font_size(False); t.set_fontsize(11.5); t.scale(1, 3.1)
for (rr, cc), cell in t.get_celld().items():
    if rr == 0: cell.set_facecolor("#eeeeee")
    if cc == 2 and rr >= 1: cell.get_text().set_color(OURS)
fig.savefig(os.path.join(OUT, "figD6-the-pair-rule-at-the-next-order.png"), dpi=110, bbox_inches="tight"); print("saved")
