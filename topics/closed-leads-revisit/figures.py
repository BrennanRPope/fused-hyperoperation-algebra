"""Pictures for closed-leads-revisit/README.md (DRAFT).

Every curve here comes from a known formula; nothing is fitted. Run: python3 figures.py
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
REOPEN, CLOSED, CORRECT, ANSWERED = "#2b8a3e", "#adb5bd", "#e67700", "#1c7ed6"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=140)
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- fig 1: overview
cats = [  # counts are the rows in README.md (some rows bundle several pages)
    ("Reopened (old closures covered by R1-R6, S1-S2)", 21, REOPEN),
    ("Answered by later work", 6, ANSWERED),
    ("Correct the record (verdict stands)", 3, CORRECT),
    ("Stays closed: known math restated", 39, CLOSED),
    ("Stays closed: failed a test the lens doesn't touch", 41, CLOSED),
    ("Stays closed: structural no-go still holds", 34, CLOSED),
]
fig, ax = plt.subplots(figsize=(10, 4.4))
y = np.arange(len(cats))[::-1]
ax.barh(y, [c[1] for c in cats], color=[c[2] for c in cats])
for yi, c in zip(y, cats):
    ax.text(c[1] + 0.6, yi, str(c[1]), va="center", fontsize=9)
ax.set_yticks(y, [c[0] for c in cats], fontsize=9)
ax.set_xlabel("number of closed leads (about 140 in all)")
ax.set_title("Fig 1. Every closed lead, re-judged with the fuller theory")
save(fig, "fig1-overview.png")

# ---------------------------------------------------------------- fig 2: negative numbers at fractional steps
x = sp.symbols("x")
N = 9
a = sp.symbols("a2:%d" % (N + 1))
h = x + sum(a[i] * x ** (i + 2) for i in range(N - 1))
target = sp.series(sp.exp(x) - 1, x, 0, N + 1).removeO()
comp = sp.series(h.subs(x, h), x, 0, N + 1).removeO()
sol = {}
for n in range(2, N + 1):
    eq = sp.expand(comp.subs(sol)).coeff(x, n) - target.coeff(x, n)
    s = sp.solve(eq, a[n - 2])[0]
    sol[a[n - 2]] = s
half = [1] + [sol[a[n - 2]] for n in range(2, N + 1)]

fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
ax = axs[0]
vals = [float(c) for c in half]
cols = [REOPEN if v >= 0 else "#c92a2a" for v in vals]
ax.bar(range(1, N + 1), vals, color=cols)
ax.set_yscale("symlog", linthresh=1e-6)
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("power of x")
ax.set_title("(a) Half rung of e^x - 1, counted around\nthe fixed point (old combinatorics test)", fontsize=9)
ax.text(6, -3e-4, "negative", fontsize=8, color="#c92a2a", ha="center")

ax = axs[1]
n_frac, c = 2.5, 0.5
z = sp.symbols("z")
pgf = sp.series((1 + c * (z - 1)) ** sp.Rational(5, 2), z, 0, 9).removeO()
pv = [float(pgf.coeff(z, k)) for k in range(9)]
ax.bar(range(9), pv, color=[REOPEN if v >= 0 else "#c92a2a" for v in pv])
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("number of claims N")
ax.set_title("(b) Actuarial: 2.5 copies of a coin flip\n(1 + c(z-1))^2.5, c = 1/2: P(N) dips below 0", fontsize=9)

ax = axs[2]
ax.axis("off")
ax.set_title("(c) One theorem family, four old closures", fontsize=9)
rows = [("whole steps", "fractional or backward steps"),
        ("counts >= 0", "half rung: negative coefficient"),
        ("P(N) >= 0", "2.5 coin flips: P(4) < 0"),
        ("Kok clock forward", "Kok clock backward: populations < 0"),
        ("Kac chain forward", "Kac chain inverse: negative entries")]
for i, (l, r) in enumerate(rows):
    w = "bold" if i == 0 else "normal"
    ax.text(0.0, 0.95 - 0.12 * i, l, fontsize=9, weight=w)
    ax.text(0.42, 0.95 - 0.12 * i, r, fontsize=9, weight=w, color="#c92a2a" if i else "k")
ax.text(0.0, 0.02, "Known names: infinite divisibility (distributions),\nKingman 1962 (Markov chains),"
        " Karlin-McGregor 1968\n(branching processes). The height flow is a branching\nprocess, so this family is native to it.",
        fontsize=8.5)
fig.suptitle("Fig 2. Reopen: the 'negative numbers at fractional steps' wall is the embeddability question", fontsize=11)
save(fig, "fig2-embeddability.png")

# ---------------------------------------------------------------- fig 3: parabolic fixed points
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
E = lambda t: math.expm1(t)
xs = np.linspace(-1.1, 0.15, 300)
ax.plot(xs, np.expm1(xs), color=REOPEN, label="E(x) = e^x - 1 (one rung)")
ax.plot(xs, xs, "k--", lw=0.8, label="y = x")
p = -1.0
for _ in range(12):
    q = E(p)
    ax.plot([p, p, q], [p, q, q], color="#495057", lw=0.7)
    p = q
ax.set_title("(a) The height flow creeps into a parabolic fixed point\n(tangent to y = x at 0)", fontsize=9)
ax.legend(fontsize=8)
ax.set_xlabel("x"); ax.set_ylabel("E(x)")
ax = axs[1]
n = np.arange(1, 2001)
vals = []
p = -1.0
for _ in n:
    p = E(p)
    vals.append(-p)
ax.loglog(n, vals, color=REOPEN, label="-E^n(-1) = chance a Poisson(1) family\nline survives n generations")
ax.loglog(n, 2 / n, "k--", lw=0.8, label="Kolmogorov 2/n")
ax.set_xlabel("generation n (rung step)")
ax.set_title("(b) Same creep, read as a critical branching process", fontsize=9)
ax.legend(fontsize=8)
ax.text(1.3, 3e-3, "Other maps with this exact shape:\n- type-I intermittency (Pomeau-Manneville)\n"
        "- diamond lattice RG at Q = 0 (cold point)\n- critical branching in general", fontsize=8)
fig.suptitle("Fig 3. Reopen: carry the height-flow fingerprint (residue, horn) to other parabolic maps", fontsize=11)
save(fig, "fig3-parabolic.png")

# ---------------------------------------------------------------- fig 4: silent layers and Beraha
fig, ax = plt.subplots(figsize=(9, 4.6))
k = np.linspace(-2, 2, 800)
for m in range(2, 8):
    # Lucas U_m with d = 1: U_0 = 0, U_1 = 1, U_{m+1} = k U_m - U_{m-1}  (= Chebyshev U_{m-1}(k/2))
    u0, u1 = np.zeros_like(k), np.ones_like(k)
    for _ in range(m - 1):
        u0, u1 = u1, k * u1 - u0
    ax.plot(k, u1, label=f"layer weight U_{m}")
    zeros = [2 * math.cos(math.pi * j / m) for j in range(1, m)]
    ax.plot(zeros, [0] * len(zeros), "o", ms=4, color="k")
ax.axhline(0, color="k", lw=0.6)
ax.set_ylim(-4, 5.6)
ax.set_xlabel("dial k (plane dial d = 1)")
ax.set_ylabel("Lucas layer weight")
ax.text(-1.95, 4.3, "black dots: silent layers at k = 2cos(pi j/m)\nthere Q = k^2 = 4cos^2(pi j/m), the Beraha numbers\n(Temperley-Lieb loop value delta = k)", fontsize=8)
ax.legend(fontsize=7, ncol=3, loc="lower right")
ax.set_title("Fig 4. Reopen: the silent layers of the Lucas weights sit exactly at the Temperley-Lieb / Beraha points", fontsize=10)
save(fig, "fig4-silent-layers-beraha.png")

# ---------------------------------------------------------------- fig 5: where k moves
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
ax = axs[0]
for pi_ in (1.0, 0.6, 0.3):
    ax.plot([0, 1], [-1 / pi_] * 2, label=f"susceptible share pi = {pi_}: k = {-1/pi_:.2f}")
ax.set_ylim(-3.6, 0.2)
ax.set_xlabel("exposure (any)")
ax.set_ylabel("dial k")
ax.set_title("(a) Two-state hidden memory: k jumps but stays flat\n(our bath = the mixture cure model)", fontsize=9)
ax.legend(fontsize=7)
ax = axs[1]
d = np.linspace(0.01, 1.5, 300)
r = lambda t: 0.5 * (1 - np.exp(-4 * t) * (1 + 2 * t))
kk = (r(2 * d) - 2 * r(d)) / r(d) ** 2
ax.plot(d, kk, color=REOPEN)
ax.axhline(-2, color="k", ls="--", lw=0.8)
ax.text(1.0, -1.93, "no interference (Haldane map): k = -2", fontsize=8)
ax.set_xlabel("interval length d (Morgans)")
ax.set_title("(b) Crossover interference (chi-square model m = 1):\nk drifts from 0 (crowded) to -2 (sparse)", fontsize=9)
ax = axs[2]
w = np.logspace(-4, 4, 400)
for g in (0.25, 0.5, 0.75):
    nocc = 1 / (w + g)
    kh = w + g - w ** g * (1 + w) ** (1 - g)
    ax.semilogx(nocc, kh, label=f"g = {g}: {2*g-1:+.2f} sparse to {g:.2f} crowded")
ax.set_xlabel("occupation n (crowding)")
ax.set_title("(c) Haldane exclusion (Wu 1994):\nk drifts with crowding", fontsize=9)
ax.legend(fontsize=7)
fig.suptitle("Fig 5. Reopen: three old 'k is not constant' closures may all be forgotten hidden state", fontsize=11)
save(fig, "fig5-where-k-moves.png")

# ---------------------------------------------------------------- fig 6: Liebig on the semiring dial
fig, ax = plt.subplots(figsize=(8, 4.4))
tA = 0.5
tB = np.linspace(0.01, 1, 300)
ax.plot(tB, tA * tB, color="#495057", label="multiplicative (the product side, dial 1)")
for T, ls in ((1.0, "-"), (0.3, "-."), (0.1, ":")):
    ax.plot(tB, (tA ** (-1 / T) + tB ** (-1 / T)) ** (-T), color=REOPEN, ls=ls,
            label=f"smooth minimum at temperature {T}")
ax.plot(tB, np.minimum(tA, tB), "k", lw=2, label="Liebig minimum (tropical stop)")
ax.set_xlabel("nutrient B sufficiency (nutrient A fixed at 0.5)")
ax.set_ylabel("growth fraction")
ax.set_title("Fig 6. Reopen: Liebig is off the product dial but is a stop on the paired sum side", fontsize=10)
ax.legend(fontsize=8)
save(fig, "fig6-liebig-tropical.png")

# ---------------------------------------------------------------- fig 7: Carnot / Curzon-Ahlborn as a knob
fig, ax = plt.subplots(figsize=(7.5, 4.4))
eta = np.linspace(0, 1, 300)
for lam, lab, col in ((1.0, "Carnot limit (lambda = 1)", "k"),
                      (0.5, "Curzon-Ahlborn = lambda = 1/2 knob", REOPEN),
                      (0.25, "lambda = 1/4", "#69db7c")):
    ax.plot(eta, 1 - (1 - eta) ** lam, color=col, label=lab)
ax.set_xlabel("Carnot efficiency eta_C")
ax.set_ylabel("efficiency")
ax.set_title("Fig 7. Small reopen: Curzon-Ahlborn is the knob x -> 1-(1-x)^lambda\nthat leaves the union rule's k = -1 alone", fontsize=10)
ax.legend(fontsize=8)
save(fig, "fig7-carnot-knob.png")

print("half-rung coefficients:", half)
