import json, math
import numpy as np, sympy as sp, mpmath as mp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 150})
INK, BLUE, RED, GREEN, GRAY, AMBER = "#1f2937", "#2563eb", "#dc2626", "#16a34a", "#9ca3af", "#d97706"
R = json.load(open("results.json"))

# F1 three singular cubics
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
X, Y = np.meshgrid(np.linspace(-1.6, 2.2, 600), np.linspace(-2.2, 2.2, 600))
specs = [("Node: two walls (mult-add)", -1, 0, [0, 1]), ("Lone dot: turning walls", 0, -1, []), ("Cusp: the border (addition)", 0, 0, [0])]
for ax, (ttl, a1, a2, slopes) in zip(axs, specs):
    ax.contour(X, Y, Y**2 + a1*X*Y - X**3 - a2*X**2, [0], colors=INK, linewidths=2)
    if a2 == -1: ax.plot([0], [0], "o", color=INK, ms=6)
    for m in slopes:
        xs = np.linspace(-1.2, 1.2, 10); ax.plot(xs, m*xs, "--", color=RED, lw=1.3)
    ax.set_title(ttl, fontsize=11); ax.set_aspect("equal"); ax.set_xlim(-1.6, 2.2); ax.set_ylim(-2.2, 2.2)
    ax.axhline(0, color=GRAY, lw=.5); ax.axvline(0, color=GRAY, lw=.5)
axs[0].text(1.25, 1.3, "branch y = x\nwall at a = -1", color=RED, fontsize=9)
axs[0].text(.15, -.55, "branch y = 0\nwall at infinity", color=RED, fontsize=9)
plt.tight_layout(); plt.savefig("figures/f1_singular_cubics.png"); plt.close()

# F2 Tate curve pinching to the node
def tate_coeffs(q, M=60):
    s = lambda k: sum(n**k*q**n/(1 - q**n) for n in range(1, M))
    return -5*s(3), -(5*s(3) + 7*s(5))/12
fig, ax = plt.subplots(figsize=(6.2, 5))
X, Y = np.meshgrid(np.linspace(-1.2, 1.5, 700), np.linspace(-1.6, 1.6, 700))
cols = [BLUE, GREEN, AMBER, INK]
for q, c in zip([0.03, 0.01, 0.002, 0.0], cols):
    a4, a6 = tate_coeffs(q) if q else (0, 0)
    ax.contour(X, Y, Y**2 + X*Y - X**3 - a4*X - a6, [0], colors=c, linewidths=2 if q == 0 else 1.5)
    ax.plot([], [], color=c, label=f"q = {q}" + ("  (the node: mult-add)" if q == 0 else ""))
ax.legend(loc="upper left", fontsize=9, frameon=False); ax.set_aspect("equal")
ax.set_title("Smoothing the node: the real Tate curve y² + xy = x³ + a₄(q)x + a₆(q)", fontsize=10)
plt.tight_layout(); plt.savefig("figures/f2_tate_pinch.png"); plt.close()

# F3 landscape with rung curve, pendulum segment, Frobenius
fig, ax = plt.subplots(figsize=(7.4, 5.4))
kk = np.linspace(-2.6, 2.6, 300); ax.plot(kk, kk**2/4, color=INK, lw=1.5, label="border (cusps): k² = 4q")
ss = np.linspace(-1.2, 2.8, 300); ax.plot(ss, ss*(ss - 1)/4, color=BLUE, lw=1.5, label="rung curve (s, s(s-1)/4)")
ax.plot([0, 0], [-0.5, -1], color=AMBER, lw=4, solid_capstyle="round", label="pendulum / sn law, m from 0 to 1")
good = [(-a/math.sqrt(p)) for nm in R["frob"] for p, a, b in R["frob"][nm] if not b]
ax.scatter(good, np.ones(len(good)) + np.random.default_rng(0).uniform(-.03, .03, len(good)), s=3, color=GREEN, alpha=.35, label="good primes: (-a_p/√p, 1), scaled by units")
ax.scatter([1, -1], [0, 0], s=80, color=RED, zorder=5, label="bad primes: k = +1 non-split, -1 split")
ax.annotate("mult-add\n(non-split)", (1, 0), (1.25, -.55), fontsize=9, arrowprops=dict(arrowstyle="-", color=GRAY))
ax.annotate("Tate curve at q=0\n(split)", (-1, 0), (-2.4, -.55), fontsize=9, arrowprops=dict(arrowstyle="-", color=GRAY))
ax.text(-2.5, 1.55, "turning country\n(no real wall)", fontsize=9, color=GRAY); ax.text(1.6, -1.2, "two-wall country", fontsize=9, color=GRAY)
ax.axhline(0, color=GRAY, lw=.5); ax.axvline(0, color=GRAY, lw=.5)
ax.set_xlabel("k: corner weight (ab)"); ax.set_ylabel("q: bend (a²b, ab²)"); ax.set_xlim(-2.7, 2.9); ax.set_ylim(-1.4, 1.9)
ax.legend(fontsize=8, frameon=False, loc="upper right"); plt.tight_layout(); plt.savefig("figures/f3_map.png"); plt.close()

# F4 Sato-Tate
st = R["sato_tate"]; e = np.array(st["edges"]); mid = (e[1:] + e[:-1])/2
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.bar(mid, st["hist"], width=(e[1]-e[0])*.9, color=GREEN, alpha=.6, label="turning angles of 4 curves, p < 4000")
ax.plot(mid, st["expected"], "o-", color=INK, label="Sato-Tate (2/π) sin²θ")
ax.set_xlabel("turning angle θ_p (a_p = 2√p cos θ_p)"); ax.legend(fontsize=9, frameon=False)
plt.tight_layout(); plt.savefig("figures/f4_sato_tate.png"); plt.close()

# F5 verdict: layer gaps vs s
s = np.linspace(-0.6, 3.0, 500)
g6 = -s*(s - 1)*(45*s**2 - 125*s + 22)/17280
g8 = s*(s - 1)*(6405*s**4 - 18116*s**3 + 4018*s**2 + 49*s + 22)/483840
fig, ax = plt.subplots(figsize=(7, 3.8))
ax.plot(s, g6*1e3, color=BLUE, label="layer 6 gap ×1000 (nothing left to tune)")
ax.plot(s, g8*1e3, color=RED, label="layer 8 gap ×1000 (after tuning a₆ at layer 7)")
ax.axhline(0, color=GRAY, lw=.7)
for v, lab in [(0, "addition\n(cusp)"), (1, "mult-add\n(node)")]:
    ax.plot([v], [0], "o", color=INK); ax.text(v, .8, lab, ha="center", fontsize=9)
ax.set_ylim(-3, 3); ax.set_xlabel("rung s"); ax.legend(fontsize=9, frameon=False, loc="lower left")
ax.set_title("Where a cubic curve can copy rung s: both gaps vanish only at s = 0 and 1", fontsize=10)
plt.tight_layout(); plt.savefig("figures/f5_verdict.png"); plt.close()

# F6 cnoidal to soliton
mp.mp.dps = 15
fig, ax = plt.subplots(figsize=(7, 3.2))
xs = np.linspace(-12, 12, 800)
for m, c in [(0.9, BLUE), (0.999, GREEN), (0.9999999, AMBER)]:
    ys = [2*m*float(mp.ellipfun('cn', x, m=m))**2 for x in xs]; ax.plot(xs, ys, color=c, lw=1.3, label=f"cnoidal wave, m = {m}")
ax.plot(xs, 2/np.cosh(xs)**2, "--", color=INK, lw=1, label="soliton 2 sech²")
ax.legend(fontsize=8, frameon=False, loc="upper right"); ax.set_yticks([])
plt.tight_layout(); plt.savefig("figures/f6_cnoidal.png"); plt.close()
print("figures done")
