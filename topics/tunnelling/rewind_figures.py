"""DRAFT · Pictures for the rewind section of tunnelling/README.md. Run: python3 rewind_figures.py (writes figures/fig7, fig8)"""
import math, os
import numpy as np
import mpmath as mp
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREEN, ORANGE, RED, GREY, PURPLE = "#2b6cb0", "#2f855a", "#dd6b20", "#c53030", "#718096", "#6b46c1"
mp.mp.dps = 40
NM = 20
tt = sp.symbols('t')
COEF = []
for n in range(NM + 1):
    poly = sp.Poly(sp.expand(sp.legendre(n, 1 + 2*tt)), tt); c = [mp.mpf(0)]*(n + 1)
    for (m,), v in zip(poly.monoms(), poly.coeffs()):
        v = sp.Rational(v); c[m] = mp.mpf(v.p)/v.q
    COEF.append(c)
def Pn(n, t):
    s = mp.mpf(0)
    for c in reversed(COEF[n]): s = s*t + c
    return s
def g(n, t): return (Pn(n, t) - 1)/(n*(n + 1))
def circ(k, a, b): return a + b + k*a*b
def uncirc(k, c, b): return (c - b)/(1 + k*b)

# ---------------------------------------------------------------- Fig 7: the divisor on each rung
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
a = ax[0]
ns = np.arange(1, 13)
for t1, col in ((0.05, GREEN), (0.25, BLUE), (1.0, PURPLE)):
    a.semilogy(ns, [float(Pn(int(n), mp.mpf(t1))) for n in ns], "o-", color=col, ms=4, label=f"barrier, t = R/T = {t1}")
a.axhline(1, color=GREY, lw=1, ls="--"); a.text(12.2, 1.08, "1", color=GREY, va="bottom")
a.set_xlabel("rung n  (fused law ∘ₖ with k = n(n + 1) = 2, 6, 12, 20, …)")
a.set_ylabel("rewind divides the error by  P_n(1 + 2t)")
a.set_title("Barriers: every rung's rewind divides by ≥ 1\nno wall, errors shrink, nothing is lost", fontsize=10)
a.legend(frameon=False, loc="upper left")
a = ax[1]
ss = np.linspace(0, 1, 801)
for l, col in ((1, GREEN), (2, BLUE), (3, ORANGE), (4, PURPLE)):
    a.plot(ss, [abs(float(mp.legendre(l, 1 - 2*s))) for s in ss], color=col, label=f"mixer rung ℓ = {l}  (k = −{l*(l + 1)})")
a.axvline(0.5, color=RED, lw=1, ls=":")
a.text(0.5, 1.03, "50/50 mixer: every odd rung's divisor is 0", color=RED, ha="center", fontsize=9)
a.set_xlabel("mixer's crossing chance s"); a.set_ylabel("rewind divides by  |P_ℓ(1 − 2s)|")
a.set_ylim(0, 1.12)
a.set_title("Mixers: the same ladder with k < 0 has walls\nwhere the divisor hits 0, that rung is erased", fontsize=10)
a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
fig.tight_layout(); fig.savefig("figures/fig7-rewind-divisors.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 8: rebuilding, noise, and the mixer forgetting
t1 = mp.mpf('0.25'); a0, b0 = mp.mpf('0.3'), mp.mpf('1.2')
raw = lambda t: (t - a0)*(b0 - t)*(mp.mpf('1.5') + mp.sin(8*t))
Z = mp.quad(raw, [a0, b0]); p0 = lambda t: raw(t)/Z
init = [mp.quad(lambda t: g(n, t)*p0(t), [a0, b0]) for n in range(1, NM + 1)]
k = [n*(n + 1) for n in range(1, NM + 1)]
fwd = [circ(k[i], init[i], g(i + 1, t1)) for i in range(NM)]   # forward law (checked against quadrature in rewind_checks.py)
LEG = []
xs = sp.symbols('x')
for j in range(NM + 1):
    poly = sp.Poly(sp.legendre(j, xs), xs); c = [mp.mpf(0)]*(j + 1)
    for (m,), v in zip(poly.monoms(), poly.coeffs()):
        v = sp.Rational(v); c[m] = mp.mpf(v.p)/v.q
    LEG.append(c)
def power_moments(gm):
    mom = [mp.mpf(1)]
    for n in range(1, len(gm) + 1):
        val = 1 + n*(n + 1)*gm[n - 1]
        mom.append((val - sum(COEF[n][j]*mom[j] for j in range(n)))/COEF[n][n])
    return mom
def rebuild(mom, J, ts):
    a, b = a0, b0
    xm = [sum(mp.binomial(m, j)*(2/(b - a))**j*mom[j]*(-(a + b)/(b - a))**(m - j) for j in range(m + 1)) for m in range(J + 1)]
    co = [(2*j + 1)/mp.mpf(2)*sum(LEG[j][m]*xm[m] for m in range(j + 1)) for j in range(J + 1)]
    return [float(sum(co[j]*mp.legendre(j, (2*t - a - b)/(b - a)) for j in range(J + 1))*2/(b - a)) for t in ts]
ts = [a0 + (b0 - a0)*i/150 for i in range(151)]
tf = [float(t) for t in ts]
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
a = ax[0]
a.plot(tf, [float(p0(t)) for t in ts], color="black", lw=2.2, label="the spread before the barrier")
mom = power_moments([uncirc(k[i], fwd[i], g(i + 1, t1)) for i in range(NM)])
for J, col in ((4, ORANGE), (8, BLUE), (14, GREEN)):
    a.plot(tf, rebuild(mom, J, ts), color=col, lw=1.2, ls="--", label=f"rewound from {J} rungs")
a.set_xlabel("t = R/T before the last barrier"); a.set_ylabel("spread")
a.set_title("Rewinding the fused rungs gets the spread back\n(one barrier, t = ¼, removed)", fontsize=10)
a.legend(frameon=False, fontsize=8, loc="upper left")
a = ax[1]
rng = np.random.default_rng(3)
tsc = ts[::5]; truth = [float(p0(t)) for t in tsc]
for eps, col in ((1e-4, ORANGE), (1e-8, BLUE), (1e-12, GREEN), (0, "black")):
    fw = [fwd[i]*(1 + mp.mpf(eps)*mp.mpf(rng.standard_normal())) for i in range(NM)]
    mm = power_moments([uncirc(k[i], fw[i], g(i + 1, t1)) for i in range(NM)])
    Js = list(range(2, NM + 1, 2))
    errs = [max(abs(r - q) for r, q in zip(rebuild(mm, J, tsc), truth)) for J in Js]
    a.semilogy(Js, errs, "o-", ms=3, color=col, label="exact data" if eps == 0 else f"data known to 1 part in {1/eps:.0e}")
a.set_xticks(range(2, NM + 1, 2)); a.set_xlabel("rungs used"); a.set_ylabel("largest error in the rebuilt spread")
a.set_title("The rewind never amplifies noise;\nthe hard step is the old moment problem", fontsize=10)
a.legend(frameon=False, fontsize=8)
a = ax[2]
N = 3_000_000
def sample(odd):
    s0 = rng.beta(2, 2, size=4*N)
    if odd:  # density 6s(1−s)(1 + 0.8(1 − 2s)) by rejection
        s0 = s0[rng.random(s0.size) < (1 + 0.8*(1 - 2*s0))/1.8]
    s0 = s0[:N]; phi = rng.uniform(0, 2*np.pi, N)
    c = 2*np.sqrt(s0*(1 - s0))*np.cos(phi)   # 1 − 2s after a 50/50 mixer
    return s0, (1 - c)/2
bins = np.linspace(0, 1, 61); mid = (bins[1:] + bins[:-1])/2
for odd, col, lab in ((False, BLUE, "spread A"), (True, ORANGE, "spread B")):
    s0, s = sample(odd)
    h0, _ = np.histogram(s0, bins, density=True); h, _ = np.histogram(s, bins, density=True)
    a.plot(mid, h0, color=col, lw=1.3, ls="--", label=f"{lab} going in")
    a.plot(mid, h, color=col, lw=2.2 if not odd else 1.2, label=f"{lab} after a 50/50 mixer")
a.set_xlabel("crossing chance s"); a.set_ylabel("spread")
a.set_title("A 50/50 mixer really forgets:\ntwo different spreads come out the same", fontsize=10)
a.legend(frameon=False, fontsize=8, loc="lower center")
fig.tight_layout(); fig.savefig("figures/fig8-rewind-rebuild.png", dpi=150); plt.close(fig)
print("wrote fig7, fig8")
