"""DRAFT · Rewinding a barrier chain with the fused operations (Nov 26, 2026).

Brennan: "Try rewinding using the fused operations" and "The goal of fusing is to keep information".
Run: python3 rewind_checks.py > rewind_checks_output.txt   (mpmath, sympy, numpy; about three minutes)

The fused ladder: for each whole n, g_n(t) = (P_n(1 + 2t) − 1)/(n(n + 1)) starts as t + ..., and its average over a
spread composes by ∘ₖ with k = n(n + 1) when a barrier is added. So the spread's information is carried by a
ladder of fused laws, one per rung n, with homotheties 2, 6, 12, 20, ...; rewinding is each law's ∘ₖ inverse.
"""
import math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
results = []
def check(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"   [{info}]" if info else ""))

tt = sp.symbols('t')
NM = 24
COEF = []   # P_n(1 + 2t) = sum_j COEF[n][j] t^j, exact rationals
for n in range(NM + 1):
    poly = sp.Poly(sp.expand(sp.legendre(n, 1 + 2*tt)), tt)
    c = [mp.mpf(0)]*(n + 1)
    for (m,), v in zip(poly.monoms(), poly.coeffs()):
        v = sp.Rational(v); c[m] = mp.mpf(v.p)/v.q
    COEF.append(c)
def Pn(n, t):  # P_n(1 + 2t) by Horner
    s = mp.mpf(0)
    for c in reversed(COEF[n]): s = s*t + c
    return s
def g(n, t): return (Pn(n, t) - 1)/(n*(n + 1))
def circ(k, a, b): return a + b + k*a*b
def uncirc(k, c, b): return (c - b)/(1 + k*b)     # the ∘ₖ inverse: the rewind

# ---------------------------------------------------------------- an arbitrary initial spread and one barrier
t1 = mp.mpf('0.25')
a0, b0 = mp.mpf('0.3'), mp.mpf('1.2')
raw = lambda t: (t - a0)*(b0 - t)*(mp.mpf('1.5') + mp.sin(8*t))
Z = mp.quad(raw, [a0, b0]); p0 = lambda t: raw(t)/Z
init = [mp.quad(lambda t: g(n, t)*p0(t), [a0, b0]) for n in range(1, NM + 1)]
M = 64; phis = [2*mp.pi*j/M for j in range(M)]
def fwd_moment(n):
    def inner(t0):
        A = t0 + t1 + 2*t0*t1; B = 2*mp.sqrt(t0*t1*(1 + t0)*(1 + t1))
        return p0(t0)*sum(g(n, A - B*mp.cos(p)) for p in phis)/M
    return mp.quad(inner, [a0, (a0 + b0)/2, b0])
print("== 1. The fused ladder carries the whole spread")
fwd = [fwd_moment(n) for n in range(1, NM + 1)]
k = [n*(n + 1) for n in range(1, NM + 1)]
err = max(abs(fwd[i] - circ(k[i], init[i], g(i + 1, t1)))/abs(fwd[i]) for i in range(NM))
check("1.1 adding a barrier to any spread: ⟨g_n⟩ → ⟨g_n⟩ ∘ₖ g_n(t₁) with k = n(n + 1), for n = 1..24 (independent 2-D quadrature)",
      err < mp.mpf(10)**-30, f"max relative error {mp.nstr(err, 3)}")
check("1.2 rung 1 is the mean law ∘₂ (g₁ = t); the homotheties run 2, 6, 12, 20, … = the hyperbolic plane's eigenvalues n(n + 1)",
      abs(g(1, mp.mpf('0.7')) - mp.mpf('0.7')) < mp.mpf(10)**-35 and k[:4] == [2, 6, 12, 20])

print("\n== 2. The rewind is each rung's ∘ₖ inverse, and a barrier never puts a wall in the way")
rew = [uncirc(k[i], fwd[i], g(i + 1, t1)) for i in range(NM)]
err = max(abs(rew[i] - init[i])/abs(init[i]) for i in range(NM))
check("2.1 rewinding every rung by (G − g)/(1 + k g) returns the initial spread's ladder exactly",
      err < mp.mpf(10)**-30, f"max relative error {mp.nstr(err, 3)}")
fac = [1 + k[i]*g(i + 1, t1) for i in range(NM)]
check("2.2 the divisor 1 + k g_n(t₁) is P_n(1 + 2t₁) ≥ 1 for every barrier and rung: no wall, ever",
      all(f >= 1 for f in fac) and all(Pn(n, x) >= 1 for n in range(1, NM + 1) for x in [mp.mpf(i)/7 for i in range(0, 60)]),
      f"t₁ = ¼: {', '.join(mp.nstr(f, 4) for f in fac[:5])}, …")
check("2.3 so each rung's rewind shrinks errors (divides by ≥ 1); information is kept at every rung", all(f >= 1 for f in fac))
D1 = mp.acosh(1 + 2*t1); r1 = mp.tanh(D1/2)
rmax = mp.tanh((mp.acosh(1 + 2*b0) + D1)/2); rmin = mp.tanh((mp.acosh(1 + 2*a0) - D1)/2)
r0max = (rmax - r1)/(1 - rmax*r1); r0min = (rmin + r1)/(1 + rmin*r1)
check("2.4 the spread's edges rewind by the (k = 0, d = −1) inverse on amplitudes (relativity's velocity subtraction)",
      abs(r0max**2/(1 - r0max**2) - b0) < mp.mpf(10)**-30 and abs(r0min**2/(1 - r0min**2) - a0) < mp.mpf(10)**-30)

print("\n== 3. Rebuilding the spread from the rewound ladder")
def power_moments(gm):
    mom = [mp.mpf(1)]
    for n in range(1, len(gm) + 1):
        val = 1 + n*(n + 1)*gm[n - 1]
        mom.append((val - sum(COEF[n][j]*mom[j] for j in range(n)))/COEF[n][n])
    return mom
LEG = []
for j in range(NM + 1):
    poly = sp.Poly(sp.legendre(j, tt), tt); c = [mp.mpf(0)]*(j + 1)
    for (m,), v in zip(poly.monoms(), poly.coeffs()):
        v = sp.Rational(v); c[m] = mp.mpf(v.p)/v.q
    LEG.append(c)
def rebuild(mom, a, b, J, ts):
    xm = [sum(mp.binomial(m, j)*(2/(b - a))**j*mom[j]*(-(a + b)/(b - a))**(m - j) for j in range(m + 1)) for m in range(J + 1)]
    co = [(2*j + 1)/mp.mpf(2)*sum(LEG[j][m]*xm[m] for m in range(j + 1)) for j in range(J + 1)]
    return [sum(co[j]*mp.legendre(j, (2*t - a - b)/(b - a)) for j in range(J + 1))*2/(b - a) for t in ts]
ts = [a0 + (b0 - a0)*i/120 for i in range(121)]
truth = [p0(t) for t in ts]
mom = power_moments(rew)
errs = {}
for J in (6, 10, 16, 24):
    rec = rebuild(mom, a0, b0, J, ts); errs[J] = max(abs(r - q) for r, q in zip(rec, truth))
check("3.1 with exact moments the rebuilt initial spread converges: max error at 6, 10, 16, 24 rungs falls",
      errs[6] > errs[10] > errs[16] > errs[24] and errs[24] < 1e-6,
      ", ".join(f"{J}: {mp.nstr(errs[J], 2)}" for J in errs))
# noise: forward ladder known to relative precision eps
rng = np.random.default_rng(11)
noisy = {}
for eps in (1e-6, 1e-10):
    fw = [fwd[i]*(1 + mp.mpf(eps)*mp.mpf(rng.standard_normal())) for i in range(NM)]
    rw = [uncirc(k[i], fw[i], g(i + 1, t1)) for i in range(NM)]
    # absolute error after/before on each rung, against 1/P_n(1 + 2t1); relative growth |G|/|G - g| for the record
    shrink = max(abs(abs(rw[i] - init[i])/abs(fw[i] - fwd[i]) - 1/Pn(i + 1, t1)) for i in range(NM))
    shrink_max = max(abs(rw[i] - init[i])/abs(fw[i] - fwd[i]) for i in range(NM))
    relg = max(abs(fwd[i])/abs(fwd[i] - g(i + 1, t1)) for i in range(NM))
    mm = power_moments(rw); best = None
    for J in range(4, NM + 1, 2):
        e = max(abs(r - q) for r, q in zip(rebuild(mm, a0, b0, J, ts[::4]), truth[::4]))
        if best is None or e < best[1]: best = (J, e)
    noisy[eps] = (shrink, best, shrink_max, relg)
check("3.2 with noisy forward data the rewind shrinks each rung's error by exactly 1/P_n(1 + 2t₁) ≤ 1 (absolute error)",
      all(noisy[e][0] < 1e-20 and noisy[e][2] <= 1 for e in noisy),
      ", ".join(f"ε={e}: largest after/before {mp.nstr(noisy[e][2], 3)}" for e in noisy))
print(f"NOTE 3.2b the relative error can grow by |G|/|G − g| (here at most ×{mp.nstr(noisy[1e-6][3], 3)}), because the barrier's own"
      "\n          part g is subtracted from the signal; the absolute error, which is what the rebuild sees, only shrinks.")
check("3.3 the hard step is the old moment problem: the usable number of rungs grows as the noise falls",
      noisy[1e-10][1][0] >= noisy[1e-6][1][0] and noisy[1e-10][1][1] < noisy[1e-6][1][1],
      "; ".join(f"ε={e}: best {noisy[e][1][0]} rungs, error {mp.nstr(noisy[e][1][1], 2)}" for e in noisy))

print("\n== 4. The oscillating spectrum: what the earlier draft got wrong")
zc = 1 + 2*t1
lam_amp = max(abs(mp.legenp(-0.5 + 1j*l, 0, zc, type=3)) for l in [0.1*i for i in range(1, 80)])
check("4.1 in the conical (oscillating) spectrum a barrier multiplies by |P₋½₊ᵢλ| ≤ 1, so a rewind there divides by ≤ 1 and blows up near zeros",
      lam_amp <= 1)
print("NOTE 4.2 (correction): a finite chain's spread has bounded support, so its Mehler–Fock transform is entire in λ and its")
print("          values at the character's zeros are fixed by the values around them. Nothing is lost; it is only ill-conditioned there.")
print("          The fused ladder (whole n) is the coordinate system where the same rewind is wall-free and error-shrinking.")

print("\n== 5. Mixers (the sphere): the same ladder has walls, and a 50/50 mixer forgets for real")
def Pl(l, x): return mp.legendre(l, x)
def h(l, s): return (1 - Pl(l, 1 - 2*s))/(l*(l + 1))
s1 = mp.mpf('0.3')
q0 = lambda s: 6*s*(1 - s)
def mix_moment(l, dens, sm):
    def inner(s0):
        A = s0 + sm - 2*s0*sm; B = 2*mp.sqrt(s0*sm*(1 - s0)*(1 - sm))
        return dens(s0)*sum(h(l, A - B*mp.cos(p)) for p in phis)/M
    return mp.quad(inner, [0, 0.5, 1])
ok = True
for l in (1, 2, 3, 5):
    lhs = mix_moment(l, q0, s1); rhs = circ(-l*(l + 1), mp.quad(lambda s: h(l, s)*q0(s), [0, 1]), h(l, s1))
    ok &= abs(lhs - rhs) < mp.mpf(10)**-25
check("5.1 mixers: ⟨h_ℓ⟩ composes by ∘ₖ with k = −ℓ(ℓ + 1), h_ℓ(s) = (1 − P_ℓ(1 − 2s))/(ℓ(ℓ + 1))", ok)
check("5.2 the divisor 1 − ℓ(ℓ + 1)h_ℓ(s₁) = P_ℓ(1 − 2s₁) has zeros: a wall on rung ℓ wherever s₁ sits on one",
      abs(Pl(3, 1 - 2*(1 - mp.sqrt(mp.mpf(3)/5))/2)) < mp.mpf(10)**-30)
qa = lambda s: 6*s*(1 - s); qb = lambda s: 6*s*(1 - s)*(1 + mp.mpf('0.8')*(1 - 2*s))
half = mp.mpf('0.5')
same = all(abs(mix_moment(l, qa, half) - mix_moment(l, qb, half)) < mp.mpf(10)**-25 for l in range(1, 9))
diff_in = abs(mp.quad(lambda s: s*(qa(s) - qb(s)), [0, 1])) > 0.01
check("5.3 two different input spreads give the identical output after a 50/50 mixer (rungs 1–8): every odd rung is erased",
      same and diff_in)
check("5.4 barriers have no such point: P_n(1 + 2t₁) ≥ 1, so two different input spreads always give different outputs",
      all(Pn(n, t1) > 1 for n in range(1, NM + 1)))

print(f"\n{sum(results)}/{len(results)} checks pass")
