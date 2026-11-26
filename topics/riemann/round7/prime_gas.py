"""Round 7, item 1: the prime gas (primon gas) with the order dial as statistics.

Single-particle levels E_p = log p.  A state of the gas is a multiset of primes
(Bose) or a word in primes (Boltzmann / infinite statistics).  Energy of n is log n.

  Bose        Z_B(b) = prod (1-p^-b)^-1      = zeta(b)
  Fermi       Z_F(b) = prod (1+p^-b)         = zeta(b)/zeta(2b)
  Boltzmann   Z_M(b) = sum_words             = 1/(1-P(b)),  P = prime zeta
  Order dial  Z_q(b) = sum_n a_q(n) n^-b,  a_q(n) = sum over orderings of the prime
              factors of n of q^inversions  (q=0 Bose, q=1 Boltzmann)

Exact identity (round 5): Z_q = (q;q)_inf sum_k q^k/(q;q)_k G_k,
G_k = prod_{j>=k} Zfug(b, q^j),  Zfug(b,z) = prod_p (1 - z p^-b)^-1  (fugacity z).
So the order-dial gas is a positive mixture of grand-canonical gases on the
fugacity ladder q^k, q^(k+1), ...
log G_k = sum_N b_k(N)/N logzeta(N b),  b_k(N) = sum_{m|N} mu(N/m) q^(k m)/(1-q^m).
"""
import json, math, itertools
import numpy as np
import mpmath as mp
from sympy import primerange, factorint, mobius, divisors

mp.mp.dps = 20
OUT = {}

def qpoch(a, q, n=None):
    if n is None:
        r = mp.mpf(1); j = 0
        while True:
            t = a * q**j
            if abs(t) < mp.mpf(10)**(-30): break
            r *= (1 - t); j += 1
            if j > 5000: break
        return r
    r = mp.mpf(1)
    for j in range(n): r *= (1 - a*q**j)
    return r

def qmultinom(es, q):
    n = sum(es)
    num = qpoch(q, q, n)
    den = mp.mpf(1)
    for e in es: den *= qpoch(q, q, e)
    return num/den

def gaussian_multinom_int(es, q):
    """exact q-multinomial as integer polynomial value for q = -1 or 1 via product of binomials."""
    from sympy import symbols, expand, Poly
    x = symbols('x')
    def qint(m): return sum(x**i for i in range(m))
    def qfact(m):
        r = 1
        for i in range(1, m+1): r *= qint(i)
        return r
    from sympy import cancel
    expr = cancel(qfact(sum(es)) / math.prod([qfact(e) for e in es]) if False else 1)
    # use iterated binomials: multinomial = prod binom(e1+..+ei, ei)
    tot = 0; val = 1
    for e in es:
        tot += e
        val *= cancel(qfact(tot)/(qfact(tot-e)*qfact(e)))
    return expand(val).subs(x, q)

def a_q_direct(n, q):
    f = factorint(n)
    es = list(f.values())
    if q in (-1, 1):
        return int(gaussian_multinom_int(es, q))
    return qmultinom(es, mp.mpf(q))

def Z_direct(beta, q, N):
    s = mp.mpf(0)
    for n in range(1, N+1):
        s += a_q_direct(n, q) * mp.mpf(n)**(-beta)
    return s

def b_k(N, k, q):
    s = mp.mpf(0)
    for m in divisors(N):
        s += mobius(N//m) * q**(k*m) / (1 - q**m)
    return s

def logZ_q_real(beta, q, Nmax=80, kmax=None):
    """exact Z_q for real beta>1 via the fugacity-ladder identity."""
    q = mp.mpf(q)
    lz = [None] + [mp.log(mp.zeta(N*beta)) for N in range(1, Nmax+1)]
    if kmax is None:
        kmax = int(60/ max(1e-9, -math.log(abs(float(q))))) + 5 if q != 0 else 0
    pref = qpoch(q, q)
    tot = mp.mpf(0)
    for k in range(kmax+1):
        lg = mp.mpf(0)
        for N in range(1, Nmax+1):
            lg += b_k(N, k, q)/N * lz[N]
        tot += (q**k/qpoch(q, q, k)) * mp.e**lg
    return mp.log(pref*tot)

def P_prime(beta, Nmax=80):
    return sum(mobius(n)/mp.mpf(n) * mp.log(mp.zeta(n*beta)) for n in range(1, Nmax+1) if mobius(n) != 0)

# ---------- (a) identity check at beta = 2 for several q ----------
print("(a) identity Z_q = (q;q)_inf sum_k q^k/(q;q)_k G_k at beta = 2")
OUT['identity'] = {}
for q in [0.3, -0.5, 0.7]:
    exact = mp.e**logZ_q_real(2, q)
    rows = []
    for N in [500, 1000, 2000, 4000]:
        d = Z_direct(2, q, N)
        rows.append((N, float(d), float(exact - d)))
    print(f"  q={q}: exact {float(exact):.10f}; direct partial sums gap", [f"{r[2]:.2e}" for r in rows])
    OUT['identity'][str(q)] = {'exact': float(exact), 'partial': rows}

# ---------- (b) the q = -1 end ----------
print("(b) q = -1 end: claim Z_-1 = (1+P(b))/(1-P(2b))")
exact_m1 = (1 + P_prime(2)) / (1 - P_prime(4))
rows = []
for N in [500, 1000, 2000, 4000]:
    d = Z_direct(2, -1, N)
    rows.append((N, float(d), float(exact_m1 - d)))
print(f"  closed form {float(exact_m1):.10f}; gaps", [f"{r[2]:.2e}" for r in rows])
# coefficient rule: nonzero iff at most one odd exponent
bad = 0
for n in range(2, 3000):
    es = list(factorint(n).values())
    odd = sum(e % 2 for e in es)
    v = a_q_direct(n, -1)
    pred = 0 if odd > 1 else math.factorial(sum(e//2 for e in es)) // math.prod(math.factorial(e//2) for e in es)
    if v != pred: bad += 1
print("  coefficient rule mismatches n<3000:", bad)
OUT['q_minus1'] = {'closed_form': float(exact_m1), 'partial': rows, 'coef_rule_mismatches': bad}
# Boltzmann end wall
bH = mp.findroot(lambda b: P_prime(b) - 1, 1.4)
print("  Boltzmann (q=1) Hagedorn wall P(b)=1 at b =", mp.nstr(bH, 8), "; q=-1 pole 1-P(2b)=0 at b =", mp.nstr(bH/2, 8))
OUT['hagedorn'] = float(bH)

# ---------- (c) critical exponent near beta = 1 ----------
print("(c) exponent: log Z_q ~ gamma log(1/(b-1)) + log A,  predicted gamma = 1/(1-q), A -> (q;q)_inf")
OUT['exponent'] = {}
for q in [-0.5, 0.0, 0.3, 0.5, 0.7]:
    pts = [mp.mpf(10)**(-e) for e in (3, 4, 5, 6)]
    vals = [logZ_q_real(1+eps, q, Nmax=60) for eps in pts]
    g = (vals[-1]-vals[-2]) / (mp.log(1/pts[-1]) - mp.log(1/pts[-2]))
    A = mp.e**(vals[-1] - g*mp.log(1/pts[-1]))
    pred = 1/(1-q)
    print(f"  q={q:+.1f}: local slope {float(g):.6f} vs 1/(1-q) = {pred:.6f};  amplitude {float(A):.5f} vs (q;q)_inf {float(qpoch(mp.mpf(q), mp.mpf(q))) if q!=0 else 1.0:.5f}")
    OUT['exponent'][str(q)] = {'slope': float(g), 'pred': pred, 'amp': float(A),
                               'qpoch': float(qpoch(mp.mpf(q), mp.mpf(q))) if q != 0 else 1.0}

# ---------- (d) thermodynamics on the real axis ----------
print("(d) energy E = -dlogZ/db and heat capacity C = b^2 dE/db... computing curves")
betas = [1.0 + 10**x for x in np.linspace(-2.5, 0.3, 40)]
thermo = {}
for q in [0.0, 0.3, 0.6]:
    lz = [float(logZ_q_real(mp.mpf(b), q, Nmax=50)) for b in betas]
    thermo[str(q)] = lz
# Boltzmann and Fermi
bz = []; bb = []
for b in np.linspace(float(bH)+0.003, 3.0, 40):
    bb.append(b); bz.append(float(-mp.log(1 - P_prime(b))))
thermo['boltz_beta'] = bb; thermo['boltz'] = bz
thermo['fermi'] = [float(mp.log(mp.zeta(b)/mp.zeta(2*b))) for b in betas]
thermo['betas'] = betas
OUT['thermo'] = thermo
json.dump(OUT, open('prime_gas.json', 'w'), indent=1)
print("saved prime_gas.json")
