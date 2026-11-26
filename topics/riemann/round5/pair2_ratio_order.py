"""Round 5, pair 2: the similarity ratio r with the order dial q.

The order dial weights each ordering of a factorisation by q^(inversions):
  Z_q(s) = sum_n W_q(n) n^-s,  W_q(n) = q-multinomial [Omega(n); e_1, ..., e_r]_q.
q = 0: one ordering per n, Z_0 = zeta (the harmonic string, lengths 1/n).
q = 1: every ordering, Z_1 = 1/(1 - P(s)) (P = prime zeta): the self-similar system with one copy
       per prime, ratio r_p = 1/p. Its similarity (Moran) dimension solves sum_p p^-D = 1, and its
       complex dimensions are the roots of P(s) = 1.
0 < q < 1: exact identity (Euler's q-exponential twice):
  W_q(n) = (q;q)_Omega * prod_i 1/(q;q)_(e_i),   1/(q^(N+1);q)_inf = sum_m q^((N+1)m)/(q;q)_m, so
  Z_q(s) = (q;q)_inf * sum_m q^m/(q;q)_m * prod_(j>=0) Z_(q^(m+j))(s)^(q^(m+j)),
  with Z_c = prod_p (1 - c p^-s)^(-1/c) the homothety cube at k = c. Each term is
  zeta(s)^(q^m/(1-q)) times a factor analytic and non-zero in Re s > 1/2.
"""
import json, math, itertools, numpy as np, mpmath as mp

res = {}
mp.mp.dps = 25

# ---- Moran / similarity dimension ----
D = mp.findroot(lambda x: mp.primezeta(x) - 1, 1.4)
res['moran_dimension'] = float(D)

# ---- complex dimensions: roots of P(s) = 1 in 1 < Re s < 1.45, 0 < t < 60 ----
mp.mp.dps = 15
sig = np.arange(1.03, 1.45, 0.02); tt = np.arange(0.0, 60.0, 0.15)
F = np.array([[abs(complex(mp.primezeta(mp.mpc(a, b))) - 1) for a in sig] for b in tt])
roots = []
for i in range(0, len(tt) - 1):
    for j in range(1, len(sig) - 1):
        lo_i = max(i - 1, 0)
        w = F[lo_i:i+2, j-1:j+2]
        if F[i, j] == w.min() and F[i, j] < 0.25:
            try:
                r = mp.findroot(lambda z: mp.primezeta(z) - 1, mp.mpc(sig[j], tt[i] + 1e-3))
                if 1.0 < r.real < 1.5 and -0.01 < r.imag < 61 and abs(mp.primezeta(r) - 1) < 1e-8:
                    rr = (round(float(r.real), 7), round(abs(float(r.imag)), 7))
                    if all(abs(rr[0]-a) + abs(rr[1]-b) > 1e-5 for a, b in roots): roots.append(rr)
            except Exception:
                pass
roots.sort(key=lambda p: p[1])
res['complex_dimensions'] = roots
res['complex_dimensions_min_re'] = min(p[0] for p in roots)
res['complex_dimensions_max_re'] = max(p[0] for p in roots)

# ---- W_q(n) as q^(inversions) over orderings, checked by brute force ----
def qint(n, q): return sum(q**i for i in range(n))
def qfact(n, q):
    out = 1.0
    for i in range(1, n + 1): out *= qint(i, q)
    return out
def qmulti(es, q):
    return qfact(sum(es), q)/math.prod(qfact(e, q) for e in es)
def brute(es, q):
    word = [i for i, e in enumerate(es) for _ in range(e)]
    tot = 0.0
    for perm in set(itertools.permutations(word)):
        inv = sum(1 for a in range(len(perm)) for b in range(a + 1, len(perm)) if perm[a] > perm[b])
        tot += q**inv
    return tot
bf = []
for es in [(1, 1), (2, 1), (2, 2), (3, 1, 1), (1, 1, 1, 1), (2, 2, 1)]:
    for q in [0.3, 0.7]:
        bf.append(abs(brute(es, q) - qmulti(es, q)))
res['qmultinomial_bruteforce_maxerr'] = max(bf)

# ---- exact identity at s = 2 ----
def primes_upto(n):
    sv = np.ones(n + 1, bool); sv[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    return np.nonzero(sv)[0]
N = 2_000_000
P = primes_upto(N)
# exponent vectors via smallest prime factor sieve
spf = np.zeros(N + 1, np.int64)
for p in P[::-1]:
    spf[p::p] = p
def evec(n):
    es = []
    while n > 1:
        p = spf[n]; e = 0
        while n % p == 0: n //= p; e += 1
        es.append(e)
    return es
ident = {}
for q in [0.3, 0.6, 0.9]:
    # direct Dirichlet sum (n up to 2e5 is plenty at s = 2 for 1e-5)
    M = 200_000
    direct = 1.0 + sum(qmulti(evec(n), q)*n**-2.0 for n in range(2, M + 1))
    # product side: log prod_p prod_j (1 - c q^j p^-2)^-1 = sum_r c^r P(2r) / (r (1 - q^r))
    PZ = [None] + [float(mp.primezeta(2*r)) for r in range(1, 80)]
    qq = 1.0
    for i in range(1, 400): qq *= (1 - q**i)
    tot = 0.0; qm = 1.0
    for m in range(0, 150):
        if m > 0: qm *= (1 - q**m)
        c = q**m
        logprod = sum(c**r*PZ[r]/(r*(1 - q**r)) for r in range(1, 80))
        tot += q**m/qm*math.exp(logprod)
    ident[str(q)] = dict(direct=direct, formula=qq*tot)
res['order_dial_identity_s2'] = ident

# ---- continuation of Z_q along the critical line through the identity (q = 0.5 and 0.3) ----
def necklace(c, Nn):
    return sum(mp.mpf(c)**d*mobius(Nn//d) for d in range(1, Nn + 1) if Nn % d == 0)/Nn
def mobius(n):
    if n == 1: return 1
    r, p, m = 1, 2, n
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
mp.mp.dps = 15
MC = {}
def mc_table(q, mmax, jmax, Nmax):
    key = (q, mmax, jmax, Nmax)
    if key not in MC:
        MC[key] = [[None, None] + [sum(necklace(mp.mpf(q)**(m + j), Nn) for j in range(jmax))
                                   for Nn in range(2, Nmax + 1)] for m in range(mmax)]
    return MC[key]
def logzeta_cont(s):
    """log zeta continued along the horizontal path from Re s = 4 (the convention behind S(t))"""
    s = mp.mpc(s)
    if mp.re(s) >= 4: return mp.log(mp.zeta(s))
    sig = np.linspace(4.0, float(mp.re(s)), 80)
    val = mp.log(mp.zeta(mp.mpc(sig[0], mp.im(s)))); im_prev = mp.im(val)
    for sg in sig[1:]:
        v = mp.log(mp.zeta(mp.mpc(sg, mp.im(s))))
        im = mp.im(v)
        im += 2*mp.pi*round(float((im_prev - im)/(2*mp.pi)))
        im_prev = im
        val = mp.mpc(mp.re(v), im)
    return val
def Zq(s, q, mmax=40, jmax=60, Nmax=12):
    s = mp.mpc(s); qq = mp.mpf(1)
    for i in range(1, 200): qq *= (1 - mp.mpf(q)**i)
    tab = mc_table(q, mmax, jmax, Nmax)
    lz = logzeta_cont(s)
    lzN = [None, None] + [logzeta_cont(Nn*s) for Nn in range(2, Nmax + 1)]
    tot = 0; qm = mp.mpf(1)
    for m in range(mmax):
        if m > 0: qm *= (1 - mp.mpf(q)**m)
        expo = mp.mpf(q)**m/(1 - q)
        g = sum(tab[m][Nn]*lzN[Nn] for Nn in range(2, Nmax + 1))
        tot += mp.mpf(q)**m/qm*mp.exp(expo*lz + g)
    return qq*tot
# sanity: formula at s = 2 equals the direct sum above
res['Zq_formula_s2_q0.6'] = float(mp.re(Zq(2, 0.6)))
line = {}
for q in [0.5, 0.3]:
    ts = np.linspace(10, 35, 126)
    line[str(q)] = dict(t=ts.tolist(), absZq=[float(abs(Zq(mp.mpc(0.5, t), q))) for t in ts],
                        absZeta=[float(abs(mp.zeta(mp.mpc(0.5, t)))) for t in ts])
res['critical_line'] = line
json.dump(res, open('pair2.json', 'w'), indent=1)
print({k: v for k, v in res.items() if k not in ('critical_line', 'complex_dimensions')})
print(len(roots), roots[:8])
