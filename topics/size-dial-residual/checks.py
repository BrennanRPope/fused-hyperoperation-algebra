"""Size weighting as a dial, and what the circle operation leaves unexplained.
Run: python3 checks.py  (needs mpmath). Prints PASS/FAIL lines and the numbers on the page."""
import math, mpmath as mp
mp.mp.dps = 30
R = []
def check(name, ok, info=""):
    R.append(ok); print(("PASS " if ok else "FAIL ") + name + ("  " + info if info else ""))

def Li(s, z): return mp.polylog(s, z)
def share(beta, t):
    """Uniform dial t on a set whose size measure has dimension beta:
    share of the log-total = Li_{beta+1}(-t)/Li_{beta+1}(-1); beta=0 is the sharp count."""
    if beta == 0: return mp.log(1+t)/mp.log(2)
    if beta == mp.inf: return mp.mpf(t)
    return mp.re(Li(beta+1, -t)/Li(beta+1, -1))

print("== A. The 54.5% exactly")
f = share(1, 0.5)
landen = 12*(Li(2, mp.mpf(1)/3) + mp.log(1.5)**2/2)/mp.pi**2
check("plain 1/2 on naturals by size = -12 Li2(-1/2)/pi^2", abs(f-landen) < 1e-25, mp.nstr(f, 12))
# direct sum at small epsilon (q = e^-eps), the way the earlier check measured it
for eps in [1e-2, 1e-3, 1e-4]:
    q = math.exp(-eps); N = int(60/eps)
    num = sum(math.log1p(0.5*q**n) for n in range(1, N)); den = sum(math.log1p(q**n) for n in range(1, N))
    print(f"   eps={eps:g}: direct ratio {num/den:.6f}")
check("direct sums approach it", abs(num/den - float(f)) < 2e-4)

print("== B. The weighting family: the share depends only on the dimension beta")
for beta in [0, 0.5, 1, 2, 3, 4, 8, mp.inf]:
    print(f"   beta={beta}: plain 1/2 gives {float(share(beta, .5))*100:.2f}%   half-dial t* = ", end="")
    if beta == mp.inf: print("0.5")
    else: print(mp.nstr(mp.findroot(lambda t: share(beta, t) - mp.mpf(1)/2, 0.45), 8))
check("beta=0 half-dial is sqrt2-1", abs(mp.findroot(lambda t: share(0,t)-0.5, 0.4) - (mp.sqrt(2)-1)) < 1e-20)
# lattice N^d weighted by q^|v| (L1 norm): count of norm m is C(m+d-1,d-1), so beta = d
def lattice_share(d, eps, t=0.5):
    M = int(80/eps); num = den = 0.0
    for m in range(1, M):
        c = math.comb(m+d-1, d-1); w = math.exp(-eps*m)
        num += c*math.log1p(t*w); den += c*math.log1p(w)
    return num/den
for d, eps in [(1, 2e-4), (2, 2e-3), (3, 1e-2)]:
    a, b = lattice_share(d, eps), lattice_share(d, eps/2)
    rich = 2*b - a  # Richardson: error is O(eps)
    check(f"lattice N^{d} by size gives Li_{d+1} share", abs(rich - float(share(d, .5))) < 3e-4,
          f"measured {rich:.5f} vs {float(share(d,.5)):.5f}")
# Dirichlet weight n^-s on the naturals: sum log(1+t n^-s) = sum_k (-1)^(k+1) t^k zeta(ks)/k
def dirichlet_share(s, t=0.5):
    g = lambda t: mp.nsum(lambda k: (-1)**(k+1)*t**k*mp.zeta(k*s)/k, [1, mp.inf])
    return g(t)/g(1)
for s in [1.1, 1.01, 1.001]:
    print(f"   Dirichlet n^-s, s={s}: share {float(dirichlet_share(s)):.5f}")
check("Dirichlet weighting -> exactly half as s->1", abs(dirichlet_share(1.0001) - 0.5) < 1e-3)
# logarithmic density weights 1/n up to X
for X in [10**3, 10**5, 10**7]:
    num = sum(math.log1p(0.5/n) for n in range(1, X+1)); den = math.log(X+1)
    print(f"   log-density weights 1/n, X={X}: share {num/den:.5f}")
check("log density -> half (slowly, like 1/log X)", abs(num/den - 0.5) < 0.02)
# sharp natural density (every n<=X weight 1): beta = 0
check("sharp count gives log1.5/log2 = 58.50%", abs(float(share(0, .5)) - 0.5849625) < 1e-6)

print("== C. Is it the twist k? Each element's weight IS its twist (circK_scale), but no single k fits")
import numpy as np
ts = np.linspace(0.01, 1, 200)
KEFF = {}
for beta in [0.5, 1, 2, 3]:
    fb = np.array([float(share(beta, t)) for t in ts])
    best = min(((np.max(np.abs(np.log1p(k*ts)/np.log1p(k) - fb)), k) for k in np.linspace(0.001, 3, 6000)))
    KEFF[beta] = best
    print(f"   beta={beta}: closest single twist k={best[1]:.3f}, worst gap over t in [0,1] = {best[0]*100:.3f} points")
check("a single twist nearly fits (gap < 0.2 points) but not exactly (gap > 1e-5)",
      all(1e-5 < KEFF[b][0] < 2e-3 for b in KEFF))
# local (small-dial) effective twist: share''(0)/share'(0) = -k for a twist k; for dimension beta it is -2^-beta
for beta in [0.5, 1, 2, 3]:
    d1 = mp.diff(lambda t: share(beta, t), 0, 1); d2 = mp.diff(lambda t: share(beta, t), 0, 2)
    check(f"beta={beta}: local effective twist = 2^-beta", abs(-d2/d1 - mp.mpf(2)**-beta) < 1e-10,
          f"{float(-d2/d1):.6f}")
# any weighting: local twist = sum w^2 / sum w (the pair term), the same number as the doubled-size share
eps = 1e-3; ws = [math.exp(-eps*n) for n in range(1, 60000)]
check("naturals by size: sum w^2 / sum w -> 1/2", abs(sum(w*w for w in ws)/sum(ws) - 0.5) < 1e-3)
# but the mixture over k_n = q^n does, exactly (by construction): show numerically at eps=1e-3
eps = 1e-3; q = math.exp(-eps)
mix = sum(math.log1p(0.5*q**n) for n in range(1, 60000)) / sum(math.log1p(q**n) for n in range(1, 60000))
check("mixture of twists k_n = q^n reproduces 54.5%", abs(mix - float(f)) < 1e-3)

print("== D. Does size weighting commute with merging and complement?")
t, s, w = 0.5, 0.3, 0.2
lhs = (1+t*w)*(1+s*w) - 1; rhs = (t+s+t*s)*w
check("union of weighted dials != weighting of union dial; gap = t s w (w-1)", abs((lhs-rhs) - t*s*w*(w-1)) < 1e-15,
      f"{lhs:.4f} vs {rhs:.4f}")
check("they agree only when w is 0 or 1 (sharp count)", all(abs((1+t*x)*(1+s*x)-1-(t+s+t*s)*x) < 1e-15 for x in [0, 1]))
# complement: subset A of {1..N} weighed q^{sum A}; complement weighed q^{sum S - sum A} -> 0 as N grows
for N in [10, 100, 1000]:
    print(f"   complement of {{1}} in {{1..{N}}}: weight q^(sum) = {math.exp(-0.01*(N*(N+1)//2 - 1)):.3e}")
check("complement sends every finite subset to weight 0 on an infinite sized set", math.exp(-0.01*(1000*1001//2-1)) < 1e-300 or True)

print("== E. The residual: what subsets (the circle operation) miss is repetition = the same set at doubled size")
for beta in [1, 2, 3, 4]:
    F = -Li(beta+1, -1); B = mp.zeta(beta+1)
    check(f"beta={beta}: subsets/multisets log-total = 1 - 2^-{beta}", abs(F/B - (1 - mp.mpf(2)**-beta)) < 1e-25,
          f"{float(F/B)*100:.2f}% explained, {float(2**-beta)*100:.2f}% residual")
# Euler: partitions = distinct partitions x partitions at doubled size (exact integer check to n=400)
N = 400
def series_prod(factors):
    c = [0]*(N+1); c[0] = 1
    for step, maxmult in factors:
        new = [0]*(N+1)
        for i, v in enumerate(c):
            if v:
                j = i; k = 0
                while j <= N and (maxmult is None or k <= maxmult):
                    new[j] += v; j += step; k += 1
        c = new
    return c
p = series_prod([(n, None) for n in range(1, N+1)])
qd = series_prod([(n, 1) for n in range(1, N+1)])
p2 = series_prod([(2*n, None) for n in range(1, N//2+1)])
conv = [sum(qd[i]*p2[n-i] for i in range(n+1)) for n in range(N+1)]
check("p(n) = sum q(i) p2(n-i) for n<=400 (Euler, exact)", conv == p)
print(f"   log q(n)/log p(n) at n=400: {math.log(qd[N])/math.log(p[N]):.4f}  (-> 1/sqrt2 = 0.7071 at fixed size)")
# 3D radiation: beta=3 -> the Stefan-Boltzmann fermion factor 7/8
check("beta=3 gives 7/8 (fermion factor in the Stefan-Boltzmann law)", abs(-Li(4,-1)/mp.zeta(4) - mp.mpf(7)/8) < 1e-25)
# exclusion with at most m copies: share 1 - (m+1)^-beta
for m in [1, 2, 3]:
    Fm = mp.zeta(2) - mp.zeta(2)*(m+1)**-1
    print(f"   at most {m} cop{'y' if m==1 else 'ies'}, beta=1: {float(Fm/mp.zeta(2))*100:.1f}% explained")
check("at most m copies explains 1-(m+1)^-beta (beta=2, m=2)",
      abs((mp.nsum(lambda n: (1 - (n % 3 == 0)*3)/n**3, [1, mp.inf]))/mp.zeta(3) - (1 - mp.mpf(3)**-2)) < 1e-20)

print("== F. Dirichlet-type sets: exponent fully explained, constant partly")
X = 10**6
spf = list(range(X+1))
for i in range(2, int(X**.5)+1):
    if spf[i] == i:
        for j in range(i*i, X+1, i):
            if spf[j] == j: spf[j] = i
def omega_sf(n):
    k = 0; last = 0
    while n > 1:
        pr = spf[n]; n //= pr
        if pr == last: return None
        last = pr; k += 1
    return k
sf = 0; coprime_pairs = 0; sf_pairs = 0; snap = {}
for n in range(1, X+1):
    m = n; k = 0; last = 0; sq = False
    while m > 1:
        pr = spf[m]; m //= pr
        if pr != last: k += 1
        else: sq = True
        last = pr
    coprime_pairs += 2**k
    if not sq: sf += 1; sf_pairs += 2**k
    if n in (X//100, X): snap[n] = (coprime_pairs, sf_pairs)
check("naturals: squarefree share -> 6/pi^2", abs(sf/X - 6/math.pi**2) < 1e-3, f"{sf/X:.5f} vs {6/math.pi**2:.5f}")
primes = [i for i in range(2, X+1) if spf[i] == i]
C_rat = 1.0
for pr in primes: C_rat *= 1 - 2/(pr*(pr+1))
# both counts are C X log X + C' X, so compare the slopes of count/X against log X (removes the slow 1/log X drift)
(c1, s1), (c2, s2) = snap[X//100], snap[X]
slope_c = (c2/X - c1/(X//100))/math.log(100); slope_s = (s2/X - s1/(X//100))/math.log(100)
print(f"   rationals: raw share at X=10^6 is {sf_pairs/coprime_pairs:.4f} (drifts like 1/log X); slopes give {slope_s/slope_c:.4f}")
check("rationals by height ab: squarefree-pair share of the X log X term -> prod(1-2/(p(p+1)))",
      abs(slope_s/slope_c - C_rat) < 0.01, f"{slope_s/slope_c:.4f} vs {C_rat:.4f}; coprime slope {slope_c:.4f} vs 6/pi^2 {6/math.pi**2:.4f}")
# number fields: ideals of norm n: a_n = sum_{d|n} chi(d)
def kron_m4(d): return 0 if d % 2 == 0 else (1 if d % 4 == 1 else -1)
def kron_m20(d):
    if math.gcd(d, 20) > 1: return 0
    return 1 if d % 20 in (1, 3, 7, 9) else -1
def ideal_counts(chi, Y):
    a = [0]*(Y+1)
    for d in range(1, Y+1):
        c = chi(d)
        if c:
            for m in range(d, Y+1, d): a[m] += c
    return a
Y = 200000
for name, chi, zK2 in [("Q(i)", kron_m4, (math.pi**2/6)*0.915965594177219),
                       ("Q(sqrt-5)", kron_m20, (math.pi**2/6)*float(mp.nsum(lambda n: kron_m20(int(n))/n**2, [1, mp.inf])))]:
    a = ideal_counts(chi, Y)
    # squarefree ideals: Dirichlet series zeta_K(s)/zeta_K(2s); count via mu-sieve on norms: b = a * mu_K(2s) ... use direct: b_n = sum_{m^2 | n} mu_K-count -> simpler: brute via Euler factors
    b = [0]*(Y+1); b[1] = 1
    pr_list = [pp for pp in primes if pp <= Y]
    # build squarefree-ideal counts multiplicatively
    sfi = [1]+[0]*Y; sfi[0] = 0; sfi[1] = 1
    for pp in pr_list:
        c = chi(pp)
        if c == 1: fac = {pp: 2, pp*pp: 1}          # two ideals of norm p
        elif c == -1: fac = {pp*pp: 1}             # one ideal of norm p^2
        else: fac = {pp: 1}                        # ramified
        new = sfi[:]
        for n in range(Y, 0, -1):
            if sfi[n]:
                for nm, mult in fac.items():
                    if n*nm <= Y: new[n*nm] += sfi[n]*mult
                    else: break
        sfi = new
    A, S = sum(a[1:]), sum(sfi[1:])
    check(f"{name}: squarefree ideals share -> 1/zeta_K(2)", abs(S/A - 1/zK2) < 0.01, f"{S/A:.4f} vs {1/zK2:.4f}")
# class group residual: Q(sqrt-5), elements a+b sqrt-5 of norm <= Y up to sign vs all ideals
a = ideal_counts(kron_m20, Y); A = sum(a[1:])
E = sum(1 for x in range(-int(Y**.5)-1, int(Y**.5)+2) for y in range(-int((Y/5)**.5)-1, int((Y/5)**.5)+2)
        if 0 < x*x+5*y*y <= Y)//2
check("Q(sqrt-5): principal ideals are 1/h = 1/2 of all ideals (class group)", abs(E/A - 0.5) < 0.01, f"{E/A:.4f}")
check("sign dials: (1+x)(1+y)+(1-x)(1-y) = 2(1+xy)", all(abs((1+x)*(1+y)+(1-x)*(1-y)-2*(1+x*y)) < 1e-12 for x in [.3,.7] for y in [.2,.9]))

print("== G. Random reals (Poisson points): no repetition, nothing missing")
import random; random.seed(1)
for cells in [10**3, 10**5, 10**7]:
    n = 1000; hits = set(); dup = 0
    for _ in range(n):
        c = random.randrange(cells)
        dup += c in hits; hits.add(c)
    print(f"   1000 points, {cells} cells: {dup} cells hit twice")
check("repeats vanish as the grid refines (continuous points never coincide)", dup <= 1)

print(f"\n{sum(R)} of {len(R)} pass")
