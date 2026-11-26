"""Round 5, pair 1: the height dial s (hypernatural log) with the homothety k on the prime cube.

Two ways to bend the logarithm:
  height s:    ln_s(X) = E^(-s)(X - 1), E(x) = e^x - 1   (regular iteration; ln_1 = log, ln_0 = X - 1)
  homothety k: l_k(X)  = (1 - X^(-k))/k                  (l_0 = log; turns x into the circle op o_(-k))
Both start u - (bend/2)u^2 with bend s and 1 + k, so they are tangent when s = 1 + k.

(A) l_k undoes the prime cube exactly: l_k(Z_k) = (1 - prod_p(1 - k p^-s))/k, a corner sum with
    coefficients mu(n) k^(omega(n)) (the generalised Mertens function M_k).
(B) l_k(zeta) is singular exactly at the zeros of zeta for every k > 0; ln_s(zeta) for s < 1 is
    singular instead where zeta = X_s = 1 + E^s(-inf) < 0 (the height slides the singular value
    from -inf at s = 0 to 0 at s = 1). Those a-points of zeta are not on the critical line.
(C) M_k(x)/sqrt(x) up to 10^7 for k = 1/2, 1, 2, 3.
"""
import json, math, numpy as np, mpmath as mp

res = {}

# ---------- regular iteration of E on the negative side (Abel function) ----------
# alpha(E(x)) = alpha(x) + 1; near 0-: alpha(x) ~ -2/x + (1/3) log(-x) + c1 x + ...
# c1 from the next order of the Abel equation (solved symbolically below, checked numerically)
import sympy as sp
xs_, c1_, c2_ = sp.symbols('x c1 c2')
Ef = sum(xs_**i/sp.factorial(i) for i in range(1, 8))
alph = lambda y: -2/y + sp.Rational(1, 3)*sp.log(-y) + c1_*y + c2_*y**2
eq = sp.series(alph(Ef) - alph(xs_) - 1, xs_, 0, 4).removeO()
sol = sp.solve([sp.expand(eq).coeff(xs_, 2), sp.expand(eq).coeff(xs_, 3)], [c1_, c2_], dict=True)[0]
C1, C2 = float(sol[c1_]), float(sol[c2_])
res['abel_coeffs'] = dict(c1=str(sol[c1_]), c2=str(sol[c2_]))

NWALK = 20000
def alpha(x):
    """Abel function of E on (-inf, 0), normalised by the asymptotic series (no constant)."""
    x = float(x)
    if x == -math.inf:
        return None
    for _ in range(NWALK):
        x = math.expm1(x)
    return -2/x + math.log(-x)/3 + C1*x + C2*x*x - NWALK

A_m1 = alpha(-1.0)
A_minf = A_m1 - 1.0              # alpha(-inf) = alpha(E(-inf)) - 1 = alpha(-1) - 1

def alpha_inv(a):
    lo, hi = -60.0, -1e-9          # alpha increasing on (-inf, 0)
    for _ in range(70):
        mid = 0.5*(lo + hi)
        if alpha(mid) < a: lo = mid
        else: hi = mid
    return 0.5*(lo + hi)

def E_s(x, s):
    return alpha_inv(alpha(x) + s)

# checks: E^1 = E, E^(1/2) twice = E
res['check_E1'] = abs(E_s(-0.7, 1.0) - math.expm1(-0.7))
h = E_s(-0.7, 0.5); res['check_half_twice'] = abs(E_s(h, 0.5) - math.expm1(-0.7))

# singular value of ln_s:  X_s = 1 + E^s(-inf),  E^s(-inf) = alpha_inv(alpha(-inf) + s)
svals = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
Xs = {}
for s in svals:
    Xs[str(s)] = 1 + alpha_inv(A_minf + s)
res['singular_value_X_s'] = Xs

# ln_s(X) near X = 0 for s = 1/2: smooth (finite derivative), unlike log
s = 0.5
ys = np.linspace(-1.2, -0.8, 9)          # y = X - 1 around -1
vals = [alpha_inv(alpha(y) - s) for y in ys]
res['ln_half_near_zero'] = dict(y=ys.tolist(), value=vals,
                                second_diff_max=float(np.max(np.abs(np.diff(vals, 2)))))

# real a-points of zeta in (1/2, 1): zeta(sigma) runs from -1.46035 (sigma = 1/2) to -inf (sigma -> 1)
mp.mp.dps = 20
zh = float(mp.zeta(0.5)); res['zeta_half'] = zh
apts = {}
for s_, X in Xs.items():
    if X < zh:
        apts[s_] = float(mp.findroot(lambda x: mp.zeta(x) - X, (0.5 + 1e-6, 1 - 1e-9), solver='anderson'))
res['real_a_points_off_line'] = apts

# complex a-points with 1/2 < sigma < 2, 0 < t < 60 for three heights (grid minima + Newton)
cpts = {}
for s_ in ['0.5', '0.8', '0.95']:
    X = Xs[s_]
    sig = np.arange(0.02, 2.0, 0.02); tt = np.arange(0.5, 60, 0.1)
    Z = np.array([[complex(mp.zeta(mp.mpc(a, b))) for a in sig] for b in tt])
    F = np.abs(Z - X)
    found = []
    for i in range(1, len(tt) - 1):
        for j in range(1, len(sig) - 1):
            w = F[i-1:i+2, j-1:j+2]
            if F[i, j] == w.min() and F[i, j] < 0.6:
                try:
                    r = mp.findroot(lambda z: mp.zeta(z) - X, mp.mpc(sig[j], tt[i]))
                    if 0.0 < r.real < 2.2 and 0 < r.imag < 61 and abs(mp.zeta(r) - X) < 1e-10:
                        rr = (round(float(r.real), 8), round(float(r.imag), 8))
                        if rr not in found: found.append(rr)
                except Exception:
                    pass
    cpts[s_] = sorted(found, key=lambda p: p[1])
res['complex_a_points'] = cpts
res['complex_a_points_side'] = {k: dict(right=sum(1 for p in v if p[0] > 0.5), left=sum(1 for p in v if p[0] < 0.5),
    closest=min(abs(p[0]-0.5) for p in v) if v else None) for k, v in cpts.items()}
# heights above 1: ln_s = log o E^-(s-1), singular where X - 1 = E^(s-1)(-1), a value in (-1, 0)
res['singular_value_above_1'] = {str(s1): 1 + E_s(-1.0, s1 - 1) for s1 in [1.25, 1.5, 2.0, 3.0]}

# ---------- (A) exact identity l_k(Z_k) = (1 - prod(1 - k p^-s))/k at s = 2 ----------
def primes_upto(n):
    sv = np.ones(n + 1, bool); sv[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    return np.nonzero(sv)[0]

P = primes_upto(2_000_000)
ident = {}
for k in [0.5, 1.0, 2.0, 3.0]:
    Zk = np.exp(np.sum(-np.log1p(-k*P.astype(float)**-2.0))/k)       # prod (1 - k p^-2)^(-1/k)
    lk = (1 - Zk**(-k))/k
    rhs = (1 - np.exp(np.sum(np.log1p(-k*P.astype(float)**-2.0))))/k
    ident[str(k)] = dict(lhs=float(lk), rhs=float(rhs))
res['lk_identity_s2'] = ident
# tangency: bends agree (s = 1 + k), third order differs
u = sp.symbols('u'); kk = sp.symbols('k')
lser = sp.series((1 - (1 + u)**(-kk))/kk, u, 0, 4).removeO()
res['lk_series'] = str(sp.expand(lser))

# ---------- (C) generalised Mertens M_k(x) = sum_{n<=x} mu(n) k^omega(n) ----------
N = 10_000_000
mu = np.ones(N + 1, np.int8); om = np.zeros(N + 1, np.int8)
small = primes_upto(N)
for p in small:
    mu[p::p] *= -1; om[p::p] += 1
    if p*p <= N: mu[p*p::p*p] = 0
mu[0] = 0
xs = np.unique(np.logspace(2, 7, 400).astype(int))
mk = {}
for k in [0.5, 1.0, 2.0, 3.0]:
    w = mu.astype(float)*np.power(k, om.astype(float)); w[0] = 0
    cs = np.cumsum(w)
    runmax = np.maximum.accumulate(np.abs(cs))
    sel = (xs >= 1e5)
    slope = float(np.polyfit(np.log(xs[sel]), np.log(runmax[xs[sel]]), 1)[0])
    entry = dict(x=xs.tolist(), M=cs[xs].tolist(), ratio=(cs[xs]/np.sqrt(xs)).tolist(), slope_runmax=slope)
    if k != int(k):
        # Selberg-Delange main term for non-integer k: G(1)/Gamma(-k) x (log x)^(-k-1), G(1) = prod (1-k/p)(1-1/p)^-k
        Pf = small.astype(float)
        G1 = float(np.exp(np.sum(np.log1p(-k/Pf) - k*np.log1p(-1/Pf))))
        main = G1/math.gamma(-k)*xs*np.log(xs)**(-k-1)
        entry['G1'] = G1; entry['main_ratio'] = (cs[xs]/main).tolist()
    mk[str(k)] = entry
res['mertens_k'] = mk
json.dump(res, open('pair1.json', 'w'), indent=1)
print({k: v for k, v in res.items() if k not in ('mertens_k', 'complex_a_points')})
print({k: v['slope_runmax'] for k, v in mk.items()}, mk['0.5']['main_ratio'][-5:])
print({k: len(v) for k, v in cpts.items()}, res['complex_a_points_side'])
