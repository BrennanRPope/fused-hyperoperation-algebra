"""Round 4: the barrier chain as composed maps, and two dials that deform it.
(A) letters: f_p(z) = p^(-s) z + c (stretch-and-shift, order-dial convention). c = 0: all share centre 0,
    they commute, the chain is the Euler product. c = 1: centre of f_p is the Euler factor 1/(1-p^-s);
    order gap f_p o f_q - f_q o f_p = c*(p^-s - q^-s) (a wedge), never 0 for 0 < Re s.
(B) the order dial q (quantum plane yx = q xy): weigh each ordering of a factorization by q^(inversions).
    q = 0: one ordering per n, the Euler product. q = 1: every ordering, 1/(1 - sum_p p^-s) (free words).
(C) the link dial (k_p, d_p): factor 1/(1 + k u + d u^2), u = p^-s. Variance of log|L| on sigma = 1/2 is
    1/2 sum_p sum_m |alpha^m+beta^m|^2 p^(-m)/m^2; the m = 1 term is k_p^2/(2p), so the chain is bounded
    on the line only if sum k_p^2/p converges (Rankin-Selberg says it never does for an L-function). -> barrier.json"""
import numpy as np, mpmath as mp, json, sympy as sp
from itertools import product
res = {}
# ---------- (A) ----------
s = mp.mpc(0.5, 14.134725)
f = lambda p, c: (lambda z: p**(-s)*z + c)
for c in [0, 1]:
    g23 = f(2, c)(f(3, c)(0)); g32 = f(3, c)(f(2, c)(0))
    res['gap_c%d' % c] = complex(g23-g32)
res['gap_formula'] = complex(2**(-s)-3**(-s))
res['centre_is_euler_factor'] = complex(1/(1-2**(-s)) - (1/(1-2**(-s))))  # centre of f_2 at c=1: z = 2^-s z + 1
# commute only where 2^-s = 3^-s: s = 2 pi i k / log(3/2)
res['commute_points_imag'] = [2*np.pi*k/np.log(1.5) for k in range(1, 4)]
# ---------- (B) order dial ----------
def qint(n, q): return sum(q**j for j in range(n))
def qfact(n, q):
    r = 1
    for j in range(1, n+1): r *= qint(j, q)
    return r
def qmult(a, q): return qfact(sum(a), q)/np.prod([qfact(x, q) for x in a])
def Zq(sv, q, P, M=40):
    tot = 0
    for a in product(range(M), repeat=len(P)):
        if sum(a) > M: continue
        tot += qmult(a, q)*np.prod([p**(-ai*sv) for p, ai in zip(P, a)])
    return tot
P = [2, 3]; sv = 2.0
res['Zq_q0'] = Zq(sv, 0, P); res['euler_23'] = 1/((1-2**-sv)*(1-3**-sv))
res['Zq_q1'] = Zq(sv, 1.0, P); res['free_23'] = 1/(1-2**-sv-3**-sv)
# closed form for two primes: sum_a x^a / prod_{j=0}^{a} (1 - q^j y)
x_, y_ = 2**-sv, 3**-sv
def closed(q, A=200):
    tot = 0
    for a in range(A):
        den = 1
        for j in range(a+1): den *= (1-q**j*y_)
        tot += x_**a/den
    return tot
res['Zq_two_prime_check'] = [(q, Zq(sv, q, P), closed(q)) for q in [0.0, 0.3, 0.7, 1.0]]
# walls: q = 1 for all primes: sigma0 with P(sigma0) = 1 (prime zeta)
primezeta = lambda x: float(mp.primezeta(x))
res['sigma0_all_primes'] = float(mp.findroot(lambda x: mp.primezeta(x)-1, 1.4))
res['sigma0_23'] = float(mp.findroot(lambda x: 2**-x+3**-x-1, 0.8))
# growth of the summed weights A_q(x) = sum_{n<=x} w_q(n), w_q(n) = q-multinomial of n's exponents
X = 2*10**6
spf = np.zeros(X+1, dtype=np.int64)
for i in range(2, X+1):
    if spf[i] == 0: spf[i::i][spf[i::i] == 0] = i
def exps(n):
    e = []
    while n > 1:
        p = spf[n]; k = 0
        while n % p == 0: n //= p; k += 1
        e.append(k)
    return e
E = [exps(n) for n in range(1, X+1)]
growth = {}
xs = np.unique(np.logspace(3, np.log10(X), 25).astype(int))
for q in [0.0, 0.25, 0.5, 0.75]:
    w = np.array([qmult(e, q) if e else 1.0 for e in E])
    A = np.cumsum(w)[xs-1]
    # fit A(x) ~ C x (log x)^beta  on the top half
    yv = np.log(A/xs); lv = np.log(np.log(xs)); h = len(xs)//2
    beta = np.polyfit(lv[h:], yv[h:], 1)[0]
    growth[str(q)] = dict(beta_fit=float(beta), beta_selberg_delange=q/(1-q), A=A.tolist())
res['growth'] = growth; res['growth_xs'] = xs.tolist()
# ---------- (C) link dial on sigma = 1/2 ----------
pr = np.array(list(sp.primerange(2, 10**6)), dtype=float)
cum = lambda v: np.cumsum(v)
checkpoints = [10**2, 10**3, 10**4, 10**5, 10**6-1]
idx = [np.searchsorted(pr, c, side='right')-1 for c in checkpoints]
def var_link(alpha, beta, sigma=0.5, Mmax=6):
    v = 0
    for m in range(1, Mmax+1):
        v = v + 0.5*np.abs(alpha**m+beta**m)**2*pr**(-2*m*sigma)/m**2
    return cum(v)[idx]
rng = np.random.default_rng(1)
# Sato-Tate angles via rejection sampling of density (2/pi) sin^2
th = []
while len(th) < len(pr):
    c = rng.uniform(0, np.pi, 2*len(pr)); u = rng.uniform(0, 1, 2*len(pr))
    th += list(c[u < np.sin(c)**2])
th = np.array(th[:len(pr)])
fam = {
    'zeta (k=-1, d=0)': (np.ones_like(pr), np.zeros_like(pr)),
    'elliptic, fixed angle pi/3 (k=-1, d=1)': (np.exp(1j*np.pi/3)*np.ones_like(pr), np.exp(-1j*np.pi/3)*np.ones_like(pr)),
    'elliptic, Sato-Tate angles (d=1)': (np.exp(1j*th), np.exp(-1j*th)),
    'k = 0, d = 1 (alpha = i, beta = -i)': (1j*np.ones_like(pr), -1j*np.ones_like(pr)),
    'hyperbolic (k=-1.2, d=0.2)': ((1.2+np.sqrt(1.44-0.8))/2*np.ones_like(pr), (1.2-np.sqrt(1.44-0.8))/2*np.ones_like(pr)),
}
res['link_dial'] = {'checkpoints': checkpoints, 'loglog': [float(np.log(np.log(c))) for c in checkpoints],
                    'var': {k: var_link(a, b).tolist() for k, (a, b) in fam.items()},
                    'k2_mean': {k: float(np.mean(np.abs(a+b)**2)) for k, (a, b) in fam.items()}}
json.dump(res, open('barrier.json', 'w'), default=str)
for k, v in res.items():
    if k not in ('growth', 'link_dial', 'growth_xs'): print(k, v)
print({q: (g['beta_fit'], g['beta_selberg_delange']) for q, g in growth.items()})
for k, v in res['link_dial']['var'].items(): print(k, [round(x, 3) for x in v])
