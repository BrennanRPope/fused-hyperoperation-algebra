"""Round 6, Euler-side pairs (u = p^-s per prime). Dials: k homothety (1 - k u)^(-1/k); t odds (1 - t u)^-1 = Z_t^t;
s height (ln_s = E^-s(X - 1)); d curvature (local factor 1 + k u + d u^2); q order (orderings weighted q^inversions);
r similarity ratio (u -> rho^s u, i.e. primes rescaled to p/rho).
Sections: t+d and d+r (root lines per prime), s+d (height detector), k+q (identity and leading power), k+r and t+r
(complex homothety), s+t/s+q/s+r (height collapse), d+q (prime-fractal dimension of quadratic fields),
k/t/q with the power dials (multiplicativity), light pass on face power (Wick rotation)."""
import json, math, numpy as np, mpmath as mp, sympy as sp
res = {}
mp.mp.dps = 25

def primes_upto(n):
    sv = np.ones(n + 1, bool); sv[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    return np.nonzero(sv)[0]

# ---------------- curve 37a: y^2 + y = x^3 - x, the genus-1 point used on the dial plane ----------------
P3 = primes_upto(3000)
def ap_37a(p):
    if p == 2:
        cnt = sum(1 for x in range(2) for y in range(2) if (y*y + y - x**3 + x) % 2 == 0)
    else:
        x = np.arange(p, dtype=np.int64); disc = (4*(x**3 - x) + 1) % p
        leg = np.array([0 if d == 0 else (1 if pow(int(d), (p - 1)//2, p) == 1 else -1) for d in disc])
        cnt = int(np.sum(1 + leg))
    return p - cnt
AP = {int(p): ap_37a(int(p)) for p in P3}
good = [p for p in AP if p != 37]
res['curve_37a_ap_first'] = {p: AP[p] for p in list(AP)[:10]}

# ---------------- t+d: odds scaling on the dial plane ----------------
kk, dd, tt, uu = sp.symbols('k d t u')
Lt = sp.expand(1 + kk*(tt*uu) + dd*(tt*uu)**2)
k_t, d_t = Lt.coeff(uu, 1), Lt.coeff(uu, 2)
res['td_disc_scaling'] = str(sp.factor(k_t**2 - 4*d_t))           # t^2 (k^2 - 4d)
res['td_parabola'] = str(sp.simplify(d_t/k_t**2 - dd/kk**2))       # 0: the point moves on d = c k^2
def root_lines(scale_fn):
    """sigma_p where the local factor 1 - a_p X + p X^2 vanishes, X = scale(s) p^-s"""
    return {p: scale_fn(p) for p in good}
td = {}
for t in [0.8, 1.25]:
    sig = np.array([0.5 + math.log(t)/math.log(p) for p in good])
    td[str(t)] = dict(min=float(sig.min()), max=float(sig.max()), at_p2=float(sig[0]), at_p2999=float(sig[-1]))
res['td_root_lines'] = td
# check numerically one prime: roots of 1 - t a_p u + t^2 p u^2 have |u| = 1/(t sqrt p)
p0 = 101; t0 = 1.25; rts = np.roots([t0**2*p0, -t0*AP[p0], 1])
res['td_check_p101'] = dict(abs_roots=[float(abs(r)) for r in rts], predicted=1/(t0*math.sqrt(p0)))

# ---------------- d+r: ratio scaling on the dial plane ----------------
dr = {}
for rho in [0.8, 1.25]:
    sig = np.array([0.5*math.log(p)/math.log(p/rho) for p in good if p > rho])
    dr[str(rho)] = dict(min=float(sig.min()), max=float(sig.max()))
res['dr_root_lines'] = dr

# ---------------- s+d: the height detector on a degree-2 local factor ----------------
# h*(u) = L L''/L'^2 for L = 1 + k u + d u^2;  with x = u + k/(2d), w^2 = Delta/(4 d^2):  h* = 1/2 - w^2/(2 x^2)
x_, w2 = sp.symbols('x w2')
Lx = dd*(x_**2 - w2)
res['sd_hstar_formula'] = str(sp.simplify(Lx*sp.diff(Lx, x_, 2)/sp.diff(Lx, x_)**2))
grid = []
for k in np.linspace(-3, 3, 41):
    for d in np.linspace(-3, 3, 41):
        if abs(d) < 1e-9: continue
        u = np.linspace(-20, 20, 40001); L = 1 + k*u + d*u*u; L1 = k + 2*d*u; L2 = 2*d
        ok = np.abs(L1) > 1e-9
        fires = bool(np.max(L[ok]*L2/L1[ok]**2) > 0.5 + 1e-9)
        grid.append((fires, k*k - 4*d < 0))
res['sd_detector_matches_delta_negative'] = dict(points=len(grid), agree=sum(a == b for a, b in grid))
hasse = [AP[p]**2 - 4*p < 0 for p in good]
res['sd_37a'] = dict(primes=len(good), detector_fires=int(sum(hasse)), note='fires at every good prime = Hasse bound |a_p| < 2 sqrt p')

# ---------------- necklace machinery: log of prod_p (1 - c p^-s)^-1 = sum_N M(c,N) log zeta(N s) ----------------
def mobius(n):
    f = sp.factorint(n)
    return 0 if any(e > 1 for e in f.values()) else (-1)**len(f)
MU = {n: mobius(n) for n in range(1, 200)}
def necklace(c, N):
    return sum(MU[N//d]*mp.power(c, d) for d in range(1, N + 1) if N % d == 0)/N
def log_odds(c, s, Nmax=90, skip_first=False):
    tot = 0
    for N in range(2 if skip_first else 1, Nmax + 1):
        tot += necklace(c, N)*mp.log(mp.zeta(N*s))
    return tot

# ---------------- k+q: order dial on homothety letters ----------------
# Z_{k,q} = [ sum_n k^Omega W_q(n) n^-s ]^(1/k) ;  inner = (q;q)_inf sum_m q^m/(q;q)_m prod_j Odds_{k q^(m+j)}
def qpoch(q, n):
    out = 1.0
    for i in range(1, n + 1): out *= 1 - q**i
    return out
def qpoch_inf(q):
    out = 1.0; i = 1
    while q**i > 1e-18: out *= 1 - q**i; i += 1
    return out
NMAX = 400000
spf = np.zeros(NMAX + 1, np.int64)
for p in primes_upto(NMAX):
    blk = spf[p::p]; blk[blk == 0] = p
def exps(n):
    e = []
    while n > 1:
        p = spf[n]; c = 0
        while n % p == 0: n //= p; c += 1
        e.append(c)
    return e
EXPS = [None, []] + [exps(n) for n in range(2, NMAX + 1)]
def direct_kq(k, q, s, N):
    tot = 0.0
    for n in range(1, N + 1):
        e = EXPS[n]; Om = sum(e)
        w = qpoch(q, Om)
        for ei in e: w /= qpoch(q, ei)
        tot += k**Om*w*n**-s
    return tot
P7 = primes_upto(10_000_000).astype(float)
def odds_direct(c, s):
    return math.exp(-np.sum(np.log1p(-c*P7**-s)))
def identity_kq(k, q, s):
    tot = 0.0
    for m in range(0, 80):
        if q**m < 1e-17: break
        prod = 1.0; j = 0
        while k*q**(m + j) > 1e-17:
            prod *= odds_direct(k*q**(m + j), s); j += 1
        tot += q**m/qpoch(q, m)*prod
    return qpoch_inf(q)*tot
kq = {}
for k, q in [(0.6, 0.3), (1.3, 0.3), (1.3, 0.6)]:
    rhs = identity_kq(k, q, 2.0)
    lhs = [direct_kq(k, q, 2.0, N) for N in (NMAX//4, NMAX//2, NMAX)]
    kq[f'k={k},q={q}'] = dict(identity=rhs, direct_N=[NMAX//4, NMAX//2, NMAX], direct=lhs, gaps=[rhs - x for x in lhs])
    print('k+q', k, q, rhs, lhs, flush=True)
res['kq_identity_s2'] = kq
# leading power near s = 1: log Z_{k,q} / log zeta -> 1/(1-q)
def inner_kq_log(k, q, s):
    tot = 0
    for m in range(0, 60):
        if q**m < 1e-15: break
        lp = 0; j = 0
        while k*q**(m + j) > 1e-15:
            lp += log_odds(k*q**(m + j), s, Nmax=60); j += 1
        tot += q**m/qpoch(q, m)*mp.exp(lp)
    return mp.log(qpoch_inf(q)*tot)
lead = {}
for k, q in [(0.6, 0.3), (1.3, 0.6)]:
    rows = []
    for eps in [1e-2, 1e-4, 1e-6, 1e-8]:
        s = 1 + mp.mpf(eps)
        rows.append([eps, float(inner_kq_log(k, q, s)/k/mp.log(mp.zeta(s)))])
    lead[f'k={k},q={q}'] = dict(rows=rows, predicted=1/(1 - q))
res['kq_leading_power'] = lead
res['kq_wall'] = {str(k): math.log(k, 2) for k in [1.3, 1.5, 2.0]}

# ---------------- k+r and t+r: homothety at the complex value k rho^s ----------------
rho1 = mp.zetazero(1)
def log_Zkr(k, rho, s):
    c = k*mp.power(rho, s)
    return mp.log(mp.zeta(s)) + log_odds(c, s, Nmax=90, skip_first=True)/c
kr = {}
for k, rho in [(0.8, 1.1), (1.2, 0.9), (1.5, 1.3)]:
    wall = math.log(k)/math.log(2/rho) if k > 1 else None
    row = dict(wall=wall)
    if wall is None or wall < 0.5:
        H = mp.exp(log_odds(k*mp.power(rho, rho1), rho1, Nmax=90, skip_first=True)/(k*mp.power(rho, rho1)))
        row['H_at_first_zero_abs'] = float(abs(H))
        # Z = zeta * H vanishes at the first zero; check by values just off it
        row['abs_Z_at_rho1_pm_1e-6'] = [float(abs(mp.exp(log_Zkr(k, rho, rho1 + d)))) for d in (1e-6, -1e-6j)]
        row['abs_zeta_prime_times_1e-6'] = float(abs(mp.zeta(rho1, derivative=1))*1e-6*abs(H))
    kr[f'k={k},rho={rho}'] = row
res['kr'] = kr
# t+r: leading exponent t rho^s is complex; at the zeros its real part sets zero (>0) or pole (<0)
gam = [float(mp.zetazero(n).imag) for n in range(1, 21)]
tr = {}
for rho in [1.1, 1.5, 2.0]:
    re_exp = [float((mp.power(rho, mp.mpc(0.5, g))).real) for g in gam]
    tr[str(rho)] = dict(zeros_kept=sum(1 for x in re_exp if x > 0), turned_to_poles=sum(1 for x in re_exp if x < 0))
res['tr_first20'] = tr

# ---------------- height collapse (s with t, q, r) ----------------
# ln_s(F) is singular where F = X_s; X_s = 0 only at s = 1 (round 5 table). Example with simple zeros: F = Z_{k rho^s}.
Xs = {'0.8': -0.14189434629030684, '0.95': -0.031641082794077846}
coll = {}
k, rho = 0.8, 1.1
F = lambda s: mp.exp(log_Zkr(k, rho, s))
for sh, X in Xs.items():
    try:
        a = mp.findroot(lambda s: F(s) - X, rho1 + 0.01)
        coll[sh] = dict(a_point=[float(a.real), float(a.imag)], off_line=float(a.real - 0.5), residual=float(abs(F(a) - X)))
    except Exception as e:
        coll[sh] = str(e)
res['height_collapse_on_kr'] = coll

# ---------------- d+q: prime-fractal dimension of quadratic fields (order dial at q = 1) ----------------
def kron(D, p):
    return int(sp.jacobi_symbol(D % p, p)) if p != 2 else (0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1))
PP = primes_upto(5_000_000)
def norms_K(D):
    out = []
    for p in PP:
        p = int(p)
        if p == 2:
            ch = 0 if D % 4 == 0 or D % 2 == 0 else (1 if D % 8 == 1 else -1)
        else:
            ch = 0 if D % p == 0 else int(sp.legendre_symbol(D % p, p))
        if ch == 1: out += [p, p]
        elif ch == 0: out.append(p)
        else: out.append(p*p)
    return np.array(out, float)
def dimension(norms, split_density):
    X = float(PP[-1])
    def f(D):
        tail = split_density*float(mp.e1((D - 1)*math.log(X)))   # sum_{p > X} p^-D ~ int x^-D/log x
        return float(np.sum(norms**-D)) + tail - 1
    lo, hi = 1.01, 3.0
    for _ in range(80):
        mid = 0.5*(lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return mid
dq = {}
dq['Q'] = dimension(PP.astype(float), 1.0)
for name, D in [('Q(i)', -4), ('Q(sqrt-3)', -3), ('Q(sqrt-5)', -20)]:
    dq[name] = dimension(norms_K(D), 1.0)
res['dq_prime_fractal_dimension'] = dq
print('dq', dq, flush=True)

# ---------------- k, t, q with the power dials: is there an Euler product? ----------------
def rep_counts(form, M, nmax):
    a, c = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1))
    v = form(a.astype(np.int64), c.astype(np.int64)).ravel()
    v = v[(v > 0) & (v <= nmax)]
    return np.bincount(v, minlength=nmax + 1)
def mult_defect(cnt, w, nmax):
    f = cnt/w; bad = tot = 0
    idx = np.nonzero(cnt)[0]; idx = idx[idx > 1]
    for i, m in enumerate(idx):
        for n in idx[i + 1:]:
            if m*n > nmax: break
            if math.gcd(int(m), int(n)) == 1:
                tot += 1; bad += abs(f[m]*f[n] - f[m*n]) > 1e-9
    return dict(coprime_pairs=tot, failures=int(bad))
NM = 200000
mult = {}
mult['square_p2 (a^2+b^2)'] = mult_defect(rep_counts(lambda a, b: a*a + b*b, 450, NM), 4, NM)
mult['square_p4 (a^4+b^4)'] = mult_defect(rep_counts(lambda a, b: a**4 + b**4, 25, NM), 4, NM)
mult['square_p6 (a^6+b^6)'] = mult_defect(rep_counts(lambda a, b: a**6 + b**6, 8, NM), 4, NM)
mult['hex_p2 (a^2+ac+c^2)'] = mult_defect(rep_counts(lambda a, c: a*a + a*c + c*c, 520, NM), 6, NM)
mult['hex_p6 ((a+c)^6+a^6+c^6)/2'] = mult_defect(rep_counts(lambda a, c: ((a + c)**6 + a**6 + c**6)//2, 8, NM), 6, NM)
res['multiplicativity'] = mult
a_, c_ = sp.symbols('a c')
res['hex_p4_is_square_of_p2'] = str(sp.expand((a_ + c_)**4 + a_**4 + c_**4 - 2*(a_**2 + a_*c_ + c_**2)**2))   # 0
# design rule: mean of cos^p over the n directions of a regular 2n-gon is constant iff p even and p <= 2n - 2
th = np.linspace(0, np.pi, 1001)
design = {}
for ngon, ndir in [('square', 2), ('hexagon', 3), ('octagon', 4)]:
    design[ngon] = [p for p in [2, 3, 4, 5, 6, 8, 10] if np.ptp(np.mean([np.abs(np.cos(th + j*np.pi/ndir))**p for j in range(ndir)], 0)) < 1e-12]
res['face_power_round_at'] = design

# ---------------- light pass on face power: Wick rotation y -> i y ----------------
xs_, ys_ = sp.symbols('x y', real=True)
ang = [sp.pi/6 + sp.pi*j/3 for j in range(6)]
light = {}
for p in [2, 4, 6, 8]:
    form = sum((sp.cos(a)*xs_ + sp.sin(a)*sp.I*ys_)**p for a in ang)
    form = sp.simplify(sp.expand(form))
    real = sp.simplify(sp.im(sp.expand(form))) == 0
    quad = sp.expand(xs_**2 - ys_**2)
    k_pow = p//2
    c0 = sp.simplify(form.subs({xs_: 1, ys_: 0}))
    is_power = sp.simplify(form - c0*quad**k_pow) == 0
    light[str(p)] = dict(form=str(sp.factor(form)), real=bool(real), equals_power_of_x2_minus_y2=bool(is_power))
res['light_wick_hexagon'] = light
# in lattice coordinates the Wick-rotated hexagon form a^2 + ac + c^2 becomes 2a^2 + 2ac - c^2 (discriminant 12, Q(sqrt 3))
res['light_lattice_form'] = dict(form='2a^2 + 2ac - c^2', discriminant=12)

json.dump(res, open('euler_pairs.json', 'w'), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k not in ('curve_37a_ap_first',)}, indent=1, default=str)[:6000])
