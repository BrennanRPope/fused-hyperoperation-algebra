"""Round 7 (DRAFT): the true size of the plane union fractal (README section 11).

Plane maps z -> c + lam (z - c) (turn-and-shrink about a centre c). Maps with the same centre commute; maps with
different centres don't. Two examples:

  ORIG  (round 2): class 1 centre 0: lam = (1+i)/4, 3/10;  class 2 centre 1: i/3;     class 3 centre i: 1/3.
  MOD   (round 7): class 1 centre 0: lam = (1+i)/4, 3/10;  class 2 centre 1: 3i/10;  class 3 centre i: 3/10.

The proof that different words give different maps has three steps (README 11.3):
  step 1  complex ping-pong: the outer class of a word is read off from where it sends a base point (certified here
          with exact rational disc covers);
  step 2  cancellation, which leaves only  a^alpha g' = b^beta h'  with alpha, beta >= 1;
  step 3  MOD only: evaluate at z = 1 3-adically. Every map keeps 3-adic integers, classes 2 and 3 land on 3-adic
          units, and lam_1^alpha lam_2^-beta has 3-adic valuation -beta, so the two sides can't match.

What this script checks (certified claims use exact fractions or interval arithmetic):
  C1  an invariant disc B, and pairwise-disjoint covers of the three class pieces X_1, X_2, X_3 (both examples)
  C2  3-adic facts for MOD: every coefficient is 3-integral; class 2/3 maps send 3-adic integers to 1 (resp. i)
      mod 3; lam_1 is a 3-adic unit and v3(lam_2) = 1
  C3  illustration of C2 on random words: the point g'(1) is a 3-adic unit for every word whose outer syllable is
      class 2 or 3
  C4  exact distinct-map counts against trace-monoid counts, by hashing the maps mod two large primes (a match
      mod p proves the maps are distinct over Q(i), because equal maps stay equal mod p)
  C5  certified D* for both examples (mpmath interval arithmetic)
  C6  certified unconditional lower bound for ORIG from the sub-system with pure a-powers in class 1
  C7  ORIG, the last step without arithmetic: for every pair (d, e) with 1 <= e <= 6 and d >= 1, the number
      rho = lam_1^d lam_2^-e is not a ratio of two points of X_2 u X_3 (exact disc test for small d, a size bound
      for the rest). Words whose class-1 syllables all have at most 6 letters b therefore never collide.
  C8  certified unconditional lower bound for ORIG from that sub-system (class-1 syllables a^alpha b^beta, beta <= 6)
Writes plane.json."""
import json, math, random, itertools
from fractions import Fraction as F
from mpmath import iv
iv.dps = 40

# ---------- Gaussian rationals as pairs of Fractions ----------
def cm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def ca(a, b): return (a[0] + b[0], a[1] + b[1])
def cs(a, b): return (a[0] - b[0], a[1] - b[1])
def n2(a): return a[0]*a[0] + a[1]*a[1]
Z0, ONE, I = (F(0), F(0)), (F(1), F(0)), (F(0), F(1))

def sqrt_up(q):
    """rational upper bound for sqrt(q), q a nonnegative Fraction"""
    q = F(q); D = 10**15
    x = F(math.isqrt(q.numerator * D * D // q.denominator) + 1, D)
    while x * x < q: x += F(1, D)
    return x

SQRT2_UP = F(99, 70)                      # 1.4142857 > sqrt 2
assert SQRT2_UP**2 > 2

EXAMPLES = {
 "ORIG": [(Z0, [(F(1,4), F(1,4)), (F(3,10), F(0))]), (ONE, [(F(0), F(1,3))]), (I, [(F(1,3), F(0))])],
 "MOD":  [(Z0, [(F(1,4), F(1,4)), (F(3,10), F(0))]), (ONE, [(F(0), F(3,10))]), (I, [(F(3,10), F(0))])],
}

def abs_up(lam):
    """rational upper bound for |lam|; exact when n2 is a square"""
    q = n2(lam); r = math.isqrt(q.numerator); s = math.isqrt(q.denominator)
    if r*r == q.numerator and s*s == q.denominator: return F(r, s)
    if q * 8 == 1: return SQRT2_UP / 4    # |(1+i)/4| = sqrt2/4
    return sqrt_up(q)

def class_elements(cls, eps):
    """non-identity elements of one class (a commutative monoid) with |lam| >= eps:
       list of (lam exact, radius-factor upper bound). Also returns the bound used for the tail (< eps)."""
    gens = cls[1]
    ups = [abs_up(l) for l in gens]
    out = []
    def rec(i, lam, up, nonid):
        if i == len(gens):
            if nonid: out.append((lam, up))
            return
        k, l, u = 0, lam, up
        while u >= eps or k == 0:
            rec(i + 1, l, u, nonid or k > 0)
            l = cm(l, gens[i]); u = u * ups[i]; k += 1
            if u < eps: break
    rec(0, ONE, F(1), False)
    return [(l, u) for l, u in out if u >= eps]

def image_disc(c, lam, lamup, disc):
    """image of disc (centre z, radius r) under z -> c + lam (z - c)"""
    z, r = disc
    return (ca(c, cm(lam, cs(z, c))), lamup * r)

def disjoint(d1, d2):
    (z1, r1), (z2, r2) = d1, d2
    return n2(cs(z1, z2)) > (r1 + r2)**2

def pingpong(name, cB, R, eps1, eps2):
    ex = EXAMPLES[name]
    gens = [(c, l) for c, ls in ex for l in ls]
    # B invariant: |g(cB) - cB| + |lam| R <= R for each generator
    inv = []
    for c, l in gens:
        gz = ca(c, cm(l, cs(cB, c)))
        inv.append(sqrt_up(n2(cs(gz, cB))) + abs_up(l) * R <= R)
    B = (cB, R)
    # B2: covers the union of h(B) over all non-identity h; tails near each centre
    B2 = []
    for cls in ex:
        c = cls[0]
        for lam, up in class_elements(cls, eps1):
            B2.append(image_disc(c, lam, up, B))
        B2.append((c, eps1 * (sqrt_up(n2(cs(cB, c))) + R)))       # tail: elements with |lam| < eps1
    # X_j: union over non-identity g in class j of g(B2)
    X = []
    for cls in ex:
        c = cls[0]; cover = []
        for lam, up in class_elements(cls, eps2):
            cover += [image_disc(c, lam, up, d) for d in B2]
        rho = max(sqrt_up(n2(cs(z, c))) + r for z, r in B2)
        cover.append((c, eps2 * rho))
        X.append(cover)
    bad = 0; mind = {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        gap = min(math.sqrt(float(n2(cs(a[0], b[0])))) - float(a[1] + b[1]) for a in X[i] for b in X[j])
        mind["%d-%d" % (i+1, j+1)] = gap
        bad += sum(1 for a in X[i] for b in X[j] if not disjoint(a, b))
    return dict(B=[[float(cB[0]), float(cB[1])], float(R)], invariant=all(inv), B2_discs=len(B2),
                X_discs=[len(x) for x in X], disjoint_pairs_failed=bad, certified=(all(inv) and bad == 0),
                min_gap_float=mind), X, B2

# ---------- maps of words; syllables outermost first ----------
def gen_maps(name):
    ex = EXAMPLES[name]
    return [[(l, cm(cs(ONE, l), c)) for l in ls] for c, ls in ex]   # z -> l z + (1-l) c

def compose(f, g):  # f after g
    return (cm(f[0], g[0]), ca(cm(f[0], g[1]), f[1]))
def apply(f, z): return ca(cm(f[0], z), f[1])

def rand_word(nsyl, rng):
    """reduced word in M = N^2 * N * N as syllables [(class, exponents)], outermost first"""
    w = []; prev = None
    for _ in range(nsyl):
        k = rng.choice([j for j in range(3) if j != prev])
        e = (rng.randint(0, 3), rng.randint(0, 3)) if k == 0 else (rng.randint(1, 4),)
        if k == 0 and e == (0, 0): e = (1, 0)
        w.append((k, e)); prev = k
    return w

def word_map(name, w):
    G = gen_maps(name); f = (ONE, Z0)
    for k, e in w:
        for gi, n in enumerate(e):
            for _ in range(n): f = compose(f, G[k][gi])
    return f

def v3(q):
    q = F(q)
    if q == 0: return math.inf
    n, d, k = q.numerator, q.denominator, 0
    while n % 3 == 0: n //= 3; k += 1
    while d % 3 == 0: d //= 3; k -= 1
    return k
def v3c(z):
    """3-adic valuation of a Gaussian rational (3 is prime in Z[i], so it is the min over parts)"""
    return min(v3(z[0]), v3(z[1]))

def pw(z, n):
    r = ONE
    for _ in range(n): r = cm(r, z)
    return r

def three_adic_constants():
    ex = EXAMPLES["MOD"]; G = gen_maps("MOD")
    integral = all(v3c(l) >= 0 and v3c(t) >= 0 for cls in G for l, t in cls)
    # class 2/3 maps: g(z) - c = lam (z - c) with v3(lam) >= 1, so g(z) = c mod 3 for 3-integral z
    contract23 = all(v3c(l) >= 1 for c, ls in ex[1:] for l in ls)
    centres_units = v3c(ONE) == 0 and v3c(I) == 0 and v3c(cs(ONE, I)) == 0
    lam1, lam2 = ex[0][1]
    return dict(all_coefficients_3_integral=integral, class23_multipliers_v3_ge_1=contract23,
                centres_1_and_i_are_3_adic_units=centres_units, v3_lam1=v3c(lam1), v3_lam2=v3c(lam2),
                v3_of_rho_is_minus_beta=all(v3c(cm(pw(lam1, al), pw((F(10,3), F(0)), be))) == -be
                                            for al in range(0, 12) for be in range(1, 12)))

def three_adic_illustration(n=3000, seed=7):
    rng = random.Random(seed); vals = {}
    for _ in range(n):
        w = rand_word(rng.randint(1, 6), rng)
        if w[0][0] == 0: w = w[1:]
        if not w: continue
        val = v3c(apply(word_map("MOD", w), ONE))
        vals[val] = vals.get(val, 0) + 1
    o = {}
    rng = random.Random(seed)
    for _ in range(n):
        w = rand_word(rng.randint(1, 6), rng)
        if w[0][0] == 0: w = w[1:]
        if not w: continue
        val = v3c(apply(word_map("ORIG", w), ONE))
        o[val] = o.get(val, 0) + 1
    return dict(MOD=dict(sorted(vals.items())), ORIG=dict(sorted(o.items())))

# ---------- C4: distinct maps vs trace counts, hashed mod two primes ----------
def trace_counts(pattern, n):
    """number of elements of length 1..n in the trace monoid N^{k1} * N^{k2} * ... via the series
       1/P = 1 + sum_j ((1-t)^{k_j} - 1)"""
    N = n + 1
    inv = [0]*N; inv[0] = 1
    for k in pattern:
        poly = [0]*N
        for j in range(min(k, n) + 1): poly[j] = math.comb(k, j) * (-1)**j
        poly[0] -= 1
        for j in range(N): inv[j] += poly[j]
    P = [0]*N; P[0] = 1
    for m in range(1, N): P[m] = -sum(inv[j] * P[m-j] for j in range(1, m+1))
    return P[1:]

PRIMES = [2147483629, 2147483549]   # both = 1 mod 4 checked below, so i exists mod p
def modp_maps(name, p):
    assert p % 4 == 1
    ii = next(pow(g, (p-1)//4, p) for g in range(2, 100) if pow(g, (p-1)//2, p) == p-1)
    def rq(x): return x.numerator * pow(x.denominator, -1, p) % p
    def red(q): return (rq(q[0]) + ii * rq(q[1])) % p
    return [(red(l), red(t)) for cls in gen_maps(name) for l, t in cls]

def distinct_counts(name, n):
    import numpy as np
    gens = [modp_maps(name, p) for p in PRIMES]
    counts = []
    # represent the level as int64 arrays per prime: (lam, t)
    L = [(np.array([1], dtype=np.int64), np.array([0], dtype=np.int64)) for _ in PRIMES]
    for step in range(n):
        newL = []
        for pi, p in enumerate(PRIMES):
            ls, ts = L[pi]; nl, nt = [], []
            for gl, gt in gens[pi]:
                nl.append(gl * ls % p); nt.append((gl * ts + gt) % p)   # g after current
            newL.append((np.concatenate(nl), np.concatenate(nt)))
        key = np.stack([newL[0][0], newL[0][1], newL[1][0], newL[1][1]], axis=1)
        _, idx = np.unique(key, axis=0, return_index=True)
        L = [(newL[pi][0][idx], newL[pi][1][idx]) for pi in range(len(PRIMES))]
        counts.append(int(len(idx)))
    return counts

# ---------- C5, C6: certified numbers ----------
R1 = iv.sqrt(iv.mpf(2)) / 4
def fD(s, ex):
    rs = {"ORIG": [[R1, iv.mpf(3)/10], [iv.mpf(1)/3], [iv.mpf(1)/3]],
          "MOD": [[R1, iv.mpf(3)/10], [iv.mpf(3)/10], [iv.mpf(3)/10]]}[ex]
    tot = iv.mpf(1)
    for sl in rs:
        pr = iv.mpf(1)
        for r in sl: pr *= 1 - r**s
        tot += pr - 1
    return tot
def certify_root(f, lo, hi, decreasing=False, tol=1e-6):
    """bisection keeping a certified bracket: sign(f(lo)) and sign(f(hi)) proved opposite"""
    sgn = (lambda x: x.b < 0) if not decreasing else (lambda x: x.a > 0)
    for _ in range(80):
        mid = (lo + hi) / 2
        val = f(iv.mpf(mid))
        if sgn(val): lo = mid
        elif (val.a > 0 if not decreasing else val.b < 0): hi = mid
        else: break
        if hi - lo < tol: break
    assert sgn(f(iv.mpf(lo))) and (f(iv.mpf(hi)).a > 0 if not decreasing else f(iv.mpf(hi)).b < 0)
    return lo, hi
def chain_f(s):
    t1 = R1**s; t3 = iv.mpf(3)**(-s)
    return (t1 / (1 - t1)) * (2*t3 / (1 - 2*t3)) - 1

def sqrt_down(q):
    q = F(q)
    return F(0) if q == 0 else q / sqrt_up(q)

def ratio_excluded(rho, rho_up, W):
    """True if rho * D_b and D_a are disjoint for all discs D_a, D_b of W (so rho is no ratio w_a / w_b)"""
    return all(n2(cs(cm(rho, zb), za)) > (rho_up * rb + ra)**2 for za, ra in W for zb, rb in W)

def last_step_orig(emax=6, eps=F(1, 60)):
    res, X, B2 = pingpong("ORIG", (F(3,5), F(3,8)), F(87,100), eps, eps)
    W = X[1] + X[2]
    m_lo = min(sqrt_down(n2(z)) - r for z, r in W); M_up = max(sqrt_up(n2(z)) + r for z, r in W)
    lam1, lam2 = EXAMPLES["ORIG"][0][1]; inv2 = (F(10, 3), F(0)); r1up = SQRT2_UP / 4
    checked = []; tails = {}
    for e in range(1, emax + 1):
        d = 1
        while True:
            rup = r1up**d * F(10, 3)**e
            if rup * M_up < m_lo:          # size bound: |rho w_b| < |w_a| for all a, b; also for every larger d
                tails[e] = d; break
            rho = ONE
            for _ in range(d): rho = cm(rho, lam1)
            for _ in range(e): rho = cm(rho, inv2)
            ok = ratio_excluded(rho, rup, W)
            checked.append((d, e, ok)); d += 1
    return dict(pingpong_certified=res["certified"], m_lo=float(m_lo), M_up=float(M_up),
                all_excluded=all(ok for _, _, ok in checked), pairs_checked_exactly=len(checked),
                size_bound_from_d=tails)

def sub_f(s, B0):
    t1 = R1**s; t2 = (iv.mpf(3)/10)**s; t3 = iv.mpf(3)**(-s)
    SU = (1/(1 - t1)) * (1 - t2**(B0 + 1)) / (1 - t2) - 1
    return SU * (2*t3 / (1 - 2*t3)) - 1

if __name__ == "__main__":
    OUT = {}
    for name, cB, R in [("ORIG", (F(3,5), F(3,8)), F(87,100)), ("MOD", (F(29,50), F(19,50)), F(17,20))]:
        res, X, B2 = pingpong(name, cB, R, F(1, 60), F(1, 60))
        OUT["C1_" + name] = res
        print("C1", name, res)
    OUT["C2"] = three_adic_constants(); print("C2", OUT["C2"])
    OUT["C3"] = three_adic_illustration(); print("C3", OUT["C3"])
    tc = trace_counts([2, 1, 1], 11)
    OUT["C4"] = {"trace": tc}
    for name in ["ORIG", "MOD"]:
        dc = distinct_counts(name, 11)
        OUT["C4"][name] = dict(distinct=dc, match=dc == tc)
        print("C4", name, dc == tc, dc[-3:], tc[-3:])
    for name in ["ORIG", "MOD"]:
        lo, hi = certify_root(lambda s: fD(s, name), 1.0, 1.4)
        OUT["C5_" + name] = [lo, hi]; print("C5", name, lo, hi)
    lo, hi = certify_root(chain_f, 0.9, 1.2, decreasing=True)
    OUT["C6_ORIG_chain"] = [lo, hi]; print("C6", lo, hi)
    OUT["C7"] = last_step_orig(); print("C7", OUT["C7"])
    lo, hi = certify_root(lambda s: sub_f(s, 6), 1.1, 1.3, decreasing=True)
    OUT["C8_ORIG_beta_le_6"] = [lo, hi]; print("C8", lo, hi)
    json.dump(OUT, open("plane.json", "w"), indent=1, default=str)
