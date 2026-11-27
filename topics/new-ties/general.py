"""Round 8 (DRAFT): a general no-collision test, and the families it covers (README section 12).

Setting: maps z -> c + lam (z - c) on the line (real coefficients) or in the plane (Gaussian-rational coefficients),
grouped into classes by centre. Words are the trace monoid N^{k_1} * ... * N^{k_N}.

Theorem P8 (no-collision test). Different words give different maps if
  (A) separation: some set Y is sent into itself by every map, and the class pieces X_i = U_{g in class i, g != id} g(Y)
      are pairwise disjoint; and
  (B) arithmetic: for every class i with two or more maps, the valuation vectors of its multipliers at the primes
      that are "good for i" have no mixed-sign integer relation. A prime p is good for i when every centre and
      multiplier has no p in its denominator, every multiplier of every OTHER class is divisible by p, and the other
      centres differ from c_i mod p.

This script
  G1  implements the test (A by exact-fraction disc covers with any number of classes; B by valuations and a kernel
      computation) and runs it on the round-7 examples: MOD passes, ORIG fails B (as it must: no prime works there);
  G2  checks the explicit sufficient condition for (A): all |lam| <= eps with eps < delta / (4 D + delta), where delta
      and D are the smallest and largest distances between centres (Lemma A');
  G3  samples three families at random and certifies (A) and (B) for each member, then confirms no collisions up to
      length 7 by hashing mod two primes, and gives each member's certified D*:
        F1 plane, 3 classes, centres 0, 1, i; class 1 = {unit-type, 3 x ...}, classes 2, 3 one multiplier 3 x ...
        F2 plane, 4 classes, centres 0, 1, i, 1+i  (all distinct mod 3)
        F3 line, 3 classes, centres 0, 1, 2; class 1 two maps, classes 2, 3 one map each, prime 3
  G4  a line member with two maps in two classes, each class with its own good prime (3 and 5)
  G5  Corollary P9: members of the dense family for 14 patterns of class sizes (line and plane), 3 each: the explicit
      separation bound, the arithmetic test, irreducibility in the plane, and no collisions to length 5 or 6
Writes general.json."""
import json, math, random, itertools
from fractions import Fraction as F
import numpy as np
import sympy
from mpmath import iv
import plane as P
iv.dps = 30
cm, ca, cs, n2, ONE, Z0 = P.cm, P.ca, P.cs, P.n2, P.ONE, P.Z0

def vq(q, p):
    q = F(q)
    if q == 0: return math.inf
    n, d, k = q.numerator, q.denominator, 0
    while n % p == 0: n //= p; k += 1
    while d % p == 0: d //= p; k -= 1
    return k

def val(z, p, field):
    """valuation of a Gaussian rational at a prime of Z[i] lying over p (inert p = 3 mod 4, or p = 2 via 1+i);
       on the line (field 'Q'), the usual p-adic valuation"""
    if field == "Q":
        assert z[1] == 0; return vq(z[0], p)
    if z == Z0: return math.inf
    if p == 2: return vq(n2(z), 2)                  # v_{1+i} = v_2(norm)
    assert p % 4 == 3
    return min(vq(z[0], p), vq(z[1], p))           # p stays prime in Z[i]

def candidate_primes(classes, field):
    ps = set()
    for c, ls in classes:
        for l in ls:
            for q in (n2(l),): ps |= set(sympy.factorint(q.numerator)) | set(sympy.factorint(q.denominator))
    for (c1, _), (c2, _) in itertools.combinations(classes, 2):
        q = n2(cs(c1, c2)); ps |= set(sympy.factorint(q.numerator)) | set(sympy.factorint(q.denominator))
    if field == "Qi": ps = {p for p in ps if p == 2 or p % 4 == 3}
    return sorted(ps)

def good_primes(classes, i, field):
    out = []
    for p in candidate_primes(classes, field):
        ok = all(val(c, p, field) >= 0 and all(val(l, p, field) >= 0 for l in ls) for c, ls in classes)
        ok &= all(val(l, p, field) >= 1 for j, (c, ls) in enumerate(classes) if j != i for l in ls)
        ok &= all(val(cs(c, classes[i][0]), p, field) == 0 for j, (c, ls) in enumerate(classes) if j != i)
        if ok: out.append(p)
    return out

def no_mixed_relation(V):
    """V: k x m integer matrix (valuations of the k multipliers at m primes). True if no integer vector w with both
       positive and negative entries has w V = 0. (Kernel of dimension 0, or 1 spanned by a one-signed vector;
       dimension >= 2 always contains a mixed-sign vector, README 12.2.)"""
    k = len(V)
    if not V or not V[0]: return k <= 1
    M = sympy.Matrix(V).T                        # m x k;  w V = 0  <=>  M w = 0
    ns = M.nullspace()
    if len(ns) == 0: return True
    if len(ns) >= 2: return False
    v = list(ns[0])
    return all(x >= 0 for x in v) or all(x <= 0 for x in v)

def test_B(classes, field):
    rep = {}
    for i, (c, ls) in enumerate(classes):
        if len(ls) < 2: rep[i] = dict(needed=False); continue
        gp = good_primes(classes, i, field)
        V = [[val(l, p, field) for p in gp] for l in ls]
        rep[i] = dict(needed=True, good_primes=gp, valuations=V, ok=no_mixed_relation(V))
    return all(r["ok"] for r in rep.values() if r["needed"]), rep

# ---------- (A): certified disc covers, any number of classes ----------
def abs_up(l): return P.sqrt_up(n2(l))
def find_disc(classes):
    from scipy.optimize import minimize
    def Rf(v):
        b = complex(*v); R = 0
        for c, ls in classes:
            cc = complex(float(c[0]), float(c[1]))
            for l in ls:
                ll = complex(float(l[0]), float(l[1])); R = max(R, abs(cc + ll*(b - cc) - b) / (1 - abs(ll)))
        return R
    cen = np.mean([[float(c[0]), float(c[1])] for c, _ in classes], axis=0)
    r = minimize(Rf, cen, method="Nelder-Mead")
    b = (F(r.x[0]).limit_denominator(200), F(r.x[1]).limit_denominator(200))
    R = F(Rf([float(b[0]), float(b[1])]) * 1.02 + 1e-6).limit_denominator(1000) + F(1, 1000)
    return b, R

def class_elements(ls, eps):
    ups = [abs_up(l) for l in ls]; out = []
    def rec(i, lam, up, nonid):
        if i == len(ls):
            if nonid: out.append((lam, up))
            return
        l, u = lam, up
        while True:
            rec(i + 1, l, u, nonid or l != lam)
            l = cm(l, ls[i]); u = u * ups[i]
            if u < eps: break
    rec(0, ONE, F(1), False)
    return out

def test_A(classes, eps=F(1, 40), level=2):
    b, R = find_disc(classes)
    inv = all(P.sqrt_up(n2(cs(ca(c, cm(l, cs(b, c))), b))) + abs_up(l) * R <= R for c, ls in classes for l in ls)
    Y = [(b, R)]
    if level == 2:                                # Y2 = B intersect (union of h(B)); covered by discs + tails
        Y2 = []
        for c, ls in classes:
            for lam, up in class_elements(ls, eps):
                Y2.append((ca(c, cm(lam, cs(b, c))), up * R))
            Y2.append((c, eps * (P.sqrt_up(n2(cs(b, c))) + R)))
        Y = Y2
    X = []
    for c, ls in classes:
        cover = []
        for lam, up in class_elements(ls, eps):
            cover += [(ca(c, cm(lam, cs(z, c))), up * r) for z, r in Y]
        cover.append((c, eps * max(P.sqrt_up(n2(cs(z, c))) + r for z, r in Y)))
        X.append(cover)
    ok = inv
    for i, j in itertools.combinations(range(len(classes)), 2):
        if not ok: break
        ok = all(n2(cs(a[0], bb[0])) > (a[1] + bb[1])**2 for a in X[i] for bb in X[j])
    return ok, dict(disc=[[float(b[0]), float(b[1])], float(R)], invariant=inv)

def lemma_A_explicit(classes):
    """Lemma A': if every |lam| <= eps and eps < delta/(4D + delta) then (A) holds with Y a disc about c_1."""
    cs_ = [complex(float(c[0]), float(c[1])) for c, _ in classes]
    ds = [abs(a - b) for a, b in itertools.combinations(cs_, 2)]
    delta, D = min(ds), max(ds)
    eps = max(abs(complex(float(l[0]), float(l[1]))) for _, ls in classes for l in ls)
    return eps < delta / (4 * D + delta), dict(eps=eps, bound=delta / (4 * D + delta))

# ---------- collisions to length n (hash mod two primes), trace counts, certified D* ----------
def distinct_counts(classes, n):
    gens_all = []
    for p in P.PRIMES:
        ii = next(pow(g, (p-1)//4, p) for g in range(2, 100) if pow(g, (p-1)//2, p) == p-1)
        rq = lambda x: x.numerator * pow(x.denominator, -1, p) % p
        red = lambda q: (rq(q[0]) + ii * rq(q[1])) % p
        gens_all.append([(red(l), red(cm(cs(ONE, l), c))) for c, ls in classes for l in ls])
    L = [(np.array([1], dtype=np.int64), np.array([0], dtype=np.int64)) for _ in P.PRIMES]; counts = []
    for _ in range(n):
        newL = []
        for pi, p in enumerate(P.PRIMES):
            ls_, ts = L[pi]
            newL.append((np.concatenate([gl * ls_ % p for gl, gt in gens_all[pi]]),
                         np.concatenate([(gl * ts + gt) % p for gl, gt in gens_all[pi]])))
        key = np.stack([newL[0][0], newL[0][1], newL[1][0], newL[1][1]], axis=1)
        _, idx = np.unique(key, axis=0, return_index=True)
        L = [(newL[pi][0][idx], newL[pi][1][idx]) for pi in range(2)]; counts.append(int(len(idx)))
    return counts

def certified_Dstar(classes):
    rs = [[iv.sqrt(iv.mpf(n2(l).numerator) / n2(l).denominator) for l in ls] for _, ls in classes]
    def f(s):
        tot = iv.mpf(1)
        for sl in rs:
            pr = iv.mpf(1)
            for r in sl: pr *= 1 - r**s
            tot += pr - 1
        return tot
    return P.certify_root(f, 0.01, 2.5)

def hochman_irreducible(classes, field):
    """plane: some first-return piece turns by an angle other than 0 or 180 degrees (enough: some multiplier is not
       real, or a product of two multipliers from different classes is not real)"""
    if field == "Q": return True
    ls = [l for _, ls in classes for l in ls]
    return any(l[1] != 0 for l in ls)

def full_test(classes, field, n_check=7):
    okA, infoA = test_A(classes)
    okB, infoB = test_B(classes, field)
    pattern = [len(ls) for _, ls in classes]
    res = dict(A=okA, B=okB, A_info=infoA, B_info={str(k): v for k, v in infoB.items()},
               irreducible=hochman_irreducible(classes, field), explicit_A=lemma_A_explicit(classes)[0])
    if okA and okB:
        dc = distinct_counts(classes, n_check); tc = P.trace_counts(pattern, n_check)
        res.update(no_collisions_to=n_check if dc == tc else None, Dstar=list(certified_Dstar(classes)))
    return res

# ---------- random families ----------
def rnd_frac(rng, lo, hi, den=60):
    while True:
        x = F(rng.randint(1, den), den) * (hi - lo) + lo
        if x != 0: return x

def fam_F1(rng):
    # class 1: one multiplier with v3 = 0 (turning), one divisible by 3; classes 2, 3: one multiplier divisible by 3
    while True:
        a = (F(rng.choice([1, 2, 4, 5, 7, 8]), rng.choice([10, 11, 13, 14, 16])) * rng.choice([1, -1]),
             F(rng.choice([1, 2, 4, 5, 7]), rng.choice([10, 11, 13, 14, 16])))
        if val(a, 3, "Qi") == 0 and math.sqrt(float(n2(a))) < 0.45: break
    b = (F(3 * rng.randint(1, 4), rng.choice([20, 25, 28, 31, 35])), F(0))
    c2 = (F(0), F(3 * rng.randint(1, 4), rng.choice([20, 25, 28, 31, 35])))
    c3 = (F(3 * rng.randint(1, 4), rng.choice([20, 25, 28, 31, 35])), F(0))
    return [(Z0, [a, b]), (ONE, [c2]), (P.I, [c3])]

def fam_F2(rng):
    base = fam_F1(rng)
    c4 = (F(3 * rng.randint(1, 3), rng.choice([25, 28, 31, 35])), F(3 * rng.randint(1, 3), rng.choice([25, 28, 31, 35])))
    return base + [((F(1), F(1)), [c4])]

def fam_F3(rng):
    a = (F(rng.choice([1, 2, 4, 5]), rng.choice([10, 11, 13, 14, 16])), F(0))
    b = (F(3 * rng.randint(1, 3), rng.choice([20, 25, 28, 31])), F(0))
    c2 = (F(3 * rng.randint(1, 3), rng.choice([20, 25, 28, 31])), F(0))
    c3 = (F(3 * rng.randint(1, 3), rng.choice([20, 25, 28, 31])), F(0))
    return [(Z0, [a, b]), (ONE, [c2]), ((F(2), F(0)), [c3])]

# ---------- Corollary P9: a dense family for every pattern of class sizes ----------
LINE_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
PLANE_PRIMES = [3, 7, 11, 19, 23, 31, 43, 47, 59, 67]          # inert in Z[i] (p = 3 mod 4)
PLANE_CENTRES = [Z0, ONE, P.I, (F(1), F(1))]

def corollary_member(pattern, field, rng):
    """Build a member of the Corollary P9 family: class i with k_i maps gets k_i - 1 primes of its own;
       lam_{i,j} (j < k_i) carries its own prime once, lam_{i,k_i} carries none of class i's primes, and every
       multiplier carries all primes of the other classes. Numerators get a random cofactor prime to all chosen primes,
       denominators are powers of 2, and every |lam| is pushed under delta/(4D + delta)."""
    N = len(pattern)
    pool = LINE_PRIMES if field == "Q" else PLANE_PRIMES
    pool = [p for p in pool if p >= N]
    it = iter(pool); own = [[next(it) for _ in range(k - 1)] for k in pattern]
    allp = [p for o in own for p in o]
    centres = [(F(i), F(0)) for i in range(N)] if field == "Q" else PLANE_CENTRES[:N]
    cc = [complex(float(c[0]), float(c[1])) for c in centres]
    ds = [abs(a - b) for a, b in itertools.combinations(cc, 2)]
    bound = min(ds) / (4 * max(ds) + min(ds))
    def cofactor():
        while True:
            if field == "Q":
                u = (rng.randint(1, 9) * rng.choice([1, -1]), 0)
            else:
                u = (rng.randint(-5, 5), rng.randint(-5, 5))
            nrm = u[0]**2 + u[1]**2
            if nrm and all(nrm % p for p in allp) and (field == "Q" or u[1] != 0 or rng.random() < 0.3): return u
    classes = []
    for i, k in enumerate(pattern):
        other = math.prod(p for jj, o in enumerate(own) if jj != i for p in o)
        ls = []
        for j in range(k):
            m = other * (own[i][j] if j < k - 1 else 1)
            u = cofactor(); num = (m * u[0], m * u[1])
            size = math.hypot(*num); n = math.ceil(math.log2(size / bound)) + rng.randint(0, 2)
            while math.hypot(*num) / 2**n >= bound: n += 1
            ls.append((F(num[0], 2**n), F(num[1], 2**n)))
        classes.append((centres[i], ls))
    return classes, own

if __name__ == "__main__":
    OUT = {}
    for name in ["MOD", "ORIG"]:
        okB, info = test_B(P.EXAMPLES[name], "Qi")
        OUT["G1_" + name] = dict(B=okB, info={str(k): v for k, v in info.items()})
        print("G1", name, "B passes:", okB, info[0])
    OUT["G2"] = {}
    ex = [(Z0, [(F(1, 10), F(1, 10)), (F(3, 25), F(0))]), (ONE, [(F(0), F(3, 25))]), (P.I, [(F(3, 25), F(0))])]
    okE, infoE = lemma_A_explicit(ex); okA, _ = test_A(ex)
    OUT["G2"] = dict(example="classes {(1+i)/10, 3/25}, {3i/25}, {3/25}", explicit=okE, info=infoE, certified_A=okA,
                     B=test_B(ex, "Qi")[0])
    print("G2", OUT["G2"])
    for fam, gen, field in [("F1", fam_F1, "Qi"), ("F2", fam_F2, "Qi"), ("F3", fam_F3, "Q")]:
        rng = random.Random({"F1": 1, "F2": 2, "F3": 3}[fam]); rows = []
        for t in range(40):
            cl = gen(rng); r = full_test(cl, field, n_check=7 if fam != "F2" else 6)
            r["classes"] = [[[str(c[0]), str(c[1])], [[str(l[0]), str(l[1])] for l in ls]] for c, ls in cl]
            rows.append(r)
        npass = sum(1 for r in rows if r["A"] and r["B"])
        nB = sum(1 for r in rows if r["B"])
        ncoll = sum(1 for r in rows if r.get("no_collisions_to"))
        OUT["G3_" + fam] = dict(members=len(rows), pass_B=nB, pass_A_and_B=npass, no_collisions_confirmed=ncoll, rows=rows)
        print("G3", fam, "members", len(rows), "B", nB, "A&B", npass, "collision-free confirmed", ncoll,
              "D* range", [round(r["Dstar"][0], 4) for r in rows if "Dstar" in r][:6])
    # G4: two maps in two classes, each class with its own good prime (3 for class 1, 5 for class 2), on the line
    ex4 = [(Z0, [(F(5, 56), F(0)), (F(15, 196), F(0))]),       # v3 = (0, 1); every multiplier divisible by 5
           (ONE, [(F(3, 56), F(0)), (F(15, 196), F(0))]),      # v5 = (0, 1); every multiplier divisible by 3
           ((F(2), F(0)), [(F(15, 196), F(0))])]               # divisible by 15
    r4 = full_test(ex4, "Q", n_check=7)
    OUT["G4"] = dict(example="line, centres 0, 1, 2; classes {5/56, 15/196}, {3/56, 15/196}, {15/196}",
                     A=r4["A"], B=r4["B"], explicit_A=r4["explicit_A"], B_info=r4["B_info"],
                     no_collisions_to=r4.get("no_collisions_to"), Dstar=r4.get("Dstar"))
    print("G4", {k: v for k, v in OUT["G4"].items() if k != "B_info"},
          {k: v.get("good_primes") for k, v in r4["B_info"].items()})
    # G5: Corollary P9 members for many patterns, line and plane
    rng = random.Random(5); rows = []
    patterns_line = [(2, 1, 1), (2, 2, 1), (3, 1, 1), (3, 2, 1), (2, 2, 2), (2, 1, 1, 1), (3, 2, 1, 1), (4, 1, 1)]
    patterns_plane = [(2, 1, 1), (2, 2, 1), (3, 1, 1), (2, 2, 2), (2, 1, 1, 1), (3, 2, 1, 1)]
    for field, pats in [("Q", patterns_line), ("Qi", patterns_plane)]:
        for pat in pats:
            for t in range(3):
                cl, own = corollary_member(pat, field, rng)
                n_check = 6 if sum(pat) <= 5 else 5
                okB, infoB = test_B(cl, field); okE, infoE = lemma_A_explicit(cl)
                dc = distinct_counts(cl, n_check); tc = P.trace_counts(list(pat), n_check)
                rows.append(dict(field=field, pattern=list(pat), own_primes=own, explicit_A=okE, eps=infoE["eps"],
                                 bound=infoE["bound"], B=okB, irreducible=hochman_irreducible(cl, field),
                                 no_collisions_to=n_check if dc == tc else None,
                                 Dstar=list(certified_Dstar(cl)),
                                 classes=[[[str(c[0]), str(c[1])], [[str(l[0]), str(l[1])] for l in ls]] for c, ls in cl]))
    OUT["G5"] = dict(members=len(rows), explicit_A=sum(r["explicit_A"] for r in rows), B=sum(r["B"] for r in rows),
                     collision_free=sum(1 for r in rows if r["no_collisions_to"]),
                     plane_irreducible=sum(1 for r in rows if r["field"] == "Qi" and r["irreducible"]),
                     plane_members=sum(1 for r in rows if r["field"] == "Qi"), rows=rows)
    print("G5", {k: v for k, v in OUT["G5"].items() if k != "rows"},
          "D* range", round(min(r["Dstar"][0] for r in rows), 3), round(max(r["Dstar"][1] for r in rows), 3))
    json.dump(OUT, open("general.json", "w"), indent=1, default=str)
