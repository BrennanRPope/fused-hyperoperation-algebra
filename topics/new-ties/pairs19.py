"""Round 19 (DRAFT): words with four N/K letters against words with one, two, three and four. Computer-assisted.

Words. A reduced word with n + 1 N/K letters is N followed by n items (r, L, M): a gap l^L m^M and then the letter
x = N (r = 19) or K (r = 10, only when M = 0); h = L + M is the gap length; the first item is never (19, 0, 0), since
N N ... = N m ... has one letter fewer. With P = 200^L 20^M = 2^(3L+2M) 5^(2L+M) and A_k = L_1 + ... + L_k,
    T_0 = 1,   T_k = 19^(A_k) + r_k P_k T_(k-1),   psi = T_n / 19^(A_n),   K (the fingerprint) = T_n with 2, 5, 19 removed.
Write T = S K, S = 2^e2 5^e5 19^g. For four letters (n = 3): g <= 3 (F3, by hand, README 23.2) and 2^e2 5^e5 <= 80 (by hand),
so S <= 80 * 19^3 = 548,720. (Two letters: S <= 19; three: S <= 14,440, Rounds 16 and 18.)

Split at an item p (1 <= p <= n): prefix = items before p, sigma = items after p (a known or unknown "suffix"),
Q_p = prod over k >= p of r_k P_k. Then exactly
    T = Q_p T(prefix) + 19^(A_p) T(sigma),   so   T = Q_p T(prefix) (1 + eta),  0 < log(1 + eta) < 0.112 * 10^-h_p.
Prefix step P(p, p'): for two words with the same K, Lambda = log(S' Q_p T(pre) / (S Q'_p' T(pre'))) is a linear form in
  log 19, log 2, log 5, log T(pre), log T(pre') with |Lambda| < 0.112 * 10^-min(h_p, h'_p') and coefficients <= 9X + 20;
  exact LLL (de Weger's lemma) for each known prefix, or Matveev's bound (Bugeaud-Mignotte-Siksek 2006, Thm 9.4), bounds it.
Suffix step S(p, p'): with sigma, sigma' known, T = 19^(A_p) T(sigma) modulo 2^v2(Q_p), so (when v2(Q_p) > e and v5(Q_p) > f,
  e, f the 2- and 5-sizes of T(sigma)) K = 19^(A_p + g_sigma - g) kappa(sigma) modulo 2^(v2(Q_p) - e). Comparing both words,
  19^m kappa(sigma)/kappa(sigma') = 1 modulo 2^min(...), |m| <= 3X + 8, and an exact 2-adic discrete logarithm bounds
  min(h_p, h'_p'). With both suffixes empty this is Round 18's D-step (also 5-adic).
Lemma N (README 23.3): Lambda != 0, and m != 0 when kappa(sigma) = kappa(sigma'), unless the two words are equal; it uses
  Theorems 15.5 to 18.1 for the parts with at most three letters, the prefix lemma and a parity argument.
Size step: 10^n 20^(sum h) <= T <= 1.12 * 19^n * 200^(sum h), and T(other) = S(other) K <= S_max(other) T(short).
Each round splits the pair into "one word is in a short box" (a tree of steps, README 23.4); the size step then caps the
other word. The first bound uses only prefix steps with Matveev's bound and the closed-form D-step. Below the final bounds
every short-box word is compared exactly with every word of the other kind.
Usage: python3 pairs19.py [all|41|42|43|44] [--skip-exact]    Writes pairs19.json.
Time: about 20 to 30 minutes in all on 4 cores, most of it the exact check of four against four. Peak memory: well under
1 GB (the exact check streams the long words and keeps only the short boxes in memory).
"""
import json, math, sys, time, resource
from decimal import Decimal
from fractions import Fraction
from multiprocessing import Pool
from baker16 import matveev_constant, largest_ok
from pairs17 import strip, v, P_of, LN, lattice_bound

THR = "0.112"
S_MAX = {1: 1, 2: 19, 3: 2 ** 3 * 5 * 19 ** 2, 4: 2 ** 4 * 5 * 19 ** 3}     # by number of N/K letters
LLL_MAX = 6000            # use LLL per prefix when there are at most this many prefixes, else Matveev's bound
SUF_MAX = 2_500_000       # use the exact suffix step when there are at most this many suffix pairs
LN10 = math.log(10)

# ---------------------------------------------------------------- words
def items_gap(h):
    return [(19, L, h - L) for L in range(h + 1)] + [(10, h, 0)]
def items_upto(hmax, first=False):
    out = []
    for h in range(hmax + 1):
        for it in items_gap(h):
            if first and it == (19, 0, 0): continue
            out.append(it)
    return out
def count_items(hmax, first=False):
    return sum(h + 2 for h in range(hmax + 1)) - (1 if first else 0)
def T_of(items):
    T = 1; A = 0
    for (r, L, M) in items:
        A += L; T = 19 ** A + r * P_of(L, M) * T
    return T
def Q_of(items):
    q = 1
    for (r, L, M) in items: q *= r * P_of(L, M)
    return q
def K_of(items):
    return strip(T_of(items))
def name(items):
    return "N" + "".join("l" * L + "m" * M + ("N" if r == 19 else "K") for (r, L, M) in items)
def hsum(items):
    return sum(L + M for (_, L, M) in items)
def logT_up(k, s):
    """an upper bound for log T of a word with k items and gap sum s: T <= 1.12 * 19^k * 200^s"""
    return math.log(1.12) + k * math.log(19) + s * math.log(200) + 1e-9

def size_other(n_short, s, n_other):
    """largest Y with 10^n_other 20^Y <= S_max(other) * 1.12 * 19^n_short * 200^s (n = number of items). Exact integers
    for s <= 5000; above that an upper bound from 700-digit logarithms, plus 1 for safety (only early rounds)."""
    if s <= 5000:
        rhs = S_MAX[n_other + 1] * Fraction(112, 100) * 19 ** n_short * 200 ** s
        return largest_ok(lambda Y: 10 ** n_other * 20 ** Y <= rhs, 0, 4 * s + 80)
    c = Decimal(S_MAX[n_other + 1] * 112 * 19 ** n_short) / Decimal(100 * 10 ** n_other)
    return int((c.ln() + Decimal(s) * Decimal(200).ln()) / Decimal(20).ln()) + 1

# ---------------------------------------------------------------- prefix steps
def P_lll(X0, Ts):
    logs = [LN[19], LN[2], LN[5]] + [Decimal(T).ln() for T in Ts]
    h, C, lam = lattice_bound(X0, logs, THR, tries=12)
    return h
def P_matveev(X0, logTs):
    C = matveev_constant(n=3 + len(logTs), alphas=(19, 2, 5)); A = 1.0
    for lt in logTs: A *= max(lt, 0.16)
    return math.floor(C * A * (1 + math.log(X0)) / LN10)

def P10(X0, a, X):
    """min(h_2, h'_1) for a W-prefix of one item with gap <= a (P(2,1)); the same bound serves P(1,2)"""
    n = count_items(a, first=True)
    if n > LLL_MAX:
        return min(P_matveev(X0, [logT_up(1, a)]), X), "Matveev", n
    worst = 0
    for it in items_upto(a, first=True):
        h = P_lll(X0, [T_of([it])])
        if h is None: h = P_matveev(X0, [math.log(T_of([it]))])
        worst = max(worst, h)
        if worst >= X: break
    return min(worst, X), "LLL", n

# ---------------------------------------------------------------- suffix steps
_PW = {}
def dlog19(Kinv, J):
    """lam mod 2^J with 19^lam = Kinv modulo 2^(J+2) (Kinv = 1 or 3 mod 8), bit by bit"""
    mod = 1 << (J + 2)
    if J not in _PW: _PW[J] = [pow(19, 1 << k, mod) for k in range(J)]
    pw = _PW[J]; cur = 1; lam = 0
    for n in range(3, J + 3):
        msk = (1 << n) - 1
        if cur & msk != Kinv & msk:
            cur = cur * pw[n - 3] % mod; lam += 1 << (n - 3)
    return lam
def jmax_ratio(a, b, M):
    """max over integers |m| <= M of v2(19^m a / b - 1), exactly (a, b odd, prime to 19, a != b)"""
    J = max(M.bit_length() + 8, 16)
    while True:
        mod = 1 << (J + 2)
        K = a * pow(b, -1, mod) % mod
        if K % 8 in (5, 7): return 2
        lam = dlog19(pow(K, -1, mod), J)
        t = 0
        while t < J:
            c = lam % (1 << (t + 1)); m = -M + ((c + M) % (1 << (t + 1)))
            if m > M: break
            t += 1
        if t < J: return 1 if t == 0 else 2 + t
        J *= 2
def flog(n, b):
    k = 0
    while b ** (k + 1) <= n: k += 1
    return k
_SD = {}
def sdata(items):
    items = tuple(items)
    if items not in _SD:
        T = T_of(items); Q = Q_of(items)
        _SD[items] = (strip(T), v(T, 2) - v(Q, 2), v(T, 5) - v(Q, 5))
    return _SD[items]
def S_bound(s1, s2, M):
    """bound for min(h_p, h'_p') from known suffixes with data s = (kappa, e - v2 Q_sigma, f - v5 Q_sigma)"""
    k1, e1, f1 = s1; k2, e2, f2 = s2
    if k1 == k2:                                         # m != 0 by Lemma N
        J2 = 2 + flog(M, 2); J5 = 1 + flog(M, 5)
        b2 = max((J2 + e1) // 2, (J2 + e2) // 2, f1, f2, 0)
        b5 = max(J5 + f1, J5 + f2, e1 // 2, e2 // 2, 0)
        return min(b2, b5)
    J = jmax_ratio(k1, k2, M)
    return max((J + e1) // 2, (J + e2) // 2, f1, f2, 0)

def S00(M, X):
    return min(S_bound(sdata(()), sdata(()), M), X)
def S10(c, M, X):
    """one known suffix item with gap <= c against an empty suffix"""
    e0 = sdata(()); worst = 0
    for it in items_upto(c):
        worst = max(worst, S_bound(sdata((it,)), e0, M))
    return min(worst, X)
def S11(c, d, M, X):
    A = items_upto(c); B = items_upto(d); worst = 0
    for s in A:
        ds = sdata((s,))
        for t in B:
            worst = max(worst, S_bound(ds, sdata((t,)), M))
    return min(worst, X)
def S02(c, d, M, X):
    """a known two-item suffix (gaps <= d, then <= c) against an empty suffix"""
    e0 = sdata(()); worst = 0
    for s in items_upto(d):
        for t in items_upto(c):
            T = T_of((s, t)); Q = Q_of((s, t))                   # not cached: there are many
            worst = max(worst, S_bound((strip(T), v(T, 2) - v(Q, 2), v(T, 5) - v(Q, 5)), e0, M))
    return min(worst, X)

def P00(X0, X):
    h = P_lll(X0, [])
    return min(h if h is not None else X, X)

# ---------------------------------------------------------------- a priori (closed forms in ell = 1 + log(9X + 20))
class Poly:
    def __init__(self, c): self.c = [float(x) for x in c]
    @staticmethod
    def of(x): return x if isinstance(x, Poly) else Poly([x])
    def __add__(self, o):
        o = Poly.of(o); n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)])
    __radd__ = __add__
    def __mul__(self, o):
        o = Poly.of(o); r = [0.0] * (len(self.c) + len(o.c) - 1)
        for i, x in enumerate(self.c):
            for j, y in enumerate(o.c): r[i + j] += x * y
        return Poly(r)
    __rmul__ = __mul__
    def at(self, ell): return sum(x * ell ** k for k, x in enumerate(self.c))
    def dominate(self, ell0):
        K = len(self.c) - 1
        return sum(x * ell0 ** (k - K) for k, x in enumerate(self.c)), K
def pmax(*ps):
    n = max(len(p.c) for p in ps)
    return Poly([max((p.c[i] if i < len(p.c) else 0) for p in ps) for i in range(n)])
ELL = Poly([0, 1])
def mv(n, *As):
    """Matveev bound for a P-step with n logarithms; As = upper bounds (Poly) for log T of the prefixes involved"""
    r = Poly([0, matveev_constant(n=n, alphas=(19, 2, 5)) / LN10])
    for A in As: r = r * A
    return r
def logT(k, s):
    return Poly([math.log(1.12) + k * math.log(19)]) + math.log(200) * s
C_D = Poly([1 - 1 / (2 * math.log(2)), 1 / (2 * math.log(2))])     # the D-step: c <= 1 + log2(3X + 8)/2
def size_poly(n_short, s, n_other):
    return (Poly([math.log(S_MAX[n_other + 1] * 1.12 * 19 ** n_short / 10 ** n_other)]) + math.log(200) * s) * (1 / math.log(20))
def solve(F, label):
    ell0 = 1 + math.log(29)
    A, K = F.dominate(ell0)
    g = lambda X: X - A * (1 + math.log(9 * X + 20)) ** K
    lo, hi = 1.0, 1e300
    for _ in range(4000):
        mid = math.sqrt(lo * hi)
        if g(mid) > 0: hi = mid
        else: lo = mid
    Xs = math.ceil(hi); ell = 1 + math.log(9 * Xs + 20)
    assert g(Xs) > 0 and A * K * ell ** (K - 1) / Xs < 1
    return dict(pair=label, A=A, K=K, X_star=Xs, check="X - A (1 + log(9X+20))^K > 0 and increasing for X >= X_star")

def a_priori_42():
    a = mv(3); c = C_D; b = mv(4, logT(1, a))
    sW = a + b + c; sZ = pmax(a, b, c)
    return solve(pmax(a, b, c, size_poly(3, sW, 1), size_poly(1, sZ, 3)), "4v2")
def a_priori_43():
    a = mv(3); c = C_D; b = mv(4, logT(1, a))                          # P(2,1) and P(1,2) alike
    fI = mv(5, logT(1, a), logT(1, b))                                 # P(2,2): prefixes x1 <= a, y1 <= b
    fII = mv(4, logT(2, a + b))                                        # P(3,1): prefix x1 x2
    fIII = mv(5, logT(1, b), logT(1, a))                               # P(2,2): prefixes x1 <= b, y1 <= a
    sW = pmax(a + b + c, a + fI + c, a + b + fII, b + fIII + c)
    sZ = pmax(b + fI, fII + c, a + b, a + fIII, a + c, b + c)
    return solve(pmax(a, b, c, fI, fII, fIII, size_poly(3, sW, 2), size_poly(2, sZ, 3)), "4v3")
def a_priori_44():
    a = mv(3); c = C_D; b = mv(4, logT(1, a))
    f = mv(5, logT(1, a), logT(1, b))                                  # P(2,2)
    g = mv(5, logT(1, a), logT(2, b + f))                              # P(2,3): prefixes x1 and y1 y2
    f1 = mv(4, logT(2, a + b))                                         # P(3,1)
    g1 = mv(5, logT(2, a + b), logT(1, f1))                            # P(3,2)
    g2 = mv(5, logT(2, a + f), logT(1, b))                             # P(3,2)
    s = pmax(a + b + c, a + f + c, a + g + c, b + f + g,               # A
             a + b + f1, a + b + g1, f1 + g1 + c,                      # B1
             b + f + c, a + f + g2, b + g2 + c)                        # B2
    return solve(pmax(a, b, c, f, g, f1, g1, g2, size_poly(3, s, 3)), "4v4")

# ---------------------------------------------------------------- rounds
def round_42(X):
    t0 = time.time(); X0 = 9 * X + 20; M = 3 * X + 8
    a = P00(X0, X); c = S00(M, X); b, meth, nb = P10(X0, a, X)
    boxW = [(a, b, c)]; boxZ = [(max(a, b, c),)]
    sW = a + b + c; sZ = max(a, b, c)
    YZ = size_other(3, sW, 1); YW = size_other(1, sZ, 3)
    Xn = max(a, b, c, min(YZ, X), min(YW, X))
    return dict(X=X, a=a, c=c, b=b, b_method=meth, b_prefixes=nb, short_W_sum=sW, short_V_sum=sZ,
                V_bound_if_W_short=YZ, W_bound_if_V_short=YW, new_X=Xn, seconds=round(time.time() - t0, 1)), Xn, \
        dict(W=boxW, Z=boxZ, YZ=YZ, YW=YW)

def round_43(X):
    t0 = time.time(); X0 = 9 * X + 20; M = 3 * X + 8
    a = P00(X0, X); c = S00(M, X); b, meth, nb = P10(X0, a, X); d = S10(c, M, X)
    boxW = [(a, b, c), (a, d, c), (a, b, d), (b, d, c)]
    boxZ = [(b, d), (b, c), (d, c), (a, b), (a, d), (a, c)]
    sW = max(sum(x) for x in boxW); sZ = max(sum(x) for x in boxZ)
    YZ = size_other(3, sW, 2); YW = size_other(2, sZ, 3)
    Xn = max(a, b, c, d, min(YZ, X), min(YW, X))
    return dict(X=X, a=a, c=c, b=b, b_method=meth, b_prefixes=nb, d=d, short_W_sum=sW, short_V_sum=sZ,
                V_bound_if_W_short=YZ, W_bound_if_V_short=YW, new_X=Xn, seconds=round(time.time() - t0, 1)), Xn, \
        dict(W=boxW, Z=boxZ, YZ=YZ, YW=YW)

def round_44(X):
    t0 = time.time(); X0 = 9 * X + 20; M = 3 * X + 8
    a = P00(X0, X); c = S00(M, X); b, meth, nb = P10(X0, a, X); d = S10(c, M, X)
    n11 = count_items(c) * count_items(d)
    if n11 <= SUF_MAX: f1 = S11(c, d, M, X); m1 = "exact S11 (%d pairs)" % n11
    else: f1 = min(P_matveev(X0, [logT_up(1, a), logT_up(1, b)]), X); m1 = "Matveev P(2,2)"
    if n11 <= SUF_MAX: f2 = S02(c, d, M, X); m2 = "exact S02 (%d suffixes)" % n11
    else: f2 = min(P_matveev(X0, [logT_up(2, a + b)]), X); m2 = "Matveev P(3,1)"
    boxes = [(a, b, c), (a, d, c), (a, f1, c), (b, f1, d),           # A
             (a, b, d), (a, b, f2), (f2, d, c),                        # B1
             (b, d, c), (a, f1, d), (b, f1, c)]                        # B2
    s = max(sum(x) for x in boxes); Y = size_other(3, s, 3)
    Xn = max(a, b, c, d, f1, f2, min(Y, X))
    return dict(X=X, a=a, c=c, b=b, b_method=meth, b_prefixes=nb, d=d, f1=f1, f1_method=m1, f2=f2, f2_method=m2,
                short_sum=s, other_bound=Y, new_X=Xn, seconds=round(time.time() - t0, 1)), Xn, dict(W=boxes, Y=Y)

def run_rounds(fn, X, log):
    rounds = []
    while True:
        rec, Xn, boxes = fn(X); rounds.append(rec); print(log, rec, flush=True)
        if Xn >= X: return rounds, X, boxes
        X = Xn

# ---------------------------------------------------------------- the exact check
def box_words(n, box, first=True):
    """all words with n items, item k having gap <= box[k]"""
    lists = [items_upto(box[0], first=True)] + [items_upto(box[k]) for k in range(1, n)]
    def rec(k, acc):
        if k == n: yield tuple(acc); return
        for it in lists[k]: yield from rec(k + 1, acc + [it])
    yield from rec(0, [])

_D = None
def _stream4(args):
    """stream all four-letter words with first item x1 and gap sum <= B; return clashes with the short-box table"""
    x1, B = args
    clashes = []; n = 0
    r1, L1, M1 = x1; h1 = L1 + M1; T1 = 19 ** L1 + r1 * P_of(L1, M1)
    for x2 in items_upto(B - h1):
        r2, L2, M2 = x2; A2 = L1 + L2; T2 = 19 ** A2 + r2 * P_of(L2, M2) * T1; rest = B - h1 - L2 - M2
        for x3 in items_upto(rest):
            r3, L3, M3 = x3; T = 19 ** (A2 + L3) + r3 * P_of(L3, M3) * T2; n += 1
            if x3 == (19, 0, 0): K = strip(T)
            else:
                K = T
                while K % 19 == 0: K //= 19
            w = _D.get(K)
            if w is not None and w != (x1, x2, x3): clashes.append((name((x1, x2, x3)), name(w)))
    return n, clashes

def exact_check(boxes43, boxes42, boxes44, B4_extra=0, procs=4):
    global _D
    t0 = time.time(); D = {}; inner = []
    def add(w):
        K = K_of(w)
        if K in D and D[K] != w: inner.append((name(w), name(D[K])))
        else: D[K] = w
    add(())                                                               # the one-letter word N (K = 1)
    nb = 0
    for bx in boxes42["W"] + boxes43["W"] + boxes44["W"]:
        for w in box_words(3, bx): add(w); nb += 1
    for bx in boxes42["Z"]:
        for w in box_words(1, bx): add(w); nb += 1
    for bx in boxes43["Z"]:
        for w in box_words(2, bx): add(w); nb += 1
    short_table = len(D)
    # planted words: a disguised m K (= l N) and a leading N N (= N m) must land on their reduced forms
    planted = found = 0; wrong = []
    for K, w in list(D.items()):
        if not w: continue
        for i, (r, L, M) in enumerate(w):
            if r == 19 and L >= 1:
                q = w[:i] + ((10, L - 1, M + 1),) + w[i + 1:]; planted += 1
                if D.get(K_of(q)) == w: found += 1
                else: wrong.append(name(q))
        r, L, M = w[0]
        if M >= 1:
            q = ((19, 0, 0), (r, L, M - 1)) + w[1:]; planted += 1
            if D.get(K_of(q)) == w: found += 1
            else: wrong.append(name(q))
    # stream the long words of each kind up to its size bound
    B2 = max(boxes42["YZ"], 0); B3 = boxes43["YZ"]; B4 = max(boxes42["YW"], boxes43["YW"], boxes44["Y"], 2) + B4_extra
    clashes = list(inner); n2 = n3 = 0
    for w1 in items_upto(B2, first=True):
        n2 += 1; K = K_of((w1,))
        if K in D and D[K] != (w1,): clashes.append((name((w1,)), name(D[K])))
    for h1 in range(B3 + 1):
        for x1 in items_gap(h1):
            if x1 == (19, 0, 0): continue
            for x2 in items_upto(B3 - h1):
                n3 += 1; w = (x1, x2); K = K_of(w)
                if K in D and D[K] != w: clashes.append((name(w), name(D[K])))
    _D = D
    jobs = [(x1, B4) for x1 in items_upto(B4, first=True)]
    jobs.sort(key=lambda j: j[0][1] + j[0][2])
    n4 = 0
    with Pool(procs) as pool:
        for n, cl in pool.imap_unordered(_stream4, jobs, chunksize=1):
            n4 += n; clashes.extend(cl)
    return dict(short_box_words=nb, short_table_distinct=short_table, planted_disguised_words=planted, planted_found=found,
                planted_wrong=wrong[:10], two_letter_streamed_up_to=B2, two_letter_words=n2, three_letter_streamed_up_to=B3,
                three_letter_words=n3, four_letter_streamed_up_to=B4, four_letter_words=n4, clashes=clashes,
                seconds=round(time.time() - t0, 1))

def four_vs_one():
    """K(W) = 1 means T = S <= 548,720, and T >= 1000 * 20^(sum h): so sum h <= 2; check those words directly"""
    Y = largest_ok(lambda Y: 1000 * 20 ** Y <= S_MAX[4], 0, 20)
    words = [w for w in box_words(3, (Y, Y, Y)) if hsum(w) <= Y]
    bad = [name(w) for w in words if K_of(w) == 1]
    return dict(sum_h_bound=Y, words_checked=len(words), with_fingerprint_1=bad)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "all"
    skip = "--skip-exact" in sys.argv
    t0 = time.time(); out = {}
    out["four_vs_one"] = four_vs_one(); print("4v1:", out["four_vs_one"], flush=True)
    res = {}
    for key, ap, fn in (("42", a_priori_42, round_42), ("43", a_priori_43, round_43), ("44", a_priori_44, round_44)):
        if which not in ("all", key): continue
        A = ap(); print("a priori", key, A, flush=True)
        rounds, Xf, boxes = run_rounds(fn, A["X_star"], key)
        res[key] = dict(a_priori=A, rounds=rounds, final_X=Xf, boxes=boxes)
    out["pairs"] = {k: {kk: vv for kk, vv in r.items()} for k, r in res.items()}
    if which == "all" and not skip:
        out["exact_check"] = exact_check(res["43"]["boxes"], res["42"]["boxes"], res["44"]["boxes"])
        print("exact:", {k: v_ for k, v_ in out["exact_check"].items()}, flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_memory_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    if which == "all":
        json.dump(out, open("pairs19.json", "w"), indent=1, default=str)
    print("done", out["seconds"], "s, peak", out["peak_memory_MB"], "MB (main process)")
