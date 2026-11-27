"""Theorem 17.1: a reduced three-letter word (not starting NN) never shares a
fingerprint with a reduced two-letter word."""
import math, time, sys
import mpmath as mp
from multiprocessing import Pool
from common import *

print('=== Theorem 17.1 ===', '(dlog excludes m=L\'-g\', i.e. D=0)' if '--exclude-D0' in sys.argv else '')
LOGB = lambda Hs: math.log(3 * Hs + 5)          # B <= 3 Hs + 5 for all forms used
C3 = matveev_C(3, [LOG2, LOG5, LOG19])
C4 = matveev_C(4, [1.0, LOG2, LOG5, LOG19])      # times log K
C5 = matveev_C(5, [1.0, 1.0, LOG2, LOG5, LOG19])  # times log K * log T_pre
print('Matveev constants: C3=%.4e  C4=%.4e (x logK)  C5=%.4e (x logK x logTpre)' % (C3, C4, C5))
LOG_TMAX0 = math.log(400.71)      # T3 <= 400.71 * 200^(H1+H2);  T' <= 20.9 * 200^H'
G3 = lambda h: 2.035 + 1.7995 * h  # g3 <= log19 T3

# ---------- lemma checks ----------
# Lambda_1 != 0 and Type-B D != 0 both reduce to: no (r2,r') with 1+v19(r2)-v19(r') = f(r2)-f(r'), f=v2-2v5
f_ = lambda r: vp(r, 2) - 2 * (1 if r == 10 else 0)
v19 = lambda r: 1 if r == 19 else 0
for r2 in (19, 10):
    for rp in (19, 10):
        assert 1 + v19(r2) - v19(rp) != f_(r2) - f_(rp)
print('r-casework for Lambda1 != 0 / type-B D != 0: ok')
# Type B generic valuations, g3 bound, eps3 bound on a sample
for w in three_letter_words(40, 6):
    T = T3(*w)
    L1, M1, x2, L2, M2, x3 = w
    r2, r3 = rval(x2), rval(x3)
    assert math.log(T, 19) <= G3(L1 + M1 + L2 + M2) + 1e-9
    eps3 = F(T, r3 * r2 * P(L2, M2) * P(L1, M1)) - 1
    assert 0 < eps3 <= F(11, 100) * F(1, 10 ** (L1 + M1))
    if not is_typeB(w):
        assert T % 2 and T % 5
print('T3 facts (type A odd & prime to 5, eps3 <= 0.11*10^-H1, g3 bound) ok on sample H1<=40,H2<=6')

def X1_matveev(Hs):
    return math.floor((C3 * (1 + LOGB(Hs)) + math.log(0.11)) / LOG10)

def X2(Hs):   # type A, 5-adic: min(H2,H') <= 1 + log5|D|, |D| <= 2.04 + 2.8 Hs
    return math.floor(1 + math.log(2.04 + 2.8 * Hs, 5))

def X2_2adic(Hs):
    return math.floor(1 + math.log2(2.04 + 2.8 * Hs) / 2)

def X2B(Hs):  # type B: min(H1,H') <= 2 + log5|D|
    return math.floor(2 + math.log(2.04 + 2.8 * Hs, 5))

def apriori(Hs):
    X1 = X1_matveev(Hs)
    x2 = X2(Hs); x2b = min(X1, X2B(Hs))
    Y = max(X1, x2, x2b)
    logK = math.log(20.9) + Y * math.log(200)
    X4 = math.floor((C4 * logK * (1 + LOGB(Hs)) + math.log(0.11)) / LOG10)
    logTpre = math.log(20.9) + X4 * math.log(200)
    X5 = math.floor((C5 * logK * logTpre * (1 + LOGB(Hs)) + math.log(0.1)) / LOG10)
    vshort = max(X4 + X5, Y)
    wshort_h = max(X1 + x2, x2b)
    wshort = max(wshort_h, math.floor((math.log(1.9) + LOG_TMAX0) / LOG10 + wshort_h * math.log10(200)))
    return max(vshort, wshort), dict(X1=X1, X2=x2, X2B=x2b, Y=Y, logK=logK, X4=X4, X5=X5)

Hs = 10 ** 200
for _ in range(100):
    nH, d = apriori(Hs)
    if nH >= Hs: break
    Hs = nH
Hs0 = Hs
x = Hs0 + 1
for _ in range(3000):
    assert apriori(x)[0] < x
    x = int(x * 1.05) + 1
print('a priori (Matveev 3/4/5 logs): Hs=max(H1+H2,H\') <= %.4e' % Hs0)
print('   pieces:', {k: ('%.4e' % v) for k, v in d.items()}, ' (2-adic version of X2 would be %d)' % X2_2adic(Hs0))

# ---------- per-K work (V short) ----------
EXCL = '--exclude-D0' in sys.argv

def per_K(args):
    K, X0, H2start, m0 = args
    fn = lambda: [mp.log(K), mp.log(2), mp.log(5), mp.log(19)]
    X4, info = reduce_bound(fn, 4, X0, math.log(0.11), LOG10)
    # H2 (type A) via exact 2-adic discrete log of K to base 19:
    # K = 19^m mod 2^{v2(r3 P2)}, m = L1+L2-g3, |m| <= 2.04+1.8(H1+H2), v2(r3 P2) >= 2 H2
    H2 = H2start
    while True:
        Mb = math.ceil(2.04 + 1.8 * (X4 + H2))
        bnd = max_v2_K19m(K, Mb, exclude=(m0 if EXCL else None))
        new = bnd // 2
        if new >= H2:
            break
        H2 = new
    return K, X4, H2, info['Cexp'], info['log10c']


def one_round(Hs, rnd, pool):
    t = time.time()
    X0 = 3 * Hs + 5
    # step 1: shared 3-log form, |log(1+Lambda1)| <= 1.1236*0.11*10^-min(H1,H')
    fn3 = lambda: [mp.log(2), mp.log(5), mp.log(19)]
    X1, info1 = reduce_bound(fn3, 3, X0, math.log(0.11 * 1.1236), LOG10)
    x2 = X2(Hs); x2b = min(X1, X2B(Hs))
    Y = max(X1, x2, x2b)
    Vs = two_letter_words(Y)
    Ks = [kappa2(*v)[0] for v in Vs]
    m0s = [v[0] - (1 if (v[2] == 'N' and v[0] >= 1) else 0) for v in Vs]   # L'-g'
    assert len(set(Ks)) == len(Ks)
    res = pool.map(per_K, [(K, X0, Hs, m0) for K, m0 in zip(Ks, m0s)], chunksize=16)
    X4max = max(r[1] for r in res); X5max = max(r[2] for r in res)
    Cexps = sorted(set(r[3] for r in res)); minc = min(r[4] for r in res)
    vshort = max(X4max + X5max, Y)
    wshort_h = max(X1 + x2, x2b)
    wshort = max(wshort_h, math.floor((math.log(1.9) + LOG_TMAX0) / LOG10 + wshort_h * math.log10(200)))
    newHs = max(vshort, wshort)
    print('round %d: Hs<=%d X0=%d' % (rnd, Hs, X0))
    print('   step1 3-log LLL: C=10^%d log10 c=%.2f -> min(H1,H\')<=%d' % (info1['Cexp'], info1['log10c'], X1))
    print('   step2 (5-adic, |D|<=%.4g): typeA min(H2,H\')<=%d ; typeB min(H1,H\')<=%d' % (2.04 + 2.8 * Hs, x2, x2b))
    print('   V-short: H\'<=%d -> %d fingerprints K, 4-log LLL each (C=10^%s, min log10 c=%.2f): H1<=%d ; 2-adic dlog: H2<=%d'
          % (Y, len(Ks), Cexps, minc, X4max, X5max))
    print('   W-short: H1<=%d H2<=%d (typeB H1<=%d) -> H\'<=%d' % (X1, x2, x2b, wshort))
    print('   new Hs<=%d   (%.1fs)' % (newHs, time.time() - t))
    return newHs, dict(X1=X1, X2=x2, X2B=x2b, Y=Y, X4max=X4max, X5max=X5max, res=res)


if __name__ == '__main__':
    pool = Pool(4)
    Hs = Hs0
    rnd = 0
    hist = []
    while True:
        rnd += 1
        nH, d = one_round(Hs, rnd, pool)
        hist.append(d)
        if nH >= Hs:
            break
        Hs = nH
    d = hist[-1]
    A = max(d['X4max'], d['X1'], d['X2B'])
    Bb = max(d['X5max'], d['X2'])
    print('final box for exact check: H1<=%d, H2<=%d' % (A, Bb))
    # exact check
    t = time.time()
    Ws = list(three_letter_words(A, Bb))
    Kw = []
    Kmax = 0
    for w in Ws:
        T = T3(*w)
        K = strip(T)
        assert K != 1
        Kw.append((K, w))
        Kmax = max(Kmax, K)
    Hcap = math.floor(math.log10(1.9 * Kmax)) + 1
    fp = {}
    for v in two_letter_words(Hcap):
        fp[kappa2(*v)[0]] = v
    hits = [(w, fp[K]) for K, w in Kw if K in fp]
    print('exact check: %d three-letter words (H1<=%d,H2<=%d) vs %d two-letter words (H\'<=%d, enough since H\'<log10(1.9K)): collisions %s  (%.1fs)'
          % (len(Ws), A, Bb, len(fp), Hcap, hits, time.time() - t))
    print('Theorem 17.1 holds by this computation:', not hits)
