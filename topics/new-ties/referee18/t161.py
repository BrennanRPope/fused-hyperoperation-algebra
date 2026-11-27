"""Theorem 16.1: no reduced three-letter word has fingerprint 1."""
import math, sys
import mpmath as mp
from common import *

print('=== Theorem 16.1 ===')
# ---------------- lemma checks ----------------
for E in range(1, 3001):
    v2 = vp(19 ** E - 1, 2); v5 = vp(19 ** E - 1, 5) if (19 ** E - 1) % 5 == 0 else 0
    assert v2 == (1 if E % 2 else 2 + vp(E, 2))
    assert v5 == (0 if E % 2 else 1 + vp(E, 5))
print('v2(19^E-1), v5(19^E-1) formulas verified for E<=3000')
# T_pre = 19^L1 + r2 P1 (u1 x2 not NN): odd, prime to 5, v19 = [x2=N and L1>=1]
for H1 in range(0, 80):
    for L1 in range(H1 + 1):
        M1 = H1 - L1
        for x2 in 'NK':
            if x2 == 'K' and M1 > 0: continue
            if x2 == 'N' and H1 == 0: continue
            Tp = 19 ** L1 + rval(x2) * P(L1, M1)
            assert Tp % 2 == 1 and Tp % 5 != 0
            assert vp(Tp, 19) == (1 if (x2 == 'N' and L1 >= 1) else 0)
print('T_pre facts verified for H1<80')

# ---------------- Case I (x3=N, u2 empty): psi = 20 + 19 r2 A1 ----------------
# claim a<=2, b<=1 except NmNN and NKN; check on all words H1<=300
exc = []
for H1 in range(0, 301):
    for L1 in range(H1 + 1):
        M1 = H1 - L1
        for x2 in 'NK':
            if x2 == 'K' and M1 > 0: continue
            if x2 == 'N' and H1 == 0: continue
            w = (L1, M1, x2, 0, 0, 'N')
            T = T3(*w)
            a, b = vp(T, 2), vp(T, 5)
            if (a, b) != (2, 1):
                exc.append((w, a, b, strip(T)))
print('Case I exceptions to (a,b)=(2,1) for H1<=300:', exc)
# elementary argument: generic case psi = 20*19^c  =>  20(19^c-1)19^L1 = 19 r2 P1
#   2-adic: 2+v2(19^c-1) = v2(r2)+3L1+2M1 ; 5-adic: 1+v5(19^c-1) = v5(r2)+2L1+M1
#   size:   19^c - 1 = (19 r2/20) A1  >= (19*10/20)*20^M1*(200/19)^L1.
# For each c, the valuations FIX (L1,M1); enumerate c until size makes it impossible.
sols = []
cmax_checked = 0
for c in range(1, 200):
    lhs2 = 2 + (1 if c % 2 else 2 + vp(c, 2))
    lhs5 = 1 + (0 if c % 2 else 1 + vp(c, 5))
    for r2 in (19, 10):
        # 3L1+2M1 = lhs2 - v2(r2), 2L1+M1 = lhs5 - v5(r2)
        A, B = lhs2 - vp(r2, 2), lhs5 - (1 if r2 == 10 else 0)
        L1 = 2 * B - A; M1 = B - 2 * L1
        # solve: 3L+2M=A, 2L+M=B -> L = 2B - A, M = B - 2L
        if L1 < 0 or M1 < 0: continue
        if r2 == 10 and M1 > 0: continue
        if 20 * (19 ** c - 1) * 19 ** L1 == 19 * r2 * P(L1, M1):
            sols.append((c, r2, L1, M1))
    cmax_checked = c
print('Case I: valuation-forced (L1,M1) checked for c<200: solutions', sols)
# size argument showing c<200 suffices: valuations give 2L1+M1 <= 2+log5 c, so
#   19^c <= 1+(361/20)*20^(2L1+M1)  ; for c>=20 this fails
for c in range(20, 2000):
    Bmax = 2 + math.log(c, 5)
    assert c * LOG19 > math.log(1 + 361 / 20 * 20 ** Bmax), c
print('Case I: for c>=20, 19^c exceeds the largest RHS allowed by the 5-adic bound -> no solutions (elementary).')

# ---------------- Case II: psi = 19^E ----------------
# derived structure (see report): E even, r3=19, L2 = 2v5(E)-v2(E) in {0,1},
# M2 = 1+2v2(E)-3v5(E), L1+L2 = 1+[x2=N,L1>=1].
# Lambda = 19 r2 A2 A1 / 19^E - 1 = -(1+19 A2)/19^E,  |Lambda| <= (1+380E^1.8614)/19^E
# exponents: b2 = v2(r2)+3s+2mu, b5 = v5(r2)+2s+mu, b19 = 1+v19(r2)-s-E, s=L1+L2<=2, mu=M1+M2<E
#            B <= 2E+7
C3 = matveev_C(3, [LOG2, LOG5, LOG19])
print('Matveev n=3 constant C3 = %.6e' % C3)
ex = math.log(20) / math.log(5)

def f(E):
    return E * LOG19 - math.log(1 + 380 * E ** ex) - C3 * (1 + math.log(2 * E + 7))

lo, hi = 1.0, 1e20
assert f(hi) > 0
for _ in range(300):
    mid = math.sqrt(lo * hi)
    if f(mid) > 0: hi = mid
    else: lo = mid
E0 = int(math.floor(hi))
print('Case II Matveev a priori bound: E <= %d  (%.4e)' % (E0, E0))

alpha_fn = lambda: [mp.log(2), mp.log(5), mp.log(19)]
Eb = E0
rounds = []
while True:
    X0 = 2 * Eb + 7
    logk = math.log(1.12 * (1 + 380 * Eb ** ex))
    newE, info = reduce_bound(alpha_fn, 3, X0, logk, LOG19)
    rounds.append((Eb, X0, info['Cexp'], info['log10c'], newE))
    print('  reduction: E<=%d  X0=%d  C=10^%d  log10(c)=%.2f  ->  E<=%d' % (Eb, X0, info['Cexp'], info['log10c'], newE))
    if newE >= Eb:
        break
    Eb = newE
Ef = Eb
print('Case II final bound: E <= %d' % Ef)

# structural exact check for E<=Ef
found = []
for E in range(2, Ef + 1, 2):
    a, b = vp(E, 2), vp(E, 5) if E % 5 == 0 else 0
    L2 = 2 * b - a; M2 = 1 + 2 * a - 3 * b
    for L1 in range(0, 3):
        for x2 in 'NK':
            for L2t in range(0, 3):
                for M2t in range(0, M2 + 3 if M2 >= 0 else 3):
                    # brute over small (L2,M2) around the forced values, any M1 by size
                    r2 = rval(x2)
                    A2 = F(200, 19) ** L2t * 20 ** M2t
                    rhs = (F(19 ** E - 1) / (19 * A2) - 1)
                    if rhs <= 0: continue
                    q = rhs / (r2 * F(200, 19) ** L1)
                    if q.denominator == 1 and q.numerator >= 1:
                        n_ = q.numerator; M1 = 0
                        while n_ % 20 == 0: n_ //= 20; M1 += 1
                        if n_ == 1 and not (x2 == 'K' and M1 > 0) and not (x2 == 'N' and L1 + M1 == 0):
                            found.append((E, L1, M1, x2, L2t, M2t))
print('Case II structural check (E<=%d, L1<=2, small u2, x3=N): solutions' % Ef, found)
# brute-force check from the closed form: every word with psi < 19^(Ef+1)
Hsum = int(math.floor((Ef + 1) * math.log10(19) - 2)) + 1
cnt = 0; bad = []
for w in three_letter_words(Hsum, Hsum):
    if w[0] + w[1] + w[3] + w[4] > Hsum: continue
    cnt += 1
    if strip(T3(*w)) == 1: bad.append(w)
print('brute: all %d words with H1+H2<=%d (covers psi<19^%d): fingerprint-1 words: %s' % (cnt, Hsum, Ef + 1, bad))
print('Theorem 16.1 holds by this computation:', not bad and not found and not sols)
