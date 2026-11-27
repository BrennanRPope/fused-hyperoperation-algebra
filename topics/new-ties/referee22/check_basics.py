"""Spot checks of the 'taken as proved' facts and of canonical-form properties,
plus Item 1 (2-adic sizes)."""
import random, sys
from fractions import Fraction as F
from itertools import product
sys.path.insert(0, '.')
from core22 import *

random.seed(22)
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); log.append(s)

# Theorem 14.1 and the recursion, on all words up to length 8
bad = 0; n = 0
for L in range(0, 9):
    for w in product('lmNK', repeat=L):
        w = ''.join(w); n += 1
        p = psi(w)
        if p != sum(psi_terms(w), F(0)): bad += 1
        if w:
            x, rest = w[0], w[1:]
            s_rest = st(rest)[0]
            d = 1 if x in NK else 0
            if p != psi(rest) + F(d) / s_rest: bad += 1
say('Thm 14.1 + recursion on all', n, 'words up to length 8: failures', bad)

# Prop 15.2 identities on random tails
bad = 0
for _ in range(3000):
    w = ''.join(random.choice('lmNK') for _ in range(random.randint(0, 9)))
    if psi('K' + w) != psi('N' + w): bad += 1
    a = psi('NN' + w)
    if not (a == psi('KN' + w) == psi('Nm' + w) == psi('Km' + w)): bad += 1
    if st('mK' + w) != st('lN' + w) or st('ml' + w) != st('lm' + w): bad += 1
    if psi('l' + w) != psi(w) or psi('m' + w) != psi(w): bad += 1
    if psi(w + 'l') != psi(w) * F(200, 19) or psi(w + 'm') != psi(w) * 20: bad += 1
say('Prop 15.2 identities on 3000 random tails: failures', bad)

# canonical form: keeps fingerprint, idempotent, local description agrees
bad_fp = bad_id = 0; n = 0
for L in range(0, 10):
    for w in product('lmNK', repeat=L):
        w = ''.join(w); n += 1
        c = canon(w)
        if fingerprint(psi(c)) != fingerprint(psi(w)): bad_fp += 1
        if canon(c) != c: bad_id += 1
say('canon on all', n, 'words up to length 9: fingerprint changed', bad_fp, '; not idempotent', bad_id)
for L in range(1, 11):
    a = set(w for w in gen_canon(L) if len(w) == L)
    b = set(w for w in gen_canon_local(L) if len(w) == L)
    assert a == b, L
say('local description of canonical words = canon(w)==w test, lengths 1..10: agree')
say('canonical words by length:', [sum(1 for w in gen_canon_local(L) if len(w) == L) for L in range(1, 13)])
# confluence of rewriting: random rewrite order gives same result
def rewrite_random(w):
    while True:
        spots = [i for i in range(len(w) - 1) if w[i:i + 2] in ('mK', 'ml')]
        if not spots: return w
        i = random.choice(spots)
        w = w[:i] + ('lN' if w[i:i + 2] == 'mK' else 'lm') + w[i + 2:]
def rewrite_det(w):
    while True:
        v = w.replace('mK', 'lN').replace('ml', 'lm')
        if v == w: return w
        w = v
bad = 0
for _ in range(20000):
    w = ''.join(random.choice('lmNK') for _ in range(random.randint(0, 14)))
    if rewrite_random(w) != rewrite_det(w): bad += 1
say('rewriting mK->lN, ml->lm: random order vs fixed order on 20000 words: differ', bad)
# no canonical word has two consecutive leading N/K terms with ratio 1/19
bad = 0
for w in gen_canon_local(13):
    tt = psi_terms(w)
    if len(tt) >= 2 and tt[1] / tt[0] == F(1, 19): bad += 1
    if any(tt[i + 1] / tt[i] > F(1, 10) for i in range(len(tt) - 1)): bad += 1
say('canonical words up to length 13: first two terms in ratio 1/19, or a ratio > 1/10:', bad)

# ---- Item 1 ----
def v2(q): return val(q, 2)
bad = 0; cnt = 0; weak = 0
for _ in range(200000):
    a = random.randint(-4, 8); b = random.randint(-6, 8)
    # kappa = 2^-a * odd/odd, mu = 2^-b * odd/odd
    odd = lambda: random.choice([1, 3, 5, 7, 9, 11, 13, 15, 19, 21, 25, 95, 125])
    kappa = F(odd(), odd()) / F(2) ** a
    mu = F(odd(), odd()) / F(2) ** b
    if kappa == 2: continue
    nu = kappa / (2 - kappa); lam = mu * (2 - kappa)
    got = (v2(1 / nu), v2(1 / kappa), v2(1 / mu), v2(1 / lam))
    if a >= 0:
        cnt += 1
        if got != (0, a, b, a + b): bad += 1
    else:
        if got != (0, a, b, a + b): weak += 1
say('Item 1: random rational kappa, mu with a>=0 (', cnt, 'cases): failures', bad)
say('Item 1: with a<0 the pattern (0,a,b,a+b) fails in', weak, 'random cases (so a>=0 is needed; a>=1 is more than needed)')
kappa, mu = F(1, 10), F(1, 20)
nu = kappa / (2 - kappa); lam = mu * (2 - kappa)
say('Item 1 example: nu, kappa, mu, lambda =', nu, kappa, mu, lam, '; 2^v2(1/.) =',
    [2 ** v2(1 / x) for x in (nu, kappa, mu, lam)], '; equal to s_N, s_K, s_m, s_l:',
    (nu, kappa, mu, lam) == (S['N'], S['K'], S['m'], S['l']))
open('check_basics_log.txt', 'w').write('\n'.join(log) + '\n')
