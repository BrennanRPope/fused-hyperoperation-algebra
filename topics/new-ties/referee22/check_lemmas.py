"""Items 3, 4, 5: exhaustive checks on canonical words, and on all words up to length 8.
Usage: python3 check_lemmas.py MAXLEN"""
import sys, time, random
from fractions import Fraction as F
from itertools import product, combinations
from collections import defaultdict
sys.path.insert(0, '.')
from core22 import *

MAXLEN = int(sys.argv[1])
t0 = time.time()
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.append(s)

D = {'l': 0, 'm': 0, 'N': 1, 'K': 1}
INV = {x: 1 / S[x] for x in S}

# ---- canonical words with psi, s, t computed incrementally (append on the right) ----
# appending y (applied first): s' = s*s_y, t' = s*t_y + t, psi' = psi/s_y + d_y
words = []   # (w, s, t)
def rec(w, s, t):
    if w[-1] in NK:
        words.append((w, s, t))
    if len(w) == MAXLEN: return
    for x in 'lmNK':
        if len(w) == 1 and x == 'N': continue
        if w[-1] == 'm' and x in 'Kl': continue
        rec(w + x, s * S[x], s * T[x] + t)
rec('N', S['N'], T['N'])
say('canonical words up to length', MAXLEN, ':', len(words), ' t=%.1fs' % (time.time() - t0))
# spot-check incremental s,t
for w, s, t in random.sample(words, min(2000, len(words))):
    assert (s, t) == st(w), w

# phi(s_w) = v2 - 2 v5 - v19 equals the number of N/K letters (hand-derived invariant)
bad = 0
for w, s, t in random.sample(words, min(20000, len(words))):
    if val(s, 2) - 2 * val(s, 5) - val(s, 19) != sum(1 for x in w if x in NK): bad += 1
say('invariant v2(s)-2v5(s)-v19(s) = #N/K on 20000 sampled canonical words: failures', bad)

# ---- Conjecture 15.3 / Theorem 21.5 consequence: distinct canonical words have distinct fingerprints ----
byfp = defaultdict(list)
for w, s, t in words:
    byfp[fingerprint(t / s)].append(w)
coll = [v for v in byfp.values() if len(v) > 1]
say('fingerprint classes:', len(byfp), '; classes with two or more canonical words:', len(coll), coll[:5])

# ---- Item 3 (Lemma 21.3) consequence on the collisions that exist: (none expected) ----
# Instead exercise the lemma's mechanics on every pair of equal-fingerprint canonical words
# found (trivially W=W') and show the canonical hypothesis is needed:
def balance_parts(Lt, Rt):
    """all minimal balancing subfamilies containing the overall largest number"""
    allv = [(0, i, v) for i, v in enumerate(Lt)] + [(1, j, v) for j, v in enumerate(Rt)]
    n = len(allv)
    big = max(range(n), key=lambda k: allv[k][2])
    bal = []
    for mask in range(1, 1 << n):
        if not mask >> big & 1: continue
        sl = sum(allv[k][2] for k in range(n) if mask >> k & 1 and allv[k][0] == 0)
        sr = sum(allv[k][2] for k in range(n) if mask >> k & 1 and allv[k][0] == 1)
        if sl == sr: bal.append(mask)
    mins = [m for m in bal if not any((o != m and o & m == o) for o in bal)]
    return [([allv[k][2] for k in range(n) if m >> k & 1 and allv[k][0] == 0],
             [allv[k][2] for k in range(n) if m >> k & 1 and allv[k][0] == 1]) for m in mins]
W, Wp = 'NN', 'Nm'      # same psi (=20), not canonical (NN)
parts = balance_parts(psi_terms(W), psi_terms(Wp))
say('non-canonical NN vs Nm: psi', psi(W), psi(Wp), '; terms', psi_terms(W), psi_terms(Wp),
    '; minimal part with largest:', parts, '-> leading terms 19 vs 20 differ: canonical hypothesis is needed')

# general fact used in 21.3: a minimal balancing subfamily containing a given number is primitive
random.seed(3)
pool = [2 ** a * 5 ** b * 19 ** c for a in range(6) for b in range(4) for c in range(3)]
tested = badprim = 0
for _ in range(4000):
    Lt = random.sample(pool, random.randint(1, 4)); Rt = random.sample(pool, random.randint(1, 4))
    d = sum(Lt) - sum(Rt)
    if d > 0: Rt.append(d)
    elif d < 0: Lt.append(-d)
    if len(Lt) + len(Rt) > 9: continue
    for pl, pr in balance_parts(Lt, Rt):
        tested += 1
        allv = [(0, v) for v in pl] + [(1, v) for v in pr]
        n = len(allv)
        for mask in range(1, (1 << n) - 1):
            sl = sum(allv[k][1] for k in range(n) if mask >> k & 1 and allv[k][0] == 0)
            sr = sum(allv[k][1] for k in range(n) if mask >> k & 1 and allv[k][0] == 1)
            if sl == sr: badprim += 1; break
say('minimal part containing the largest is primitive: tested', tested, 'parts, non-primitive', badprim)

# ---- Item 4 (Lemma 21.4): u determined by (canon(u), t_u) over canonical W = N u ----
groups = defaultdict(set)
for w, s, t in words:
    u = w[1:]
    su, tu = st(u) if len(u) < 4 else (s / S['N'], (t - T['N']) / S['N'])
    groups[(canon(u), tu)].add(u)
badu = [sorted(v) for v in groups.values() if len(v) > 1]
say('Lemma 21.4: groups (canon(u), t_u):', len(groups), '; groups with two different u:', len(badu), badu[:5])
# how much does canon(u) alone lose? (to see the lemma is not vacuous)
g2 = defaultdict(set)
for w, s, t in words: g2[canon(w[1:])].add(w[1:])
say('  for comparison: canon(u) classes', len(g2), 'holding', sum(len(v) for v in g2.values()), 'distinct u; largest class size',
    max(len(v) for v in g2.values()))
# check of (s_u, t_u) formulas
for w, s, t in random.sample(words, min(2000, len(words))):
    assert (s / S['N'], (t - T['N']) / S['N']) == st(w[1:])
# drop the canonical hypothesis on W: u ranging over ALL words up to length 8 (not just suffixes of canonical words)
g3 = defaultdict(set)
for L in range(0, 9):
    for u in product('lmNK', repeat=L):
        u = ''.join(u)
        if u and u[-1] not in NK: continue      # keep only u with no trailing l/m
        g3[(canon(u), st(u)[1])].add(u)
ex = sorted((sorted(v, key=lambda z: (len(z), z)) for v in g3.values() if len(v) > 1), key=lambda g: (len(g[0]), g))
say('  without the canonical hypothesis (all u up to length 8 ending in N/K):', len(ex), 'groups with several u; shortest:',
    [g[:3] for g in ex[:6]])

# ---- Theorem 21.5 consequence for ALL words up to length 8: same fingerprint => same canonical form ----
g4 = defaultdict(set)
for L in range(0, 9):
    for w in product('lmNK', repeat=L):
        w = ''.join(w); g4[fingerprint(psi(w))].add(canon(w))
b4 = [v for v in g4.values() if len(v) > 1]
say('15.3 on all 87381 words up to length 8: fingerprint classes', len(g4), '; classes with two canonical forms', len(b4))

# ---- The induction step of 21.5 run mechanically on canonical pairs with equal fingerprint
#      (only trivial pairs exist), and the reduction W=Nu -> u preserving equal fingerprints ----
# check: for canonical W, tau_0 = 1/s_u and psi(W) = tau_0 + psi(u)
bad = 0
for w, s, t in random.sample(words, min(5000, len(words))):
    u = w[1:]; su, tu = st(u)
    if psi_terms(w)[0] != 1 / su or psi(w) != 1 / su + (tu / su if True else 0): bad += 1
say('tau_0 = 1/s_u and psi(Nu) = tau_0 + psi(u) on 5000 canonical words: failures', bad)
say('done t=%.1fs' % (time.time() - t0))
open('check_lemmas_log_%d.txt' % MAXLEN, 'w').write('\n'.join(log) + '\n')
