"""Round 14 (DRAFT): every psi-collision between reduced words of up to n letters, found exhaustively.

psi(w) = t_w / s_w for the word's affine map w(y) = s_w y + t_w (letters as in rdist13: l=(19/200,0), m=(1/20,0),
N=(1/19,1/19), K=(1/10,1/10)); equivalently psi(w) = -w^{-1}(0). Theorem 14.1: a second collision for the bad pair (a, b)
is the same as two words h (first letter N) and g with psi(h) = psi(g) and s_h / s_g = lambda^a / mu^b.

All reduced words (after m only m or N) of 1..n letters whose first letter is N or K are generated level by level with
psi kept modulo two 31-bit primes; words whose two residues agree are grouped, rebuilt from parent pointers and compared in
exact rational arithmetic. A reported collision is therefore exact, and none is missed, because an exact collision agrees
modulo both primes (neither prime divides a denominator). Pairs are kept when h starts with N and the
last letters differ (a common last letter cancels: psi(wx) = r_x psi(w) + d_x).
Usage: python3 psi14.py [n=14]   (n = 14: 29 million words, about 15 s, a few GB of memory). Writes psi14.json.
"""
import sys, json, numpy as np
from fractions import Fraction as F
from collections import Counter, defaultdict
n = int(sys.argv[1]) if len(sys.argv) > 1 else 14
P = [2147483647, 2147483629]
LET = "lmNK"
def inv(a, p): return pow(a, p - 2, p)
r = [[200 * inv(19, p) % p, 20, 19, 10] for p in P]     # psi(wx) = r_x psi(w) + d_x, r_x = 1/s_x
d = [0, 0, 1, 1]
allowed = {0: [0, 1, 2, 3], 1: [1, 2], 2: [0, 1, 2, 3], 3: [0, 1, 2, 3]}
k1 = np.array([1, 1], dtype=np.int64); k2 = np.array([1, 1], dtype=np.int64)
last = np.array([2, 3], dtype=np.int8)
par = [np.array([-1, -1])]; lets = [last.copy()]
store = [(k1, k2, last)]; level_sizes = [2]
for lev in range(1, n):
    nk1, nk2, nl, npar = [], [], [], []
    for x in range(4):
        idx = np.nonzero(np.isin(last, [y for y in range(4) if x in allowed[y]]))[0]
        nk1.append((k1[idx] * r[0][x] + d[x]) % P[0]); nk2.append((k2[idx] * r[1][x] + d[x]) % P[1])
        nl.append(np.full(len(idx), x, dtype=np.int8)); npar.append(idx)
    k1 = np.concatenate(nk1); k2 = np.concatenate(nk2); last = np.concatenate(nl)
    par.append(np.concatenate(npar)); lets.append(last.copy()); store.append((k1, k2, last)); level_sizes.append(len(k1))
lvl = np.concatenate([np.full(len(s[0]), i, dtype=np.int8) for i, s in enumerate(store)])
pos = np.concatenate([np.arange(len(s[0])) for s in store])
K = np.concatenate([s[0] for s in store]) * P[1] + np.concatenate([s[1] for s in store])
LA = np.concatenate([s[2] for s in store])
o = np.argsort(K, kind="stable"); Ks = K[o]; LAs = LA[o]
br = np.nonzero(np.diff(Ks))[0] + 1; starts = np.r_[0, br]; ends = np.r_[br, len(Ks)]
def word(L, p):
    w = []
    while L >= 0: w.append(LET[lets[L][p]]); p = par[L][p]; L -= 1
    return "".join(reversed(w))
RV = {"l": F(200, 19), "m": F(20), "N": F(19), "K": F(10)}
def exact(w):
    ps = F(1); Rw = RV[w[0]]
    for x in w[1:]: ps = RV[x] * ps + (1 if x in "NK" else 0); Rw *= RV[x]
    return ps, Rw            # psi(w), R_w = 1/s_w
rat = Counter(); by_len = defaultdict(Counter); ex = {}; groups = 0; false = 0
for s, e in zip(starts, ends):
    if e - s < 2 or len(set(LAs[s:e].tolist())) < 2: continue
    groups += 1
    mem = [word(int(lvl[o[i]]), int(pos[o[i]])) for i in range(s, e)]
    vals = [exact(w) for w in mem]
    for i, h in enumerate(mem):
        if h[0] != "N": continue
        for j, g in enumerate(mem):
            if g[-1] == h[-1]: continue
            if vals[i][0] != vals[j][0]: false += 1; continue
            q = vals[j][1] / vals[i][1]                  # s_h / s_g
            rat[q] += 1; by_len[str(q)][max(len(h), len(g))] += 1; ex.setdefault(q, (h, g))
def v(q, p):
    k = 0; a, b = q.numerator, q.denominator
    while a % p == 0: a //= p; k += 1
    while b % p == 0: b //= p; k -= 1
    return k
out = dict(n=n, words=int(len(K)), level_sizes=level_sizes, groups=groups, hash_false_pairs=false,
           ratios={str(q): dict(count=c, example=ex[q], v19=v(q, 19), v2=v(q, 2), v5=v(q, 5),
                                by_max_length=dict(sorted(by_len[str(q)].items()))) for q, c in rat.items()})
json.dump(out, open("psi14.json", "w"), indent=1)
print("n", n, "words", len(K), "groups", groups, "false hash pairs", false)
for q, c in sorted(rat.items(), key=lambda t: -t[1]): print(" ratio s_h/s_g =", q, " pairs", c, " example", ex[q], " v19", v(q, 19))
