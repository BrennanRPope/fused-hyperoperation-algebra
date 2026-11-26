"""Helpers for checks_trotter.py: word polynomials, descent weights, piecewise Magnus terms."""
import itertools, math
from fractions import Fraction as F

# ---------- noncommutative polynomials: dict word(tuple) -> Fraction
def padd(*ps):
    out = {}
    for p in ps:
        for w, c in p.items():
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c != 0}
def pscale(p, s): return {w: c*s for w, c in p.items() if c*s != 0}
def pmul(p, q, maxdeg=None):
    out = {}
    for w1, c1 in p.items():
        for w2, c2 in q.items():
            if maxdeg is not None and len(w1) + len(w2) > maxdeg: continue
            w = w1 + w2; out[w] = out.get(w, 0) + c1*c2
    return {w: c for w, c in out.items() if c != 0}
def comm(p, q): return padd(pmul(p, q), pscale(pmul(q, p), -1))
def L(*letters): return {tuple(letters): F(1)}

def descents(s): return sum(1 for i in range(len(s)-1) if s[i] > s[i+1])
def weight(n, d): return F((-1)**d, n*math.comb(n-1, d))

def magnus_piecewise(blocks, n):
    """degree-n term of log(e^{b1} e^{b2} ...), blocks = list of (letter, length) in time order, piecewise constant.
    Omega_n = sum over orderings sigma of weight * integral over the time simplex of H(t_sigma)."""
    T = sum(l for _, l in blocks)
    out = {}
    # an increasing time tuple t_1 < ... < t_n lands j_1 points in block 1, j_2 in block 2 ...: volume prod l^j / j!
    for js in itertools.product(range(n+1), repeat=len(blocks)):
        if sum(js) != n: continue
        vol = F(1)
        letters = []
        for (let, l), j in zip(blocks, js):
            vol *= F(l)**j / math.factorial(j); letters += [let]*j
        for sgm in itertools.permutations(range(n)):
            w = tuple(letters[sgm[i]] for i in range(n))   # factor i of the word is H at time t_{sigma(i)}
            # convention: the product is e^{b1} e^{b2} ..., earlier block on the left (checked numerically in T4, T5)
            c = weight(n, descents(sgm)) * vol
            out[w] = out.get(w, 0) + c
    return {w: c for w, c in out.items() if c != 0}
