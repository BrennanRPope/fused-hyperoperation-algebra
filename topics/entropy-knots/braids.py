"""Braids of periodic orbits of the thickened logistic map (the horseshoe), and braid entropy
from the growth of free-group words under the braid's Artin action."""
import numpy as np, math, json
from scipy.optimize import brentq

def orbit_braid(xs, c=0.5):
    """xs: the orbit points (any order). Phase 1: every point right of c swings over to the left
    of c along its own semicircle in the upper half plane, with slightly different speeds so no
    two crossings happen at once. Phase 2: an order-preserving slide to the targets (no
    crossings). Returns the generator word as a list of signed indices (1-based)."""
    xs = sorted(xs); n = len(xs)
    rho = [x - c for x in xs]
    a = [1 + 0.05 * k / n for k in range(n)]
    def pos(j, s):
        if rho[j] <= 0: return complex(xs[j], 0.0)
        th = math.pi * s ** a[j]
        return complex(c + rho[j] * math.cos(th), rho[j] * math.sin(th))
    events = []
    for j in range(n):
        for k in range(j + 1, n):
            f = lambda s: pos(j, s).real - pos(k, s).real
            S = np.linspace(0, 1, 4001)
            v = [f(s) for s in S]
            for i in range(len(S) - 1):
                if v[i] == 0 or v[i] * v[i + 1] < 0:
                    t = brentq(f, S[i], S[i + 1]) if v[i] != 0 else S[i]
                    if 0 < t < 1: events.append((t, j, k))
    events.sort()
    order = list(range(n))            # order[p] = point at position p (left to right)
    word = []
    for t, j, k in events:
        pj, pk = order.index(j), order.index(k)
        assert abs(pj - pk) == 1, 'non-adjacent crossing'
        p = min(pj, pk); left = order[p]; right = order[p + 1]
        over_left = pos(left, t).imag > pos(right, t).imag
        word.append((p + 1) if over_left else -(p + 1))
        order[p], order[p + 1] = order[p + 1], order[p]
    return word, n

# ---- free group words: generators 1..n, inverses negative ----
def reduce(w):
    out = []
    for x in w:
        if out and out[-1] == -x: out.pop()
        else: out.append(x)
    return out

def artin(sig, n):
    """Automorphism of F_n for sigma_i^{+-1}: returns dict gen -> word."""
    i = abs(sig); m = {g: [g] for g in range(1, n + 1)}
    if sig > 0:
        m[i] = [i, i + 1, -i]; m[i + 1] = [i]
    else:
        m[i] = [i + 1]; m[i + 1] = [-(i + 1), i, i + 1]
    return m

def apply(m, w):
    out = []
    for x in w:
        out.extend(m[x] if x > 0 else [-y for y in reversed(m[-x])])
    return reduce(out)

def compose_word(word, n):
    """The automorphism of the whole braid word, as images of generators."""
    imgs = {g: [g] for g in range(1, n + 1)}
    for sig in word:            # apply sigma's in order: image under phi_word = phi_1(phi_2(...))
        m = artin(sig, n)
        imgs = {g: apply(m, imgs[g]) for g in imgs}
    return imgs

def cyc_reduce(w):
    w = reduce(w)
    while len(w) > 1 and w[0] == -w[-1]: w = w[1:-1]
    return w

def braid_entropy(word, n, iters=400, cap=2000000):
    if not word: return 0.0
    phi = compose_word(word, n)
    w = list(range(1, n + 1))            # the word x1 x2 ... xn (boundary) is fixed; use a test loop
    w = [1, -2] if n >= 2 else [1]
    w = [1, 2, -1, 3, -2] if n >= 3 else w
    lens = []
    for _ in range(iters):
        w = cyc_reduce(apply(phi, w))
        lens.append(len(w))
        if len(w) > cap: break
    if len(lens) < 6: return float('nan')
    L = np.log(np.array(lens, float))
    k = len(L); a = k // 2
    x = np.arange(a, k)
    return float(max(0.0, np.polyfit(x, L[a:], 1)[0]))

if __name__ == '__main__':
    # tests: sigma1 sigma2^-1 on 3 strands has dilatation phi^2
    print('s1 s2^-1:', braid_entropy([1, -2], 3), 2 * math.log((1 + 5 ** .5) / 2))
    print('s1 s2 (periodic):', braid_entropy([1, 2], 3))
    print('s1^2 s2^-1:', braid_entropy([1, 1, -2], 3), math.log((3 + 5 ** .5) / 2 + 0) )
