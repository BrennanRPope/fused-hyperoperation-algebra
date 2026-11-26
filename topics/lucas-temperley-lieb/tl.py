"""Temperley-Lieb diagram algebra TL_n(delta), built from scratch (planar matchings, stacking, loop counting).
Used to test the Lucas weights U_m(k, d) against Temperley-Lieb objects without assuming Chebyshev facts."""
from functools import lru_cache
import itertools

def matchings(n):
    """All non-crossing perfect matchings of 2n boundary points 0..2n-1 placed on a circle
    (0..n-1 along the top left to right, then n..2n-1 along the bottom right to left)."""
    @lru_cache(None)
    def nc(pts):
        if not pts:
            return [()]
        out = []
        a = pts[0]
        for i in range(1, len(pts), 2):
            b = pts[i]
            for left in nc(pts[1:i]):
                for right in nc(pts[i + 1:]):
                    out.append(((a, b),) + left + right)
        return out
    res = []
    for m in nc(tuple(range(2 * n))):
        d = {}
        for a, b in m:
            d[a] = b; d[b] = a
        res.append(tuple(d[i] for i in range(2 * n)))
    return res

def top(n, i): return i            # top point i (left to right)
def bot(n, i): return 2 * n - 1 - i  # bottom point i (left to right)

def compose(n, D1, D2):
    """D1 on top of D2. Returns (matching, loops)."""
    # nodes: ('a', p) for D1 points, ('b', p) for D2 points; D1 bottom i glued to D2 top i
    def partner(side, p):
        D = D1 if side == 'a' else D2
        return D[p]
    def glue(side, p):
        # if p is an internal point, return the glued point on the other diagram
        if side == 'a' and p >= n:          # D1 bottom
            i = 2 * n - 1 - p
            return ('b', top(n, i))
        if side == 'b' and p < n:           # D2 top
            return ('a', bot(n, p))
        return None
    result = [None] * (2 * n)
    seen = set()
    def outer_label(side, p):  # map to the composite's boundary numbering
        return p if side == 'a' else p   # D1 top stays top; D2 bottom stays bottom (same numbering)
    for side0, rng in (('a', range(n)), ('b', range(n, 2 * n))):
        for p0 in rng:
            if (side0, p0) in seen:
                continue
            side, p = side0, p0
            seen.add((side, p))
            while True:
                q = partner(side, p); seen.add((side, q))
                g = glue(side, q)
                if g is None:
                    result[p0] = q; result[q] = p0
                    break
                side, p = g; seen.add((side, p))
    loops = 0
    for i in range(n):  # remaining internal points form loops
        node = ('b', top(n, i))
        if node in seen:
            continue
        side, p = node
        while (side, p) not in seen:
            seen.add((side, p))
            q = partner(side, p); seen.add((side, q))
            side, p = glue(side, q)
        loops += 1
    return tuple(result), loops

class TL:
    def __init__(self, n, delta):
        self.n, self.delta = n, delta
        self.basis = matchings(n)
        self.index = {b: i for i, b in enumerate(self.basis)}
        self.one = tuple(bot(n, i) if j < n else None for j, i in []) if False else self._identity()
    def _identity(self):
        n = self.n; m = [None] * (2 * n)
        for i in range(n):
            m[top(n, i)] = bot(n, i); m[bot(n, i)] = top(n, i)
        return tuple(m)
    def e(self, i):  # cup-cap at strands i, i+1 (0-based)
        n = self.n; m = list(self._identity())
        m[top(n, i)] = top(n, i + 1); m[top(n, i + 1)] = top(n, i)
        m[bot(n, i)] = bot(n, i + 1); m[bot(n, i + 1)] = bot(n, i)
        return tuple(m)
    def elem(self, terms):  # dict matching -> coeff
        return {k: v for k, v in terms.items() if v != 0}
    def mul(self, X, Y):
        out = {}
        for a, ca in X.items():
            for b, cb in Y.items():
                m, l = compose(self.n, a, b)
                out[m] = out.get(m, 0) + ca * cb * self.delta ** l
        return {k: v for k, v in out.items() if v != 0}
    def add(self, X, Y, s=1):
        out = dict(X)
        for k, v in Y.items():
            out[k] = out.get(k, 0) + s * v
        return {k: v for k, v in out.items() if v != 0}
    def scal(self, c, X): return {k: c * v for k, v in X.items() if c * v != 0}
    def closure_loops(self, D):
        """Markov closure: join top i to bottom i around the right side; count loops."""
        n = self.n
        seen = set(); loops = 0
        for start in range(2 * n):
            if start in seen:
                continue
            p = start
            while p not in seen:
                seen.add(p); q = D[p]; seen.add(q)
                p = (2 * n - 1 - q) if q < n else (2 * n - 1 - q)  # top i <-> bottom i
            loops += 1
        return loops
    def trace(self, X):
        return sum(c * self.delta ** self.closure_loops(D) for D, c in X.items())
    def gram(self):
        """Meander (Gram) matrix <a, b> = delta^(loops of a* stacked on b, then closed)."""
        n = self.n
        def star(D):  # flip top and bottom
            m = [None] * (2 * n)
            for p in range(2 * n):
                q = D[p]
                fp = (2 * n - 1 - p); fq = (2 * n - 1 - q)
                m[fp] = fq
            return tuple(m)
        B = self.basis
        G = [[None] * len(B) for _ in B]
        for i, a in enumerate(B):
            for j, b in enumerate(B):
                m, l = compose(n, star(a), b)
                G[i][j] = self.delta ** (l + self.closure_loops(m))
        return G
