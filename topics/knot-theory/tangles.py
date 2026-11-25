"""Tangle diagrams and the Kauffman bracket, two independent ways.

1. diagrams: a tangle is a set of crossing tiles plus joins; the bracket of a closed diagram is
   the state sum over the 2^n corners of the cube (loops counted by union-find).
2. algebra: a tangle's bracket is a pair (f, g) with <T> = f<[0]> + g<[infinity]>, and
   sums are tsum (side by side) and vsum (stacked).
Positions: NW=(-1,1), NE=(1,1), SW=(-1,-1), SE=(1,-1). Tile type +1 has its over strand NW-SE,
type -1 has its over strand SW-NE. Strands go straight through: NW<->SE, SW<->NE.
"""
import sympy as sp, itertools
A = sp.symbols('A')
d = -A**2 - A**-2
POS = {'NW': (-1, 1), 'NE': (1, 1), 'SW': (-1, -1), 'SE': (1, -1)}
THRU = {'NW': 'SE', 'SE': 'NW', 'SW': 'NE', 'NE': 'SW'}

class Tangle:
    def __init__(self):
        self.tiles = []      # list of (type, {pos: label})
        self.joins = []      # pairs of labels joined by an arc
        self.ends = {}       # pos -> label
        self.n = 0
    def fresh(self):
        self.n += 1; return self.n

def zero():
    t = Tangle(); l = {p: t.fresh() for p in POS}; t.ends = l
    t.joins += [(l['NW'], l['NE']), (l['SW'], l['SE'])]; return t
def inf():
    t = Tangle(); l = {p: t.fresh() for p in POS}; t.ends = l
    t.joins += [(l['NW'], l['SW']), (l['NE'], l['SE'])]; return t
def cross(s):
    t = Tangle(); l = {p: t.fresh() for p in POS}; t.ends = dict(l); t.tiles.append((s, l)); return t

def _merge(S, T):
    """Copy T's labels after S's."""
    R = Tangle(); R.n = S.n + T.n; off = S.n
    R.tiles = list(S.tiles) + [(s, {p: v + off for p, v in l.items()}) for s, l in T.tiles]
    R.joins = list(S.joins) + [(a + off, b + off) for a, b in T.joins]
    return R, {p: v + off for p, v in T.ends.items()}

def hsum(S, T):
    R, Te = _merge(S, T)
    R.joins += [(S.ends['NE'], Te['NW']), (S.ends['SE'], Te['SW'])]
    R.ends = {'NW': S.ends['NW'], 'SW': S.ends['SW'], 'NE': Te['NE'], 'SE': Te['SE']}; return R
def vsum(S, T):
    R, Te = _merge(S, T)
    R.joins += [(S.ends['SW'], Te['NW']), (S.ends['SE'], Te['NE'])]
    R.ends = {'NW': S.ends['NW'], 'NE': S.ends['NE'], 'SW': Te['SW'], 'SE': Te['SE']}; return R

def closeN(T):
    R, _ = _merge(T, Tangle()); R.ends = {}
    R.joins += [(T.ends['NW'], T.ends['NE']), (T.ends['SW'], T.ends['SE'])]; return R
def closeD(T):
    R, _ = _merge(T, Tangle()); R.ends = {}
    R.joins += [(T.ends['NW'], T.ends['SW']), (T.ends['NE'], T.ends['SE'])]; return R

def _components(n, pairs):
    par = list(range(n + 1))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in pairs:
        ra, rb = f(a), f(b)
        if ra != rb: par[ra] = rb
    return len({f(x) for x in range(1, n + 1)})

def bracket(K):
    """State sum over the cube of smoothings: sum A^(#A-#B) d^(loops-1)."""
    tot = 0
    for st in itertools.product((0, 1), repeat=len(K.tiles)):
        pairs = list(K.joins); na = 0
        for (s, l), b in zip(K.tiles, st):
            horiz = [(l['NW'], l['NE']), (l['SW'], l['SE'])]
            vert = [(l['NW'], l['SW']), (l['NE'], l['SE'])]
            Asm = horiz if s == 1 else vert          # A-smoothing of type +1 gives [0]
            if b == 0: pairs += Asm; na += 1
            else: pairs += (vert if Asm is horiz else horiz)
        L = _components(K.n, pairs)
        tot += A**(na - (len(K.tiles) - na)) * d**(L - 1)
    return sp.expand(tot)

def writhe(K):
    """Orient every component by walking it; sign = sign of cross(over, under)."""
    inc = {}
    for i, (s, l) in enumerate(K.tiles):
        for p, v in l.items(): inc.setdefault(v, []).append((i, p))
    # arcs: labels joined; build adjacency label->label
    adj = {}
    for a, b in K.joins: adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    dirs = {}  # (tile, strand) -> (entry, exit)
    seen = set()
    def other_end(label, came_from_tile_end):
        # walk along joined labels to the next tile end
        prev, cur = None, label
        while True:
            ends = [e for e in inc.get(cur, []) if e != came_from_tile_end]
            if ends: return ends[0]
            nxt = [x for x in adj.get(cur, []) if x != prev]
            prev, cur = cur, nxt[0]
    for i, (s, l) in enumerate(K.tiles):
        for start in ('NW', 'SW'):
            if (i, start) in seen: continue
            ti, p = i, start
            while (ti, p) not in seen:
                q = THRU[p]; seen.add((ti, p)); seen.add((ti, q))
                dirs[(ti, frozenset((p, q)))] = (p, q)
                lab = K.tiles[ti][1][q]
                ti, p = other_end(lab, (ti, q))
    w = 0
    for i, (s, l) in enumerate(K.tiles):
        over = frozenset(('NW', 'SE')) if s == 1 else frozenset(('SW', 'NE'))
        under = frozenset(('SW', 'NE')) if s == 1 else frozenset(('NW', 'SE'))
        def vec(pr):
            a, b = dirs[(i, pr)]; return (POS[b][0] - POS[a][0], POS[b][1] - POS[a][1])
        o, u = vec(over), vec(under)
        w += 1 if o[0] * u[1] - o[1] * u[0] > 0 else -1
    return w

t = sp.symbols('t')
def jones(K):
    """V(t) = (-A^3)^(-w) <K>, with A = t^(-1/4): returned as a Laurent poly in s = t^(1/4)... we
    return the polynomial in A and let callers substitute."""
    return sp.expand((-A**3)**(-writhe(K)) * bracket(K))

def A_to_t(fA):
    """Substitute A = t^(-1/4) for knots (powers of A multiples of 4)."""
    fA = sp.expand(fA); out = 0
    for (e,), c in sp.Poly(sp.expand(fA * A**200), A).terms():
        e -= 200; assert e % 4 == 0, e
        out += c * t**(-e // 4)
    return sp.expand(out)

# ---------------- algebra -----------------
Z0 = (sp.Integer(1), sp.Integer(0)); ZI = (sp.Integer(0), sp.Integer(1))
XP = (A, A**-1); XM = (A**-1, A)
def ts(S, T, dd=d):
    return (sp.expand(S[0]*T[0]), sp.expand(S[0]*T[1] + S[1]*T[0] + dd*S[1]*T[1]))
def vs(S, T, dd=d):
    return (sp.expand(S[0]*T[1] + S[1]*T[0] + dd*S[0]*T[0]), sp.expand(S[1]*T[1]))
def cN(T, dd=d): return sp.expand(dd*T[0] + T[1])
def cD(T, dd=d): return sp.expand(T[0] + dd*T[1])

def build(word):
    """word: list of 'h+','h-','v+','v-' applied to [0]; returns (diagram, algebra pair)."""
    D, P = zero(), Z0
    for w in word:
        s = 1 if w[1] == '+' else -1
        X = XP if s == 1 else XM
        if w[0] == 'h': D, P = hsum(D, cross(s)), ts(P, X)
        else: D, P = vsum(D, cross(s)), vs(P, X)
    return D, P

def cf_word(cf):
    """[a1, a2, ..., an] -> h^a1 v^a2 h^a3 ... from [0] (signs from the a's)."""
    w = []
    for i, a in enumerate(cf):
        k = 'h' if i % 2 == 0 else 'v'
        w += [k + ('+' if a > 0 else '-')] * abs(a)
    return w

def cf_value(cf):
    """Fraction of h^a1 v^a2 ...: h adds to F, v adds to 1/F (start F = 0)."""
    F = sp.Integer(0)
    for i, a in enumerate(cf):
        if i % 2 == 0: F = F + a
        else: F = 1 / (1 / F + a) if F != 0 else sp.Rational(1, a)
    return F
