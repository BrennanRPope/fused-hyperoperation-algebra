"""Round 14 (DRAFT): neighbour maps of two reduced words, read from the left, in exact rational arithmetic.

For prefixes P of h and Q of g the neighbour map is T = Q^{-1} o P, T(y) = a y + b. Peel the next letter of h when a >= 1,
else of g; then a stays in [1/20, 20). Every letter maps J = [0, 1/9] into J, so T(J) must meet J (real cut). T = identity
again means the prefixes are equal as maps. States are deduplicated with the last letters (for the grammar), so the
search is complete for all pairs whose combined length is at most the depth.
Facts this records:
 - linear states (b = 0) are the stretches lambda^i mu^-j; a nonlinear state can only be born from a linear one whose a lies
   in a window: [1, 19/9] (h takes N), [1, 10/9] (h takes K), [9/19, 1) (g takes N), [9/10, 1) (g takes K);
 - every return to the identity is trivial (both words take the same letter);
 - from each linear state found, every trip off the line, followed to its end with no length limit (excursions()), is
   finite and never comes back to the line. A bad pair's trip starts at the linear state mu^b/lambda^a (Lemma 14.3), which
   the search reaches when a + b <= depth, so this settles every bad pair with a + b <= depth.
Usage: python3 nb14.py [depth=300]. Writes nb14.json.
"""
import sys, json
from fractions import Fraction as F
from collections import Counter
L = {'l': (F(19, 200), F(0)), 'm': (F(1, 20), F(0)), 'N': (F(1, 19), F(1, 19)), 'K': (F(1, 10), F(1, 10))}
NXT = {None: 'lmNK', 'l': 'lmNK', 'm': 'mN', 'N': 'lmNK', 'K': 'lmNK'}
J1 = F(1, 9)
def ok(a, b): return b <= J1 and b + a * J1 >= 0
def run(depth):
    start = (F(1), F(0), None, None)
    seen = {start: 0}; births = []; frontier = [start]; hits = []; per = []
    born = {start: 0}; returns_to_linear = []; life = 0
    for dep in range(1, depth + 1):
        new = []
        for (a, b, lh, lg) in frontier:
            if a >= 1: moves = [((a * L[x][0], a * L[x][1] + b), x, lg, 'h', x) for x in NXT[lh]]
            else: moves = [((a / L[y][0], (b - L[y][1]) / L[y][0]), lh, y, 'g', y) for y in NXT[lg]]
            for (na, nb), nlh, nlg, side, let in moves:
                if not ok(na, nb): continue
                if b == 0 and nb != 0: births.append((a, side, let))
                if b != 0 and nb == 0: returns_to_linear.append((dep, na))
                st = (na, nb, nlh, nlg)
                if na == 1 and nb == 0: hits.append((dep, nlh, nlg)); continue
                if st in seen: continue
                seen[st] = dep; new.append(st)
                born[st] = dep if (nb == 0 or b == 0) else born[(a, b, lh, lg)]
                if nb != 0: life = max(life, dep - born[st])
        frontier = new
        per.append(dict(depth=dep, new=len(new), linear=sum(1 for s in new if s[1] == 0)))
    alive = [s for s in frontier if s[1] != 0]
    return seen, births, hits, per, returns_to_linear, life, alive, born
def excursions(a0, lh, lg, cap=10000):
    """All trips off the line from the linear state T(y) = a0*y (last letters lh, lg), followed to the end, deduplicated
    within this one tree. Returns (number of states, longest trip in letters, returns to the line, finished)."""
    seen = {}; stack = []; rets = []; longest = 0
    moves = ([(a0 * L[x][0], a0 * L[x][1], x, lg) for x in NXT[lh] if x in 'NK'] if a0 >= 1 else
             [(a0 / L[y][0], -L[y][1] / L[y][0], lh, y) for y in NXT[lg] if y in 'NK'])
    for na, nb, x, y in moves:
        if ok(na, nb): stack.append((na, nb, x, y, 1))
    while stack:
        a, b, lh2, lg2, k = stack.pop()
        if (a, b, lh2, lg2) in seen and seen[(a, b, lh2, lg2)] <= k: continue
        seen[(a, b, lh2, lg2)] = k; longest = max(longest, k)
        if len(seen) > cap: return len(seen), longest, rets, False
        if a >= 1: nxt = [(a * L[x][0], a * L[x][1] + b, x, lg2) for x in NXT[lh2]]
        else: nxt = [(a / L[y][0], (b - L[y][1]) / L[y][0], lh2, y) for y in NXT[lg2]]
        for na, nb, x, y in nxt:
            if not ok(na, nb): continue
            if nb == 0: rets.append((na, k + 1)); continue
            stack.append((na, nb, x, y, k + 1))
    return len(seen), longest, rets, True

if __name__ == "__main__":
    depth = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seen, births, hits, per, rets, life, alive, born = run(depth)
    win = {('h', 'N'): (F(1), F(19, 9)), ('h', 'K'): (F(1), F(10, 9)), ('g', 'N'): (F(9, 19), F(1)), ('g', 'K'): (F(9, 10), F(1))}
    inwin = all(win[(s, x)][0] <= a <= win[(s, x)][1] for a, s, x in births)
    nontriv = [h for h in hits if h[1] != h[2]]
    lin = sum(1 for s in seen if s[1] == 0)
    sample = [dict(a=float(s[0]), b=float(s[1]), depth=d) for s, d in seen.items()]
    out = dict(depth=depth, states=len(seen), linear_states=lin, nonlinear_states=len(seen) - lin, births=len(births),
               births_inside_windows=inwin, identity_returns=len(hits), nontrivial_identity_returns=len(nontriv),
               returns_to_linear_from_nonlinear=len(rets), longest_excursion=life,
               nonlinear_alive_at_end=len(alive), latest_birth_alive=max([born[s] for s in alive], default=None),
               per_depth=per, states_ab=sample)
    tot = lng = nret = unf = 0
    for (a0, b0, lh, lg) in seen:
        if b0 != 0 or lh is None: continue
        ns, l, r, fin = excursions(a0, lh, lg); tot += ns; lng = max(lng, l); nret += len(r); unf += (not fin)
    out.update(per_stretch_trip_states=tot, per_stretch_longest_trip=lng, per_stretch_returns=nret, per_stretch_unfinished=unf)
    json.dump(out, open("nb14.json", "w"))
    print("trips from each of the %d stretches, each followed to its end: %d states, longest %d letters, %d returns to the line, %d unfinished"
          % (lin - 1, tot, lng, nret, unf))
    print("depth", depth, "states", len(seen), "linear", lin, "nonlinear", len(seen) - lin)
    print("births of nonlinear states", len(births), "all inside the windows:", inwin)
    print("returns to the identity", len(hits), "nontrivial", len(nontriv))
    print("nonlinear -> linear returns", len(rets), " longest excursion (letters)", life,
          " nonlinear alive at the end", len(alive), " their birth depths", sorted(set(born[s] for s in alive)))
