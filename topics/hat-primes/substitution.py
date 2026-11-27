"""Part A and B: the hat's metatile substitution matrix, its eigenvalues, frequencies,
and what its counts do modulo each prime. Exact arithmetic (sympy)."""
import json, sympy as sp
import hatgen as h

KINDS = ['H', 'T', 'P', 'F']
# children of each new metatile, read off the construction (hatgen NEW_*_CH)
def row(chs): return [sum(1 for c in chs if h.PATCH_KINDS[c] == k) for k in KINDS]
A = sp.Matrix([row(h.NEW_H_CH), row(h.NEW_T_CH), row(h.NEW_P_CH), row(h.NEW_F_CH)])  # A[new][old]
x = sp.symbols('x')
cp = sp.factor(A.charpoly(x).as_expr())
phi = (1 + sp.sqrt(5)) / 2
lam = sp.nsimplify(phi**4)
# Perron vector: frequencies of H,T,P,F metatiles (left/right: counts_{n+1} = A counts_n)
# frequencies of metatile types = left Perron vector (row u with u A = lam u)
ev = [v for (val, mult, vecs) in A.T.eigenvects() if sp.simplify(val - lam) == 0 for v in vecs][0]
ev = sp.simplify(ev / sum(ev))
hats_per = sp.Matrix([[4, 1, 2, 2]])  # hats in a level-0 H,T,P,F
refl_per = sp.Matrix([[1, 0, 0, 0]])
hats_total = sp.simplify((hats_per * ev)[0])
refl = sp.simplify((refl_per * ev)[0])
ratio_unrefl_refl = sp.nsimplify(sp.simplify((hats_total - refl) / refl))

# counts by level from the matrix, compared with the generator (levels 0..4)
v = sp.Matrix([1, 0, 0, 0])   # one H metatile
seq = []
M = sp.eye(4)
for n in range(9):
    tot = (M * hats_per.T)[0, 0]  # hats in a level-n H: (A^n hats)_H
    rr = (M * refl_per.T)[0, 0]
    seq.append((int(tot), int(rr)))
    M = A * M
gen = [len(h.hats_of(q[0])) for q in h.levels(4)]

# B: periods of the H-count sequence mod p
def period_mod(p, nmax=200000):
    Ap = [[int(a) % p for a in A.row(i)] for i in range(4)]
    def mm(X, Y): return [[sum(X[i][k]*Y[k][j] for k in range(4)) % p for j in range(4)] for i in range(4)]
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    X = Ap
    for n in range(1, nmax):
        if X == I: return n
        X = mm(X, Ap)
    return None

primes = list(sp.primerange(2, 400))
per = {}
for p in primes:
    T = period_mod(p)
    if p == 5: cls = 'ramified'
    elif p % 5 in (1, 4): cls = 'split'
    else: cls = 'inert'
    per[p] = dict(period=T, cls=cls, divides_pminus1=(T is not None and (p - 1) % T == 0),
                  divides_pplus1=(T is not None and (p + 1) % T == 0))

out = dict(A=[[int(a) for a in A.row(i)] for i in range(4)], charpoly=str(cp),
           perron=str(lam), perron_float=float(lam),
           freq_HTPF=[str(sp.radsimp(e)) for e in ev], freq_float=[float(e) for e in ev],
           hats_per_metatile=str(sp.radsimp(hats_total)),
           ratio_unreflected_to_reflected=str(ratio_unrefl_refl), ratio_float=float(ratio_unrefl_refl),
           level_counts_matrix=seq, level_counts_generator=gen, periods=per)
json.dump(out, open('substitution.json', 'w'), indent=1)
print('A =', out['A']); print('charpoly', cp); print('Perron', lam, float(lam))
print('freq', out['freq_HTPF']); print('unreflected:reflected =', ratio_unrefl_refl, float(ratio_unrefl_refl))
print('counts (matrix)', seq[:6]); print('counts (generator)', gen)
bad = [p for p, d in per.items() if d['cls'] == 'split' and not d['divides_pminus1']]
bad2 = [p for p, d in per.items() if d['cls'] == 'inert' and not d['divides_pplus1']]
print('split primes whose period does not divide p-1:', bad)
print('inert primes whose period does not divide p+1:', bad2)
print({p: (d['period'], d['cls']) for p, d in list(per.items())[:16]})
