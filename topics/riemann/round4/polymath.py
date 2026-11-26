"""Round 4: a toy of the positive quantity Polymath 15 used at large x.
Their effective approximation (Theorem 1.3) writes H_t/B_t = sum_{n<=N} b_n^t / n^{s*} + gamma*(...) + errors,
b_n^t = exp((t/4) log^2 n),  Re s* ~ (1+y)/2 + (t/4) log(x/4pi),  N ~ sqrt(x/4pi).
Toy margin (B-sum only, triangle inequality, no mollifier, no A-sum, no errors):
  m(t) = 1 - sum_{n=2}^{N} n^{-(1+y)/2} exp(-(t/4) log n * log(N^2/n)).
It is a dominance margin: the first Dirichlet term must beat all the others.  -> polymath.json"""
import numpy as np, json
def margin(t, x, y):
    N = int(np.sqrt(x/(4*np.pi))); n = np.arange(2, N+1, dtype=float); ln = np.log(n)
    return 1 - np.sum(n**(-(1+y)/2)*np.exp(-(t/4)*ln*(2*np.log(N)-ln)))
xs = [1e6, 1e8, 1e10, 6e10, 1e12, 1e14]
ts = np.linspace(0, 0.5, 101)
res = {'y': 0.2, 'xs': xs, 'ts': ts.tolist(), 'm': {str(x): [margin(t, x, 0.2) for t in ts] for x in xs}}
# t where the toy margin turns positive, for each x
res['t_cross'] = {}
for x in xs:
    m = np.array(res['m'][str(x)]); k = np.argmax(m > 0)
    res['t_cross'][str(x)] = float(ts[k]) if m[k] > 0 else None
res['m_at_t0'] = {str(x): margin(0.0, x, 0.2) for x in xs}
res['m_at_0.2'] = {str(x): margin(0.2, x, 0.2) for x in xs}
print(res['t_cross']); print(res['m_at_0.2']); print(res['m_at_t0'])
json.dump(res, open('polymath.json', 'w'))

# Euler-mollified toy: multiply the B-sum by prod_{p in P}(1 - b_p p^{-s}) and redo the triangle bound
from itertools import combinations
def margin_moll(t, x, y, P=(2, 3, 5, 7)):
    N = int(np.sqrt(x/(4*np.pi))); M = N*int(np.prod(P))
    b = lambda n: np.exp((t/4)*np.log(n)**2)
    c = np.zeros(M+1); n = np.arange(1, N+1, dtype=float); bn = b(n)
    for r in range(len(P)+1):
        for S in combinations(P, r):
            d = int(np.prod(S)) if S else 1
            c[d*np.arange(1, N+1)] += (-1)**r*np.prod([b(p) for p in S])*bn
    m = np.arange(2, M+1, dtype=float)
    w = m**(-(1+y)/2)*np.exp(-(t/2)*np.log(N)*np.log(m))
    return c[1] - np.sum(np.abs(c[2:])*w)
res['moll'] = {str(x): [margin_moll(t, x, 0.2) for t in ts[::5]] for x in [1e8, 6e10]}
res['moll_ts'] = ts[::5].tolist()
for x in [1e8, 6e10]:
    m = np.array(res['moll'][str(x)]); k = np.argmax(m > 0)
    res['t_cross_moll_' + str(x)] = float(ts[::5][k]) if m[k] > 0 else None
print({k: v for k, v in res.items() if k.startswith('t_cross_moll')})
json.dump(res, open('polymath.json', 'w'))
