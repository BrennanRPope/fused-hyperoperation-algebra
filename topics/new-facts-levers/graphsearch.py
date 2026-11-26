"""Random holographic graph models: does any min-cut entropy vector make the corner sum negative at a k above the
best flower graph's root?  Boundary vertices = n named regions + O; up to 5 internal vertices; brute-force min cut."""
import numpy as np, itertools, json, sys
rng = np.random.default_rng(1)
def entropies(W, n, nin):
    # W: symmetric weights on V = n+1 boundary + nin internal
    B = n + 1; V = B + nin
    ins = list(itertools.product([0, 1], repeat=nin))
    S = {}
    for m in range(1, 1 << n):
        side = np.zeros(V, int)
        for i in range(n): side[i] = m >> i & 1
        best = 1e18
        for a in ins:
            side[B:] = a
            cut = W[np.ix_(side == 1, side == 0)].sum()
            best = min(best, cut)
        S[m] = best
    return S
def f(k, S, n):
    return sum(k ** (bin(m).count('1') - 1) * v for m, v in S.items())
def search(n, k, trials):
    worst = (1e9, None)
    for t in range(trials):
        nin = rng.integers(1, 6); V = n + 1 + nin
        W = rng.exponential(1, (V, V)) * (rng.random((V, V)) < rng.uniform(0.3, 1))
        W = np.triu(W, 1); W = W + W.T
        S = entropies(W, n, nin); s = sum(v for m, v in S.items() if bin(m).count('1') == 1)
        if s <= 1e-9: continue
        val = f(k, S, n) / s
        if val < worst[0]: worst = (val, None)
    return worst[0]
out = {}
for n, ks in [(5, [-0.3819, -0.37, -0.36]), (6, [-0.3293, -0.32, -0.30])]:
    for k in ks:
        out[f'{n},{k}'] = search(n, k, int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
        print(n, k, 'min corner sum / s over random graphs:', out[f'{n},{k}'])
json.dump(out, open('graphsearch.json', 'w'), indent=1)
