"""Four set-operations facts as levers (Nov 26). Numbers used by checks_setops.py and fig7."""
import itertools, json, math
import numpy as np
from scipy.optimize import linprog
from scipy.special import gamma, erf

def E_k(k):
    """Mean number of uniform draws until x1 o_k x2 o_k ... passes 1 (closed form)."""
    if abs(k) < 1e-12: return math.e
    return (1 + k) ** (1 / k) + k / (1 + k)

def E_k_mc(k, runs=200000, seed=1):
    rng = np.random.default_rng(seed)
    t = math.log1p(k) if k != 0 else 1.0
    tot = np.zeros(runs); n = np.zeros(runs, int); alive = np.ones(runs, bool)
    while alive.any():
        x = rng.random(alive.sum())
        tot[alive] += np.log1p(k * x) if k != 0 else x
        n[alive] += 1
        alive[alive] = (tot[alive] <= t) if k >= 0 else (tot[alive] >= t)
    return n.mean()

def draws_p(p, N=400):
    """Sum over n of the cube-corner volume of the l_p ball = mean draws until sum x^p passes 1."""
    g = gamma(1 + 1 / p)
    return sum(math.exp(n * math.log(g) - math.lgamma(1 + n / p)) for n in range(N))

def neg_mass(S):
    """Least total negative weight on the 16 deterministic strategies reproducing
    isotropic CHSH correlations of value S (zero marginals)."""
    lams = list(itertools.product([1, -1], repeat=4))  # a0 a1 b0 b1
    E = {(0, 0): S / 4, (0, 1): S / 4, (1, 0): S / 4, (1, 1): -S / 4}
    A, b = [], []
    for x, y in E:
        for a in (1, -1):
            for bb in (1, -1):
                row = [1.0 if (l[x] == a and l[2 + y] == bb) else 0.0 for l in lams]
                A.append(row + [-r for r in row]); b.append((1 + a * bb * E[(x, y)]) / 4)
    c = [0] * 16 + [1] * 16
    r = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)] * 32, method="highs")
    return r.fun

if __name__ == "__main__":
    out = {"E_k": {str(k): E_k(k) for k in [-0.5, 0, 0.5, 1, math.e - 1, 3]},
           "E_k_mc": {str(k): E_k_mc(k) for k in [-0.5, 0.5, 1, math.e - 1]},
           "draws_p": {str(p): draws_p(p) for p in [1, 2, 3]},
           "neg_mass": {f"{S:.4f}": neg_mass(S) for S in [2.0, 2.4, 2 * math.sqrt(2), 3.2, 3.6, 4.0]}}
    json.dump(out, open("setops_levers.json", "w"), indent=1); print(json.dumps(out, indent=1))
