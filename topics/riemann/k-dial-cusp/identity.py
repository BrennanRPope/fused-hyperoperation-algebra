"""The k dial on the chain step: (u, v) -> (u o_k v, v), (u, u o_k v), with a o_k b = a + b + k ab.
Claim: since 1 + k(a o_k b) = (1 + ka)(1 + kb), the tree's sums u o_k v are exactly the q-integers
[m]_q = (q^m - 1)/(q - 1), q = 1 + k, m = p + q over coprime pairs. Check on the tree. DRAFT."""
import numpy as np
from math import gcd
def tree(k, depth):
    out = []; stack = [(1.0, 1.0, 1, 1, 0)]
    while stack:
        u, v, a, b, d = stack.pop(); w = u + v + k*u*v; out.append((w, a + b, a, b))
        if d < depth: stack += [(w, v, a + b, b, d + 1), (u, w, a, a + b, d + 1)]
    return out
if __name__ == "__main__":
    for k in (0.0, 0.05, 0.3):
        T = tree(k, 14); q = 1 + k
        qint = lambda m: m if k == 0 else (q**m - 1)/k
        err = max(abs(w - qint(m))/qint(m) for w, m, a, b in T)
        cop = all(gcd(a, b) == 1 for _, _, a, b in T); uniq = len({(a, b) for _, _, a, b in T}) == len(T)
        print("k = %.2f: %d nodes, max relative gap to [p+q]_q = %.2e, all pairs coprime %s, each pair once %s" % (k, len(T), err, cop, uniq))
