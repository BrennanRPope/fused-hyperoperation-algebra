"""Small helpers shared by figures.py: Eisenstein prime and square-free tests for
points z = u + v*w6 of the hexagon-centre lattice (w6 = e^{i pi/3}, N = u^2 + u v + v^2)."""
import numpy as np, sympy as sp

def eis_prime_mask(P):
    out = []
    for u, v in P.tolist():
        n = u*u + u*v + v*v
        if n < 2: out.append(False); continue
        if sp.isprime(n): out.append(True); continue
        q = int(round(n**0.5))
        out.append(q*q == n and sp.isprime(q) and q % 3 == 2 and u % q == 0 and v % q == 0)
    return np.array(out)

def squarefree_mask(P):
    """True where no Eisenstein prime squared divides z (z != 0)."""
    out = []
    for u, v in P.tolist():
        n = u*u + u*v + v*v
        ok = n > 0
        for p, e in sp.factorint(n).items():
            if not ok: break
            if p == 3:      # (1+w6)^2 ~ 3: square divides iff 3 | u and 3 | v
                ok = not (u % 3 == 0 and v % 3 == 0)
            elif p % 3 == 2:  # inert, norm p^2: p^2 | z iff p^2 | u and p^2 | v
                ok = not (u % (p*p) == 0 and v % (p*p) == 0)
            else:           # split: pi^2 | z for one of the two primes over p
                for r in [t for t in range(p) if (t*t - t + 1) % p == 0]:
                    q = p*p; r2 = (r - (r*r - r + 1)*pow(2*r - 1, -1, q)) % q
                    if (u + v*r2) % q == 0: ok = False
        out.append(ok)
    return np.array(out)
