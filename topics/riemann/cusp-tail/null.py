"""Null test: the same fit with the wrong cusp shape should not help. DRAFT."""
import numpy as np, json, correct as C
shapes = {"right n^(1-b)": lambda n, b: n**(1-b), "no phase n^(-1/2)": lambda n, b: n**(-0.5+0j),
          "wrong phase n^(1-b)*n^(0.5i)": lambda n, b: n**(1-b+0.5j), "wrong phase n^(1-b)*n^(3i)": lambda n, b: n**(1-b+3j),
          "wrong phase n^(1-b)*n^(10i)": lambda n, b: n**(1-b+10j), "mirrored phase n^(1-conj b)": lambda n, b: n**(1-np.conj(b)), "geometric 0.9^n": lambda n, b: 0.9**n+0j}
out = {}
for name, f in shapes.items():
    def F(ns, b, model, f=f):
        M = np.array([[1, f(n, b)] for n in ns], complex); y = np.array([C.Z(n, b) for n in ns])
        return np.linalg.lstsq(M, y, rcond=None)[0][0]
    C.F = F; res = []
    for N in range(27, 34):
        ns = list(range(N-N//2, N+1)); z = [C.zero(ns, g, None) for g in C.gam]
        res.append(float(np.sqrt(np.mean([(v.real-1.5)**2 for v in z]))))
    out[name] = res; print("%-30s mean miss N=27..33: %.5f  %s" % (name, np.mean(res), np.round(res, 5)), flush=True)
json.dump(out, open("null.json", "w"), indent=1)
