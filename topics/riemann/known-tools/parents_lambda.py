"""Lambda is monotone under Polya universal factors: if K = P * U with U(u) = phi(iu), phi in the Laguerre-Polya
class (a universal factor), then e^{tu^2} P real-rooted => e^{tu^2} K real-rooted, so Lambda_K <= Lambda_P.
With Rodgers-Tao (Lambda_Phi >= 0): every parent of Phi has Lambda_P >= 0.
Here: a local version in the window z < 68 for the bath parents P = Phi / cosh(J u) (U = cosh(Ju) = cos(iJu)):
t*(J) = the largest heat time at which H_t of the parent still has a non-real zero in the window."""
import json, sys
import numpy as np
sys.path.insert(0, "../unified-hypothesis")
import kern
K = kern.riemann(500)
res = {}
for J in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
    wt = (lambda u, J=J: 1 / np.cosh(J * u))
    tstar = None
    for t in np.round(np.arange(4.0, -6.01, -0.1), 2):
        cz = kern.complex_zeros(K, 74.0, 16.0, t=float(t), weight=wt)
        if any((c.imag > 1e-4) and (c.real < 68) for c in cz):
            tstar = float(t); break
    res[str(J)] = tstar
    print(J, tstar, flush=True)
json.dump(res, open("parents_lambda.json", "w"), indent=1)
