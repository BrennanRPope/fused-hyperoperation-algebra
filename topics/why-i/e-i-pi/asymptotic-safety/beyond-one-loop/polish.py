"""Polish for checks.py L8 (Nov 26): 'the UV exponents turn real at a finite N_V' is read off a grid
(N_V = 0..59 in steps of 1, then 70 log-spaced points to 1e6, about 15% apart): checks.py reports the first grid point
where both exponents are real. Here the crossing itself is found by bisection on N_V of the discriminant
Delta = (tr J/2)^2 - det J of the stability matrix at the UV fixed point (exponents real iff Delta >= 0), with the
fixed point re-converged by Newton at every trial N_V, the same way checks.py follows it. numpy only (flows.py)."""
import time
import numpy as np
from flows import eh_beta, jac, newton
t0 = time.process_time()
NV = np.concatenate([np.arange(0, 60, 1.0), np.logspace(np.log10(60), 6, 70)])
v = newton(eh_beta(), [0.7, 0.2]); prev = None
Dl = lambda J: (np.trace(J)/2)**2 - np.linalg.det(J)
for N in NV:                                          # same continuation as checks.py, stopped at the first real grid point
    bt = eh_beta(-4*N, 2*N); v = newton(bt, v, 40); th = np.linalg.eigvals(-jac(bt, v))
    if np.all(np.abs(th.imag) < 1e-9): break
    prev = (N, v.copy())
grid_N = N; lo, vlo = prev; hi = N
def disc(n, vstart):
    bt = eh_beta(-4*n, 2*n); vv = newton(bt, vstart, 40); return Dl(jac(bt, vv)), vv
dlo, _ = disc(lo, vlo); dhi, _ = disc(hi, vlo)
for _ in range(60):
    mid = 0.5*(lo + hi); dm, vm = disc(mid, vlo)
    if np.sign(dm) == np.sign(dlo): lo, dlo, vlo = mid, dm, vm
    else: hi = mid
    if hi - lo < 1e-9*hi: break
Nstar = 0.5*(lo + hi)
# sensitivity to the finite-difference step inside jac (h = 1e-7): redo the last bisection with h = 1e-5
def disc_h(n, vstart, h):
    bt = eh_beta(-4*n, 2*n); vv = newton(bt, vstart, 40); return Dl(jac(bt, vv, h))
a, b = Nstar*0.99, Nstar*1.01; da = disc_h(a, vlo, 1e-5)
for _ in range(60):
    m = 0.5*(a + b); dm = disc_h(m, vlo, 1e-5)
    if np.sign(dm) == np.sign(da): a, da = m, dm
    else: b = m
Nstar_h = 0.5*(a + b)
lines = ["Polish of checks.py L8: where the UV exponents turn real (hybrid Litim graviton + DEP matter, N_V gauge fields)", "",
         f"grid answer (first grid point with real exponents): N_V ~ {grid_N:.0f}   (previous grid point {prev[0]:.1f}, still complex)",
         f"bisection on the discriminant: N_V* = {Nstar:.4f}",
         f"same with finite-difference step 1e-5 instead of 1e-7 in jac: {Nstar_h:.4f}  (agree to {abs(Nstar-Nstar_h):.1e})",
         f"so 'turn real near N_V ~ {grid_N:.0f}' -> 'near N_V ~ {Nstar:.0f}' (polished)",
         "Verdict unchanged: they still turn real at a finite N_V and still head to {2, 4}.",
         f"CPU time: {time.process_time() - t0:.2f} s"]
out = "\n".join(lines); print(out); open('polish_output.txt', 'w').write(out + "\n")
