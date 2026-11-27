"""Pólya bumps e^{-u^p} as Laguerre-multiplier parents of Phi.
r_k(p) = m_{2k}(Phi) / m_{2k}(e^{-u^p}),  m_{2k}(e^{-u^p}) = Gamma((2k+1)/p)/p.
p = 2 is the classical case: r_k(2) multiplier <=> RH (Pólya 1927; Turán = Csordas-Norfolk-Varga 1986;
Jensen = GORZ 2019). Even p >= 4 gives a different sufficient condition (Pólya: e^{-u^p} transform real-rooted).
Turán and Jensen are invariant under rescaling u -> lambda u (r_k -> q^k r_k)."""
import json, mpmath as mp
from multiplier import Phi, jensen_check
mp.mp.dps = 50
KMAX = 120
mP = [mp.quad(lambda u: u**(2*k) * Phi(u), [0, 0.25, 0.5, 0.75, 1, 1.25, 1.6]) for k in range(KMAX + 1)]
out = {}
for p in [2, 4, 6, 8, 12]:
    r = [mP[k] / (mp.gamma(mp.mpf(2*k + 1) / p) / p) for k in range(KMAX + 1)]
    tur = [k for k in range(1, KMAX) if r[k]**2 - r[k-1]*r[k+1] < 0]
    jb = jensen_check(r, nmax=40)
    out[str(p)] = {'turan_fail': tur, 'jensen_fail': jb}
    print(f'p={p}: Turán fails at {len(tur)} of {KMAX-1} k; first fails {tur[:8]} last {tur[-3:] if tur else []}; non-real Jensen degrees {[d for d,_ in jb][:10]}', flush=True)
json.dump(out, open('multiplier2.json', 'w'), indent=1)
