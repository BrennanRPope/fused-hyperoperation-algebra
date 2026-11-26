"""Probe-spin (precession) picture of zeta: a probe spin coupled to the primon gas precesses at the
frequencies log n with weights n^(-beta); its coherence is L_beta(t) = zeta(beta + i t)/zeta(beta)
(the Wei–Liu / Peng et al. reading: coherence = partition function at an imaginary field).
The bath only exists for beta > 1 (sum n^-beta converges). -> probe.json"""
import json, mpmath as mp
mp.mp.dps = 20
ts = [i * 0.05 for i in range(0, 1001)]          # t in [0, 50]
out = {"t": ts, "abs_zeta": {}}
for beta in (2.0, 1.5, 1.1, 0.75, 0.5):
    out["abs_zeta"][str(beta)] = [float(abs(mp.zeta(mp.mpc(beta, t)))) for t in ts]
# coherence minimum over t in (0, 50] for beta > 1 (never zero: Euler product converges, no factor vanishes)
mins = {}
for beta in (2.0, 1.5, 1.1):
    v = out["abs_zeta"][str(beta)][1:]
    mins[str(beta)] = min(v) / float(mp.zeta(beta))
out["min_coherence"] = mins
# lower bound for beta > 1: |zeta(beta+it)| >= zeta(2 beta)/zeta(beta)
out["lower_bound"] = {str(b): float(mp.zeta(2 * b) / mp.zeta(b) / mp.zeta(b)) for b in (2.0, 1.5, 1.1)}
json.dump(out, open("probe.json", "w"))
print(mins, out["lower_bound"])
