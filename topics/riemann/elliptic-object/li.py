"""Li / Keiper-Li coefficients lambda_n, n = 1..NMAX, straight from zeta (no zeros used).

lambda_n = n * [z^n] log xi(1/(1-z)),  xi(s) = s(s-1)/2 * pi^(-s/2) * Gamma(s/2) * zeta(s).
Coefficients by the Cauchy integral on |z| = R (trapezoid rule, M points); log xi is
analytic for |z| < 1 as long as xi has no zero with Re s > 1/2 inside the image disc,
which here is Re s in [0.556, 5], |Im s| <= 2.3 (no zeros there). -> li.json
"""
import json
from multiprocessing import Pool
import mpmath as mp

DPS, R, M, NMAX = 140, mp.mpf("0.8"), 4096, 1000

def xi_at(j):
    mp.mp.dps = DPS
    z = R * mp.expjpi(mp.mpf(2 * j) / M)
    s = 1 / (1 - z)
    v = s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)
    return (mp.nstr(v.real, DPS), mp.nstr(v.imag, DPS))

if __name__ == "__main__":
    with Pool(4) as P:
        raw = P.map(xi_at, range(M), chunksize=64)
    mp.mp.dps = DPS
    vals = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in raw]
    # continuous log: log|xi| + unwrapped phase
    logs, prev = [], None
    for v in vals:
        ph = mp.arg(v)
        if prev is not None:
            while ph - prev > mp.pi: ph -= 2 * mp.pi
            while ph - prev < -mp.pi: ph += 2 * mp.pi
        prev = ph
        logs.append(mp.log(abs(v)) + 1j * ph)
    wind = (logs[-1].imag - logs[0].imag)
    roots = [mp.expjpi(-mp.mpf(2 * j) / M) for j in range(M)]
    lam = []
    for n in range(1, NMAX + 1):
        acc = mp.fsum(logs[j] * roots[(j * n) % M] for j in range(M))
        a_n = acc / M / R ** n
        lam.append(n * a_n)
    out = {"R": str(R), "M": M, "dps": DPS, "phase_change_last_minus_first": mp.nstr(wind, 10),
           "lambda": [mp.nstr(l.real, 30) for l in lam],
           "max_abs_imag": mp.nstr(max(abs(l.imag) for l in lam), 5)}
    json.dump(out, open("li.json", "w"), indent=0)
    print("lambda_1..3", [mp.nstr(l.real, 20) for l in lam[:3]])
    print("lambda_100,500,1000", [mp.nstr(lam[i].real, 15) for i in (99, 499, 999)])
    print("min", mp.nstr(min(l.real for l in lam), 10), "max imag", out["max_abs_imag"])
