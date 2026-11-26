"""Split Li's lambda_n into the 'prime at infinity' part and the finite-primes part.
log xi(s) = [log(s/2) - (s/2) log pi + log Gamma(s/2)]   (archimedean, 'prime at infinity')
          + [log((s - 1) zeta(s))]                        (finite primes: Euler product, continued)
Each part's coefficients by the same Cauchy integral as li.py (|z| = 0.8, 4096 points, 140 digits).
lambda_n = arch_n + prime_n. -> li_split.json"""
import json
from multiprocessing import Pool
import mpmath as mp
DPS, R, M, NMAX = 140, mp.mpf("0.8"), 4096, 1000
def parts(j):
    mp.mp.dps = DPS
    z = R * mp.expjpi(mp.mpf(2 * j) / M); s = 1 / (1 - z)
    a = s / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)
    b = (s - 1) * mp.zeta(s)
    return [(mp.nstr(v.real, DPS), mp.nstr(v.imag, DPS)) for v in (a, b)]
def unwrap_log(vals):
    out, prev = [], None
    for v in vals:
        ph = mp.arg(v)
        if prev is not None:
            while ph - prev > mp.pi: ph -= 2 * mp.pi
            while ph - prev < -mp.pi: ph += 2 * mp.pi
        prev = ph; out.append(mp.log(abs(v)) + 1j * ph)
    return out
if __name__ == "__main__":
    with Pool(4) as P: raw = P.map(parts, range(M), chunksize=64)
    mp.mp.dps = DPS
    A = unwrap_log([mp.mpc(mp.mpf(r[0][0]), mp.mpf(r[0][1])) for r in raw])
    Bv = unwrap_log([mp.mpc(mp.mpf(r[1][0]), mp.mpf(r[1][1])) for r in raw])
    roots = [mp.expjpi(-mp.mpf(2 * j) / M) for j in range(M)]
    arch, prim = [], []
    for n in range(1, NMAX + 1):
        fa = mp.fsum(A[j] * roots[(j * n) % M] for j in range(M)) / M / R ** n
        fb = mp.fsum(Bv[j] * roots[(j * n) % M] for j in range(M)) / M / R ** n
        arch.append(mp.nstr((n * fa).real, 25)); prim.append(mp.nstr((n * fb).real, 25))
    json.dump({"arch": arch, "prime": prim}, open("li_split.json", "w"))
    print("n=1", arch[0], prim[0]); print("n=1000", arch[-1], prim[-1])
