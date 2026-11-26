"""Far coefficients (up to k ~ M/2) of F_t with a radix-2 FFT in mpmath. Usage: far2.py t M rho [dps]"""
import json, sys
import mpmath as mp
sys.path.insert(0, '.')
import flow
from flow import F
mp.mp.dps = int(sys.argv[4]) if len(sys.argv) > 4 else 50
flow.NSTEP = 1500
# rebuild the Fatou constants at the working precision (flow.py builds them at 40 digits)
flow.C = mp.mpf(flow.C_FR.numerator) / flow.C_FR.denominator
flow.A = [mp.mpf(v.numerator) / v.denominator for v in flow.A_FR]


def fft(a):
    n = len(a)
    if n == 1:
        return a
    ev, od = fft(a[0::2]), fft(a[1::2])
    out = [0] * n
    for k in range(n // 2):
        w = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]
        out[k] = ev[k] + w
        out[k + n // 2] = ev[k] - w
    return out


t, M, rho = mp.mpf(sys.argv[1]), int(sys.argv[2]), mp.mpf(sys.argv[3])
vals = [F(t, rho * mp.expjpi(2 * mp.mpf(j) / M)) if sys.argv[1] != 'f' else 0 for j in range(M)]
co = fft(vals)
re = [co[k].real / M / rho ** k for k in range(M // 2)]
negs = [k for k, v in enumerate(re) if v < 0]
print(sys.argv[1], 'M', M, 'first negatives:', negs[:8], 'count:', len(negs), flush=True)
json.dump({"t": sys.argv[1], "dps": mp.mp.dps, "coeffs": [mp.nstr(v, 12) for v in re], "neg": negs},
          open("far2_t%s.json" % sys.argv[1], "w"))
