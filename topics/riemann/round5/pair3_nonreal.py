"""Round 5, pair 3 (second part): the non-real zeros that the height detector points at.
Newton on the series of F_p at 90 digits, seeded at the detector's peak (x = where h* > 1) and
heights y = 1 .. 5, plus the next seeds along the row. Writes into pair3.json."""
import json, mpmath as mp
d = json.load(open('pair3.json'))
def Fc(p, z, N=500):
    return 2/mp.mpf(p)*mp.fsum((-1)**k*mp.gamma((2*k+1)/mp.mpf(p))/mp.factorial(2*k)*z**(2*k) for k in range(N))
mp.mp.dps = 90
out = {}
for p in ['2.5', '3.0']:
    peaks = [g['at'] for g in d['kernels'][p]['gaps'] if g['hmax'] > 1]
    seeds = [mp.mpc(x0, y0) for x0 in peaks for y0 in [1, 1.5, 2, 2.5, 3.5]] + [mp.mpc(7, 2), mp.mpc(9.5, 3.75), mp.mpc(11.5, 5)]
    found = []
    for z0 in seeds:
        try:
            r = mp.findroot(lambda z: Fc(float(p), z), z0)
        except Exception:
            continue
        rr = (round(float(r.real), 6), round(abs(float(r.imag)), 6))
        if rr[1] > 1e-6 and rr[0] < 14 and all(abs(rr[0]-a)+abs(rr[1]-b) > 1e-4 for a, b in found):
            found.append(rr)
    out[p] = sorted(found)
d['nonreal_zeros'] = out
d['detector_peaks'] = {p: [g['at'] for g in d['kernels'][p]['gaps'] if g['hmax'] > 1] for p in ['2.5', '3.0']}
json.dump(d, open('pair3.json', 'w'), indent=1)
print(out, d['detector_peaks'])
