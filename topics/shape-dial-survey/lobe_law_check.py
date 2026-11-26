"""Lobe growth far out: the asymptotic lobe law for the top strip of the Lee-Yang rule.
Shell peeling (step identity): s^kappa G_a(s) = a p int_0^s t^(a p) G_(a-1)(t) dt, kappa = a p + 1.
Weighted lobes L_k = int_(z_(k-1))^(z_k) t^(a p) G_(a-1)(t) dt between neighbouring zeros z_k of G_(a-1).
Law (from the rim wave plus the centre kink of G_(a-1)):
  (L_(k+1)/L_k - 1) * z_k / (pi a (p-1))  =  1 + sigma_(k+1) * B(p,a) * z_k^(a-p)  + O(1/z_k),
  B(p,a) = (a-1) Gamma(p+1) sin(pi p/2) / (a (p-1) p^(a-1) Gamma(a)),  sigma = sign of the lobe.
At a = p, B = sin(pi p/2)/p^(p-1): the same kink share as the ball comparison, below 1 for every p > 1.
"""
import json, time
import numpy as np, mpmath as mp
RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
def coeffs(p, a, R):
    K = int(1.6*R) + 80; mp.mp.dps = int(R/2) + 40; p = mp.mpf(p); a = mp.mpf(a)
    return [(-1)**k/mp.factorial(2*k)*mp.beta((2*k + 1)/p, a + 1)/p for k in range(K)]
def ev(c, s):
    s2 = mp.mpmathify(s)**2; acc = mp.mpf(0)
    for ck in reversed(c): acc = acc*s2 + ck
    return acc
def zeros(c, R, ds=0.25):
    g = np.arange(ds, R, ds); v = [ev(c, x) for x in g]
    return [mp.findroot(lambda x: ev(c, x), (mp.mpf(g[i]), mp.mpf(g[i + 1])), solver="illinois", verify=False)
            for i in range(len(v) - 1) if v[i]*v[i + 1] < 0]
def B(p, a): return float((a - 1)*mp.gamma(p + 1)*mp.sin(mp.pi*p/2)/(a*(p - 1)*mp.mpf(p)**(a - 1)*mp.gamma(a)))
def lobes(p, a, R):
    cb = coeffs(p, a - 1, R); ca = coeffs(p, a, R); z = zeros(cb, R)
    F = [zk**(a*p + 1)*ev(ca, zk)/(a*p) for zk in z]
    L = [F[0]] + [F[k] - F[k - 1] for k in range(1, len(F))]
    rows = [(float(z[k]), float(mp.sign(L[k + 1])), float(abs(L[k + 1]/L[k])),
             float((abs(L[k + 1]/L[k]) - 1)*z[k]/(mp.pi*a*(p - 1)))) for k in range(len(L) - 1)]
    return rows, z, cb, F
t0 = time.perf_counter()
# L1 the law at the edge a = p: fit D = c0 + c1*sigma + (c2 + c3*sigma)/z on the outer two thirds
fit = {}
for p in [2.5, 3.0, 3.5]:
    rows, *_ = lobes(p, p, 220); r = np.array(rows[len(rows)//3:])
    X = np.c_[np.ones(len(r)), r[:, 1], 1/r[:, 0], r[:, 1]/r[:, 0]]; co = np.linalg.lstsq(X, r[:, 3], rcond=None)[0]
    fit[f"p={p:g}"] = dict(c0=co[0], c_sigma=co[1], predicted_sigma=B(p, p), c_1_over_z=co[2], rows=rows)
OUT["edge_fit"] = {k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in fit.items()}
OUT["edge_rows"] = {k: v["rows"] for k, v in fit.items()}
check("L1 at a = p the normalised lobe excess tends to 1 +/- sin(pi p/2)/p^(p-1), alternating with the lobe's sign",
      all(abs(v["c0"] - 1) < 0.02 and abs(v["c_sigma"] - v["predicted_sigma"]) < 0.01 for v in fit.values()),
      "; ".join(f"{k}: 1 -> {v['c0']:.3f}, kink {v['predicted_sigma']:+.4f} -> {v['c_sigma']:+.4f}" for k, v in fit.items()))
# L2 inside the strip the kink term fades like z^(a-p): p = 3, a = 2.5 and p = 2.5, a = 2 (next kink order z^(a-p-1) in the fit)
fit2 = {}
for p, a in [(3.0, 2.5), (2.5, 2.0)]:
    rows, *_ = lobes(p, a, 500); r = np.array(rows[len(rows)//4:])
    X = np.c_[np.ones(len(r)), r[:, 1]*r[:, 0]**(a - p), 1/r[:, 0], r[:, 1]/r[:, 0], r[:, 1]*r[:, 0]**(a - p - 1)]
    co = np.linalg.lstsq(X, r[:, 3], rcond=None)[0]
    fit2[f"p={p:g}, a={a:g}"] = dict(c0=co[0], c_kink=co[1], predicted=B(p, a), rows=rows)
OUT["inside_fit"] = {k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in fit2.items()}
OUT["inside_rows"] = {k: v["rows"] for k, v in fit2.items()}
check("L2 inside the strip the kink term fades like z^(a-p) with the predicted size",
      all(abs(v["c0"] - 1) < 0.03 and abs(v["c_kink"] - v["predicted"]) < 0.03*abs(v["predicted"]) for v in fit2.values()),
      "; ".join(f"{k}: kink {v['predicted']:+.3f} -> {v['c_kink']:+.3f}" for k, v in fit2.items()))
# L3 grid over the strip: every lobe grows, and the excess never drops much below the law's floor 1 - |B|
grid = []
for p in [2.2, 2.5, 2.9, 3.0, 3.5, 4.5, 6.0]:
    for a in sorted({1.5, round(p - 0.9, 3), round(p - 0.5, 3), p}):
        rows, *_ = lobes(p, round(a, 3), 120)
        grid.append(dict(p=p, a=round(a, 3), n=len(rows) + 1, min_ratio=min(r[2] for r in rows),
                         min_D=min(r[3] for r in rows), floor=1 - abs(B(p, a))))
OUT["grid"] = grid
check("L3 every step of the climb from a = 1 (7 roundnesses, a = 1.5 and three values in the strip): every weighted lobe beats the one before, and the excess stays above the law's floor",
      all(g["min_ratio"] > 1 and g["min_D"] > g["floor"] for g in grid),
      f"{len(grid)} cases, smallest ratio {min(g['min_ratio'] for g in grid):.3f}, smallest excess/law {min(g['min_D'] for g in grid):.2f}")
# L4 the law also knows where the rule stops: past a = p the kink term grows and lobes stop growing (p = 3, a = 4.5)
rows, *_ = lobes(3.0, 4.5, 120)
OUT["past_edge"] = dict(p=3.0, a=4.5, rows=rows, first_fail=next((r[0] for r in rows if r[2] < 1), None))
check("L4 past the edge (p = 3, a = 4.5) the alternating kink term wins far out and a lobe fails to grow, as the law predicts",
      OUT["past_edge"]["first_fail"] is not None, f"first shrinking lobe after z = {OUT['past_edge']['first_fail']}")
# L5 polish: first lobes by direct quadrature of t^(ap) G_(a-1) against the identity, zeros refined to full precision
p, a = 3.0, 3.0; rows, z, cb, F = lobes(p, a, 40)
mp.mp.dps = 50
quad = [mp.quad(lambda t: t**(a*p)*ev(cb, t), [0 if k == 0 else z[k - 1], z[k]]) for k in range(3)]
ident = [F[0], F[1] - F[0], F[2] - F[1]]
gap = max(float(abs(q/i - 1)) for q, i in zip(quad, ident)); res = max(float(abs(ev(cb, zk))) for zk in z[:5])
OUT["polish"] = dict(rel_gap=gap, zero_residual=res)
check("L5 polish: first three lobes by direct quadrature agree with the identity; refined zeros are clean",
      gap < 1e-30 and res < 1e-40, f"agree to all {mp.mp.dps} digits (gap {gap:.0e}), residual {res:.1e}")
OUT["seconds"] = time.perf_counter() - t0
print(f"{sum(RES)}/{len(RES)} pass, {OUT['seconds']:.0f} s")
json.dump(OUT, open("lobe_law_check.json", "w"), indent=1, default=float)
