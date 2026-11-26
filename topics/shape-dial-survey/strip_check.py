"""The top strip p-1 < a <= p of the Lee-Yang rule, seen through the shape dial.
Tool 1, the roundness dial at fixed a: slide p -> 2.  The rim (x near 1) looks the same at every roundness, so far out
G at roundness p matches (p/2)^a times the ball's G.  Where |G_p - (p/2)^a G_2| < |(p/2)^a G_2| on a circle (Rouche),
G_p has exactly as many zeros inside as the ball, and the ball's are all real.  Beyond that radius Lee-Yang is automatic.
Tool 2, shell peeling (the step identity): s^kappa G_a = a p * int_0^s t^(kappa-1) G_(a-1)(t) dt, kappa = a p + 1.
Inside the radius, G_a gets a zero between neighbouring zeros of G_(a-1) when the weighted lobes
L_k = int_(z_(k-1))^(z_k) t^(kappa-1) G_(a-1) grow in size.
"""
import json, time, sys
import numpy as np
import mpmath as mp
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
def real_zeros(p, a, R, ds=0.05, polish=False):
    c = coeffs(p, a, R); g = np.arange(ds, R + 1e-9, ds); v = [ev(c, x) for x in g]
    br = [(g[i], g[i + 1]) for i in range(len(v) - 1) if v[i]*v[i + 1] < 0]
    if not polish: return [float((l + r)/2) for l, r in br], c
    return [float(mp.findroot(lambda x: ev(c, x), (mp.mpf(l), mp.mpf(r)), solver="anderson")) for l, r in br], c
def crest(a, target):
    m = round((target - np.pi*(a/2 + 1))/np.pi); return float(np.pi*(a/2 + 1) + m*np.pi - np.pi/2)
def margin(p, a, r, n=360):
    cp = coeffs(p, a, r); c2 = coeffs(2.0, a, r); f = mp.mpf(p/2)**a; w = 0
    for k in range(n + 1):
        s = r*mp.expj(mp.pi*k/(2*n)); b = f*ev(c2, s); w = max(w, float(abs(ev(cp, s) - b)/abs(b)))
    return w
t0 = time.perf_counter()
# S1 ball comparison: margin on crest circles
PTS = [(2.5, 2.0), (2.5, 2.5), (3.0, 2.5), (3.0, 3.0), (5.0, 4.5), (5.0, 5.0), (7.0, 7.0)]
TARGETS = {2.5: [20, 40, 70], 3.0: [20, 40, 70], 5.0: [40, 70, 110], 7.0: [110, 160, 220]}
mg = {}
for p, a in PTS:
    mg[f"p={p:g}, a={a:g}"] = {f"{crest(a, t):.1f}": margin(p, a, crest(a, t)) for t in TARGETS[p]}
    print("    margin", p, a, {k: round(v, 3) for k, v in mg[f"p={p:g}, a={a:g}"].items()}, flush=True)
OUT["margin"] = mg
def far_radius(d):
    ks = sorted(d, key=float); good = [k for i, k in enumerate(ks) if all(d[kk] < 0.9 for kk in ks[i:])]
    return float(good[0]) if good else None
OUT["far_radius"] = {k: far_radius(v) for k, v in mg.items()}
check("S1 far zone: past a radius, the roundness-p function stays within 90% of the ball's size on every tested circle (p = 2.5, 3, 5, 7)",
      all(OUT["far_radius"][k] is not None for k in mg),
      "; ".join(f"{k}: from r = {v}" for k, v in OUT["far_radius"].items()))
# S2 near zone: inside that radius, G_p has as many real zeros as the ball (so all its zeros there are real)
near = {}
for p, a in PTS:
    R = OUT["far_radius"][f"p={p:g}, a={a:g}"]; zp, _ = real_zeros(p, a, R); z2, _ = real_zeros(2.0, a, R)
    near[f"p={p:g}, a={a:g}"] = dict(R=R, real_p=len(zp), ball=len(z2))
OUT["near"] = near
check("S2 near zone: inside the far-zone radius, G at roundness p has exactly as many real zeros as the ball",
      all(v["real_p"] == v["ball"] for v in near.values()),
      "; ".join(f"{k}: {v['real_p']} = {v['ball']}" for k, v in near.items()))
# S3 shell peeling: lobes of G_(a-1) weighted by t^(kappa-1) grow, from the first lobe, out to the far-zone radius
lob = {}
for p, a in PTS:
    R = OUT["far_radius"][f"p={p:g}, a={a:g}"]; kap = a*p + 1
    za, _ = real_zeros(p, a - 1, R, polish=True); cA = coeffs(p, a, R)
    Fz = [mp.mpf(z)**kap*ev(cA, z)/(a*p) for z in za]
    L = [Fz[0]] + [Fz[k] - Fz[k - 1] for k in range(1, len(Fz))]
    ratios = [float(abs(L[k + 1]/L[k])) for k in range(len(L) - 1)]
    lob[f"p={p:g}, a={a:g}"] = dict(n=len(L), min_ratio=min(ratios), log10_abs_L=[float(mp.log10(abs(x))) for x in L],
                                   alternates=all(Fz[k]*Fz[k + 1] < 0 for k in range(len(Fz) - 1)))
OUT["lobes"] = lob
check("S3 shell peeling: the weighted lobes grow from the first one out to the far-zone radius, so G at a changes sign "
      "between every pair of zeros of G at a - 1",
      all(v["min_ratio"] > 1 and v["alternates"] for v in lob.values()),
      "smallest growth ratio " + f"{min(v['min_ratio'] for v in lob.values()):.2f}")
# S4 the whole strip at p = 3 on a grid of a: margin on the circle r = crest near 40, and the near-zone count
strip = []
for a in np.round(np.arange(2.05, 3.0001, 0.05), 3):
    r = crest(a, 40); m = margin(3.0, a, r, n=180); zp, _ = real_zeros(3.0, a, r); z2, _ = real_zeros(2.0, a, r)
    strip.append(dict(a=float(a), r=r, margin=m, real_p=len(zp), ball=len(z2)))
OUT["strip_p3"] = strip
check("S4 whole strip at p = 3 (a from 2.05 to 3 in steps of 0.05): margin under 0.9 on a circle near r = 40 and equal counts inside",
      all(d["margin"] < 0.9 and d["real_p"] == d["ball"] for d in strip),
      f"{len(strip)} values of a, largest margin {max(d['margin'] for d in strip):.3f}")
# S5 polish: zeros refined with findroot, and the first two lobes by direct quadrature against the identity
tp = time.perf_counter(); worst = 0
for p, a in [(3.0, 2.5), (5.0, 4.5)]:
    kap = a*p + 1; za, c1 = real_zeros(p, a - 1, 20, polish=True); cA = coeffs(p, a, 20); mp.mp.dps = 40
    for k in range(2):
        lo = 0 if k == 0 else za[k - 1]
        d = mp.quad(lambda t: t**(kap - 1)*ev(c1, t), [lo, za[k]])
        i_ = (mp.mpf(za[k])**kap*ev(cA, za[k]) - (0 if k == 0 else mp.mpf(lo)**kap*ev(cA, lo)))/(a*p)
        worst = max(worst, float(abs(d - i_)/abs(i_)))
OUT["polish_seconds"] = round(time.perf_counter() - tp, 2)
check("S5 polish: lobes by direct quadrature agree with the identity at refined zeros", worst < 1e-10,
      f"rel gap {worst:.1e}, {OUT['polish_seconds']} s")
OUT["seconds"] = round(time.perf_counter() - t0, 1)
json.dump(OUT, open("strip_check.json", "w"), indent=1, default=float)
print(f"\n{sum(RES)} of {len(RES)} checks pass  ({OUT['seconds']} s)")
