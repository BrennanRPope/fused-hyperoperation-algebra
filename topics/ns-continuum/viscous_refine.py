"""Refined check near the dyadic transition at alpha = 1: does the cascade's deepest shell
saturate as N grows (regular) or follow N (blow-up)? Compare with the K41 estimate
j_d = -ln(nu) / ((2 gam - 2/3) ln lam)."""
import numpy as np, json, time
from viscous_scan import outcome
nu = 1e-2; res = []
for n in [3.0, 3.5, 3.75, 4.0, 4.25, 4.5]:
    theta = 1 + n / 2; lam = 2 ** theta; gam = 1 / theta
    jd = -np.log(nu) / ((2 * gam - 2 / 3) * np.log(lam)) if gam > 1/3 + 1e-12 else float('inf')
    row = []
    for N in [24, 32, 40, 48, 56]:
        t0 = time.time(); o = outcome(n, 1.0, N, nu, Tmax=10.0)
        o.update(n=n, N=N, nu=nu, gam=gam, jd=jd, secs=time.time() - t0); res.append(o); row.append(o)
    print(f"n={n} gam={gam:.4f} K41 j_d~{jd:.1f}: " + " | ".join(
        f"N={o['N']} {'REACH' if o['reached'] else 'stop'} deep={o['deepest']} t={o['t']:.4f}" for o in row), flush=True)
json.dump(res, open('viscous_refine.json', 'w'))
