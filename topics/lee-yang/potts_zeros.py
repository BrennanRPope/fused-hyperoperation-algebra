from potts_strips import *
import mpmath as mp, time, pickle, numpy as np
mp.mp.dps = 80
res = {}
for L, ns in ((1, [12, 24, 40]), (2, [8, 14, 20]), (3, [6, 10, 14]), (4, [5, 8, 10])):
    for n in ns:
        t = time.time()
        p = cyclic_strip_chromatic(L, n)
        while p[-1] == 0: p.pop()
        roots = mp.polyroots([mp.mpf(c) for c in p[::-1]], maxsteps=2000, extraprec=2000)
        zs = [complex(r) for r in roots]
        res[(L, n)] = zs
        real = sorted(round(z.real, 4) for z in zs if abs(z.imag) < 1e-9)
        # value at Q=-1 (parity) and sign checks
        print(f"L={L} n={n} deg={len(p)-1} t={time.time()-t:.1f}s real roots {real}  P(-1)={evalp(p,-1)}  min|Q-x| x=0,1,2,3: "
              + " ".join(f"{min(abs(z-x) for z in zs):.3g}" for x in (0, 1, 2, 3)), flush=True)
pickle.dump(res, open('potts_zeros.pkl', 'wb'))
