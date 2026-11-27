"""The law of the swing: at each zeta zero rho, the chain's error Z_n(rho) turns by a steady angle per length.
Measure that angle, its dependence on the height gamma, and compare the swing period with the best window. DRAFT."""
import numpy as np, json, sys, mpmath as mp
from wide import lev
gam = [float(mp.zetazero(k).imag) for k in range(1, 9)]
if __name__ == "__main__":
    L = {n: lev(n) for n in range(10, 31)}; out = dict(gam=gam, turn=[], turn_spread=[], mag=[])
    for g in gam:
        E = np.array([(c*np.exp(-(1.5+1j*g)*np.log(q))).sum() for q, c in (L[n] for n in range(10, 31))])
        r = E[1:]/E[:-1]; ang = np.angle(r[6:])                   # n = 16..30
        m = np.angle(np.mean(np.exp(1j*ang))); R = np.abs(np.mean(np.exp(1j*ang)))
        out["turn"].append(float(m)); out["turn_spread"].append(float(R)); out["mag"].append(float(np.exp(np.mean(np.log(np.abs(r[6:]))))))
        print("gamma %6.2f  turn per length %+.3f rad (steadiness %.2f), period %.1f lengths, |ratio| %.3f" % (g, m, R, 2*np.pi/abs(m), out["mag"][-1]), flush=True)
    sl = np.polyfit(gam, out["turn"], 1); out["fit"] = sl.tolist(); print("fit turn = %.4f*gamma %+.3f" % tuple(sl))
    w = json.load(open("windows.json")); best = []
    for j in range(5):
        ks = list(range(0, 13, 2)); m = [np.mean([abs(w["rows"]["tri k=%d N=%d" % (k, N)][j][0]) for N in range(16, 25)]) for k in ks]
        best.append(ks[int(np.argmin(m))] + 1)
    out["best_window_per_zero"] = best; print("best triangle window per zero (lengths):", best, " swing periods:", [round(2*np.pi/abs(t), 1) for t in out["turn"][:5]])
    json.dump(out, open("swing.json", "w"), indent=1)
