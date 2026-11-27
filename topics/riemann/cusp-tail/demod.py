"""Is the swing the cusp? (1) per-length turn of Z_n(rho) against -gamma*log((n+1)/n);
(2) demodulate Z_n(rho) by the cusp phase n^(rho-1): if the swing is the cusp, the result sits still. DRAFT."""
import numpy as np, sys, json, mpmath as mp
sys.path.insert(0, "../length-smoothing"); from wide import lev
gam = [float(mp.zetazero(k).imag) for k in range(1, 7)]
if __name__ == "__main__":
    ns = np.arange(8, 34); L = {n: lev(n) for n in ns}; out = dict(gam=gam, ns=ns.tolist(), E=[], raw=[], demod=[])
    for g in gam:
        rho = 1.5+1j*g; E = np.array([(c*np.exp(-rho*np.log(q))).sum() for q, c in (L[n] for n in ns)])
        out["E"].append([[float(e.real), float(e.imag)] for e in E])
        sel = ns >= 18; D = E*ns**(rho-1)
        out["raw"].append(float(abs(np.mean(E[sel]/abs(E[sel]))))); out["demod"].append(float(abs(np.mean(D[sel]/abs(D[sel])))))
        print("gamma %5.2f  phase steadiness, lengths 18..33: raw %.2f   cusp-demodulated %.2f" % (g, out["raw"][-1], out["demod"][-1]))
    json.dump(out, open("demod.json", "w"), indent=1)
