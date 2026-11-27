"""Rewind dial: rotate the field from the magnetic one (theta = 0, where the fixed chain passes Lee-Yang) to the
height field that p's Wick turn makes (theta = 90 deg, where the zeros are Knauf's temperature zeros). DRAFT.
F_theta(z) = sum over states of e^{-sigma0 H} e^{z X_theta},  X_theta = cos(theta) a M + sin(theta) H,
on the doubled (compass-spin) chain. Re z = 0 is the circle at theta = 0 and Knauf's line Re beta = 3/2 at 90 deg."""
import numpy as np, json, sys
sys.path.insert(0, "../knauf-field"); from knauf import chain
def doubled_states(n):
    phi = chain(n).astype(float); H = np.log(phi); M = n - 2*np.array([bin(k).count("1") for k in range(2**n)])
    H2 = np.concatenate([H, H[::-1]]); M2 = np.concatenate([M+1, M-1]); return H2, M2.astype(float)
def setup(n, sigma0=1.5):
    H2, M2 = doubled_states(n); W = np.exp(-sigma0*H2)
    # group identical (M, H) pairs
    key = np.round(H2, 12) + 1000*M2; u, inv = np.unique(key, return_inverse=True)
    w = np.bincount(inv, weights=W); Hu = np.bincount(inv, weights=H2)/np.bincount(inv); Mu = np.bincount(inv, weights=M2)/np.bincount(inv)
    p = W/W.sum(); a = np.sqrt((p*(H2-(p*H2).sum())**2).sum()/(p*M2**2).sum())
    return dict(w=w, H=Hu, M=Mu, a=a, n=n)
def F(S, z, th):
    X = np.cos(th)*S["a"]*S["M"] + np.sin(th)*S["H"]; e = S["w"]*np.exp(np.multiply.outer(z, X))
    return e.sum(-1), (e*X).sum(-1)
def start_zeros(S, T):
    n = S["n"]; a = S["a"]; c = np.zeros(n+2)
    j = np.rint((S["M"] + n + 1)/2).astype(int); np.add.at(c, j, S["w"])
    r = np.roots(c[::-1]); ph = np.angle(r) % (2*np.pi)
    zs = []
    for m in range(0, int(T*a/np.pi)+2):
        for p in ph:
            z = 1j*(p + 2*np.pi*m)/(2*a)
            if 0 < z.imag <= T: zs.append(z + np.log(abs(r[list(ph).index(p)]))/(2*a))
    return np.array(sorted(zs, key=lambda z: z.imag))
def newton(S, z, th, it=40):
    for _ in range(it):
        f, d = F(S, z, th); dz = f/d; z = z - dz
        if np.all(np.abs(dz) < 1e-13): break
    return z
def track(n, T=36.0, steps=900):
    S = setup(n); z = start_zeros(S, T); ths = np.linspace(0, np.pi/2, steps+1); path = [z.copy()]
    for k in range(1, len(ths)):
        zp = newton(S, z, ths[k])
        jump = np.abs(zp - z) > 0.5            # sub-step if any zero jumped
        if jump.any():
            zz = z.copy()
            for t in np.linspace(ths[k-1], ths[k], 21)[1:]: zz = newton(S, zz, t)
            zp = zz
        z = zp; path.append(z.copy())
    return S, ths, np.array(path)
if __name__ == "__main__":
    out = {}
    for n in (10, 12, 14):
        S, ths, P = track(n)
        f, _ = F(S, P[-1], np.pi/2); res = float(np.abs(f).max()/S["w"].sum())
        dev = np.abs(P.real).max(axis=1)
        print("n=%d a=%.4f zeros tracked %d | max|Re z| at 0deg %.1e, 10deg %.3f, 45deg %.3f, 90deg %.3f | resid %.1e" %
              (n, S["a"], P.shape[1], dev[0], dev[100], dev[450], dev[-1], res), flush=True)
        end = P[-1]; low = end[np.abs(end.imag) < 34]
        print("   theta=90 zeros as beta = 1.5 - z (first ones):", ["%.3f%+.3fi" % ((1.5-z).real, (1.5-z).imag) for z in sorted(low, key=lambda z: abs(z.imag))[:12]])
        out[n] = dict(a=S["a"], theta=ths.tolist(), re=P.real.tolist(), im=P.imag.tolist())
    json.dump(out, open("rotate.json", "w"))
