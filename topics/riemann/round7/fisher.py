"""Round 7, item 1b: Fisher zeros of the order-dial prime gas Z_q(beta), complex beta.

Continuation: Z_q = (q;q)_inf sum_k q^k/(q;q)_k exp( sum_N b_k(N)/N logzeta(N beta) ),
with logzeta continued along horizontal lines from Re = 4 (the usual convention; under the
numerically verified RH up to these heights, Re beta > 1/2 is free of cuts).
Fisher zeros = zeros of Z_q in Re beta > 1/2.  At q = 0 they are the Riemann zeros (on the edge).
"""
import json, math
import numpy as np
import mpmath as mp
from fastl import hurwitz
from sympy import primerange, divisors, mobius

mp.mp.dps = 15
PR = np.array(list(primerange(2, 3000)), dtype=float)
LOGP = np.log(PR)

def bk_table(q, K, NMAX):
    B = np.zeros((K+1, NMAX+1))
    for k in range(K+1):
        for N in range(1, NMAX+1):
            B[k, N] = sum(mobius(N//m) * q**(k*m)/(1-q**m) for m in divisors(N))
    return B

def qp(q, n):
    r = 1.0
    for j in range(1, n+1): r *= (1-q**j)
    return r

def row_logzeta(t, sig_grid, scale=1, start=4.0, step=0.01):
    """continued log zeta(scale*(sig+i t)) for sig in sig_grid (descending fine path from start)."""
    smin = sig_grid.min()
    path = np.arange(start, smin - 1e-12, -step)
    path = np.append(path, smin) if path[-1] > smin + 1e-12 else path
    vals = [complex(mp.zeta(mp.mpc(scale*s, scale*t))) for s in path]
    lg = np.log(np.abs(vals)); ph = np.unwrap(np.angle(vals))
    # principal at start
    L = lg + 1j*ph
    return np.interp(sig_grid, path[::-1], L.real[::-1]) + 1j*np.interp(sig_grid, path[::-1], L.imag[::-1])

def logzeta_point(beta, scale=1, start=4.0, nstep=400):
    s0 = mp.mpf(start)
    sig = mp.re(beta); t = mp.im(beta)
    path = [s0 + (sig - s0)*j/nstep for j in range(nstep+1)]
    vals = [mp.zeta(scale*mp.mpc(s, t)) for s in path]
    ph = np.unwrap([float(mp.arg(v)) for v in vals])
    return mp.log(abs(vals[-1])) + 1j*mp.mpf(ph[-1])

def Zq_point(beta, q, B, K, NMAX):
    L = {1: logzeta_point(beta, 1), 2: logzeta_point(beta, 2)}
    for N in range(3, NMAX+1):
        L[N] = mp.log(mp.zeta(N*beta))
    tot = 0
    for k in range(K+1):
        e = sum(B[k, N]/N * L[N] for N in range(1, NMAX+1))
        tot += (q**k/qp(q, k)) * mp.e**e
    return qp(q, 200) * tot

def grid_Z(q, sig, ts, NMAX=40):
    """vectorised: continued log zeta along horizontal rows (fine path from Re 4), Hurwitz-EM zeta."""
    K = 0 if q == 0 else int(36/(-math.log(abs(q)))) + 2
    B = bk_table(q, K, NMAX)
    W = np.array([q**k/qp(q, k) for k in range(K+1)])
    path = np.arange(4.0, sig.min()-1e-9, -0.005)
    def contlog(scale):
        S = scale*(path[None, :] + 1j*ts[:, None])
        v = hurwitz(S, 1.0)
        L = np.log(np.abs(v)) + 1j*np.unwrap(np.angle(v), axis=1)
        out = np.zeros((len(ts), len(sig)), dtype=complex)
        for i in range(len(ts)):
            out[i] = np.interp(sig, path[::-1], L[i, ::-1].real) + 1j*np.interp(sig, path[::-1], L[i, ::-1].imag)
        return out
    Ls = [contlog(1), contlog(2)]
    G = sig[None, :] + 1j*ts[:, None]
    for N in range(3, NMAX+1):
        Ls.append(np.log(hurwitz(N*G, 1.0)))
    Lmat = np.array(Ls)                          # (NMAX, nt, ns)
    E = np.tensordot(B[:, 1:]/np.arange(1, NMAX+1), Lmat, axes=(1, 0))   # (K+1, nt, ns)
    Z = qp(q, 200)*(W[:, None, None]*np.exp(E)).sum(axis=0)
    return Z, B, K

def cell_zeros(Z, sig, ts):
    ph = np.angle(Z)
    out = []
    for i in range(len(ts)-1):
        for j in range(len(sig)-1):
            loop = [ph[i, j], ph[i, j+1], ph[i+1, j+1], ph[i+1, j], ph[i, j]]
            w = sum(((loop[a+1]-loop[a]+np.pi) % (2*np.pi)) - np.pi for a in range(4))
            if abs(w) > np.pi:
                out.append((sig[j]+sig[j+1])/2 + 1j*(ts[i]+ts[i+1])/2)
    return out

if __name__ == '__main__':
    sig = np.round(np.arange(0.505, 1.605, 0.01), 4)
    ts = np.round(np.arange(1.0, 50.0, 0.05), 4)
    res = {'sig': sig.tolist(), 'ts': ts.tolist()}
    NMAX = 40
    for q in [0.2, 0.5, 0.8, -0.5, -0.8]:
        Z, B, K = grid_Z(q, sig, ts, NMAX)
        cz = cell_zeros(Z, sig, ts)
        refined = []
        for z0 in cz:
            try:
                r = mp.findroot(lambda b: Zq_point(b, q, B, K, NMAX), mp.mpc(z0.real, z0.imag), tol=1e-18)
                refined.append([float(mp.re(r)), float(mp.im(r))])
            except Exception as e:
                refined.append([z0.real, z0.imag, 'unrefined'])
        print(f"q = {q}: {len(cz)} Fisher zeros in 0.51<Re<1.61, 1<Im<50")
        for r in refined: print("   ", r)
        res[str(q)] = {'zeros': refined, 'absZ': np.log10(np.abs(Z)+1e-300).tolist()}
    rz = [float(mp.im(mp.zetazero(n))) for n in range(1, 11)]
    res['riemann'] = rz
    json.dump(res, open('fisher.json', 'w'))
    print("saved")
