"""Round 8, item 3a: a twist on a flat surface.  Square torus T vs Klein bottle K.

Torus R^2/Z^2: Laplace eigenvalues 4 pi^2 (m^2+n^2); spectral zeta (2pi)^{-2s} * 4 zeta(s) L(s, chi_-4).
Klein bottle = T / glide g(x,y) = (x+1/2, -y)  (double cover lattice Z x Z, glide squares to x+1).
On e_{m,n} = exp 2pi i (m x + n y):  g* e_{m,n} = (-1)^m e_{m,-n}.
Invariant (even) sector  = Klein bottle:  pairs {(m,n),(m,-n)} n>0, all m;  n=0 needs m even.
Anti-invariant (odd) sector: pairs n>0, all m;  n=0 needs m odd.
Dropping (2pi)^{-2s}:  Z_even = 2 zeta L - (1 - 2^{1-2s}) zeta(2s),  Z_odd = 2 zeta L + (1 - 2^{1-2s}) zeta(2s).
Z_even + Z_odd = torus.  Twist strength eps: Z_eps = 2 zeta L - eps (1-2^{1-2s}) zeta(2s).
"""
import json, numpy as np, mpmath as mp
from fastl import hurwitz, dirichlet_series_periodic

mp.mp.dps = 20
def Zeps_np(s, eps):
    zl = hurwitz(s, 1.0)*dirichlet_series_periodic(s, [1, 0, -1, 0])
    return 2*zl - eps*(1-np.exp((1-2*s)*np.log(2)))*hurwitz(2*s, 1.0)

def L4(s): return (mp.zeta(s, 0.25) - mp.zeta(s, 0.75))*mp.power(4, -s)
def Zeps(s, eps): return 2*mp.zeta(s)*L4(s) - eps*(1-mp.power(2, 1-2*s))*mp.zeta(2*s)

# direct check of the Klein bottle spectrum formula at s = 2
def direct(s, sector, R=400):
    tot = 0.0
    for m in range(-R, R+1):
        for n in range(0, R+1):
            if m == 0 and n == 0: continue
            lam = m*m + n*n
            if n > 0: tot += lam**(-s)
            elif (m % 2 == 0) == (sector == 'even'): tot += lam**(-s)
    return tot
for sec, eps in (('even', 1), ('odd', -1)):
    d = direct(2.0, sec); f = float(mp.re(Zeps(2, eps)))
    print(f"{sec} sector at s=2: direct lattice sum {d:.6f}, formula {f:.6f}")

def cell_zeros(Z, sig, ts):
    ph = np.angle(Z); out = []
    for i in range(len(ts)-1):
        for j in range(len(sig)-1):
            loop = [ph[i, j], ph[i, j+1], ph[i+1, j+1], ph[i+1, j], ph[i, j]]
            w = sum(((loop[a+1]-loop[a]+np.pi) % (2*np.pi)) - np.pi for a in range(4))
            if abs(w) > np.pi: out.append(complex((sig[j]+sig[j+1])/2, (ts[i]+ts[i+1])/2))
    return out

sig = np.arange(-0.495, 1.6, 0.01); ts = np.arange(2.0, 80.0, 0.02)
S = sig[None, :] + 1j*ts[:, None]
res = {}
for name, eps in (('torus', 0.0), ('klein_even', 1.0), ('klein_odd', -1.0)):
    Z = Zeps_np(S, eps)
    cz = cell_zeros(Z, sig, ts)
    ref = []
    for z0 in cz:
        try:
            r = mp.findroot(lambda s: Zeps(s, eps), mp.mpc(z0.real, z0.imag))
            ref.append([float(mp.re(r)), float(mp.im(r))])
        except Exception:
            ref.append([z0.real, z0.imag])
    ref = sorted({(round(a, 9), round(b, 9)) for a, b in ref}, key=lambda x: x[1])
    on = sum(1 for a, b in ref if abs(a-0.5) < 1e-7)
    print(f"{name}: {len(ref)} zeros with -0.5<Re<1.6, 2<Im<80; on Re=1/2: {on}; Re range {min(a for a,b in ref):.4f}..{max(a for a,b in ref):.4f}")
    res[name] = ref

# zero trajectories as the twist is turned on, eps 0 -> 1 and 0 -> -1, from torus zeros below Im 40
tor = [z for z in res['torus'] if z[1] < 40]
traj = {}
for sgn in (1, -1):
    paths = []
    for a, b in tor:
        s0 = mp.mpc(a, b); path = [[a, b]]
        for e in np.linspace(0, sgn, 41)[1:]:
            try:
                s0 = mp.findroot(lambda s: Zeps(s, e), s0)
            except Exception:
                break
            path.append([float(mp.re(s0)), float(mp.im(s0))])
        paths.append(path)
    traj[str(sgn)] = paths
res['traj'] = traj
json.dump(res, open('klein.json', 'w'))
print('saved')
