"""Round 4: follow the local fused discriminant D_e through the rewind.
(1) toy polynomial, exact heat flow, through the merger of the close pair;
(2) Lehmer's pair window (round-3 ODE), D_e from the zeros at each step.
Conventions: z units, dH/dt = -H_zz, rewind = t decreasing. -> rewind4.json"""
import json, numpy as np, mpmath as mp
from numpy.polynomial import polynomial as Pn
out = {}
# ---------- (1) toy ----------
mp.mp.dps = 50
roots = [1, 2, 2.3, 3.5, 5]
P = [mp.mpf(1)]
for r in roots:
    P = Pn.polymul(P, [-mp.mpf(r)**2, 0, 1])
P = [mp.mpf(c) for c in P]
def deriv(c, k=1):
    for _ in range(k): c = [i*c[i] for i in range(1, len(c))]
    return c
def heat(c, t):  # e^{-t D^2} c, exact for polynomials
    res = [mp.mpf(0)]*len(c); dk = c; k = 0; fac = mp.mpf(1)
    while any(dk):
        term = [(-t)**k/fac*v for v in dk]
        for i, v in enumerate(term): res[i] += v
        dk = deriv(dk, 2); k += 1; fac *= k
    return res
ev = lambda c, x: mp.polyval(c[::-1], x)
def De_toy(t):
    c = heat(P, t); c1 = deriv(c); c2 = deriv(c, 2); c3 = deriv(c, 3); c4 = deriv(c, 4)
    xe = mp.findroot(lambda x: ev(c1, x), (mp.mpf(1.9), mp.mpf(2.45)), solver='anderson')
    H0, H2, H3, H4 = ev(c, xe), ev(c2, xe), ev(c3, xe), ev(c4, xe)
    De = -8*H0/H2; kap = H4/H2 - (H3/H2)**2
    rts = sorted([r for r in mp.polyroots(c[::-1], maxsteps=200, extraprec=200)], key=lambda r: abs(mp.re(r)-2.15))[:2]
    return float(xe), float(De), float(8+De*kap), [complex(r) for r in rts]
ts = np.linspace(0.004, -0.030, 69)
toy = []
for t in ts:
    xe, De, rate, rts = De_toy(mp.mpf(t))
    toy.append(dict(t=float(t), xe=xe, De=De, rate=rate, im=max(abs(r.imag) for r in rts), re=[r.real for r in rts]))
# numerical derivative check of the rate identity at t=0
h = mp.mpf('1e-12')
num_rate = (De_toy(h)[1]-De_toy(-h)[1])/(2*float(h))
# merger time: where the two roots meet (D_e = 0)
tm = mp.findroot(lambda t: De_toy(t)[1], (mp.mpf(-0.005), mp.mpf(-0.02)), solver='anderson')
out['toy'] = dict(roots=roots, traj=toy, De0=De_toy(0)[1], rate0=De_toy(0)[2], num_rate0=num_rate,
                  t_merge=float(tm), lone_pair=-(0.3**2)/8, De0_over_8=-De_toy(0)[1]/8)
print('toy', out['toy']['De0'], out['toy']['rate0'], num_rate, float(tm), out['toy']['lone_pair'])
# ---------- (2) Lehmer window ----------
G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
L = json.load(open('../round3/zeros_lehmer.json')); GL = np.array([float(x) for x in L['gamma']]); NL = L['n']
G2 = np.concatenate([G, GL]); z_all = 2*G2
dgl = np.diff(GL); il = int(np.argmin(dgl[152:-152])) + 152; j = len(G) + il
W = 150; lo, hi = j-W, j+W+2
zw = z_all[lo:hi].copy(); outz = np.concatenate([z_all[:lo], z_all[hi:]])
zmax = 2*GL[-1]
wq = np.linspace(zmax, zmax*400, 400001); rho = np.log(wq/(4*np.pi))/(4*np.pi)
def bg(x, k):
    """k=1: background part of (log H)'(x); k=2: of (log H)''(x) (frozen zeros, mirrors, tail)"""
    if k == 1:
        s = np.sum(1/(x-outz)) + np.sum(1/(x+z_all)) + np.trapezoid((1/(x-wq)+1/(x+wq))*rho, wq)
    else:
        s = -np.sum(1/(x-outz)**2) - np.sum(1/(x+z_all)**2) - np.trapezoid((1/(x-wq)**2+1/(x+wq)**2)*rho, wq)
    return s
bgv = np.array([2*bg(zi, 1) for zi in zw])
def vel(z):
    D = z[:, None]-z[None, :]; np.fill_diagonal(D, np.inf)
    return 2*np.sum(1/D, axis=1) + bgv
def De_win(z):
    a, b = z[W], z[W+1]
    xm = 0.5*(a+b); B1 = bg(xm, 1); B2 = bg(xm, 2)   # background is smooth across the pair
    f1 = lambda x: np.sum(1/(x-z)) + B1
    lo_, hi_ = a+1e-12*(b-a), b-1e-12*(b-a)
    for _ in range(200):
        m = 0.5*(lo_+hi_)
        if np.sign(f1(m)) == np.sign(f1(lo_)): lo_ = m
        else: hi_ = m
    x = 0.5*(lo_+hi_)
    l2 = -np.sum(1/(x-z)**2) + B2
    return x, -8/l2
t = 0.0; traj = []; g0 = zw[W+1]-zw[W]
x0, D0 = De_win(zw); traj.append((t, D0, zw[W+1]-zw[W]))
while True:
    g = zw[W+1]-zw[W]
    if g < 2e-3*g0: break
    dt = -0.002*g*g
    k1 = vel(zw); k2 = vel(zw+dt/2*k1); k3 = vel(zw+dt/2*k2); k4 = vel(zw+dt*k3)
    zw = zw + dt/6*(k1+2*k2+2*k3+k4); t += dt
    x, De = De_win(zw); traj.append((t, De, g))
traj = np.array(traj)
# extrapolate D_e linearly to 0 from last points
p = np.polyfit(traj[-20:, 0], traj[-20:, 1], 1); t_zero = -p[1]/p[0]
rate0 = (traj[1, 1]-traj[0, 1])/(traj[1, 0]-traj[0, 0])
out['lehmer'] = dict(n=NL[il], gamma=float(GL[il]), dgamma=float(dgl[il]), De0_z=float(D0), De0_gamma=float(D0/4),
                     rate0_z=float(rate0), t_zero_De=float(t_zero), traj=traj[::max(1, len(traj)//300)].tolist())
print('lehmer', out['lehmer']['De0_gamma'], dgl[il]**2, rate0, t_zero)
json.dump(out, open('rewind4.json', 'w'))
