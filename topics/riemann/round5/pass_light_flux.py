"""Round 5, light and flux-tube pass.

(L1) Light rotates the rewind. The de Bruijn-Newman flow H_t(z) = int_0^inf e^(t u^2) Phi(u) cos(z u) du
     solves dH/dt = -H''. A zero moves with dz/dt = H''/H' (real: along the line). Wick-rotate
     t = i tau (heat -> Schroedinger, the light-cone / wave side): dz/dtau = i H''/H', the same speed
     turned 90 degrees, so every zero leaves the line at once. Checked against zeros of H_(i tau)
     computed straight from Phi.
(L2) Berry-Keating's H = xp: its orbits xp = E are boost orbits (light-cone coordinates), the phase-space
     area under the hyperbola with cutoffs x >= lx, p >= lp is E*eta - E + lx*lp with rapidity span
     eta = log(E/(lx lp)). With lx lp = 2 pi: N(E) = (E/2pi)(log(E/2pi) - 1) + 1 (+ Maslov -1/8), the smooth
     zero count. Compared with the 5,000 zeros.
(F)  Flux-tube reading: the effective dimension of the zero count D_eff = d log N / d log T
     = log(T/2pi)/(log(T/2pi) - 1) falls to 1 (a tube) from above; measured from the zeros.
"""
import json, math, numpy as np, mpmath as mp, sympy as sp

res = {}
G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])

# ---------- (L1) zero velocities from the zeros (gamma units, then z = 2 gamma) ----------
def vel_gamma(n):
    g = G[n]
    others = np.delete(G, n)
    s = np.sum(2/(g - others)) + np.sum(2/(g + G))
    # tail beyond the last zero: density log(x/2pi)/2pi, pairs give 4g/(g^2 - x^2)
    w = np.linspace(G[-1], 400*G[-1], 400001); rho = np.log(w/(2*np.pi))/(2*np.pi)
    s += np.trapezoid(4*g/(g*g - w*w)*rho, w)
    return s          # Xi''/Xi' at gamma_n (gamma variable)
v = np.array([vel_gamma(n) for n in range(40)])
res['v_gamma_first40'] = v.tolist()          # heat: d gamma/dt' ; light: i times it

# ---------- direct check from Phi ----------
mp.mp.dps = 30
def Phi(u):
    tot = mp.mpf(0)
    for n in range(1, 12):
        e = mp.exp(4*u); a = mp.pi*n*n
        tot += (2*a*a*mp.exp(9*u) - 3*a*mp.exp(5*u))*mp.exp(-a*e)
    return tot
def H(z, t):   # t may be complex (Wick-rotated)
    return mp.quad(lambda u: mp.exp(t*u*u)*Phi(u)*mp.cos(z*u), [0, 0.25, 0.5, 1, 2])
def dH(z, t):
    return mp.quad(lambda u: -u*mp.exp(t*u*u)*Phi(u)*mp.sin(z*u), [0, 0.25, 0.5, 1, 2])
def d2H(z, t):
    return mp.quad(lambda u: -u*u*mp.exp(t*u*u)*Phi(u)*mp.cos(z*u), [0, 0.25, 0.5, 1, 2])
z1 = mp.findroot(lambda z: H(z, 0), 2*G[0])
res['first_zero_z'] = float(z1); res['first_zero_gamma_check'] = float(z1/2 - G[0])
vz = d2H(z1, 0)/dH(z1, 0)
res['v_z_from_Phi'] = float(mp.re(vz))
res['v_z_from_zeros'] = float(v[0]/2)          # H''/H' in z = (1/2) Xi''/Xi' in gamma
wick = []
for tau in [0.002, 0.005, 0.01, 0.02]:
    zt = mp.findroot(lambda z: H(z, mp.mpc(0, tau)), z1)
    pred = z1 + mp.mpc(0, 1)*tau*vz
    wick.append(dict(tau=tau, zero=[float(mp.re(zt)), float(mp.im(zt))],
                     predicted=[float(mp.re(pred)), float(mp.im(pred))]))
heat = []
for t in [0.002, 0.005, 0.01, 0.02]:
    zt = mp.findroot(lambda z: H(z, t), z1)
    heat.append(dict(t=t, zero=float(mp.re(zt)), predicted=float(z1 + t*vz), imag=float(mp.im(zt))))
res['wick_check'] = wick; res['heat_check'] = heat

# ---------- (L2) Berry-Keating area under the boost hyperbola ----------
E, lx, lp, x = sp.symbols('E l_x l_p x', positive=True)
area = sp.integrate(E/x - lp, (x, lx, E/lp))
res['bk_area'] = str(sp.simplify(area))
res['bk_area_is_E_eta'] = bool(sp.simplify(area - (E*sp.log(E/(lx*lp)) - E + lx*lp)) == 0)
T = np.linspace(20, G[-1] - 1, 4000)
Ncount = np.searchsorted(G, T)
Nbk = T/(2*np.pi)*(np.log(T/(2*np.pi)) - 1) + 7/8
res['bk_residual_mean'] = float(np.mean(Ncount - Nbk))
res['bk_residual_maxabs'] = float(np.max(np.abs(Ncount - Nbk)))
res['bk_curve'] = dict(T=T[::20].tolist(), N=Ncount[::20].tolist(), Nbk=Nbk[::20].tolist())

# ---------- (F) effective dimension of the count ----------
Tm = np.exp(np.linspace(np.log(60), np.log(G[-1] - 1), 60))
Deff_meas = []
for a, b in zip(Tm[:-1], Tm[1:]):
    Na, Nb = np.searchsorted(G, a), np.searchsorted(G, b)
    Deff_meas.append((math.log(Nb) - math.log(Na))/(math.log(b) - math.log(a)))
Tc = np.sqrt(Tm[:-1]*Tm[1:])
Deff_smooth = np.log(Tc/(2*np.pi))/(np.log(Tc/(2*np.pi)) - 1)
res['Deff'] = dict(T=Tc.tolist(), measured=Deff_meas, smooth=Deff_smooth.tolist())
res['Deff_max_dev'] = float(np.max(np.abs(np.array(Deff_meas) - Deff_smooth)))
json.dump(res, open('pass_light_flux.json', 'w'), indent=1)
print({k: v for k, v in res.items() if k not in ('bk_curve', 'Deff', 'v_gamma_first40')})
