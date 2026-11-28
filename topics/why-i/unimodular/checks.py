"""Checks for 'Unimodular gravity and why quantum needs i' (DRAFT).
Model: flat FLRW, unimodular gauge N a^3 = 1 (Unruh 1989; Fabris et al. 2025, Eqs. 18-19), z = a^3.
Run: python3 checks.py -> PASS/FAIL lines, results.json."""
import json, os, numpy as np, sympy as sp
from scipy.linalg import expm, eigh
HERE = os.path.dirname(os.path.abspath(__file__)); res, out = {}, []
def check(n, ok, msg): out.append(f"[{'PASS' if ok else 'FAIL'}] {n}: {msg}"); res[n] = {'pass': bool(ok), 'msg': msg}

# ---------- classical ----------
# M1 unimodular time is the four-volume to the past, per unit comoving volume: dT = a^3 dt
t, H, T0 = sp.symbols('t H T0', positive=True)
a = sp.exp(H*t); T = sp.integrate(a**3, (t, 0, t))
dt_from_N = sp.simplify(sp.diff(T, t) - a**3)
check('M1 unimodular time = four-volume', dt_from_N == 0 and sp.simplify(T - (sp.exp(3*H*t) - 1)/(3*H)) == 0,
      'with N a^3 = 1, dT = a^3 dt; for de Sitter T = (e^{3Ht} - 1)/(3H): each cell keeps its four-volume while space stretches')

# M2 the conserved energy of unimodular time is (minus) the cosmological constant
A, P, E = sp.symbols('a pi_a E', real=True); rho = sp.Function('rho')
Hu = -P**2/(24*A**4) + rho(A)                      # Fabris et al. Eq. 18 without the scalar field; matter density rho(a)
Psq = sp.solve(sp.Eq(Hu, E), P**2)[0]              # on a solution of energy E
dadT = sp.diff(Hu, P)                              # Hamilton: da/dT
hub2 = sp.simplify(((A**3*dadT)/A)**2).subs(P**2, Psq)   # (da/dt / a)^2 with dt = dT / a^3
check('M2 Lambda is the unimodular energy', sp.simplify(hub2 - (rho(A) - E)/6) == 0,
      '(adot/a)^2 = (rho - E)/6: Friedmann with 8piG/3 -> 1/6 and the cosmological term set by -E, a constant of the motion')

# ---------- quantum: discretize H_u = (3/8)(1/z) d/dz (z d/dz) + rho_dust(z) + c, inner product z dz ----------
# with Psi = u / sqrt(z) this is (3/8)(u'' + u/(4 z^2)) + ..., a symmetric matrix on u
Ngrid, L = 600, 30.0; h = L/Ngrid; z = h*np.arange(1, Ngrid + 1)
lap = (np.diag(-2*np.ones(Ngrid)) + np.diag(np.ones(Ngrid - 1), 1) + np.diag(np.ones(Ngrid - 1), -1))/h**2
def Hu_mat(c, m=1.0, robin=0.0):
    K = 3/8*(lap + np.diag(1/(4*z**2)))
    K = K.copy(); K[0, 0] += 3/8*robin/h**2           # a different self-adjoint extension at the big bang (z -> 0)
    return K + np.diag(m/z) + c*np.eye(Ngrid)
u0 = np.exp(-(z - 8)**2/2)*np.exp(-1j*0.6*z); u0 /= np.linalg.norm(u0)
def zexp(Hm, TT, u0=u0): u = expm(-1j*Hm*TT) @ u0; return float(np.real(np.vdot(u, z*u)))

# M3 vacuum energy is a common phase in unimodular time: expansion history unchanged
Ts = [0.0, 2.0, 5.0]
d = max(abs(zexp(Hu_mat(cv), TT) - zexp(Hu_mat(0.0), TT)) for cv in (1.0, 5.0, 40.0) for TT in Ts)
moves = [zexp(Hu_mat(0.0), TT) for TT in Ts]
check('M3 vacuum energy changes nothing in unimodular time', d < 1e-9 and np.ptp(moves) > 0.1,
      f'<a^3>(T) = {[round(x,3) for x in moves]} at T = 0, 2, 5 (the universe evolves); change for vacuum energy c = 1, 5, 40: {d:.1e}')

# M4 the spectrum shifts rigidly and the eigenstates don't move: exactly why-i's assumption A4
w0, v0 = eigh(Hu_mat(0.0)); w5, v5 = eigh(Hu_mat(5.0))
shift = np.max(np.abs(w5 - w0 - 5.0)); ovl = np.min(np.abs(np.sum(v0*v5, axis=0)))
check('M4 A4 holds: A(H + c) = A(H) + c A(1)', shift < 1e-9 and ovl > 1 - 1e-9,
      f'every eigenvalue moves by exactly c (error {shift:.1e}); eigenstates unchanged (min overlap {ovl:.12f})')

# M5 A4 holds for every choice of big-bang boundary condition (self-adjoint extension), which is a real free parameter
un = np.exp(-(z - 1.0)**2/0.5); un = un/np.linalg.norm(un)      # a packet near the big bang feels the boundary condition
dr = max(abs(zexp(Hu_mat(cv, robin=0.7), 3.0, un) - zexp(Hu_mat(0.0, robin=0.7), 3.0, un)) for cv in (1.0, 5.0))
dext = abs(zexp(Hu_mat(0.0, robin=0.7), 3.0, un) - zexp(Hu_mat(0.0), 3.0, un))
check('M5 boundary condition at the big bang is free, A4 holds for each', dr < 1e-9 and dext > 1e-3,
      f'changing the extension changes <a^3>(3) by {dext:.2e}; vacuum energy then changes it by {dr:.1e}')

# M6 contrast: proper-time (Wheeler-DeWitt) form. The constraint is z * H_u, so vacuum energy enters as c * z,
# whose generator z does not commute with the kinetic term: the physical (zero-mode) state changes shape
Zm = np.diag(z)
def wdw_state(c):
    C = Zm @ Hu_mat(0.0) ; C = (C + C.T)/2 + c*Zm        # symmetrized z H_u plus c z
    w, v = eigh(C); i = np.argmin(np.abs(w)); return v[:, i], w[i]
s0, e0 = wdw_state(0.0); s1, e1 = wdw_state(0.05)
comm = np.linalg.norm(Zm @ Hu_mat(0.0) - Hu_mat(0.0) @ Zm)
ov = abs(np.dot(s0, s1))
check('M6 without fixed volume the vacuum energy is physical', comm > 1 and ov < 0.99 and np.isrealobj(s0),
      f'||[z, H_u]|| = {comm:.1f}; zero-mode overlap for c = 0 vs 0.05 is {ov:.3f}; the zero mode is a real vector (no time, no i needed)')

# M7 the unimodular dynamics is not real: e^{-i H T} has no real form unless H = 0 (so a shared J is needed)
U = expm(-1j*Hu_mat(0.0)*1.0)
check('M7 unimodular evolution needs i', np.linalg.norm(U.imag) > 1, f'||Im e^(-i H_u T)|| = {np.linalg.norm(U.imag):.1f} at T = 1; real form uses J with J^2 = -1')

# M8 Sorkin's estimate: Lambda ~ 1/sqrt(V), with V the four-volume of the observable universe in Planck units
H0 = 67.4e3/3.0857e22; tP = 5.391e-44; OmL = 0.685
R = 1/(H0*tP); V = R**4; lam_sorkin = V**-0.5; lam_obs = 3*OmL*(H0*tP)**2
ratio = lam_obs/lam_sorkin
check('M8 Sorkin: Lambda ~ 1/sqrt(V) has the observed order', 0.1 < ratio < 10,
      f'1/sqrt(V) = {lam_sorkin:.2e}, observed Lambda = {lam_obs:.2e} (Planck units), ratio {ratio:.2f}; naive vacuum energy ~ 1')

res['numbers'] = {'zexp_T': moves, 'lam_sorkin': lam_sorkin, 'lam_obs': lam_obs}
n = sum(v['pass'] for v in res.values() if 'pass' in v); tot = sum(1 for v in res.values() if 'pass' in v)
out.append(f'\n{n}/{tot} checks pass'); print('\n'.join(out))
json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
