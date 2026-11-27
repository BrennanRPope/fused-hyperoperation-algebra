"""DRAFT, round 19 (Nov 27, 2026): do Kerr's orbits over the poles follow the shifted source's exact polar law?

Round 18 found that Newton's orbit over the poles of the shifted source has an exact law: at fixed distance r,
Omega_theta^2 r^3/M = (pi/2)^2/(2K(k) - E(k))^2 with k = a/r. Does Kerr's polar orbit (no spin about the axis) follow
it, in the same Boyer-Lindquist radius r?

1. Kerr's spherical polar orbits have exact constants (from Carter's R(r) = R'(r) = 0 with L = 0):
     E^2 = r Delta^2 / [(r^2+a^2)(r^3 + a^2 r - 3M r^2 + M a^2)],
     Q   = M r^2 [(r^2+a^2)^2 - 4M a^2 r] / [(r^2+a^2)(r^3 + a^2 r - 3M r^2 + M a^2)],
   and the up-and-down rate in elliptic integrals with modulus
     m = (a^2/r^2) [(r^2+a^2)^2 - 4M r^3] / [(r^2+a^2)^2 - 4M a^2 r]   (m = a^2/r^2 when M -> 0).
2. Answer, part one: at Newtonian strength Kerr follows the shifted source exactly, at every order in a/r: the u^0 row
   of Omega_theta^2 r^3/M is 1 - 3k^2/2 + 33k^4/32 - 41k^6/64 + ..., and the same holds for the ratio (Omega_r/Omega_theta)^2,
   whose Newtonian part is (1 + k^2)(2K(k)/pi)^2.
3. Answer, part two: beyond Newton, over the poles mass and shape mix (-4 u k^2 in Omega_theta^2 r^3/M, +15 u k^2 in the
   ratio), while in the equator the ratio (Omega_r/Omega_phi)^2 = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 has no mixing at all.
4. An independent check: polar orbits flown in Kerr-Schild Cartesian coordinates (no Carter separation used) keep r
   fixed and beat at the exact rate.
Run: python3 polar.py  (writes polar.json, about 1 minute)
"""
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
r, M, a = sp.symbols("r M a", positive=True)
u, k, mm = sp.symbols("u k m", positive=True)


# ---------------------------------------------------------------- 1. exact constants and rates
def constants():
    E2, Q = sp.symbols("E2 Q")
    D = r**2 - 2*M*r + a**2
    R = E2*((r**2 + a**2)**2 - a**2*D) - D*(r**2 + Q)            # Carter's R(r) with L = 0 (K = Q + a^2 E^2)
    sol = sp.solve([R, sp.diff(R, r)], [E2, Q], dict=True)[0]
    E2s, Qs = sp.factor(sol[E2]), sp.factor(sol[Q])
    m = sp.factor(a**2*(1 - E2s)/Qs)
    x = sp.symbols("x", positive=True)
    Rx = E2s*(x**2 + a**2)**2 - (x**2 - 2*M*x + a**2)*(x**2 + a**2*E2s + Qs)
    Ur2 = sp.factor(-sp.diff(Rx, x, 2).subs(x, r)/2)              # radial Mino rate^2 of a nearby orbit
    return dict(E2=E2s, Q=Qs, m=m, Ups_r2=Ur2, pref=sp.factor(Ur2/Qs), D=D)


def omega_theta_exact(Mv, av, rv):
    """Up-and-down rate in Boyer-Lindquist (= Kerr-Schild) coordinate time, exact."""
    Mv, av, rv = mp.mpf(Mv), mp.mpf(av), mp.mpf(rv)
    D = rv**2 - 2*Mv*rv + av**2
    den = (av**2 + rv**2)*(Mv*av**2 - 3*Mv*rv**2 + av**2*rv + rv**3)
    E = mp.sqrt(rv*D**2/den); Q = Mv*rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv)/den
    m = av**2*((rv**2 + av**2)**2 - 4*Mv*rv**3)/(rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv))
    K, Ee = mp.ellipk(m), mp.ellipe(m)
    T = 4*E/mp.sqrt(Q)*(((rv**2 + av**2)**2/D)*K - av**2*(K - (K - Ee)/m))
    return 2*mp.pi/T


def newton_polar(kv):
    kv = mp.mpf(kv)
    return (mp.pi/2)**2/(2*mp.ellipk(kv**2) - mp.ellipe(kv**2))**2


def double_series(N=10, U=4):
    """Omega_theta^2 r^3/M and (Omega_r/Omega_theta)^2 as double series in k = a/r and u = M/r (r = 1)."""
    c = constants()
    sub = {r: 1, M: u, a: k}
    E2, Q, m, D = (c[x].subs(sub) for x in ("E2", "Q", "m", "D"))
    Ks = sp.series(sp.elliptic_k(mm), mm, 0, N//2 + 1).removeO()
    Es = sp.series(sp.elliptic_e(mm), mm, 0, N//2 + 1).removeO()
    KmE = sp.expand(sp.cancel((Ks - Es)/mm))
    T = 4*sp.sqrt(E2)/sp.sqrt(Q)*((1 + k**2)**2/D*Ks.subs(mm, m) - k**2*(Ks.subs(mm, m) - KmE.subs(mm, m)))
    F = sp.series(sp.simplify((2*sp.pi/T)**2/u), k, 0, N).removeO()
    ratio = sp.series(sp.simplify(c["pref"].subs(sub)*(2*Ks.subs(mm, m)/sp.pi)**2), k, 0, N).removeO()
    rows = {}
    for name, ser in (("omega_theta2", F), ("ratio", ratio)):
        rows[name] = [[sp.nsimplify(sp.series(sp.simplify(ser.coeff(k, j)), u, 0, U).removeO().coeff(u, i))
                       for i in range(U)] for j in range(0, N, 2)]
    m_ = sp.symbols("m", positive=True)
    rows["newton_omega_theta2"] = [sp.nsimplify(c_) for c_ in sp.Poly(sp.series(
        (sp.pi/2)**2/(2*sp.elliptic_k(k**2) - sp.elliptic_e(k**2))**2, k, 0, N).removeO(), k).all_coeffs()[::-1][::2]]
    rows["newton_ratio"] = [sp.nsimplify(c_) for c_ in sp.Poly(sp.series(
        (1 + k**2)*(2*sp.elliptic_k(k**2)/sp.pi)**2, k, 0, N).removeO(), k).all_coeffs()[::-1][::2]]
    return rows


# ---------------------------------------------------------------- 2. independent flight in Kerr-Schild coordinates
def ks_hamiltonian():
    X, Y, Z, Ms, As = sp.symbols("x y z M a", real=True)
    rho2 = X**2 + Y**2 + Z**2
    rr = sp.sqrt((rho2 - As**2 + sp.sqrt((rho2 - As**2)**2 + 4*As**2*Z**2))/2)
    f = 2*Ms*rr**3/(rr**4 + As**2*Z**2)
    lu = [-1, (rr*X + As*Y)/(rr**2 + As**2), (rr*Y - As*X)/(rr**2 + As**2), Z/rr]   # l^mu, ingoing
    P = sp.symbols("pt px py pz", real=True)
    eta = [-1, 1, 1, 1]
    H = sp.Rational(1, 2)*(sum(eta[i]*P[i]**2 for i in range(4)) - f*sum(lu[i]*P[i] for i in range(4))**2)
    grad = [sp.diff(H, p) for p in P] + [sp.diff(H, c) for c in (X, Y, Z)]
    return (sp.lambdify([X, Y, Z, *P, Ms, As], grad, "numpy"), sp.lambdify([X, Y, Z, As], rr, "numpy"),
            sp.lambdify([X, Y, Z, *P, Ms, As], H, "numpy"))


def fly(Mv, av, rv, laps=8):
    """Fly the polar spherical orbit from the equator with Hamilton's equations of the exact Kerr-Schild metric."""
    fn, rfn, Hfn = ks_hamiltonian()
    D = rv**2 - 2*Mv*rv + av**2
    den = (av**2 + rv**2)*(Mv*av**2 - 3*Mv*rv**2 + av**2*rv + rv**3)
    E = np.sqrt(rv*D**2/den); Q = Mv*rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv)/den
    pBL = np.array([-E, 0.0, -np.sqrt(Q), 0.0])                   # (p_t, p_r, p_theta, p_phi) at theta = pi/2
    J = np.zeros((4, 4))                                          # d(t_KS, x, y, z)/d(t, r, theta, phi) there
    J[0, 0] = 1; J[0, 1] = 2*Mv*rv/D
    J[1, 1] = 1 - av*av/D; J[1, 3] = -av
    J[2, 1] = rv*av/D;     J[2, 3] = rv
    J[3, 2] = -rv
    pKS = np.linalg.solve(J.T, pBL)
    s0 = [0.0, rv, av, 0.0, *pKS]

    def flow(tau, s):
        o = fn(s[1], s[2], s[3], *s[4:], Mv, av)
        return [o[0], o[1], o[2], o[3], 0.0, -o[4], -o[5], -o[6]]

    def up(tau, s):
        return s[3]
    up.direction = 1
    Om = float(omega_theta_exact(Mv, av, rv))
    sol = solve_ivp(flow, [0, 1.6*laps*2*np.pi/Om], s0, rtol=1e-12, atol=1e-12, events=up, dense_output=True)
    te = sol.y_events[0][:, 0]
    per = (te[-1] - te[0])/(len(te) - 1)
    ts = np.linspace(0, sol.t[-1], 3000); Yv = sol.sol(ts)
    return dict(M=Mv, a=av, r=rv, flown=2*np.pi/per, exact=Om, laps=len(te) - 1,
                max_dr=float(np.max(np.abs(rfn(Yv[1], Yv[2], Yv[3], av) - rv))),
                H0=float(Hfn(*s0[1:], Mv, av)), path=[Yv[1].tolist()[::6], Yv[2].tolist()[::6], Yv[3].tolist()[::6]])


def main():
    out = {}
    c = constants()
    out["constants"] = {x: str(c[x]) for x in ("E2", "Q", "m", "pref")}
    out["newton_limit"] = {"m": str(sp.simplify(c["m"].subs(M, 0))), "pref": str(sp.simplify(c["pref"].subs(M, 0))),
                           "schwarzschild_pref": str(sp.simplify(c["pref"].subs(a, 0)))}
    rows = double_series()
    out["series"] = {kk: [[str(x) for x in row] if isinstance(row, list) else str(row) for row in v]
                     for kk, v in rows.items()}
    print("series", out["series"], flush=True)
    flights = [fly(1.0, 0.9, 6.0), fly(1.0, 0.5, 8.0), fly(0.05, 0.8, 3.0)]
    out["flights"] = flights
    print("flights", [(f["flown"], f["exact"], f["max_dr"]) for f in flights], flush=True)
    # curves: Omega_theta^2 r^3/M against k = a/r at several field strengths u, and the shifted source's law
    curves = {}
    ks = [0.02*i for i in range(1, 46)]
    curves["shifted source (Newton)"] = [float(newton_polar(x)) for x in ks]
    for uv in (0.01, 0.03, 0.1):
        curves[f"Kerr, u = {uv}"] = [float(omega_theta_exact(uv, x, 1.0)**2/uv) for x in ks]
    out["curves"] = {"k": ks, **curves}
    json.dump(out, open(os.path.join(HERE, "polar.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
