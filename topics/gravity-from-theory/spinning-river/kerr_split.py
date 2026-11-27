"""DRAFT, round 16 (Nov 27, 2026): Kerr's orbit factor split into its pieces, two ways.

A slightly eccentric orbit in the equator of a spinning black hole turns by 2 pi (1/sqrt(T) - 1) per lap, with
    T = (Omega_r / Omega_phi)^2 = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2,   u = M/r, chi = a/M (prograde chi > 0)
in Boyer-Lindquist r. This is exact (the known Kerr epicyclic frequencies). Which piece of the field makes each number?

Way 1, the double copy. Kerr's Kerr-Schild form is exact: g = flat + (2M/r) k k, so a particle's mass shell is
    F = -E^2 + ((r^2+a^2)/r^2) p_r^2 + L^2/(r^2+a^2) + 1 - (2M/r) (k.p)^2,   k.p = E + P + S,
with P = p_r and S = -a L/(r^2+a^2) (the twist of the rays). (k.p)^2 has six pieces, EE, PP, SS, 2EP, 2ES, 2PS. Each
gets a dial; T is computed exactly (to 60 digits) and its u, chi u^(3/2), u^2 and chi^2 u^2 numbers read off.
The single copy (electromagnetism's version) is -(2M/r)(E + P + S), linear.

Way 2, the rules. A test body around a spinning body in slow-motion form with dials,
    g_tt = -1 + 2U - 2 beta U^2,  g_tphi = -2 kappa M a / r,  g_rr = 1 + 2 gamma U,  g_phph = r^2 (1 + 2 gamma U),
    U = M/r + q M a^2 / (2 r^3)  (q = 1: the body's shape, Kerr's quadrupole Q = -M a^2),
and rule A (fields sent at light speed, so a moving source's field is its rest field boosted) fixes kappa.

Also: does the orbit page's product form 1 - 6u = (1 - 2u)(1 - (2v)^2) carry over? (No.)
Run: python3 kerr_split.py  (writes kerr_split.json, about 2 minutes)
"""
import itertools
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60
M, a, r, E, L, p = sp.symbols("M a r E L p", real=True)
W = dict(zip(("EE", "PP", "SS", "EP", "ES", "PS"), sp.symbols("wEE wPP wSS wEP wES wPS")))
RHO2 = r**2 + a**2
P_, S_ = p, -a*L/RHO2


def ks_shell(w):
    """Kerr-Schild mass shell with a dial on each piece of (k.p)^2. w: dict name -> 0/1 (missing = 1)."""
    d = {W[k]: w.get(k, 1) for k in W}
    sq = (d[W["EE"]]*E**2 + d[W["PP"]]*P_**2 + d[W["SS"]]*S_**2 + 2*d[W["EP"]]*E*P_ + 2*d[W["ES"]]*E*S_
          + 2*d[W["PS"]]*P_*S_)
    return -E**2 + (RHO2/r**2)*p**2 + L**2/RHO2 + 1 - (2*M/r)*sq


def single_shell(wP=1, wS=1):
    return -E**2 + (RHO2/r**2)*p**2 + L**2/RHO2 + 1 - (2*M/r)*(E + wP*P_ + wS*S_)


def pn_shell(beta=1, kappa=1, gamma=1, q=1):
    U = M/r + q*M*a**2/(2*r**3)
    G = sp.Matrix([[-1 + 2*U - 2*beta*U**2, -2*kappa*M*a/r], [-2*kappa*M*a/r, r**2*(1 + 2*gamma*U)]])
    Gi = G.inv()
    return Gi[0, 0]*E**2 - 2*Gi[0, 1]*E*L + Gi[1, 1]*L**2 + p**2/(1 + 2*gamma*U) + 1


def orbit_factor(Fs):
    """T(u, chi) = (Omega_r/Omega_phi)^2 for the circular orbit at r = 1 with M = u, a = chi u, for a mass shell
    F = A p^2 + 2 B p + C = 0. Radial motion about the circle sees F_eff = C - B^2/A; T = 2 A F_eff'' / F_L^2."""
    Fs = sp.expand(Fs)
    A = Fs.coeff(p, 2); B = Fs.coeff(p, 1)/2; C = Fs.coeff(p, 0)
    Feff = C - B**2/A
    args = (r, M, a, E, L)
    f0, f1, f2 = (sp.lambdify(args, sp.diff(Feff, r, k), "mpmath") for k in (0, 1, 2))
    fA = sp.lambdify(args, A, "mpmath")
    FL = sp.lambdify(args, sp.diff(Fs, L).subs(p, -B/A), "mpmath")

    def T(u, chi):
        Mv = mp.mpf(u); av = mp.mpf(chi)*Mv; rv = mp.mpf(1)
        e, l = mp.findroot(lambda e, l: (f0(rv, Mv, av, e, l), f1(rv, Mv, av, e, l)), (1 - Mv/2, mp.sqrt(Mv)))
        return 2*fA(rv, Mv, av, e, l)*f2(rv, Mv, av, e, l)/FL(rv, Mv, av, e, l)**2
    return T


def series(T, K=8):
    """The numbers in T - 1 = c_u u + c_s chi u^(3/2) + (c_uu + c_cc chi^2) u^2 + ..., by an exact fit in s = sqrt(u)
    at small u for chi = 0 and +-1/2. Returns rationals (as strings) for u, chi u^1.5, u^2, chi^2 u^2."""
    def fit(chi):
        ss = [mp.mpf(10)**(-3 - mp.mpf(j)/4) for j in range(K)]
        rows = [[s**k for k in range(2, K + 2)] for s in ss]
        c = mp.lu_solve(mp.matrix(rows), mp.matrix([T(s**2, chi) - 1 for s in ss]))
        return [c[i] for i in range(K)]
    c0, cp, cm = fit(0), fit(mp.mpf("0.5")), fit(mp.mpf("-0.5"))
    nice = lambda x: sp.nsimplify(float(x), rational=True, tolerance=1e-9)
    return {"u": nice(c0[0]), "chi_u15": nice(cp[1] - cm[1]), "u2": nice(c0[2]),
            "chi2_u2": nice(((cp[2] + cm[2])/2 - c0[2])/mp.mpf("0.25"))}


def exact_T(u, chi):
    return 1 - 6*u + 8*chi*u**mp.mpf(1.5) - 3*chi**2*u**2


def isco(chi):
    """Kerr's innermost stable circular orbit (Bardeen, Press & Teukolsky 1972), M = 1; chi < 0 retrograde."""
    z1 = 1 + (1 - chi**2)**(mp.mpf(1)/3)*((1 + chi)**(mp.mpf(1)/3) + (1 - chi)**(mp.mpf(1)/3))
    z2 = mp.sqrt(3*chi**2 + z1**2)
    return 3 + z2 - mp.sign(chi)*mp.sqrt((3 - z1)*(3 + z1 + 2*z2)) if chi != 0 else mp.mpf(6)


def zamo_speed(rr, chi):
    """Speed of a circular equatorial orbit seen by the local frame that turns with the hole (ZAMO); M = 1.
    Bardeen, Press & Teukolsky (1972); chi < 0 is retrograde."""
    a_ = abs(chi); s = 1 if chi >= 0 else -1
    D = rr**2 - 2*rr + a_**2
    return s*(rr**2 - s*2*a_*mp.sqrt(rr) + a_**2)/(mp.sqrt(D)*(rr**mp.mpf(1.5) + s*a_))


def product_form_gap(rr, chi):
    """T - (Delta/r^2)(1 - 4 v^2): zero for Schwarzschild, not for Kerr."""
    u = 1/mp.mpf(rr); v = zamo_speed(mp.mpf(rr), chi)
    return exact_T(u, chi) - (rr**2 - 2*rr + chi**2)/rr**2*(1 - 4*v**2)


def simulate(w, u, chi, laps=12):
    """Periapsis advance per lap of a slightly eccentric orbit, integrating Hamilton's equations of the exact
    Kerr-Schild mass shell (with piece dials w); compared with 2 pi (1/sqrt(T) - 1)."""
    from scipy.integrate import solve_ivp
    Fw = ks_shell(w).subs({M: u, a: chi*u})
    dFp, dFr, dFL = (sp.lambdify((r, p, E, L), sp.diff(Fw, x)) for x in (p, r, L))
    Fe = sp.expand(Fw); A = Fe.coeff(p, 2); B = Fe.coeff(p, 1)/2; C = Fe.coeff(p, 0); Feff = C - B**2/A
    f0 = sp.lambdify((r, E, L), Feff, "mpmath"); f1 = sp.lambdify((r, E, L), sp.diff(Feff, r), "mpmath")
    e0, l0 = (float(x) for x in mp.findroot(lambda e, l: (f0(1, e, l), f1(1, e, l)), (1 - u/2, mp.sqrt(u))))
    pst = lambda rr, ee: float((-B/A).subs({r: rr, E: ee, L: l0}))
    r0 = 1.0 + 1e-4
    En = sp.symbols("En")
    e1 = float(sp.nsolve(Fw.subs({L: l0, r: r0, p: pst(r0, e0)}).subs(E, En), En, e0))
    p1 = pst(r0, e1)
    per = 2*np.pi/abs(dFL(1, pst(1, e0), e0, l0))
    sol = solve_ivp(lambda t, y: [dFp(y[0], y[2], e1, l0), dFL(y[0], y[2], e1, l0), -dFr(y[0], y[2], e1, l0)],
                    [0, per*laps], [r0, 0.0, p1], rtol=1e-12, atol=1e-14, dense_output=True, max_step=per/400)
    t = np.linspace(0, per*laps, 400000); rr = sol.sol(t)[0]; dt = t[1] - t[0]
    phs = []
    for i in range(1, len(t) - 1):
        if rr[i] < rr[i - 1] and rr[i] <= rr[i + 1]:
            den = rr[i - 1] - 2*rr[i] + rr[i + 1]
            phs.append(sol.sol(t[i] - 0.5*dt*(rr[i + 1] - rr[i - 1])/den)[1])
    sim = float(np.mean(np.diff(phs) - 2*np.pi))
    pred = 2*np.pi*(1/float(mp.sqrt(orbit_factor(ks_shell(w))(u, chi))) - 1)
    return sim, pred


def main():
    out = {}
    # exactness of the Kerr-Schild orbit factor
    Tk = orbit_factor(ks_shell({}))
    out["exact"] = {f"u={u},chi={c}": mp.nstr(Tk(mp.mpf(u), mp.mpf(c)) - exact_T(mp.mpf(u), mp.mpf(c)), 5)
                    for u, c in (("0.05", "0.7"), ("0.05", "-0.7"), ("0.1", "0.3"))}
    # way 1: all 32 on/off settings of the five non-Newtonian pieces (EE always on)
    names = ["PP", "EP", "ES", "PS", "SS"]
    grid = {}
    for combo in itertools.product([0, 1], repeat=5):
        w = {"EE": 1, **dict(zip(names, combo))}
        key = "+".join(n for n, v in zip(names, combo) if v) or "E2 only"
        grid[key] = {k: str(v) for k, v in series(orbit_factor(ks_shell(w))).items()}
        print(key, grid[key], flush=True)
    out["double_copy"] = grid
    out["single_copy"] = {f"P{wP}S{wS}": {k: str(v) for k, v in series(orbit_factor(single_shell(wP, wS))).items()}
                          for wP in (0, 1) for wS in (0, 1)}
    # way 2: the dials
    dials = {"Kerr (beta=gamma=kappa=q=1)": dict(beta=1, kappa=1, gamma=1, q=1),
             "no push (kappa=0)": dict(beta=1, kappa=0, gamma=1, q=1),
             "half push (kappa=1/2)": dict(beta=1, kappa=sp.Rational(1, 2), gamma=1, q=1),
             "no shape (q=0)": dict(beta=1, kappa=1, gamma=1, q=0),
             "double shape (q=2)": dict(beta=1, kappa=1, gamma=1, q=2),
             "no space warp, no squared clock (beta=gamma=0)": dict(beta=0, kappa=1, gamma=0, q=1)}
    out["pn"] = {k: {kk: str(vv) for kk, vv in series(orbit_factor(pn_shell(**v))).items()} for k, v in dials.items()}
    # the push from rule A: boost the rest field h_00 = 2U, h_ij = 2U delta_ij to speed v: h_0i = -(2U + 2U) v
    out["rule_A_boost"] = "h_0i = -(h_00 + h_xx) v = -(2U [time part] + 2U [space part]) v, so kappa = (t + s)/2 = 1"
    # landmarks
    out["isco_zamo_speed"] = {str(c): mp.nstr(abs(zamo_speed(isco(mp.mpf(c)), mp.mpf(c))), 4)
                              for c in ("0", "0.3", "0.6", "0.9", "0.99", "0.999999", "-0.5", "-0.9")}
    out["isco_T"] = {str(c): mp.nstr(exact_T(1/isco(mp.mpf(c)), mp.mpf(c)), 5) for c in ("0", "0.5", "0.9", "-0.5")}
    out["product_form_gap"] = {f"r={rr},chi={c}": mp.nstr(product_form_gap(rr, mp.mpf(c)), 6)
                               for rr, c in ((10, "0"), (10, "0.5"), (30, "0.5"), (10, "-0.5"))}
    # simulated orbits
    sims = {}
    for nm, w in (("Kerr", {}), ("E2+ES only", dict(EE=1, PP=0, SS=0, EP=0, ES=1, PS=0)),
                  ("E2 only", dict(EE=1, PP=0, SS=0, EP=0, ES=0, PS=0))):
        for c in (0.7, -0.7):
            s_, p_ = simulate(w, 0.05, c)
            sims[f"{nm}, chi={c}"] = {"sim": s_, "from_T": p_, "rel_diff": abs(s_ - p_)/abs(p_)}
            print(nm, c, sims[f"{nm}, chi={c}"], flush=True)
    out["simulated"] = sims
    json.dump(out, open(os.path.join(HERE, "kerr_split.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
