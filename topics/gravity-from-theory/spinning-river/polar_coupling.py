"""DRAFT, round 23 (Nov 27, 2026): is the mass's coupling to the turning the source of the polar a^2 M/r^3 terms?

Round 22 wrote Carter's radial equation in Mino time as R = r^4 (E^2-1) + 2Mr^3 - K0 r^2 + 2Mr x^2 and found one law for
both wobbles of an equatorial orbit. Round 23 asks whether the same coupling, 2Mr x^2 with x = -aE over the poles, makes
round 19's mixed terms (+15 (a/r)^2 (M/r) in (Omega_r/Omega_theta)^2, -4 (a/r)^2 (M/r) in Omega_theta^2 r^3/M).

Findings (all exact unless marked):
1. Kerr's orbit equations in Mino time are Newton's equations for the ring-shifted mass (Euler's problem with the two
   centres at z = +-ia, in oblate spheroidal coordinates), with E^2 - 1 in place of twice Newton's energy, plus one term
   in each equation that carries the mass:
     radial   R_Kerr   = R_Newton + 2Mr (K - a^2)          (K = Carter's constant, Q + (L - aE)^2)
     polar    Theta_Kerr = Theta_Newton                     (no extra term at all)
     clock    dt/dlam  = E Sigma + 2Mr P/Delta              (Newton: Sigma; P = E(r^2+a^2) - aL)
     turning  dphi/dlam = Newton's + 2Mr a P/((r^2+a^2) Delta)
   In Schwarzschild the radial term is 2Mr L^2, the textbook term behind perihelion precession.
2. The stiffness law holds for every spherical orbit, not only the equator: with W the r^2 coefficient of R,
     omega_r^2 = M(r - 3K/r) + 4a^2 Q/r^2,   W = M(r + 3K/r) - 2a^2 Q/r^2      (Mino time).
   In the equator K = x^2 and Q = 0: round 22's law.
3. Over the poles (L = 0): (Omega_r/Omega_theta)^2 = (1+m)(2K(m)/pi)^2 * B, with B = omega_r^2/W, which is exactly 1 at
   Newtonian strength. At order a^2, to every order in u = M/r:
     ratio = 1 - 6u + k^2 [ (3/2)(1-4u)(1-6u) + 6u(5-8u) ] + O(k^4),   k = a/r.
4. Answer: the coupling's x^2 = a^2 E^2 is only part of the source. Its a^2 part (2M a^2 r) is Newtonian: it is the ring's
   oblateness and it is in Euler's problem too; switching it off breaks the Newtonian limit. Its binding part
   2M a^2 r (E^2-1) is relativistic and gives 6 of the 15 and -24 of the -12. The rest (9 and +12) comes from the other
   half of the same one term, 2MrQ (Schwarzschild's 2Mr L^2 with L^2 -> Q), acting on the ring's oblate orbit. With the
   whole term 2Mr(K - a^2) off, every mixed term vanishes and Newton's polar law returns exactly.
5. Omega_theta^2 r^3/M also needs the clock: its -4 is -21/2 from the radial term and +13/2 from the clock terms, and
   at order a^2 it is exact: Omega_theta^2 r^3/M = 1 - k^2 (3/2 + 4u) + O(k^4), with no higher powers of u.
6. Independent check: perturbed polar orbits flown in Cartesian coordinates (Mino time) in Newton's ring field beat at
   (1+k^2)(2K(k^2)/pi)^2, which is R_Newton's prediction (round 18's law); the gap shrinks in step with the kick
   (about 9 delta at k = 0.9, so 1e-5 at delta = 1e-6).
Run: python3 polar_coupling.py  (writes polar_coupling.json, about 40 seconds)
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
r, M, a, u, k, mm = sp.symbols("r M a u k m", positive=True)
E, L, Q, th = sp.symbols("E L Q theta", real=True)
N, U = 8, 4


# ---------------------------------------------------------------- 1. Kerr = Newton's ring + one mass term
def newton_plus_one():
    e = E**2 - 1
    D = r**2 - 2*M*r + a**2
    P = E*(r**2 + a**2) - a*L
    K = Q + (L - a*E)**2
    C = Q + L**2
    S = r**2 + a**2*sp.cos(th)**2
    R_K = P**2 - D*(r**2 + K)
    R_N = (r**2 + a**2)*(e*r**2 + 2*M*r - C) + a**2*L**2           # Newton, Mino time, oblate coordinates
    Th_K = Q - sp.cos(th)**2*(a**2*(1 - E**2) + L**2/sp.sin(th)**2)
    Th_N = C + e*a**2*sp.cos(th)**2 - L**2/sp.sin(th)**2
    T_K = (r**2 + a**2)*P/D - a*(a*E*sp.sin(th)**2 - L)
    F_K = a*P/D - (a*E - L/sp.sin(th)**2)
    F_N = L/sp.sin(th)**2 - a**2*L/(r**2 + a**2)
    return {
        "radial_R_Kerr_minus_Newton_minus_2Mr(K-a^2)": sp.expand(R_K - R_N - 2*M*r*(K - a**2)) == 0,
        "polar_Theta_Kerr_minus_Newton": sp.simplify(Th_K - Th_N) == 0,
        "clock_minus_(E_Sigma+2MrP/Delta)": sp.simplify(T_K - E*S - 2*M*r*P/D) == 0,
        "turning_minus_(Newton+2Mr_aP/((r^2+a^2)Delta))": sp.simplify(F_K - F_N - 2*M*r*a*P/((r**2 + a**2)*D)) == 0,
        "schwarzschild_extra_term_is_2Mr(Q+L^2)": sp.expand((2*M*r*(K - a**2)).subs(a, 0) - 2*M*r*(Q + L**2)) == 0,
        "newton_ring_potential_is_Mr/Sigma": True,   # see ring_field(): Re 1/sqrt(rho^2+(z-ia)^2) = r/Sigma, checked there
    }


def ring_field():
    """Re 1/sqrt(rho^2 + (z - ia)^2) = r/(r^2 + a^2 cos^2 theta) in oblate coordinates: x = sqrt(r^2+a^2) sin(th)."""
    c = sp.symbols("c", real=True)
    s = r - sp.I*a*c                                                  # candidate square root
    w = (r**2 + a**2)*(1 - c**2) + (r*c - sp.I*a)**2
    return {"square_root": sp.expand(s**2 - w) == 0,
            "real_part": sp.simplify(sp.re(sp.expand_complex(1/s)) - r/(r**2 + a**2*c**2)) == 0}


# ---------------------------------------------------------------- 2. the stiffness law for every spherical orbit
def stiffness_law():
    e, W, Kc = sp.symbols("e W K")
    R = r**4*e + 2*M*r**3 - W*r**2 + 2*M*r*Kc - a**2*Q
    D = r**2 - 2*M*r + a**2
    P = E*(r**2 + a**2) - a*L
    R_K = P**2 - D*(r**2 + Q + (L - a*E)**2)
    form = sp.expand(R.subs({e: E**2 - 1, W: Q + L**2 - a**2*(E**2 - 1), Kc: Q + (L - a*E)**2}) - R_K) == 0
    sol = sp.solve([R, sp.diff(R, r)], [e, W], dict=True)[0]
    wr2 = sp.simplify((-sp.diff(R, r, 2)/2).subs(sol))
    return {"R_in_mino_form_all_orbits": form,
            "omega_r2_law": sp.simplify(wr2 - (M*r - 3*M*Kc/r + 4*a**2*Q/r**2)) == 0,
            "W_law": sp.simplify(sol[W] - (M*r + 3*M*Kc/r - 2*a**2*Q/r**2)) == 0,
            "sum_2Mr_plus_2a2Q/r2": sp.simplify(wr2 + sol[W] - 2*M*r - 2*a**2*Q/r**2) == 0,
            "equator_reduces_to_round_22": sp.simplify((wr2 - (M*r - 3*M*Kc/r)).subs(Q, 0)) == 0}


# ---------------------------------------------------------------- 3. polar orbits with dials on each mass term
Ks = sp.series(sp.elliptic_k(mm), mm, 0, N//2 + 1).removeO()
Es = sp.series(sp.elliptic_e(mm), mm, 0, N//2 + 1).removeO()
COS2 = sp.expand(sp.series(sp.cancel((Ks - Es)/(mm*Ks)), mm, 0, N//2).removeO())   # Mino-time mean of cos^2 theta


def polar_orbit(eq=1, mu=1):
    """Spherical L = 0 orbit of R = Newton + eq*2MrQ + mu*2M a^2 r (E^2-1); Kerr at eq = mu = 1."""
    e = sp.symbols("e")
    R = r**4*e + 2*M*r**3 - (Q - a**2*e)*r**2 + 2*M*a**2*r + eq*2*M*r*Q + mu*2*M*a**2*r*e - a**2*Q
    s = sp.solve([R, sp.diff(R, r)], [e, Q], dict=True)[0]
    wr2 = sp.factor((-sp.diff(R, r, 2)/2).subs(s))
    W = sp.factor((Q - a**2*e).subs(s))
    return dict(e=sp.factor(s[e]), Q=sp.factor(s[Q]), wr2=wr2, W=W, m=sp.factor(-a**2*s[e]/s[Q]))


def ser2(expr):
    """Rational function of (u, k) -> polynomial truncated below k^N and u^U."""
    ser = sp.series(sp.cancel(expr), k, 0, N).removeO()
    return sp.expand(sum(sp.series(sp.cancel(ser.coeff(k, j)), u, 0, U).removeO()*k**j for j in range(N)))


def trunc(p):
    p = sp.expand(p)
    if p == 0:
        return p
    return sum(c*k**i*u**j for (i, j), c in sp.Poly(p, k, u).terms() if i < N and j < U)


def tpow(p, n):
    out = sp.Integer(1)
    for _ in range(n):
        out = trunc(out*p)
    return out


def tsub(poly_m, m_s):
    """A polynomial in m with m -> m_s (m_s starts at k^2)."""
    cs = sp.Poly(poly_m, mm).all_coeffs()[::-1]
    return trunc(sum(c*tpow(m_s, i) for i, c in enumerate(cs)))


def tinv_sq(p):
    """1/p^2 for p = 1 + (terms that vanish at k = u = 0)."""
    x = trunc(p - 1)
    return trunc(sum((n + 1)*(-1)**n*tpow(x, n) for n in range(N//2 + U + 1)))


def rows_of(p):
    return [[sp.nsimplify(sp.expand(p).coeff(k, j).coeff(u, i)) for i in range(U)] for j in range(0, N, 2)]


def polar_quantities(eq=1, mu=1, c1=1, c2=1):
    """Double series (k^0..k^6, u^0..u^3) of the ratio and of Omega_theta^2 r^3/M, and their exact k^2 parts."""
    o = polar_orbit(eq, mu)
    sub = {r: 1, M: u, a: k}
    e, Qv, wr2, m = (o[x].subs(sub) for x in ("e", "Q", "wr2", "m"))
    m_s, e_s = ser2(m), ser2(e)
    g = tsub(sp.expand(2*Ks/sp.pi), m_s)                                   # 2K(m)/pi
    ratio_s = trunc(ser2(wr2/Qv)*trunc(g*g))
    Ev_s = trunc(sum(sp.binomial(sp.Rational(1, 2), n)*tpow(e_s, n) for n in range(U + 1)))
    clock_s = trunc((1 + c1*(Ev_s - 1))*(1 + k**2*tsub(COS2, m_s)) + c2*2*u*Ev_s*ser2((1 + k**2)/(1 - 2*u + k**2)))
    om2_s = trunc(ser2(Qv/u)*tinv_sq(g)*tinv_sq(clock_s))
    # exact k^2 parts: put k^2 = t and take d/dt at t = 0
    t = sp.symbols("t", positive=True)
    st = {k: sp.sqrt(t)}
    m_t, e_t = sp.cancel(m.subs(st)), sp.cancel(e.subs(st))
    ratio_t = sp.cancel((wr2/Qv).subs(st))*(2*Ks.subs(mm, m_t)/sp.pi)**2
    Ev = sp.sqrt(1 + e_t)
    clock_t = (1 + c1*(Ev - 1))*(1 + t*COS2.subs(mm, m_t)) + c2*2*u*Ev*(1 + t)/(1 - 2*u + t)
    om2_t = (sp.pi**2/4)*sp.cancel((Qv/u).subs(st))/(Ks.subs(mm, m_t)**2*clock_t**2)
    k2 = lambda f: sp.factor(sp.simplify(sp.diff(f, t).subs(t, 0)))
    return dict(ratio=rows_of(ratio_s), omega_theta2=rows_of(om2_s), ratio_k2_exact=k2(ratio_t),
                omega_theta2_k2_exact=k2(om2_t))



def dials():
    out = {}
    cases = {"Kerr": (1, 1, 1, 1), "binding part of x^2 off": (1, 0, 1, 1), "2MrQ off": (0, 1, 1, 1),
             "whole term 2Mr(K-a^2) off": (0, 0, 1, 1), "clock terms off": (1, 1, 0, 0),
             "everything off (Newton)": (0, 0, 0, 0)}
    for name, cs in cases.items():
        out[name] = polar_quantities(*cs)
        print(name, out[name]["ratio"][1], out[name]["omega_theta2"][1], flush=True)
    return out


def polar_checks(tab):
    o = polar_orbit()
    pref19 = ((r**2 + a**2)**3 - 2*M*r*(3*r**4 - 2*a**2*r**2 + 3*a**4))/(r**2*((r**2 + a**2)**2 - 4*M*a**2*r))
    Kc = o["Q"] + a**2*(1 + o["e"])
    B = sp.cancel((o["wr2"]/o["W"]).subs({r: 1, M: u, a: k}))
    newton_ratio = sp.series((1 + k**2)*(2*sp.elliptic_k(k**2)/sp.pi)**2, k, 0, N).removeO()
    newton_om2 = sp.series((sp.pi/2)**2/(2*sp.elliptic_k(k**2) - sp.elliptic_e(k**2))**2, k, 0, N).removeO()
    nt = tab["everything off (Newton)"]
    return {
        "pref_equals_omega_r2_over_Q_round19": sp.simplify(o["wr2"]/o["Q"] - pref19) == 0,
        "omega_r2_law_holds_on_polar_orbit": sp.simplify(o["wr2"] - (M*r - 3*M*Kc/r + 4*a**2*o["Q"]/r**2)) == 0,
        "W_equals_Q(1+m)": sp.simplify(o["W"] - o["Q"]*(1 + o["m"])) == 0,
        "B_is_1_at_Newton": sp.simplify(B.subs(u, 0)) == 1,
        "kerr_rows_match_round19": [[str(x) for x in row] for row in tab["Kerr"]["ratio"]] ==
            [["1", "-6", "0", "0"], ["3/2", "15", "-12", "0"], ["27/32", "-429/16", "54", "-33"],
             ["39/64", "1203/32", "-1179/8", "429/2"]],
        "kerr_omega_theta2_k2u_is_-4": tab["Kerr"]["omega_theta2"][1][1] == -4,
        "k2_law": sp.simplify(tab["Kerr"]["ratio_k2_exact"] - (1 - 6*u)*0 -
                              (sp.Rational(3, 2)*(1 - 4*u)*(1 - 6*u) + 6*u*(5 - 8*u))) == 0,
        "newton_returns_ratio": all(nt["ratio"][j][0] == sp.nsimplify(newton_ratio.coeff(k, 2*j)) and
                                    all(x == 0 for x in nt["ratio"][j][1:]) for j in range(N//2)),
        "newton_returns_omega_theta2": all(nt["omega_theta2"][j][0] == sp.nsimplify(newton_om2.coeff(k, 2*j)) and
                                           all(x == 0 for x in nt["omega_theta2"][j][1:]) for j in range(N//2)),
        "whole_term_off_kills_mixed_ratio_terms": all(all(x == 0 for x in row[1:])
                                                      for row in tab["whole term 2Mr(K-a^2) off"]["ratio"]),
        "split_k2u_ratio_15=9+6": (tab["Kerr"]["ratio"][1][1], tab["binding part of x^2 off"]["ratio"][1][1],
                                   tab["2MrQ off"]["ratio"][1][1]) == (15, 9, 6),
        "split_k2u2_ratio": [str(tab[c]["ratio"][1][2]) for c in ("Kerr", "binding part of x^2 off", "2MrQ off")],
        "split_k2u_omega_theta2_-4=-21/2+13/2": (tab["clock terms off"]["omega_theta2"][1][1],
                                                 tab["whole term 2Mr(K-a^2) off"]["omega_theta2"][1][1]) ==
                                                (sp.Rational(-21, 2), sp.Rational(13, 2)),
    }


# ---------------------------------------------------------------- 4. independent flight in Newton's ring field
def ring_flight(kv, laps=40, delta=1e-6):
    Mv, r0, av = 1.0, 1.0, float(kv)
    # Newton's spherical polar orbit from R_Newton = R_Newton' = 0 (L = 0): e and C
    Rn = lambda e, C, x: (x**2 + av**2)*(e*x**2 + 2*Mv*x - C)
    A = np.array([[r0**2*(r0**2 + av**2), -(r0**2 + av**2)], [4*r0**3 + 2*av**2*r0, -2*r0]])
    b = -np.array([2*Mv*r0*(r0**2 + av**2), 2*Mv*(3*r0**2 + av**2)])
    e, C = np.linalg.solve(A, b)
    v0 = np.sqrt(C)/r0*(1 + delta)
    s0 = [np.sqrt(r0**2 + av**2), 0.0, 0.0, 0.0, 0.0, -v0]

    def flow(lam, s):                                                  # Mino time: dt/dlam = Sigma
        x, y, z = s[:3]
        w = np.sqrt(complex(x*x + y*y + (z - 1j*av)**2))
        f = -Mv*np.real(np.array([x, y, z - 1j*av])/w**3)
        S2 = x*x + y*y + z*z - av*av
        rr2 = (S2 + np.sqrt(S2*S2 + 4*av*av*z*z))/2
        Sig = rr2 + av*av*z*z/rr2
        return [Sig*s[3], Sig*s[4], Sig*s[5], *(Sig*f)]

    def r2dot(t, s):
        x, y, z, vx, vy, vz = s
        S = x*x + y*y + z*z - av*av; Sd = 2*(x*vx + y*vy + z*vz)
        return Sd + (S*Sd + 4*av*av*z*vz)/np.sqrt(S*S + 4*av*av*z*z)
    r2dot.direction = -1

    def up(t, s):
        return s[2]
    up.direction = 1
    sol = solve_ivp(flow, [0, laps*2*np.pi/np.sqrt(C)], s0, rtol=1e-12, atol=1e-13, events=[up, r2dot])
    tz, tr = sol.t_events
    Om_th = 2*np.pi*(len(tz) - 1)/(tz[-1] - tz[0])
    Om_r = 2*np.pi*(len(tr) - 1)/(tr[-1] - tr[0])
    law = (1 + kv**2)*(2*ellipk(kv**2)/np.pi)**2
    return dict(k=kv, delta=delta, newton_e=float(e), newton_C=float(C), flown_ratio2=float((Om_r/Om_th)**2),
                law=float(law), rel_err=float(abs((Om_r/Om_th)**2/law - 1)), theta_laps=len(tz) - 1,
                radial_beats=len(tr) - 1)


def x2_whole_off():
    """Switch off all of x^2 = a^2 E^2 in the coupling (its Newtonian a^2 part too): the Newtonian column moves."""
    e = sp.symbols("e")
    R = r**4*e + 2*M*r**3 - (Q - a**2*e)*r**2 + 2*M*r*Q - a**2*Q
    s = sp.solve([R, sp.diff(R, r)], [e, Q], dict=True)[0]
    sub = {r: 1, M: u, a: k}
    wr2 = (-sp.diff(R, r, 2)/2).subs(s).subs(sub)
    m = (-a**2*s[e]/s[Q]).subs(sub)
    ratio = sp.cancel(wr2/s[Q].subs(sub))*(2*Ks.subs(mm, m)/sp.pi)**2
    return {"ratio_k2_newton_coefficient": str(sp.limit(sp.diff(ratio, k, 2).subs(k, 0)/2, u, 0)),
            "kerr_value": "3/2"}


def curves(tab):
    us = [0.005*i for i in range(0, 31)]
    out = {"u": us}
    for name in ("Kerr", "binding part of x^2 off", "2MrQ off", "whole term 2Mr(K-a^2) off"):
        f = sp.lambdify(u, tab[name]["ratio_k2_exact"])
        out["ratio_k2: " + name] = [float(f(x)) for x in us]
    for name in ("Kerr", "clock terms off", "whole term 2Mr(K-a^2) off", "everything off (Newton)"):
        f = sp.lambdify(u, tab[name]["omega_theta2_k2_exact"])
        out["omega_theta2_k2: " + name] = [float(f(x)) for x in us]
    # the radial function of one Kerr polar orbit split into Newton's part and the one mass term (M = 1, a = 0.9, r0 = 6)
    o = polar_orbit()
    vals = {M: 1, a: sp.Rational(9, 10), r: 6}
    ev, Qv = float(o["e"].subs(vals)), float(o["Q"].subs(vals))
    Kv = Qv + 0.81*(1 + ev)
    xs = np.linspace(3.0, 12.0, 181)
    RN = (xs**2 + 0.81)*(ev*xs**2 + 2*xs - Qv)
    extra = 2*xs*(Kv - 0.81)
    out["R_split"] = {"r": xs.tolist(), "newton": RN.tolist(), "extra": extra.tolist(), "kerr": (RN + extra).tolist(),
                      "r0": 6, "a": 0.9}
    return out


def main():
    out = {"newton_plus_one": newton_plus_one(), "ring_field": ring_field(), "stiffness_law": stiffness_law()}
    print(out, flush=True)
    tab = dials()
    out["polar_checks"] = {kk: (v if isinstance(v, (bool, list)) else bool(v)) for kk, v in polar_checks(tab).items()}
    print(out["polar_checks"], flush=True)
    out["dials"] = {n: {"ratio": [[str(x) for x in row] for row in d["ratio"]],
                        "omega_theta2": [[str(x) for x in row] for row in d["omega_theta2"]],
                        "ratio_k2_exact": str(d["ratio_k2_exact"]),
                        "omega_theta2_k2_exact": str(d["omega_theta2_k2_exact"])} for n, d in tab.items()}
    out["x2_whole_off"] = x2_whole_off()                      # not a counted check: shows the a^2 part is Newton's
    print(out["x2_whole_off"], flush=True)
    out["flights"] = [ring_flight(x, delta=d) for d in (1e-5, 1e-6) for x in (0.3, 0.6, 0.9)]   # error ~ delta
    print([(f["k"], f["flown_ratio2"], f["law"], f["rel_err"]) for f in out["flights"]], flush=True)
    out["curves"] = curves(tab)
    json.dump(out, open(os.path.join(HERE, "polar_coupling.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
