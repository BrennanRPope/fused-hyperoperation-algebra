"""DRAFT, round 25 (Nov 27, 2026): is round 23's one mass term exactly what the river's push adds?

Round 23 found Kerr's orbit equations in Mino time are Newton's for the ring-shifted mass plus one term carrying the
mass in each (radial 2Mr(K - a^2), polar none, clock and turning 2Mr P/Delta). Round 20 found the river's time part
(E.E) alone gives the ring's polar law and the push makes the extra polar terms. Is the one term the push?

The Kerr-Schild shell with a dial on each piece of (k.p)^2, k.p = E + P + S (round 20):
  Sigma(2H + 1) = A p^2 + 2 B p + C_r + Theta-part,  A = r^2+a^2 - 2Mr wPP,  B = -2Mr (wEP E + wPS S),
  C_r = -a^2 L^2/(r^2+a^2) - (E^2-1) r^2 + C - 2Mr (wEE E^2 + wSS S^2 + 2 wES E S),  S = -aL/(r^2+a^2),
and R = (dr/dlambda)^2 = B^2 - A C_r for every setting.

Findings (all exact):
1. The time part alone is Newton's ring with mass M E^2 instead of M: R_EE = R_Newton(M -> M E^2). The orbit's energy,
   not its rest mass, is what the river's time part pulls on.
2. The one term is shared by three pieces, each a single 2Mr times something:
     2Mr(K - a^2) = 2Mr (E^2-1)(r^2+a^2)                    time part (M -> M E^2)
                  + 2Mr [C - (E^2-1) r^2 - a^2 L^2/(r^2+a^2)] push, P.P
                  + 2Mr [a^2 L^2/(r^2+a^2) - 2aEL]           twist, S.S and 2E.S
   E.P and P.S give only M^2 terms, and those cancel P.P's own M^2 terms exactly. That needs the pieces to form one
   square: the M^2 part vanishes exactly when wEE wPP = wEP^2, wES wPP = wEP wPS and wSS wPP = wPS^2.
3. So, read with round 23's Newton energy E^2 - 1, the term is not the push alone. But measured per unit of the
   orbit's energy (momenta divided by E, Newton's energy read as 1 - 1/E^2, Carter's C -> C/E^2, L -> L/E), Kerr's
   equations are exactly E^2 times Newton's ring (mass M) plus the push plus the twist (with E -> 1 in 2E.S). The time
   part is then exactly Newton's pull and nothing else in the path; its only other effect is the clock's 2Mr. Since
   every Mino rate scales by the same E, all of the orbit's shape (radial, up-and-down and turning rates against
   each other) is Newton's ring plus push and twist. Over the poles L = 0, so the push alone: Schwarzschild's 1 - 6u
   and round 23's (3/2)(1 + 10u - 8u^2) come from it, at every order and every eccentricity. On spherical polar orbits
   this shows as R and Theta being of degree one in (E^2 - 1, Q, Newton's mass): the time part's share only rescales.
4. The coordinate-time rate Omega_theta^2 needs the time part's clock too. Clock: time 2MrE, twist -2MraL/(r^2+a^2),
   push 4M^2 r^2 E/Delta, push-twist 4M^2 r^2 S/Delta. The turning rate comes from the twist pieces only.
5. Flight check: eccentric polar orbits (a = 0.9, r from 6 to 10; a = 0.5, r from 7 to 14) flown in the exact
   Kerr-Schild metric in Cartesian coordinates, in Mino time, beat at the ratio that Newton's ring + push predicts, to
   about 1e-10. Newton's ring alone is far off (1.011 against 0.493).
Run: python3 river_term.py  (writes river_term.json, about a minute)
"""
import json
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
r, M, a, E, L, C, p, th = sp.symbols("r M a E L C p theta", real=True)
u, k = sp.symbols("u k", positive=True)
NAMES = ("EE", "PP", "SS", "EP", "ES", "PS")
WS = dict(zip(NAMES, sp.symbols("wEE wPP wSS wEP wES wPS")))
D = r**2 + a**2
Dl = r**2 - 2*M*r + a**2
Qs = C - L**2
Kc = Qs + (L - a*E)**2
Pc = E*D - a*L


def shell(w):
    d = {n: w.get(n, WS[n]) for n in NAMES}
    S = -a*L/D
    A = D - 2*M*r*d["PP"]
    B = -2*M*r*(d["EP"]*E + d["PS"]*S)
    Cr = -a**2*L**2/D - (E**2 - 1)*r**2 + C - 2*M*r*(d["EE"]*E**2 + d["SS"]*S**2 + 2*d["ES"]*E*S)
    return A, B, Cr


def R_of(w):
    A, B, Cr = shell(w)
    return sp.cancel(B**2 - A*Cr)


def R_newton(m, e=None, c=None, l2=None):
    e = E**2 - 1 if e is None else e
    c = C if c is None else c
    l2 = L**2 if l2 is None else l2
    return D*(e*r**2 + 2*m*r - c) + a**2*l2


ALL = {n: 1 for n in NAMES}
ONLY_EE = {"EE": 1, "PP": 0, "SS": 0, "EP": 0, "ES": 0, "PS": 0}


def piece_identities():
    zero = lambda x: sp.simplify(x) == 0
    RK = Pc**2 - Dl*(r**2 + Kc)
    Rg = sp.expand(R_of({}))
    m2 = sp.factor(Rg.coeff(M, 2))
    al, be, ga = sp.symbols("alpha beta gamma")
    sq = {WS["EE"]: al**2, WS["PP"]: be**2, WS["SS"]: ga**2, WS["EP"]: al*be, WS["ES"]: al*ga, WS["PS"]: be*ga}
    first = {n: sp.factor(sp.diff(Rg, WS[n]).subs({WS[x]: 0 for x in NAMES}).coeff(M, 1)*M) for n in NAMES}
    time_share = 2*M*r*(E**2 - 1)*D
    push_share = 2*M*r*(C - (E**2 - 1)*r**2 - a**2*L**2/D)
    twist_share = 2*M*r*(a**2*L**2/D - 2*a*E*L)
    # clock and turning in Mino time, the part that does not average away (p -> -B/A on the orbit)
    A, B, Cr = shell({})
    pBL = -B/A
    clock = E*a**2*sp.cos(th)**2 - sp.Rational(1, 2)*sp.diff(2*B*p + Cr, E).subs(p, pBL)
    turn = L/sp.sin(th)**2 + sp.Rational(1, 2)*sp.diff(2*B*p + Cr, L).subs(p, pBL)
    Sg = r**2 + a**2*sp.cos(th)**2
    clock_full = sp.simplify(clock.subs({WS[n]: 1 for n in NAMES}))
    turn_full = sp.simplify(turn.subs({WS[n]: 1 for n in NAMES}))
    S = -a*L/D
    clock_pieces = {"time (EE)": 2*M*r*E, "twist (ES)": 2*M*r*S, "push (EP, with PP)": 4*M**2*r**2*E/Dl,
                    "push-twist (EP.PS, with PP)": 4*M**2*r**2*S/Dl}
    turn_pieces = {"twist (SS, ES)": 2*M*r*a*(E + S)/D, "push-twist (PS with EP, PP)": 4*M**2*r**2*a*(E + S)/(D*Dl)}
    # per unit energy: Kerr = E^2 x (Newton(M) + push + twist with E -> 1) at e~ = 1 - 1/E^2, C~ = C/E^2, L~ = L/E
    et, Ct, Lt = 1 - 1/E**2, C/E**2, L/E
    per_energy = E**2*(R_newton(M, et, Ct, Lt**2) + 2*M*r*(Ct - et*r**2 - a**2*Lt**2/D)
                       + 2*M*r*(a**2*Lt**2/D - 2*a*Lt))
    return {
        "shell_with_all_pieces_is_carter": zero(R_of(ALL) - RK),
        "kerr_is_newton_plus_one_term": zero(RK - R_newton(M) - 2*M*r*(Kc - a**2)),
        "time_part_is_newton_with_mass_M_E2": zero(R_of(ONLY_EE) - R_newton(M*E**2)),
        "theta_equation_has_no_dial": True,
        "first_order_pieces_from_flat": {n: str(v) for n, v in first.items()},
        "M2_part_general": str(m2),
        "M2_part_vanishes_with_all_pieces": zero(m2.subs({WS[n]: 1 for n in NAMES})),
        "M2_part_vanishes_for_any_square": zero(m2.subs(sq)),
        "M2_part_nonzero_for_EP_alone": not zero(m2.subs({**{WS[n]: 0 for n in NAMES}, WS["EP"]: 1})),
        "one_term_split_time_push_twist": zero(time_share + push_share + twist_share - 2*M*r*(Kc - a**2)),
        "first_order_PP_is_push_share_SS_ES_twist": zero(first["PP"] - push_share)
            and zero(first["SS"] + first["ES"] - twist_share) and zero(first["EE"] - 2*M*r*D - time_share),
        "shares": {"time": str(time_share), "push": str(sp.factor(push_share)), "twist": str(sp.factor(twist_share))},
        "clock_full_is_round23": zero(clock_full - (E*Sg + 2*M*r*Pc/Dl)),
        "clock_pieces_add_up": zero(sum(clock_pieces.values()) - 2*M*r*Pc/Dl),
        "clock_pieces": {n: str(v) for n, v in clock_pieces.items()},
        "turning_full_is_round23": zero(turn_full - (L/sp.sin(th)**2 - a**2*L/D + a*2*M*r*Pc/(D*Dl))),
        "turning_pieces_add_up": zero(sum(turn_pieces.values()) - a*2*M*r*Pc/(D*Dl)),
        "turning_has_no_time_or_push_alone_piece": zero(sp.diff(turn, WS["EE"])) and zero(
            turn.subs({WS["SS"]: 0, WS["ES"]: 0, WS["PS"]: 0}) - (L/sp.sin(th)**2 - a**2*L/D)),
        "kerr_per_unit_energy_is_newton_plus_push_plus_twist": zero(per_energy - RK),
        "clock_per_unit_energy": zero(clock_full - E*(Sg + 2*M*r*(D - a*Lt)/Dl)),
        "turning_per_unit_energy": zero(turn_full - E*(Lt/sp.sin(th)**2 - a**2*Lt/D + 2*M*r*a*(D - a*Lt)/(D*Dl))),
        "time_part_clock_per_unit_energy": "2Mr (redshift); push 4M^2 r^2/Delta; twist -2MraL~/(r^2+a^2) and "
                                           "-4M^2 r^2 a L~/((r^2+a^2) Delta)",
    }


# ---------------------------------------------------------------- over the poles: the shape is the push's
def polar_homogeneity():
    """R = (r^2+a^2)(e r^2 + 2 mu r - Q) + 2Mr(Q - e r^2) on spherical orbits: m and omega_r^2/Q don't depend on mu."""
    mu, e, Q = sp.symbols("mu e Q")
    R = D*(e*r**2 + 2*mu*r - Q) + 2*M*r*(Q - e*r**2)
    s = sp.solve([R, sp.diff(R, r)], [e, Q], dict=True)[0]
    wr2Q = sp.factor((-sp.diff(R, r, 2)/2).subs(s)/s[Q])
    m = sp.factor(-a**2*s[e]/s[Q])
    sc = sp.symbols("lambda")
    return {"degree_one": sp.expand(R.subs({e: sc*e, Q: sc*Q, mu: sc*mu}) - sc*R) == 0,
            "m_free_of_mu": sp.diff(m, mu) == 0, "omega_r2_over_Q_free_of_mu": sp.diff(wr2Q, mu) == 0,
            "m": str(m), "omega_r2_over_Q": str(wr2Q)}


def polar_k2_table():
    """Exact (a/r)^2 parts with the one term's time and push shares on dials (twist is absent over the poles)."""
    import polar_coupling as pc
    mm, Q = pc.mm, pc.Q

    def polar(tau, pi):
        e = sp.symbols("e")
        R = D*(e*r**2 + 2*M*r - Q) + tau*2*M*r*e*D + pi*2*M*r*(Q - e*r**2)
        s = sp.solve([R, sp.diff(R, r)], [e, Q], dict=True)[0]
        return [x.subs({r: 1, M: u, a: k}) for x in (sp.factor(s[e]), sp.factor(s[Q]),
                sp.factor((-sp.diff(R, r, 2)/2).subs(s)), sp.factor(-a**2*s[e]/s[Q]))]

    def parts(tau, pi, c1, c2):
        e, Qv, wr2, m = polar(tau, pi)
        t = sp.symbols("t", positive=True); st = {k: sp.sqrt(t)}
        m_t, e_t = sp.cancel(m.subs(st)), sp.cancel(e.subs(st))
        ratio_t = sp.cancel((wr2/Qv).subs(st))*(2*pc.Ks.subs(mm, m_t)/sp.pi)**2
        Ev = sp.sqrt(1 + e_t)
        clock_t = (1 + c1*(Ev - 1))*(1 + t*pc.COS2.subs(mm, m_t)) + c1*2*u*Ev + c2*4*u**2*Ev/(1 - 2*u + t)
        om2_t = (sp.pi**2/4)*sp.cancel((Qv/u).subs(st))/(pc.Ks.subs(mm, m_t)**2*clock_t**2)
        k0 = lambda f: sp.factor(sp.simplify(f.subs(t, 0)))
        k2 = lambda f: sp.factor(sp.simplify(sp.diff(f, t).subs(t, 0)))
        return [k0(ratio_t), k2(ratio_t), k0(om2_t), k2(om2_t)]

    cases = {"Newton": (0, 0, 0, 0), "time part only": (1, 0, 1, 0), "push only": (0, 1, 0, 1),
             "time + push, radial only": (1, 1, 0, 0), "Kerr (time + push)": (1, 1, 1, 1)}
    out = {}
    for name, cs in cases.items():
        res = parts(*cs)
        out[name] = {"ratio_k0": str(res[0]), "ratio_k2": str(res[1]),
                     "ratio_k2_series": str(sp.series(res[1], u, 0, 4).removeO()),
                     "omega_theta2_k0_series": str(sp.series(res[2], u, 0, 4).removeO()),
                     "omega_theta2_k2_series": str(sp.series(res[3], u, 0, 4).removeO()),
                     "_exact": [str(x) for x in res]}
        print(name, out[name]["ratio_k2_series"], out[name]["omega_theta2_k2_series"], flush=True)
    return out


# ---------------------------------------------------------------- eccentric polar orbit: quadrature and flight
def kerr_polar_constants(Mv, av, r1, r2):
    """E, Q of the L = 0 Kerr orbit with turning points r1, r2."""
    Ev, Qv = sp.symbols("Ev Qv", positive=True)
    RK = lambda rr: (Ev*(rr**2 + av**2))**2 - (rr**2 - 2*Mv*rr + av**2)*(rr**2 + Qv + av**2*Ev**2)
    sol = sp.nsolve([RK(r1), RK(r2)], [Ev, Qv], [0.96, 12.0], prec=30)
    return float(sol[0]), float(sol[1])


def mino_ratio(Rpoly, e_th, Q_th, av, r1, r2):
    """Upsilon_r / Upsilon_theta for a polynomial R(r) with turning points r1 < r2 and Theta = Q + e a^2 cos^2."""
    rr = sp.symbols("rr")
    Rr = sp.Poly(sp.expand(Rpoly(rr)), rr)
    # R = (r - r1)(r2 - r) G(r), G > 0 between: substitute r = mid + half cos(chi)
    G = sp.lambdify(rr, sp.cancel(Rr.as_expr()/((rr - r1)*(r2 - rr))), "numpy")
    mid, half = (r1 + r2)/2, (r2 - r1)/2
    Ir = quad(lambda c: 1/np.sqrt(G(mid + half*np.cos(c))), 0, np.pi, epsabs=1e-14, epsrel=1e-14)[0]
    Ith = quad(lambda ps: 1/np.sqrt(Q_th + e_th*av**2*np.cos(ps)**2), 0, 2*np.pi, epsabs=1e-14, epsrel=1e-14)[0]
    return (np.pi/Ir)/(2*np.pi/Ith)          # radial: half period Ir in Mino time; theta: full circuit Ith


def eccentric_polar(Mv=1.0, av=0.9, r1=6.0, r2=10.0, laps=60):
    Ev, Qv = kerr_polar_constants(Mv, av, r1, r2)
    e = Ev**2 - 1
    Dn = lambda rr: rr**2 + av**2
    kerr = lambda rr: Dn(rr)*(e*rr**2 + 2*Mv*rr - Qv) + 2*Mv*rr*(Qv + av**2*e)
    et, Qt = e/Ev**2, Qv/Ev**2
    push = lambda rr: Dn(rr)*(et*rr**2 + 2*Mv*rr - Qt) + 2*Mv*rr*(Qt - et*rr**2)
    rat_kerr = mino_ratio(kerr, e, Qv, av, r1, r2)
    rat_push = mino_ratio(push, et, Qt, av, r1, r2)
    # Newton's ring on its own with the same turning points (constants refitted)
    es, Qn = sp.symbols("es Qn")
    RN = lambda rr: Dn(rr)*(es*rr**2 + 2*Mv*rr - Qn)
    sN = sp.solve([RN(r1), RN(r2)], [es, Qn], dict=True)[0]
    rat_newton = mino_ratio(lambda rr: Dn(rr)*(float(sN[es])*rr**2 + 2*Mv*rr - float(sN[Qn])),
                            float(sN[es]), float(sN[Qn]), av, r1, r2)
    # flight in the exact Kerr-Schild metric, from the pericentre at the equator
    import polar as pl
    fn, rfn, Hfn = pl.ks_hamiltonian()
    X, Y, Z, As = sp.symbols("x y z a", real=True)
    rho2 = X**2 + Y**2 + Z**2
    rr_ = sp.sqrt((rho2 - As**2 + sp.sqrt((rho2 - As**2)**2 + 4*As**2*Z**2))/2)
    grad_r = sp.lambdify([X, Y, Z, As], [sp.diff(rr_, c) for c in (X, Y, Z)], "numpy")
    Dv = r1**2 - 2*Mv*r1 + av**2
    J = np.zeros((4, 4)); J[0, 0] = 1; J[0, 1] = 2*Mv*r1/Dv
    J[1, 1] = 1 - av*av/Dv; J[1, 3] = -av; J[2, 1] = r1*av/Dv; J[2, 3] = r1; J[3, 2] = -r1
    pKS = np.linalg.solve(J.T, np.array([-Ev, 0.0, -np.sqrt(Qv), 0.0]))

    def flow(lam, s):                      # Mino time: d/dlambda = Sigma d/dtau, so both motions are periodic
        o = fn(s[1], s[2], s[3], *s[4:], Mv, av)
        rv = rfn(s[1], s[2], s[3], av)
        Sg = rv**2 + av**2*s[3]**2/rv**2
        return [Sg*o[0], Sg*o[1], Sg*o[2], Sg*o[3], 0.0, -Sg*o[4], -Sg*o[5], -Sg*o[6]]

    def up(tau, s):
        return s[3]
    up.direction = 1

    def peri(tau, s):
        o = fn(s[1], s[2], s[3], *s[4:], Mv, av)
        g = grad_r(s[1], s[2], s[3], av)
        return g[0]*o[1] + g[1]*o[2] + g[2]*o[3]
    peri.direction = 1
    Tr = 2*np.pi*np.sqrt(((r1 + r2)/2)**3/Mv)/((r1 + r2)/2)**2   # rough radial period in Mino time
    sol = solve_ivp(flow, [0, laps*Tr], [0.0, r1, av, 0.0, *pKS], rtol=1e-12, atol=1e-12, events=[up, peri],
                    method="DOP853", dense_output=True)
    tu, tp = sol.t_events
    ts = np.linspace(0, min(sol.t[-1], 8*(tp[-1] - tp[0])/max(len(tp) - 1, 1)), 900)
    Yp = sol.sol(ts)
    flown = ((len(tp) - 1)/(tp[-1] - tp[0]))/((len(tu) - 1)/(tu[-1] - tu[0]))
    return dict(M=Mv, a=av, r1=r1, r2=r2, E=Ev, Q=Qv, ratio_kerr_quadrature=float(rat_kerr),
                ratio_newton_ring_plus_push=float(rat_push), ratio_newton_ring_alone=float(rat_newton),
                ratio_flown=float(flown), radial_laps=len(tp) - 1, theta_circuits=len(tu) - 1,
                path=[Yp[1].tolist(), Yp[2].tolist(), Yp[3].tolist()])


def main():
    out = {"pieces": piece_identities()}
    print(json.dumps({k_: v for k_, v in out["pieces"].items() if isinstance(v, bool)}, indent=1), flush=True)
    out["polar_homogeneity"] = polar_homogeneity()
    print(out["polar_homogeneity"], flush=True)
    out["polar_k2"] = polar_k2_table()
    out["eccentric"] = [eccentric_polar(1.0, 0.9, 6.0, 10.0), eccentric_polar(1.0, 0.5, 7.0, 14.0)]
    for x in out["eccentric"]:
        print({k_: v for k_, v in x.items() if k_.startswith("ratio") or k_ in ("radial_laps", "theta_circuits")},
              flush=True)
    json.dump(out, open(os.path.join(HERE, "river_term.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
