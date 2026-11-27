"""DRAFT, round 20 (Nov 27, 2026): Kerr's polar orbit split into the field's pieces.

Round 19 found that Kerr's orbit over the poles follows the shifted source's exact polar law at Newtonian strength,
but adds terms in a^2 M / r^3 that the equator doesn't have. Are they frame dragging acting twice (the rays' twist S,
squared), or the shape meeting the mass?

The Kerr-Schild mass shell (round 16, now off the equator): in the flat space's own oblate coordinates,
  Sigma (2H + 1) = (r^2+a^2) p_r^2 + p_theta^2 + L^2 Sigma/((r^2+a^2) sin^2) - (E^2-1) Sigma - 2 M r (k.p)^2,
  k.p = E + P + S,  P = p_r (the river's radial push),  S = -a L/(r^2+a^2) (the rays' twist).
The factor 2Mr/Sigma = v^2 is the river's speed squared (Part 4's velocity rule), so each piece of (k.p)^2 is the
river's flow meeting the particle's energy (E), its radial motion (P) or its turning (S). The shell separates in r and
theta for every setting of the six dials (EE, PP, SS, 2EP, 2ES, 2PS), because the 2Mr factor depends on r alone.

Findings:
1. For an orbit with no spin about the axis, S = 0, so the twist pieces (SS, ES, PS) drop out of the radial and
   up-and-down rates and of the clock entirely: frame dragging twice is ruled out. The twist only turns the orbit's
   plane, at the exact rate 2 M r a E/Delta per unit Mino time, from ES and PS (first power of the twist).
2. E.E alone gives exactly the shifted source's polar ratio (Omega_r/Omega_theta)^2 = (1 + a^2/r^2)(2K(a/r)/pi)^2, at every
   M, as in the equator, where E.E alone gives exactly 1 - 3a^2/r^2.
3. The extra a^2 M/r^3 terms come from the river's radial push: in the ratio, 2E.P gives (-4 + 16k^2 - 209k^4/8) u and
   P.P gives (-2 - k^2 - 11k^4/16) u; they add to Kerr's (-6 + 15k^2 - 429k^4/16) u. 2E.P also moves the up-and-down
   motion's elliptic modulus off a^2/r^2.
4. The equator, by the same split: the push brings the same kind of term there too (k^2 u = chi^2 u^3: P.P +8, 2E.P +20),
   and the twist cancels it exactly (S.S -16, 2E.S -12, P.S 0; the pieces add in all 32 settings). So the equator is
   clean because frame dragging cancels the push's shape terms; over the poles there is no twist to cancel them.
Run: python3 polar_split.py  (writes polar_split.json, about 3 minutes)
"""
import itertools
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
r, M, a, E, L, x = sp.symbols("r M a E L x", positive=True)
E2, Kc = sp.symbols("E2 Kc", positive=True)
u, k, mm = sp.symbols("u k m", positive=True)
NAMES = ("EE", "PP", "SS", "EP", "ES", "PS")
WS = dict(zip(NAMES, sp.symbols("wEE wPP wSS wEP wES wPS")))


def radial_quadratic(w, Lv=L, rr=x):
    """The radial part of the shell as A p^2 + 2 B p + C (theta part moved to the separation constant Kc).
    w: dict piece -> weight (missing = symbol)."""
    d = {n: w.get(n, WS[n]) for n in NAMES}
    S = -a*Lv/(rr**2 + a**2)
    A = rr**2 + a**2 - 2*M*rr*d["PP"]
    B = -2*M*rr*(d["EP"]*E + d["PS"]*S)
    C = (-a**2*Lv**2/(rr**2 + a**2) - (E**2 - 1)*rr**2 + Kc
         - 2*M*rr*(d["EE"]*E**2 + d["SS"]*S**2 + 2*d["ES"]*E*S))
    return A, B, C


def twist_drops_out():
    """With L = 0 every r, theta and t quantity is free of the three twist dials; with L != 0 it is not."""
    A, B, C = radial_quadratic({}, Lv=0)
    R0 = sp.expand(B**2 - A*C)
    A1, B1, C1 = radial_quadratic({})
    R1 = sp.expand(B1**2 - A1*C1)
    tw = (WS["SS"], WS["ES"], WS["PS"])
    return all(sp.diff(R0, t) == 0 for t in tw), any(sp.diff(R1, t) != 0 for t in tw)


def polar_constants(w):
    """Spherical polar orbit (L = 0) for dial setting w: E^2, Kc, radial Mino rate^2, modulus, clock term."""
    A, B, C = radial_quadratic(w, Lv=0)
    R = sp.expand((B**2 - A*C).subs(E**2, E2).subs(E, sp.sqrt(E2)))
    sol = sp.solve([R.subs(x, r), sp.diff(R, x).subs(x, r)], [E2, Kc], dict=True)[0]
    E2s, Ks = sp.factor(sol[E2]), sp.factor(sol[Kc])
    Ur2 = sp.factor((-sp.diff(R, x, 2)/2).subs(x, r).subs({E2: E2s, Kc: Ks}))
    m = sp.factor((1 - E2s)*a**2/Ks)
    d = {n: w.get(n, 1) for n in NAMES}
    Ar = A.subs(x, r)
    clock = 2*M*r*(d["EE"] + 2*M*r*d["EP"]**2/Ar)     # dt/dlambda = E [Sigma + clock], with p0 = -B/A on the orbit
    return dict(E2=E2s, Kc=Ks, Ur2=Ur2, m=m, clock=clock)


def ee_alone_exact():
    c = polar_constants({"EE": 1, "PP": 0, "EP": 0})
    return (sp.simplify(c["Ur2"]/c["Kc"] - (1 + a**2/r**2)), sp.simplify(c["m"] - a**2/r**2), c)


def series_rows(w, N=6, U=3):
    c = polar_constants(w)
    sub = {r: 1, M: u, a: k}
    E2v, Kv, Uv, mv, cl = (sp.simplify(c[z].subs(sub)) for z in ("E2", "Kc", "Ur2", "m", "clock"))
    Ks = sp.series(sp.elliptic_k(mm), mm, 0, N//2 + 1).removeO()
    Es = sp.series(sp.elliptic_e(mm), mm, 0, N//2 + 1).removeO()
    KmE = sp.expand(sp.cancel((Ks - Es)/mm))
    Kq = Ks.subs(mm, mv)
    T = 4*sp.sqrt(E2v)/sp.sqrt(Kv)*(Kq + k**2*KmE.subs(mm, mv) + cl*Kq)
    F = sp.series((2*sp.pi/T)**2/u, k, 0, N).removeO()
    rat = sp.series(Uv/Kv*(2*Kq/sp.pi)**2, k, 0, N).removeO()
    rows = {}
    for name, ser in (("omega_theta2", F), ("ratio", rat)):
        rows[name] = [[sp.nsimplify(sp.series(sp.simplify(ser.coeff(k, j)), u, 0, U).removeO().coeff(u, i))
                       for i in range(U)] for j in range(0, N, 2)]
    return rows


def node_rate_split():
    """dphi/dlambda at L = 0 from the shell: only ES and PS contribute; together they give 2 M r a E / Delta."""
    A, B, C = radial_quadratic({}, Lv=L, rr=r)
    p0 = -B/A
    # dphi/dlambda = (1/2) dG/dL with G = A p^2 + 2 B p + C (+ theta part), at the orbit's p = p0, then L -> 0
    G = A*sp.Symbol("p")**2 + 2*B*sp.Symbol("p") + C
    dphi = sp.simplify((sp.diff(G, L)/2).subs(sp.Symbol("p"), p0).subs(L, 0))
    full = sp.simplify(dphi.subs({w_: 1 for w_ in WS.values()}))
    D = r**2 - 2*M*r + a**2
    return dphi, sp.simplify(full - 2*M*r*a*E/D)


def node_rate_exact(Mv, av, rv):
    """Turning rate of the polar orbit's plane in coordinate time (exact), and its up-and-down rate."""
    Mv, av, rv = mp.mpf(Mv), mp.mpf(av), mp.mpf(rv)
    D = rv**2 - 2*Mv*rv + av**2
    den = (av**2 + rv**2)*(Mv*av**2 - 3*Mv*rv**2 + av**2*rv + rv**3)
    E_ = mp.sqrt(rv*D**2/den); Q = Mv*rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv)/den
    m = av**2*((rv**2 + av**2)**2 - 4*Mv*rv**3)/(rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv))
    K, Ee = mp.ellipk(m), mp.ellipe(m)
    Tt = 4*E_/mp.sqrt(Q)*(((rv**2 + av**2)**2/D)*K - av**2*(K - (K - Ee)/m))
    Lam = 4*K/mp.sqrt(Q)
    return (2*Mv*rv*av*E_/D)*Lam/Tt, 2*mp.pi/Tt


def node_from_flight(Mv, av, rv, laps=6):
    """Fly the polar orbit in Kerr-Schild Cartesian coordinates and measure how far its equator crossing turns per lap."""
    import polar as pl
    fn, rfn, _ = pl.ks_hamiltonian()
    D = rv**2 - 2*Mv*rv + av**2
    den = (av**2 + rv**2)*(Mv*av**2 - 3*Mv*rv**2 + av**2*rv + rv**3)
    E_ = np.sqrt(rv*D**2/den); Q = Mv*rv**2*((rv**2 + av**2)**2 - 4*Mv*av**2*rv)/den
    J = np.zeros((4, 4)); J[0, 0] = 1; J[0, 1] = 2*Mv*rv/D
    J[1, 1] = 1 - av*av/D; J[1, 3] = -av; J[2, 1] = rv*av/D; J[2, 3] = rv; J[3, 2] = -rv
    pKS = np.linalg.solve(J.T, np.array([-E_, 0.0, -np.sqrt(Q), 0.0]))

    def flow(tau, s):
        o = fn(s[1], s[2], s[3], *s[4:], Mv, av)
        return [o[0], o[1], o[2], o[3], 0.0, -o[4], -o[5], -o[6]]

    def up(tau, s):
        return s[3]
    up.direction = 1
    nr, Om = node_rate_exact(Mv, av, rv)
    sol = solve_ivp(flow, [0, 1.6*laps*2*np.pi/float(Om)], [0.0, rv, av, 0.0, *pKS], rtol=1e-12, atol=1e-12,
                    events=up)
    ev = sol.y_events[0]
    ang = np.unwrap(np.arctan2(ev[:, 2], ev[:, 1]))
    rate = (ang[-1] - ang[0])/(ev[-1, 0] - ev[0, 0])
    return float(rate), float(nr)


def equator_k2u(w):
    """Coefficient of k^2 u (= chi^2 u^3) and of k^2 in the equatorial orbit factor T for dial setting w (round 16's
    mass shell, EE always on). Odd powers of the twist are removed by averaging prograde and retrograde."""
    import kerr_split as ks
    mp.mp.dps = 60
    T = ks.orbit_factor(ks.ks_shell(w))
    kk = mp.mpf("1e-12")

    def g(uv):
        return (T(uv, kk/uv) + T(uv, -kk/uv) - 2*T(uv, mp.mpf(0)))/(2*kk**2)
    us = [mp.mpf("1e-6")*j for j in (1, 2, 3, 4)]
    c = mp.lu_solve(mp.matrix([[1, v, v**2, v**3] for v in us]), mp.matrix([g(v) for v in us]))
    mp.mp.dps = 30
    return c[0], c[1]


def main():
    out = {}
    eq = {}
    for bits in itertools.product((0, 1), repeat=5):
        w = dict(zip(("PP", "EP", "SS", "ES", "PS"), bits))
        name = "+".join(["EE"] + [n for n, b in w.items() if b])
        c0, c1 = equator_k2u(w)
        eq[name] = [mp.nstr(c0, 15), mp.nstr(c1, 15)]
    out["equator_k2_and_k2u"] = eq
    print("equator", eq, flush=True)
    z, nz = twist_drops_out()
    out["twist_drops_out_at_L0"] = bool(z); out["twist_matters_at_L_nonzero"] = bool(nz)
    d1, d2, c = ee_alone_exact()
    out["ee_alone"] = {"ratio_prefactor_minus_shifted": str(d1), "modulus_minus_shifted": str(d2),
                       "E2": str(c["E2"]), "Kc": str(c["Kc"])}
    tab = {}
    for pp, ep in itertools.product((0, 1), (0, 1)):
        key = "EE" + ("+PP" if pp else "") + ("+EP" if ep else "")
        tab[key] = {kk: [[str(v) for v in row] for row in vv] for kk, vv in series_rows({"EE": 1, "PP": pp, "EP": ep}).items()}
        print(key, tab[key], flush=True)
    out["settings"] = tab
    dphi, chk = node_rate_split()
    out["node_rate_mino"] = str(dphi); out["node_rate_minus_2MraE_over_Delta"] = str(chk)
    fl = [node_from_flight(1.0, 0.9, 6.0), node_from_flight(1.0, 0.5, 8.0)]
    out["node_flights"] = fl
    lt = []
    for rv in (6, 10, 20, 50, 100, 200):
        nr, Om = node_rate_exact(1.0, 0.9, rv)
        lt.append([rv, float(nr), float(2*0.9/mp.mpf(rv)**3)])
    out["node_vs_lense_thirring"] = lt
    print("node", fl, lt, flush=True)
    json.dump(out, open(os.path.join(HERE, "polar_split.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
