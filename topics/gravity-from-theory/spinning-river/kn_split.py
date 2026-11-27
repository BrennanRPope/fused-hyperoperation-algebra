"""DRAFT, round 27 (Nov 27, 2026): charge on the spinning black hole, piece by piece.

Rounds 25 and 26 found that, counted per unit energy, Kerr's orbits are Newton's spinning ring (the river's time part)
plus the push plus the twist, with the time part exactly Newton's pull. Brennan asked whether a charged spinning black
hole's extra terms come from the push alone.

Findings (all exact):
1. Charge enters the river's shell in one place: the field 2Mr (k.p)^2 becomes (2Mr - Q^2)(k.p)^2 (Kerr-Schild form of
   Kerr-Newman; Q is the charge, Carter's constant is called K here). That is mass M replaced by the mass function
   M - Q^2/(2r), the mass less the electric field's energy outside r (round 10's C1). So charge sits in the same slot
   as mass in every piece: time part, push and twist.
2. The time part alone is Newton's ring around that mass function (with M E^2 as before): per unit energy, Newton's
   pull of the hole plus its field's energy, nothing else. Kerr-Newman = Newton's ring (mass function) +
   (2Mr - Q^2)(K - a^2), round 23's one term with 2Mr -> 2Mr - Q^2. The pieces still cancel at second order because
   the field is one square.
3. At the equator, per unit energy, round 26's split carries over with 2Mr -> f = 2Mr - Q^2:
     Newton's ring  e r^4 + f r^2 - (L^2 - e a^2) r^2 + f a^2
     push           -f e r^2 + f L^2 - f a^2 L^2/(r^2+a^2)
     twist          -2 f a L + f a^2 L^2/(r^2+a^2)
     sum            e r^4 + f(1 - e) r^2 - (L^2 - e a^2) r^2 + f (L - a)^2
   i.e. R = e r^4 + 2 alpha r^3 - W r^2 + 2 beta r - gamma with alpha = M(1-e), beta = M x^2, gamma = Q^2 x^2.
   The radial stiffness is alpha r - 3 beta/r + 4 gamma/r^2: round 23's spherical-orbit law with a^2 K -> Q^2 x^2.
4. The rate ratio (Omega_r/Omega_phi)^2 of circular equatorial orbits, with v^2 = M/r, q = a/r, s = Q/r and
   V^2 = v^2 - s^2 (Newton's circular speed around the mass function), is exactly
     [1 - 3q^2 + s^2 (1 + q^2)/V^2]  +  [-6v^2 + 3s^2 - s^4/V^2]  +  [8 q V]
       Newton's ring (time part)          push                        twist
   The first bracket is exactly what the time part alone gives; at a = 0 the first two brackets are exactly
   Reissner-Nordstrom's. So charge's leading term s^2/V^2 = Q^2/(Mr - Q^2) is Newton's pull of the field's energy
   (time part); its relativistic terms +3s^2 - s^4/V^2 are the push's; and the twist's 8qv becomes 8qV. With the push
   or the twist alone the half-squares (+-28 q^2 v^2) appear again and cancel together.
   Equivalently: Kerr's law with v -> V, D_V - 3(V - q)^2, plus (s^2/V^2)(1 + q^2 - s^2 - 3V^2).
5. Independent check: the same ratio from Carter's Kerr-Newman equations in Boyer-Lindquist form, solved numerically
   for the circular orbit at 40 digits, with no river pieces.
Run: python3 kn_split.py  (writes kn_split.json, about a minute)
"""
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
r, M, a, L, p, et, Qe, E, C = sp.symbols("r M a L p e Q E C", real=True)
v, q, s, chi, lam = sp.symbols("v q s chi lambda", positive=True)
NAMES = ("EE", "PP", "SS", "EP", "ES", "PS")
CASES = {"Newton's ring (time part)": ("EE",), "ring + push": ("EE", "PP", "EP"),
         "ring + twist": ("EE", "SS", "ES"), "Kerr-Newman (all)": NAMES}
D = r**2 + a**2


def shell(on, energy=None, carter=None, w=None):
    """The Kerr-Schild shell with field (2Mr - Q^2)(k.p)^2. energy=None: per unit energy at the equator (C = L^2)."""
    d = {n: (1 if n in on else 0) for n in NAMES} if w is None else w
    f = 2*M*r - Qe**2
    S = -a*L/D
    E_ = 1 if energy is None else energy
    Cc = L**2 if carter is None else carter
    ee = et if energy is None else energy**2 - 1
    A = D - f*d["PP"]
    B = -f*(d["EP"]*E_ + d["PS"]*S)
    Cr = -a**2*L**2/D - ee*r**2 + Cc - f*(d["EE"]*E_**2 + d["SS"]*S**2 + 2*d["ES"]*E_*S)
    return A, B, Cr


def R_ring(m2r, e, c, l2):
    """Newton's spinning ring in Mino time, with 2 m r given (a mass function allowed)."""
    return D*(e*r**2 + m2r - c) + a**2*l2


def split():
    zero = lambda x: sp.simplify(x) == 0
    K = C - L**2 + (L - a*E)**2
    RKN = (E*D - a*L)**2 - (r**2 - 2*M*r + a**2 + Qe**2)*(r**2 + K)
    A, B, Cr = shell(NAMES, E, C)
    A0, B0, C0 = shell(("EE",), E, C)
    f = 2*M*r - Qe**2
    # second-order (f^2) part with a dial on each piece: vanishes for any square
    WS = dict(zip(NAMES, sp.symbols("wEE wPP wSS wEP wES wPS")))
    Ag, Bg, Cg = shell((), E, C, w=WS)
    Fs = sp.Symbol("F")
    Rg = sp.expand(sp.cancel((Bg**2 - Ag*Cg).subs(Qe**2, 2*M*r - Fs)))
    f2 = sp.factor(Rg.coeff(Fs, 2))
    al, be, ga = sp.symbols("alpha beta gamma")
    sq = {WS["EE"]: al**2, WS["PP"]: be**2, WS["SS"]: ga**2, WS["EP"]: al*be, WS["ES"]: al*ga, WS["PS"]: be*ga}
    # equator, per unit energy
    ring = et*r**4 + f*r**2 - (L**2 - et*a**2)*r**2 + f*a**2
    push = -f*et*r**2 + f*L**2 - f*a**2*L**2/D
    twist = -2*f*a*L + f*a**2*L**2/D
    Ae, Be, Ce = shell(NAMES)
    Ar, Br, Crr = shell(("EE",))
    Ap, Bp, Cp = shell(("EE", "PP", "EP"))
    At, Bt, Ct = shell(("EE", "SS", "ES"))
    Req = sp.expand(sp.cancel(Be**2 - Ae*Ce))
    x2 = (L - a)**2
    return {
        "shell_with_all_pieces_is_kerr_newman": zero(B**2 - A*Cr - RKN),
        "time_part_is_newton_ring_around_mass_function": zero(B0**2 - A0*C0 - R_ring(2*M*E**2*r - Qe**2*E**2,
                                                                                       E**2 - 1, C, L**2)),
        "kerr_newman_is_ring_plus_one_term": zero(RKN - R_ring(f, E**2 - 1, C, L**2) - f*(K - a**2)),
        "second_order_vanishes_with_all_pieces": zero(f2.subs({WS[n]: 1 for n in NAMES})),
        "second_order_vanishes_for_any_square": zero(f2.subs(sq)),
        "equator_ring_is_time_part": zero(Br**2 - Ar*Crr - ring),
        "equator_ring_plus_push": zero(Bp**2 - Ap*Cp - ring - push),
        "equator_ring_plus_twist": zero(Bt**2 - At*Ct - ring - twist),
        "equator_all_three": zero(Req - ring - push - twist),
        "equator_square_form": zero(Req - (et*r**4 + f*(1 - et)*r**2 - (L**2 - et*a**2)*r**2 + f*x2)),
        "equator_coefficients_alpha_beta_gamma": zero(Req.coeff(r, 3) - 2*M*(1 - et)) and zero(Req.coeff(r, 1)
            - 2*M*x2) and zero(Req.coeff(r, 0) + Qe**2*x2),
    }


def generic_stiffness():
    """R = e r^4 + 2 alpha r^3 - W r^2 + 2 beta r - gamma, circular: -R''/2 = alpha r - 3beta/r + 4gamma/r^2."""
    al, be, ga, W, e = sp.symbols("alpha beta gamma W e")
    R = e*r**4 + 2*al*r**3 - W*r**2 + 2*be*r - ga
    sol = sp.solve([R, sp.diff(R, r)], [e, W], dict=True)[0]
    w2 = sp.simplify((-sp.diff(R, r, 2)/2).subs(sol))
    return {"radial": sp.simplify(w2 - (al*r - 3*be/r + 4*ga/r**2)) == 0,
            "W": sp.simplify(sol[W] - (al*r + 3*be/r - 2*ga/r**2)) == 0}


V = sp.sqrt(v**2 - s**2)
LAW = {"ring": 1 - 3*q**2 + s**2*(1 + q**2)/V**2, "push": -6*v**2 + 3*s**2 - s**4/V**2, "twist": 8*q*V}


def ratios():
    sub = {M: v**2, a: q, Qe: s}
    out, res = {}, {}
    for name, on in CASES.items():
        A, B, Cr = shell(on)
        Rt = sp.together(B**2 - A*Cr)
        Rn, Rd = sp.numer(Rt), sp.denom(Rt)
        ets = sp.solve(Rn, et)[0]
        dR = sp.numer(sp.together(sp.diff(Rn, r).subs(et, ets).subs(r, 1)))
        pick = None
        for R0 in sp.solve(dR.subs(sub), L):                     # the prograde circular orbit: L ~ +v
            lead = sp.series(R0.subs({q: chi*v**2, s: lam*v**2}), v, 0, 2).removeO()
            if sp.simplify(lead - v) == 0:
                pick = R0
        turn = sp.Rational(1, 2)*sp.diff(2*B*p + Cr, L).subs(p, -B/A)
        wr2 = (-sp.diff(Rn, r, 2)/(2*Rd)).subs(et, ets).subs(sub).subs(r, 1)
        ratio = (wr2/turn.subs(sub).subs(r, 1)**2).subs(L, pick)
        ser = sp.expand(sp.series(ratio.subs({q: chi*v**2, s: lam*v**2}), v, 0, 9).removeO())
        out[name] = ratio
        res[name] = {"L_circular": str(pick), "series_q=chi v^2, s=lambda v^2": str(sp.collect(ser, v))}
        print(name, sp.collect(ser, v), flush=True)
    # exactness checks at 40 digits on random points (closed forms too large for simplify)
    mp.mp.dps = 40
    rng = np.random.default_rng(27)
    pts = [tuple(mp.mpf(x) for x in (rng.uniform(.2, .45), rng.uniform(0, .3), rng.uniform(0, .15)))
           for _ in range(12)]
    fn = {n: sp.lambdify((v, q, s), out[n], "mpmath") for n in CASES}
    law = {k: sp.lambdify((v, q, s), x, "mpmath") for k, x in LAW.items()}
    rn = sp.lambdify((v, s), (v**2 - 6*v**4 + 9*s**2*v**2 - 4*s**4)/(v**2 - s**2), "mpmath")
    err = lambda xs: float(max(abs(x) for x in xs))
    res["checks"] = {
        "ring_exact": err([fn["Newton's ring (time part)"](*t) - law["ring"](*t) for t in pts]) < 1e-30,
        "kerr_newman_is_ring_plus_push_plus_twist": err([fn["Kerr-Newman (all)"](*t) - law["ring"](*t)
                                                         - law["push"](*t) - law["twist"](*t) for t in pts]) < 1e-30,
        "reissner_nordstrom_at_a0": err([fn["Kerr-Newman (all)"](t[0], 0, t[2]) - rn(t[0], t[2]) for t in pts]) < 1e-30
        and err([fn["ring + push"](t[0], 0, t[2]) - rn(t[0], t[2]) for t in pts]) < 1e-30,
        "kerr_form_plus_charge_term": sp.simplify(sum(LAW.values()) - ((1 - 3*V**2 + 2*q*V) - 3*(V - q)**2
                                                  + s**2/V**2*(1 + q**2 - s**2 - 3*V**2))) == 0,
    }
    # partial cases: the half-squares
    ser_p = sp.expand(sp.series(out["ring + push"].subs({q: chi*v**2, s: lam*v**2}), v, 0, 7).removeO())
    ser_t = sp.expand(sp.series(out["ring + twist"].subs({q: chi*v**2, s: lam*v**2}), v, 0, 7).removeO())
    res["checks"]["half_squares_pm28"] = (sp.expand(ser_p.coeff(v, 6)).coeff(chi, 2) - 28).subs(lam, 0) == 0 and \
        (sp.expand(ser_t.coeff(v, 6)).coeff(chi, 2) + 28).subs(lam, 0) == 0
    res["exact_law"] = {k: str(x) for k, x in LAW.items()}
    # curves for the figure: a = 0.5M, Q = 0.6M, r from 4M to 20M; and Kerr (Q = 0) for comparison
    rs = np.linspace(4.0, 20.0, 160)
    av, Qv = 0.5, 0.6
    cur = {}
    for n in CASES:
        g = sp.lambdify((v, q, s), out[n], "mpmath")
        cur[n] = [float(mp.re(g(mp.sqrt(1/x), av/x, Qv/x))) for x in rs]
    g = sp.lambdify((v, q, s), out["Kerr-Newman (all)"], "mpmath")
    cur["Kerr (Q = 0)"] = [float(mp.re(g(mp.sqrt(1/x), av/x, 0))) for x in rs]
    res["curves"] = {"r": rs.tolist(), "a": av, "Q": Qv, **cur}
    return res


def carter_numeric(av=0.5, Qv=0.6, radii=(5.0, 8.0, 13.0)):
    """Kerr-Newman's own Carter equations (Boyer-Lindquist), no river pieces: circular prograde orbit found by
    Newton's method on R = R' = 0, then (Omega_r/Omega_phi)^2 = (-R''/2)/(dphi/dlambda)^2."""
    mp.mp.dps = 40
    out = []
    for rv in radii:
        rv = mp.mpf(rv)
        Dl = lambda x: x**2 - 2*x + av**2 + Qv**2
        R = lambda x, E_, L_: (E_*(x**2 + av**2) - av*L_)**2 - Dl(x)*(x**2 + (L_ - av*E_)**2)
        F = lambda E_, L_: [R(rv, E_, L_), mp.diff(lambda x: R(x, E_, L_), rv)]
        E0, L0 = mp.findroot(F, (mp.mpf("0.95"), mp.sqrt(rv)))
        wr2 = -mp.diff(lambda x: R(x, E0, L0), rv, 2)/2
        phid = -(av*E0 - L0) + av*(E0*(rv**2 + av**2) - av*L0)/Dl(rv)
        ratio = wr2/phid**2
        vv, qq, ss = mp.sqrt(1/rv), av/rv, Qv/rv
        law = sum(sp.lambdify((v, q, s), x, "mpmath")(vv, qq, ss) for x in LAW.values())
        out.append({"r": float(rv), "E": float(E0), "L": float(L0), "ratio_carter": float(ratio),
                    "ratio_three_pieces": float(law), "diff": float(abs(ratio - law))})
    return out


def main():
    res = {"split": split(), "generic_stiffness": generic_stiffness()}
    print(res, flush=True)
    res["ratios"] = ratios()
    print(res["ratios"]["checks"], flush=True)
    res["carter_numeric"] = carter_numeric()
    print(res["carter_numeric"], flush=True)
    json.dump(res, open(os.path.join(HERE, "kn_split.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
