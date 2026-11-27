"""DRAFT, round 26 (Nov 27, 2026): the equator's push-twist cancellation, piece by piece.

Round 25 found that per unit energy Kerr's orbits are Newton's spinning ring plus the river's push plus its twist, and
that over the poles (no twist) the push alone makes every relativistic change to the orbit's shape. In the equator the
twist is there, and round 20 found the push's +28 (a/r)^2 (M/r) cancelled by the twist's 28. Is that the same rule?

Findings (exact unless marked):
1. Same rule. Per unit energy (E = 1, mass-shell constant e = 1 - 1/E^2) the equatorial radial function splits as
     Newton's ring   e r^4 + 2M r^3 - (L^2 - e a^2) r^2 + 2M r a^2
     push            -2M e r^3 + 2M r L^2 - 2M r a^2 L^2/(r^2 + a^2)
     twist           -4M r a L + 2M r a^2 L^2/(r^2 + a^2)
   and the sum is e r^4 + 2M(1 - e) r^3 - (L^2 - e a^2) r^2 + 2M r (L - a)^2: the a^2 L^2/(r^2+a^2) parts cancel, and
   the push's L^2, the twist's -2aL and the ring's a^2 complete the square (L - a)^2 = x^2, round 22's coupling.
2. For any R = e r^4 + 2 alpha r^3 - W r^2 + 2 beta r on a circular orbit, the stiffnesses are exactly
   omega_r^2 = alpha r - 3 beta/r and W = alpha r + 3 beta/r. Kerr per unit energy has alpha = M(1 - e) (Newton's M,
   less the push's M e) and beta = M x^2 (ring a^2 + push L^2 + twist -2aL). This is round 22's law
   D -+ 3(v - q)^2 read piece by piece: in (Omega_r/Omega_phi)^2 = 1 - 6v^2 + 8qv - 3q^2 the -3q^2 is Newton's ring,
   -6v^2 is the push, +8qv is the twist (v^2 = M/r, q = a/r).
3. With only one of push or twist the square is broken and R is not even a polynomial (the a^2 L^2/(r^2+a^2) piece
   survives). Ring + push gives 1 - 6u - 3q^2 + 28 q^2 u - 24 q^2 u^2 + ...; ring + twist gives
   1 + 8qv - 3q^2 - 28 q^2 u + 16 q^3 v^3 + ... (u = v^2). The +28 and -28 are the half-squares; together they leave
   exactly D - 3(v - q)^2 at every strength. Round 20's +28 = 28 is the square closing.
Run: python3 equator_split.py  (writes equator_split.json, about 10 seconds)
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
r, M, a, L, p, et = sp.symbols("r M a L p e", real=True)
v, q, chi = sp.symbols("v q chi", positive=True)
NAMES = ("EE", "PP", "SS", "EP", "ES", "PS")
CASES = {"Newton's ring (time part)": ("EE",), "ring + push": ("EE", "PP", "EP"),
         "ring + twist": ("EE", "SS", "ES"), "Kerr (all)": NAMES}
LSOL = {"Newton's ring (time part)": v*(q**2 + 1),
        "ring + push": v*(q**2 + 1)*sp.sqrt(q**4 - 4*q**2*v**2 + 2*q**2 + 1)/(q**2 - 2*v**2 + 1),
        "ring + twist": v*(q**2 + 1)*(q*v*(q**2 + 3) - sp.sqrt(q**4 + 4*q**2*v**2 + 2*q**2 + 1))
        / (q**4*v**2 + 5*q**2*v**2 - q**2 - 1),
        "Kerr (all)": v*(q**2 - 2*q*v + 1)/(q*v - 2*v**2 + 1)}


def shell(on):
    """Per unit energy (E = 1), equator (Q = 0, so Carter's C = L^2)."""
    d = {n: (1 if n in on else 0) for n in NAMES}
    D = r**2 + a**2
    S = -a*L/D
    A = D - 2*M*r*d["PP"]
    B = -2*M*r*(d["EP"] + d["PS"]*S)
    Cr = -a**2*L**2/D - et*r**2 + L**2 - 2*M*r*(d["EE"] + d["SS"]*S**2 + 2*d["ES"]*S)
    return A, B, Cr


def split():
    D = r**2 + a**2
    ring = et*r**4 + 2*M*r**3 - (L**2 - et*a**2)*r**2 + 2*M*r*a**2
    push = -2*M*et*r**3 + 2*M*r*L**2 - 2*M*r*a**2*L**2/D
    twist = -4*M*r*a*L + 2*M*r*a**2*L**2/D
    A, B, Cr = shell(NAMES)
    Rk = B**2 - A*Cr
    A0, B0, C0 = shell(("EE",))
    Ap, Bp, Cp = shell(("EE", "PP", "EP"))
    At, Bt, Ct = shell(("EE", "SS", "ES"))
    zero = lambda x: sp.simplify(x) == 0
    return {
        "ring_is_time_part": zero(B0**2 - A0*C0 - ring),
        "ring_plus_push_setting": zero(Bp**2 - Ap*Cp - ring - push),
        "ring_plus_twist_setting": zero(Bt**2 - At*Ct - ring - twist),
        "all_three_is_kerr": zero(Rk - ring - push - twist),
        "kerr_per_unit_energy_square": zero(Rk - (et*r**4 + 2*M*(1 - et)*r**3 - (L**2 - et*a**2)*r**2
                                                 + 2*M*r*(L - a)**2)),
        "square_from_three_sources": zero(a**2 + L**2 - 2*a*L - (L - a)**2),
        "push_alone_not_polynomial": sp.denom(sp.together(ring + push)) != 1,
        "twist_alone_not_polynomial": sp.denom(sp.together(ring + twist)) != 1,
    }


def generic_stiffness():
    """R = e r^4 + 2 alpha r^3 - W r^2 + 2 beta r, circular: omega_r^2 = alpha r - 3 beta/r, W = alpha r + 3 beta/r."""
    e, al, W, be = sp.symbols("e alpha W beta")
    R = e*r**4 + 2*al*r**3 - W*r**2 + 2*be*r
    s = sp.solve([R/r, sp.diff(R, r)], [e, W], dict=True)[0]
    w2 = sp.simplify((-sp.diff(R, r, 2)/2).subs(s))
    return {"omega_r2": sp.simplify(w2 - (al*r - 3*be/r)) == 0,
            "W": sp.simplify(s[W] - (al*r + 3*be/r)) == 0}


def ratios():
    out = {}
    for name, on in CASES.items():
        A, B, Cr = shell(on)
        Rt = sp.together(B**2 - A*Cr)
        Rn, Rd = sp.numer(Rt), sp.denom(Rt)
        ets = sp.solve(Rn, et)[0]
        turn = sp.Rational(1, 2)*sp.diff(2*B*p + Cr, L).subs(p, -B/A)       # whole turning rate (Cr holds C = L^2)
        sub = {M: v**2, a: q}
        wr2 = (-sp.diff(Rn, r, 2)/(2*Rd)).subs(et, ets).subs(sub).subs(r, 1)
        # the circular orbit really is one: Rn' = 0 there
        circ = sp.simplify(sp.diff(Rn, r).subs(et, ets).subs(sub).subs(r, 1).subs(L, LSOL[name]))
        ratio = sp.simplify((wr2/turn.subs(sub).subs(r, 1)**2).subs(L, LSOL[name]))
        ser = sp.expand(sp.series(ratio.subs(q, chi*v**2), v, 0, 11).removeO())
        out[name] = {"ratio": ratio, "series_chi": ser, "circular_ok": circ == 0}
        print(name, ser, flush=True)
    Dw = 1 - 3*v**2 + 2*q*v
    res = {n: {"series_in_chi_v": str(d["series_chi"]), "circular_ok": d["circular_ok"]} for n, d in out.items()}
    res["checks"] = {
        "ring_exact_1_minus_3q2": sp.simplify(out["Newton's ring (time part)"]["ratio"] - (1 - 3*q**2)) == 0,
        "kerr_exact_D_minus_3(v-q)^2": sp.simplify(out["Kerr (all)"]["ratio"] - (Dw - 3*(v - q)**2)) == 0,
        "push_series_head": sp.expand(out["ring + push"]["series_chi"]
                                      - (1 - 6*v**2 - 3*chi**2*v**4 + 28*chi**2*v**6 - 24*chi**2*v**8)).coeff(v, 6) == 0
        and out["ring + push"]["series_chi"].coeff(v, 2) == -6 and out["ring + push"]["series_chi"].coeff(v, 3) == 0
        and out["ring + push"]["series_chi"].coeff(v, 6) == 28*chi**2,
        "twist_series_head": out["ring + twist"]["series_chi"].coeff(v, 3) == 8*chi
        and out["ring + twist"]["series_chi"].coeff(v, 2) == 0 and out["ring + twist"]["series_chi"].coeff(v, 6) == -28*chi**2,
    }
    # curves for the figure: a = 0.9 M, r from 3.5M to 20M
    curves = {}
    rs = np.linspace(3.5, 20, 160)
    for name in CASES:
        f = sp.lambdify((v, q), out[name]["ratio"], "numpy")
        curves[name] = [float(np.real(f(np.sqrt(1/x), 0.9/x))) for x in rs]
    res["curves"] = {"r": rs.tolist(), "a": 0.9, **curves}
    return res


def main():
    out = {"split": split(), "generic_stiffness": generic_stiffness()}
    print(out, flush=True)
    out["ratios"] = ratios()
    print(out["ratios"]["checks"], flush=True)
    json.dump(out, open(os.path.join(HERE, "equator_split.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
