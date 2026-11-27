"""DRAFT, round 22 (Nov 27, 2026): why the push's +28 equals the twist's 28 in the equator.

Round 20 split the (a/r)^2 (M/r) term of the equatorial orbit factor T = (Omega_r/Omega_phi)^2 by Kerr-Schild piece:
the push brings +20 + 8, the twist takes -16 - 12. Round 21 showed it isn't the one jet (P.S adds 0 there). This round
finds the reason: one exact law for both wobbles of a circular equatorial orbit, in which spin squared can enter only
one way.

Findings (v = sqrt(M/r), q = a/r, prograde q > 0; x = L - aE, the orbit's angular momentum about the spin):
1. Carter's radial function R = [E(r^2+a^2) - aL]^2 - Delta (r^2 + x^2) gives, in Mino time,
       (dr/dlambda)^2 = r^4 (E^2 - 1) + 2 M r^3 - K0 r^2 + 2 M r x^2,   K0 = L^2 + a^2 (1 - E^2).
   The mass meets the orbit twice: 2 M r^3 (Newton) and 2 M r x^2 (the relativistic term). The second is one power of r,
   so it adds nothing to the radial stiffness. On every circular orbit the Mino-time stiffnesses are
       radial  omega_r^2 = M (r - 3 x^2/r),   up-and-down  omega_theta^2 = K0 = M (r + 3 x^2/r),
   so they always add to 2 M r, with no spin in it.
2. On the orbit the river's jet energy is exact and simple: k.p = E + P + S (energy + push + twist) = 1/sqrt(D),
   D = 1 - 3v^2 + 2qv. The orbit's Mino turning rate is r v k.p, and x = r (v - q) k.p.
3. Dividing: T_r = D - 3 (v - q)^2 and T_theta = D + 3 (v - q)^2. These equal Kerr's known ratios
   1 - 6v^2 + 8qv - 3q^2 and 1 - 4qv + 3q^2 (Okazaki, Kato & Fukue 1987; Kato 1990).
4. So spin squared enters only through (v - q)^2, the angular momentum about the spin, squared; D has the spin only to
   the first power. The q^2 term is therefore exactly -3 (radial) and +3 (up-and-down) at every strength: no
   (a/r)^2 (M/r)^n term for any n >= 1. Round 20's +28 = 28 is the n = 1 case of this.
5. The piece split is bookkeeping, not two effects. The same jet energy on the exact orbit, cut by piece, gives
   a^2 M parts E alone -1, push of E +4, twist -3, push of the twist 0 (exact functions below); round 20's
   cut of the field gives +28 and -28. Different cuts, different numbers, same zero.
6. Checked in flight: a slightly eccentric orbit integrated in the exact Kerr-Schild field turns per lap by
   2 pi (1/sqrt(D - 3(v - q)^2) - 1), prograde and retrograde.
Run: python3 wobble_law.py  (writes wobble_law.json, about 1 minute; peak memory well under 1 GB)
"""
import json
import os

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
v, q, r = sp.symbols("v q r", positive=True)
eps = sp.symbols("epsilon", positive=True)
D = 1 - 3*v**2 + 2*q*v
E_C = (1 - 2*v**2 + q*v)/sp.sqrt(D)                 # Bardeen, Press & Teukolsky 1972, in v and q
L_C = r*v*(1 - 2*q*v + q**2)/sp.sqrt(D)
X_C = L_C - q*r*E_C                                  # x = L - aE
T_R_KNOWN = 1 - 6*v**2 + 8*q*v - 3*q**2
T_TH_KNOWN = 1 - 4*q*v + 3*q**2


def carter():
    """1: the circular orbit solves R = R' = 0; Mino stiffnesses M(r -+ 3x^2/r); their sum 2Mr."""
    rr, E, L, M, a = sp.symbols("rr E L M a")
    R = (E*(rr**2 + a**2) - a*L)**2 - (rr**2 - 2*M*rr + a**2)*(rr**2 + (L - a*E)**2)
    R_mino = rr**4*(E**2 - 1) + 2*M*rr**3 - (L**2 + a**2*(1 - E**2))*rr**2 + 2*M*rr*(L - a*E)**2
    s = {E: E_C, L: L_C, rr: r, M: v**2*r, a: q*r}
    x2 = (L - a*E)**2
    wr2 = -sp.diff(R, rr, 2)/2; wth2 = L**2 + a**2*(1 - E**2)
    return {
        "R_equals_mino_form": sp.simplify(sp.expand(R - R_mino)) == 0,
        "R_on_orbit": sp.simplify(R.subs(s)) == 0,
        "dR_on_orbit": sp.simplify(sp.diff(R, rr).subs(s)) == 0,
        "mass_x2_term_stiffness": sp.diff(2*M*rr*(L - a*E)**2, rr, 2) == 0,
        "radial_stiffness": sp.simplify((wr2 - M*(rr - 3*x2/rr)).subs(s)) == 0,
        "vertical_stiffness": sp.simplify((wth2 - M*(rr + 3*x2/rr)).subs(s)) == 0,
        "sum_2Mr": sp.simplify((wr2 + wth2 - 2*M*rr).subs(s)) == 0,
    }


def jet_on_orbit():
    """2: E + push + twist = 1/sqrt(D); Mino turning rate r v k.p; x = r (v - q) k.p."""
    a, M = q*r, v**2*r
    Delta = r**2 - 2*M*r + a**2
    S = -a*L_C/(r**2 + a**2); push = 2*M*r*(E_C + S)/Delta
    kp = E_C + push + S
    phi_mino = X_C + a*(E_C*(r**2 + a**2) - a*L_C)/Delta
    return {
        "kp_times_sqrtD": sp.simplify(kp*sp.sqrt(D)),
        "kp_equals_P_over_Delta": sp.simplify(kp - (E_C*(r**2 + a**2) - a*L_C)/Delta) == 0,
        "turning_rate_over_r_v_kp": sp.simplify(phi_mino/(r*v*kp)),
        "x_over_r_kp": sp.simplify(X_C/(r*kp)),
    }


def law():
    """3: T_r = D - 3(v-q)^2, T_theta = D + 3(v-q)^2, from the stiffnesses over the turning rate squared."""
    M = v**2*r
    turn2 = (r*v/sp.sqrt(D))**2
    Tr = M*(r - 3*X_C**2/r)/turn2; Tth = M*(r + 3*X_C**2/r)/turn2
    return {
        "T_r_from_stiffness_minus_law": sp.simplify(Tr - (D - 3*(v - q)**2)),
        "T_theta_from_stiffness_minus_law": sp.simplify(Tth - (D + 3*(v - q)**2)),
        "T_r_law_minus_known": sp.expand(D - 3*(v - q)**2 - T_R_KNOWN),
        "T_theta_law_minus_known": sp.expand(D + 3*(v - q)**2 - T_TH_KNOWN),
        "q2_coefficient_radial": sp.Poly(sp.expand(D - 3*(v - q)**2), q).coeff_monomial(q**2),
        "q2_coefficient_vertical": sp.Poly(sp.expand(D + 3*(v - q)**2), q).coeff_monomial(q**2),
        "D_q2_coefficient": sp.Poly(D, q).coeff_monomial(q**2),
    }


def jet_pieces():
    """5: spin-squared part of 1/(k.p)^2 on the exact orbit, built up piece by piece (r = 1, a = eps)."""
    qq = eps
    Dq = 1 - 3*v**2 + 2*qq*v
    Ek = (1 - 2*v**2 + qq*v)/sp.sqrt(Dq); Lk = v*(1 - 2*qq*v + qq**2)/sp.sqrt(Dq)
    Delta = 1 - 2*v**2 + eps**2
    S = -eps*Lk/(1 + eps**2); pE = 2*v**2*Ek/Delta; pS = 2*v**2*S/Delta
    parts = {"E alone": Ek, "E + push of E": Ek + pE, "E + twist": Ek + S,
             "E + push of E + twist": Ek + pE + S, "all (the one jet)": Ek + pE + S + pS}
    out = {}
    for name, e in parts.items():
        c = sp.factor(sp.simplify(sp.diff(e**-2, eps, 2).subs(eps, 0)/2))
        lin = sp.series(c, v, 0, 3).removeO().coeff(v, 2)
        out[name] = {"exact": str(c), "a2M": str(lin)}
    return out


def flight(u, chi, laps=10):
    """6: periapsis advance per lap from Hamilton's equations of the exact Kerr-Schild field, vs the law."""
    import kerr_split as ks
    sim, _ = ks.simulate({}, u, chi, laps=laps)
    vv = mp.sqrt(u); qq = chi*u
    law_T = (1 - 3*vv**2 + 2*qq*vv) - 3*(vv - qq)**2
    return sim, float(2*mp.pi*(1/mp.sqrt(law_T) - 1))


def main():
    out = {"carter": carter(), "jet_on_orbit": {k: str(x) for k, x in jet_on_orbit().items()},
           "law": {k: str(x) for k, x in law().items()}, "jet_pieces": jet_pieces()}
    out["carter"] = {k: bool(x) for k, x in out["carter"].items()}
    out["flight"] = {f"u={u}, chi={c}": flight(u, c) for u, c in ((0.1, 0.9), (0.1, -0.9))}
    mp.mp.dps = 30
    grid = []
    for vv in (0.1, 0.2, 0.3, 0.4):
        grid.append(float(sp.Poly(sp.expand((D - 3*(v - q)**2).subs(v, sp.Rational(str(vv)))), q).coeff_monomial(q**2)))
    out["q2_term_at_v_0.1_to_0.4"] = grid
    with open(os.path.join(HERE, "wobble_law.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
