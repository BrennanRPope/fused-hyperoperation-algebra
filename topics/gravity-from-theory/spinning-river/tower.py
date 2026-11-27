"""DRAFT, round 17 (Nov 27, 2026): why the shifted source's whole tower of shapes leaves only one -3.

The source shifted by i*a (Part 4) is not a point: seen from far away it has a whole tower of shapes, a mass
quadrupole, hexadecapole, ... (Kerr's "no-hair" tower M (ia)^l). Yet Kerr's exact orbit factor
T = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 has a single shape term. Why?

1. Newton first. In the equator the shifted source's pull is Phi = -M/sqrt(R^2 - a^2) in flat radius R, i.e.
   exactly -M/r in the shift's own radius r = sqrt(R^2 - a^2) (the real distance to the shifted source). Then
   Kepler's law holds exactly in r (Omega^2 = M/r^3) and T = 1 - 3a^2/r^2 exactly: the tower sums to one term.
   In flat radius the same T is 1 - 3(a/R)^2 - 3(a/R)^4 - ..., every shape contributing.
2. The staircase: each higher shape of the shifted source cancels exactly the leftover of the ones below it.
   A body with Kerr's quadrupole but other higher shapes (quadrupole only; a ring of radius sqrt(2) a) leaves an
   endless series.
3. Kepler in r is the whole story: T = 4 + dln(Omega^2)/dln R, so T = 1 - 3a^2/r^2 for every orbit holds exactly
   when Omega^2 = M/r^3, and that fixes the shifted source.
4. Gravity: Kerr's radial motion in the same radius is Carter's quartic R(r), whose only spin terms are a^2 (shape)
   and (L - aE) (twist), and Kerr's orbital period is Kepler's in r plus the clock gap: 1/Omega = r^(3/2)/sqrt(M) + a.
   The exact T follows from the quartic. Written against the measured period with the clock gap taken off,
   T is a finite polynomial in a measurable frequency.
Run: python3 tower.py  (writes tower.json, about 10 seconds)
"""
import json
import os

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
R, a, M, r, x = sp.symbols("R a M r x", positive=True)
NMAX = 5                                     # shapes up to (a/R)^10


def T_newton(Phi):
    """Orbit factor (Omega_r/Omega_phi)^2 of near-circular orbits in an equatorial potential Phi(R)."""
    return sp.simplify(3 + R*sp.diff(Phi, R, 2)/sp.diff(Phi, R))


def omega2(Phi):
    return sp.simplify(sp.diff(Phi, R)/R)


SHIFTED = -M/sp.sqrt(R**2 - a**2)            # Re M/(distance to the point i*a on the axis), in the equator


def shifted_tower(n_max=NMAX):
    """Equatorial coefficients of the shifted source's shapes: -M/sqrt(R^2-a^2) = -sum c_n M a^(2n)/R^(2n+1)."""
    return [sp.binomial(2*n, n)/sp.Integer(4)**n for n in range(n_max + 1)]


def ring_tower(radius, n_max=NMAX):
    """Equatorial coefficients of a thin ring of mass M and the given radius: P_2n(0)^2 radius^(2n)."""
    return [sp.legendre(2*n, 0)**2*(radius/a)**(2*n) for n in range(n_max + 1)]


def potential(coeffs):
    return -sum(c*M*a**(2*n)/R**(2*n + 1) for n, c in enumerate(coeffs))


def in_r(expr, order=2*NMAX + 1):
    """Series in a/r of an expression in R, with R = sqrt(r^2 + a^2)."""
    return sp.expand(sp.series(sp.simplify(expr.subs(R, sp.sqrt(r**2 + a**2))), a, 0, order).removeO())


def coeff_list(ser, var):
    return [sp.nsimplify(sp.expand(ser).coeff(a, 2*n)*var**(2*n)) for n in range(1, NMAX + 1)]  # var: r or R


def in_x(Phi, order=2*NMAX + 1):
    """T as a series in x = (M Omega)^(2/3) (measurable), with chi = a/M: returns coefficients of chi^(2n) x^(2n)."""
    w2 = omega2(Phi); T = T_newton(Phi)
    # x^3 = M^2 Omega^2 at radius R; invert R(x) as a series: start from R = M/x and iterate
    chi = sp.symbols("chi", positive=True)
    Tsub = sp.simplify(T.subs(a, chi*M)); wsub = sp.simplify(w2.subs(a, chi*M))
    Rx = M/x
    for _ in range(NMAX + 1):
        # x^3 = M^2 w2(R)  ->  R = M/x * (M w2 R^3 / M^2 ... ) fixed point
        corr = sp.series(sp.simplify((M**2*wsub*R**3/M**3)).subs(R, Rx), chi, 0, order).removeO()
        Rx = sp.expand(M/x*corr**sp.Rational(1, 3))
        Rx = sp.series(Rx, chi, 0, order).removeO()
    Tx = sp.expand(sp.series(Tsub.subs(R, Rx), chi, 0, order).removeO())
    return [sp.nsimplify(sp.simplify(Tx.coeff(chi, 2*n)/x**(2*n))) for n in range(1, NMAX + 1)], Tx


def carter():
    """Kerr in the shift's radius: Carter's radial quartic and the exact orbit factor from it."""
    E, L = sp.symbols("E L", positive=True)
    D = r**2 - 2*M*r + a**2; P = E*(r**2 + a**2) - a*L
    Rq = sp.expand(P**2 - D*(r**2 + (L - a*E)**2))
    quartic = (E**2 - 1)*r**4 + 2*M*r**3 + (a**2*(E**2 - 1) - L**2)*r**2 + 2*M*(L - a*E)**2*r
    s, chi = sp.symbols("s chi", positive=True)
    den = sp.sqrt(1 - 3*s**2 + 2*chi*s**3)
    sub = {r: M/s**2, a: chi*M, E: (1 - 2*s**2 + chi*s**3)/den,
           L: sp.sqrt(M*M/s**2)*(1 - 2*chi*s**3 + chi**2*s**4)/den}        # Bardeen, Press & Teukolsky (1972)
    circ = [sp.simplify(Rq.subs(sub)), sp.simplify(sp.diff(Rq, r).subs(sub))]
    phid = -(a*E - L) + a*P/D                                              # d phi / d(Mino time)
    T = sp.simplify((-sp.diff(Rq, r, 2)/2/phid**2).subs(sub))             # (omega_r/omega_phi)^2
    # period: dt/dlambda = a(L - aE) + (r^2+a^2) P/Delta ; 1/Omega = (dt/dl)/(dphi/dl)
    tdot = a*(L - a*E) + (r**2 + a**2)*P/D
    inv_om = sp.simplify((tdot/phid).subs(sub))
    return {"quartic_ok": sp.simplify(Rq - quartic) == 0, "circular_ok": all(c == 0 for c in circ),
            "T": sp.expand(T), "T_minus_formula": sp.simplify(T - (1 - 6*s**2 + 8*chi*s**3 - 3*chi**2*s**4)),
            "period_minus_kepler_plus_gap": sp.simplify(inv_om - (M/s**3 + chi*M))}


def main():
    out = {}
    Ts = T_newton(SHIFTED)
    out["shifted_T_flat"] = str(sp.factor(Ts))
    out["shifted_T_r"] = str(sp.simplify(Ts.subs(R, sp.sqrt(r**2 + a**2))))
    out["shifted_omega2_r"] = str(sp.simplify(omega2(SHIFTED).subs(R, sp.sqrt(r**2 + a**2))))
    out["shifted_T_flat_series"] = [str(c) for c in coeff_list(sp.series(Ts, a, 0, 2*NMAX + 1).removeO(), R)]
    # each shape alone (leading order), and the staircase: the tower cut after shape n
    c = shifted_tower()
    alone = []
    for n in range(1, NMAX + 1):
        Tn = T_newton(-M/R - c[n]*M*a**(2*n)/R**(2*n + 1))
        alone.append(str(sp.nsimplify(sp.series(Tn, a, 0, 2*n + 1).removeO().coeff(a, 2*n)*R**(2*n))))
    out["shape_alone_leading"] = alone
    stair = {}
    for n in range(1, NMAX + 1):
        stair[f"up to (a/R)^{2*n}"] = [str(v) for v in coeff_list(in_r(T_newton(potential(c[:n + 1]))), r)]
        print("staircase", n, stair[f"up to (a/R)^{2*n}"], flush=True)
    out["staircase_in_r"] = stair
    # rivals with Kerr's quadrupole
    rivals = {"quadrupole only": [1, sp.Rational(1, 2)] + [0]*(NMAX - 1),
              "ring of radius sqrt(2) a": ring_tower(sp.sqrt(2)*a)}
    out["rivals_in_r"] = {k: [str(v) for v in coeff_list(in_r(T_newton(potential(v))), r)] for k, v in rivals.items()}
    print("rivals", out["rivals_in_r"], flush=True)
    # measurable version: against orbital frequency x = (M Omega)^(2/3)
    xs = {}
    for k, Phi in (("shifted source", SHIFTED), ("quadrupole only", potential(rivals["quadrupole only"])),
                   ("ring of radius sqrt(2) a", potential(rivals["ring of radius sqrt(2) a"]))):
        xs[k] = [str(v) for v in in_x(Phi)[0]]
        print("in x", k, xs[k], flush=True)
    out["in_frequency"] = xs
    # Kepler in r fixes the source: Omega^2 = M/r^3  ->  Phi' = R M/(R^2 - a^2)^(3/2)  ->  the shifted source
    Phi_from_kepler = sp.integrate(R*M/(R**2 - a**2)**sp.Rational(3, 2), R)
    out["kepler_fixes_source"] = str(sp.simplify(Phi_from_kepler - SHIFTED))
    out["T_from_omega"] = str(sp.simplify(T_newton(SHIFTED) - (4 + R*sp.diff(sp.log(omega2(SHIFTED)), R))))
    ca = carter()
    out["carter"] = {k: str(v) for k, v in ca.items()}
    print("carter", out["carter"], flush=True)
    json.dump(out, open(os.path.join(HERE, "tower.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
