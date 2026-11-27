"""DRAFT, round 21 (Nov 27, 2026): the river's push and twist as one jet, or as two velocities?

Brennan's jet rule: wherever velocities combine, use the jet. Each river cell's map pushes its radial leg along the
time leg (the push, speed v) and along the azimuth leg (the twist, -v*alpha). Both are shears into the same leg, so
as jets they compose into ONE shear, I + aN_P + bN_S (N_P N_S = N_S N_P = 0), and the field is that one ray squared:
    (k.p)^2 = (E + P + S)^2 = EE + PP + SS + 2EP + 2ES + 2PS.
Adding the push and the twist as two separate velocities instead (superposing their fields) gives
    (E + P)^2 + (E + S)^2 - E^2,
the same except for the cross piece 2PS. So the cross piece P.S is exactly what the jet rule adds.

Findings:
1. Field equations: one jet is Kerr (Ricci ~ 1e-50 at 50 digits); two velocities are not vacuum (Ricci ~ 1e-2 near the
   hole), whether the pieces are added in the metric or in its inverse.
2. Equator, spin term: Kerr's T = (Omega_r/Omega_phi)^2 = 1 - 6u + 8 chi u^1.5 - 3 chi^2 u^2 exactly. With two
   velocities the spin term becomes 8 chi u^1.5 (1 - 2u)(1 - 3u) exactly (to first order in the spin): an extra
   -40 chi u^2.5 + 48 chi u^3.5. P.S needs the push E.P to act; with it, it takes those terms out at every order.
3. Equator, shape: two velocities give -3k^2 + 108 k^2 u^2 - 240 k^2 u^3 + 192 k^2 u^4 + ... (k = a/r); one jet gives
   exactly -3k^2.
4. Clock gap: one jet gives Kerr's 2 pi a per lap between prograde and retrograde at every radius (Part 1's cell
   clock gap); two velocities give 2 pi a (1 - 4M/r) at small spin, which vanishes at r = 4M.
5. Round 20's chi^2 u^3 cancellation does NOT come from the one jet: P.S adds 0 to it. Grouped by the jet's layers it is
   the first layer 2E.J = 2EP + 2ES = +20 - 12 = +8 against the second J.J = PP + SS + 2PS = 8 - 16 + 0 = -8.
6. Where the combined boost lives: on a circular or spherical orbit the push answers P = 2Mr(E + S)/Delta with one jet,
   i.e. the orbit's energy E - Omega_s L as seen by the shell turning at Omega_s = a/(r^2+a^2) (Part 3's shell rate),
   so k.p = (r^2+a^2)(E - Omega_s L)/Delta. With two velocities P = 2MrE/Delta: the push sees E alone.
7. Polar orbits (L = 0): S = 0, so their r, theta and clock motion are the same either way. Only the plane's turning
   differs: one jet 2MraE/Delta (Mino time; Kerr's), two velocities 2MraE/(r^2+a^2), slower by Delta/(r^2+a^2)
   (0.674 at a = 0.9M, r = 6M). The jet's cross piece is the strong-field part of polar frame dragging.
Run: python3 jet_boost.py (writes jet_boost.json, a few minutes)
"""
import itertools
import json
import os

import mpmath as mp
import sympy as sp

import kerr_split as ks
from numgr import prepare, ricci_max

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ("PP", "EP", "SS", "ES", "PS")


def composition():
    """Push and twist shears into the radial leg (legs ordered t, r, phi) compose additively: one shear."""
    a_, b_ = sp.symbols("a b")
    NP = sp.zeros(3, 3); NP[1, 0] = 1
    NS = sp.zeros(3, 3); NS[1, 2] = 1
    one = (sp.eye(3) + a_*NP)*(sp.eye(3) + b_*NS)
    return (NP*NS == sp.zeros(3, 3) and NS*NP == sp.zeros(3, 3)
            and one == sp.eye(3) + a_*NP + b_*NS and one == (sp.eye(3) + b_*NS)*(sp.eye(3) + a_*NP))


def jet_minus_two_velocities():
    E, P, S = sp.symbols("E P S")
    return sp.expand((E + P + S)**2 - ((E + P)**2 + (E + S)**2 - E**2))


# ---------------------------------------------------------------- 1. vacuum
t, r, th, ph, M, a = sp.symbols("t r theta phi M a", real=True)
X = (t, r, th, ph)
S2 = r**2 + a**2*sp.cos(th)**2
R2 = r**2 + a**2
V2 = 2*M*r/S2
ETA = sp.diag(-1, S2/R2, S2, R2*sp.sin(th)**2)
ETAI = sp.diag(-1, R2/S2, 1/S2, 1/(R2*sp.sin(th)**2))


def _outer(u, v):
    return sp.Matrix(4, 4, lambda i, j: u[i]*v[j])


def _pieces(E_, P_, S_, w):
    m = sp.zeros(4, 4)
    for n, u, v in (("EE", E_, E_), ("PP", P_, P_), ("SS", S_, S_)):
        m += w.get(n, 1)*_outer(u, v)
    for n, u, v in (("EP", E_, P_), ("ES", E_, S_), ("PS", P_, S_)):
        m += w.get(n, 1)*(_outer(u, v) + _outer(v, u))
    return m


def metric(w, where):
    """Kerr-Schild in oblate coordinates (t, r, theta, phi), dials on the pieces of k k, in the metric ('lower') or
    its inverse ('upper')."""
    if where == "lower":
        return ETA + V2*_pieces([-1, 0, 0, 0], [0, -S2/R2, 0, 0], [0, 0, 0, a*sp.sin(th)**2], w)
    return (ETAI - V2*_pieces([1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 0, a/R2], w)).inv(method="LU")


PTS = [(0, 3.1, 0.7, 0.3), (0, 5.3, 1.2, 1.0), (0, 8.0, 0.4, 2.0)]


def vacuum():
    mp.mp.dps = 50  # the check reads Ricci at 50 digits, whatever precision an earlier step left set
    out = {}
    for where in ("lower", "upper"):
        for name, w in (("one jet", {}), ("two velocities", {"PS": 0})):
            f = prepare(metric(w, where).subs({M: 1, a: sp.Rational(9, 10)}), X)
            out[f"{name} ({where})"] = [float(ricci_max(f, P)) for P in PTS]
    return out


# ---------------------------------------------------------------- 2-4. equator
def spin_term_exact(w):
    """d T / d a at a = 0 for the circular equatorial orbit at r = 1, M = u (so the spin term is this times a)."""
    Ms, as_, rs, Es, Ls, ps = ks.M, ks.a, ks.r, ks.E, ks.L, ks.p
    u, eps, E1, L1 = sp.symbols("u epsilon E1 L1")
    F = sp.expand(ks.ks_shell(w)); A = F.coeff(ps, 2); B = F.coeff(ps, 1)/2; C = F.coeff(ps, 0)
    Feff = C - B**2/A
    FL = sp.diff(F, Ls).subs(ps, -B/A)
    E0 = (1 - 2*u)/sp.sqrt(1 - 3*u); L0 = sp.sqrt(u)/sp.sqrt(1 - 3*u)
    sub = {Ms: u, as_: eps, Es: E0 + eps*E1, Ls: L0 + eps*L1}
    f0 = Feff.subs(sub).subs(rs, 1); f1 = sp.diff(Feff, rs).subs(sub).subs(rs, 1)
    e0 = sp.series(f0, eps, 0, 2).removeO().coeff(eps, 1); e1 = sp.series(f1, eps, 0, 2).removeO().coeff(eps, 1)
    sol = sp.solve([e0, e1], [E1, L1], dict=True)[0]
    Tx = (2*A*sp.diff(Feff, rs, 2)/FL**2).subs(sub).subs(rs, 1)
    return sp.factor(sp.simplify(sp.series(Tx, eps, 0, 2).removeO().coeff(eps, 1).subs(sol))), u


def spin_tower(w):
    """Coefficients of chi u^1.5, chi u^2.5, chi u^3.5 in T (linear in the spin), by an exact fit at small u."""
    mp.mp.dps = 60
    T = ks.orbit_factor(ks.ks_shell(w)); h = mp.mpf("1e-15")
    us = [mp.mpf("1e-5")*j for j in (1, 2, 3, 4, 5)]
    A = mp.matrix([[v**k for k in range(5)] for v in us])
    b = mp.matrix([(T(v, h) - T(v, -h))/(2*h)/v**mp.mpf(1.5) for v in us])
    s = mp.lu_solve(A, b)
    return [float(s[0]), float(s[1]), float(s[2])]


def shape_tower(w):
    """Coefficients of k^2, k^2 u, ..., k^2 u^4 in T (even in the spin), k = a/r."""
    mp.mp.dps = 80
    T = ks.orbit_factor(ks.ks_shell(w)); kk = mp.mpf("1e-20")
    us = [mp.mpf("1e-4")*j for j in range(1, 8)]
    g = [(T(v, kk/v) + T(v, -kk/v) - 2*T(v, mp.mpf(0)))/(2*kk**2) for v in us]
    c = mp.lu_solve(mp.matrix([[v**k for k in range(7)] for v in us]), mp.matrix(g))
    mp.mp.dps = 30
    return [str(sp.nsimplify(float(x), rational=True, tolerance=1e-6)) for x in c[:5]]


def clock_gap_fn(w):
    """dt/dphi of the circular equatorial orbit (period per radian); prograde minus retrograde is the clock gap / 2 pi."""
    mp.mp.dps = 50
    Ms, as_, rs, Es, Ls, ps = ks.M, ks.a, ks.r, ks.E, ks.L, ks.p
    F = sp.expand(ks.ks_shell(w)); A = F.coeff(ps, 2); B = F.coeff(ps, 1)/2; C = F.coeff(ps, 0); Feff = C - B**2/A
    args = (rs, Ms, as_, Es, Ls)
    f0, f1 = (sp.lambdify(args, sp.diff(Feff, rs, k), "mpmath") for k in (0, 1))
    FL = sp.lambdify(args, sp.diff(F, Ls).subs(ps, -B/A), "mpmath")
    FE = sp.lambdify(args, sp.diff(F, Es).subs(ps, -B/A), "mpmath")

    def dt_dphi(rv, Mv, av):
        e, l = mp.findroot(lambda e, l: (f0(rv, Mv, av, e, l), f1(rv, Mv, av, e, l)),
                           (1 - Mv/(2*rv), mp.sqrt(Mv*rv)*(1 + av/rv)))
        return -FE(rv, Mv, av, e, l)/FL(rv, Mv, av, e, l)
    return dt_dphi


def clock_gaps():
    out = {}
    for name, w in (("one jet", {}), ("two velocities", {"PS": 0})):
        f = clock_gap_fn(w)
        h = mp.mpf("1e-12")
        out[name] = {"a=0.5": [[rv, float(f(mp.mpf(rv), 1, mp.mpf("0.5")) - f(mp.mpf(rv), 1, mp.mpf("-0.5")))]
                               for rv in (6, 10, 30, 100)],
                     "per unit a, small spin": [[rv, float((f(mp.mpf(rv), 1, h) - f(mp.mpf(rv), 1, -h))/(2*h))]
                                                for rv in (4, 6, 10, 30, 100)]}
    return out


def push_on_orbit():
    """The push on a circular orbit, P = -B/A, for one jet and for two velocities (equatorial shell)."""
    out = {}
    for name, w in (("one jet", {}), ("two velocities", {"PS": 0})):
        F = sp.expand(ks.ks_shell(w)); A = F.coeff(ks.p, 2); B = F.coeff(ks.p, 1)/2
        out[name] = sp.factor(sp.simplify(-B/A))
    D = ks.r**2 - 2*ks.M*ks.r + ks.a**2
    Ss = -ks.a*ks.L/(ks.r**2 + ks.a**2)
    kp = sp.simplify(ks.E + Ss + out["one jet"] - (ks.r**2 + ks.a**2)*(ks.E - ks.a*ks.L/(ks.r**2 + ks.a**2))/D)
    return out, kp


# ---------------------------------------------------------------- 7. polar
def polar_node_rates():
    import polar_split as px
    E_, M_, a_, r_ = px.E, px.M, px.a, px.r
    full = px.node_rate_split()[0]
    sub = lambda d: sp.simplify(full.subs({px.WS[k]: v for k, v in d.items()}))
    one = sub({"EP": 1, "PP": 1, "ES": 1, "PS": 1})
    two = sub({"EP": 1, "PP": 1, "ES": 1, "PS": 0})
    D = r_**2 - 2*M_*r_ + a_**2
    ratio = sp.simplify(two/one)
    # r, theta and clock motion at L = 0 don't see P.S: S = 0
    A, B, C = px.radial_quadratic({}, Lv=0)
    no_ps = all(sp.diff(q, px.WS["PS"]) == 0 for q in (A, B, C))
    nums = [[rv, float(ratio.subs({M_: 1, a_: sp.Rational(9, 10), r_: rv}))] for rv in (3, 4, 6, 10, 20, 50, 100)]
    return (str(one), str(two), str(sp.simplify(one - 2*M_*r_*a_*E_/D)),
            str(sp.simplify(two - 2*M_*r_*a_*E_/(r_**2 + a_**2))), str(sp.factor(ratio)), no_ps, nums)


def main():
    out = {"composition_one_shear": bool(composition()), "jet_minus_two_velocities": str(jet_minus_two_velocities())}
    print(out, flush=True)
    out["vacuum_ricci_max"] = vacuum(); print(out["vacuum_ricci_max"], flush=True)
    s1, u = spin_term_exact({}); s2, _ = spin_term_exact({"PS": 0})
    out["spin_term_per_a"] = {"one jet": str(s1), "two velocities": str(s2)}
    print(out["spin_term_per_a"], flush=True)
    tow = {}
    for bits in itertools.product((0, 1), repeat=5):
        w = dict(zip(NAMES, bits))
        tow["+".join(["EE"] + [n for n, b in w.items() if b])] = spin_tower(w)
    out["spin_tower_chi_u15_u25_u35"] = tow
    out["shape_tower_k2"] = {"one jet": shape_tower({}), "two velocities": shape_tower({"PS": 0})}
    print(out["shape_tower_k2"], flush=True)
    out["clock_gap"] = clock_gaps(); print(out["clock_gap"], flush=True)
    po, kp = push_on_orbit()
    out["push_on_orbit"] = {k: str(v) for k, v in po.items()}; out["kp_minus_turning_shell_form"] = str(kp)
    one, two, c1, c2, ratio, no_ps, nums = polar_node_rates()
    out["polar"] = {"node_rate_one_jet": one, "node_rate_two_velocities": two, "one_minus_2MraE_over_Delta": c1,
                    "two_minus_2MraE_over_r2a2": c2, "ratio_two_over_one": ratio, "r_theta_clock_free_of_PS": no_ps,
                    "ratio_numbers_a0.9": nums}
    print(out["polar"], flush=True)
    eq = json.load(open(os.path.join(HERE, "polar_split.json")))["equator_k2_and_k2u"]
    out["round20_by_jet_layer"] = {"2E.J = 2EP + 2ES": 20 - 12, "J.J = PP + SS + 2PS": 8 - 16 + 0,
                                   "full": float(eq["EE+PP+EP+SS+ES+PS"][1]),
                                   "full without PS": float(eq["EE+PP+EP+SS+ES"][1])}
    json.dump(out, open(os.path.join(HERE, "jet_boost.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
