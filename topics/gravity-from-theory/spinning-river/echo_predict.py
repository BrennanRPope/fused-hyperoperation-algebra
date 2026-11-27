"""DRAFT, round 29 (Nov 27, 2026): the radar echo delay past the Sun, PREDICTED piece by piece before comparing.

Brennan asked for the Shapiro delay (the extra time a radar echo takes when its path grazes the Sun) split into the
river's time part, push and twist, and, before looking at Cassini's measurement, for the dial settings to be fixed from
earlier results and the prediction written down with a timestamp. This script makes the prediction and writes
echo_prediction.json. It contains no measured delay, no measured gamma and no Cassini number; echo_compare.py does the
comparison afterwards.

The dials. The shell's field is f (k.p)^2 with k.p = E + P + S (time part, push, twist). A dial that changes how strong
the push or the twist is, while keeping the field one square (the Kerr-Schild form, and round 21's one-jet rule), is
k.p = E + lambda P + tau S, so the push enters as 2 lambda E P + lambda^2 P P. Rounds 7-8 derived the river with
lambda = tau = 1 (that is exact Kerr), and rounds 25-28 used those settings unchanged. They were never fitted to any
delay. The shadow (round 28) gives a second, data-only setting for lambda.

What the script finds (no Cassini input):
1. Per leg, to first order in M, the light time from the closest approach r0 out to r is
     sqrt(r^2 - r0^2) + M [ (1 + lambda^2) ln((r + sqrt(r^2 - r0^2))/r0) + sqrt((r - r0)/(r + r0)) ] (+- 2 lambda M ln r,
   a clock relabelling that cancels out of anything that depends on the path).
   The time part gives the 1 in (1 + lambda^2): its half of the delay is clocks running slow near the Sun. The push gives
   lambda^2. Newton's corpuscle (light that speeds up as it falls) would give -1: it arrives early.
2. In the usual bookkeeping (Eddington's gamma, delay ~ (1 + gamma)), the river predicts gamma = lambda^2.
3. The same dial sets light bending, 2(1 + lambda^2) M/b, and the perihelion advance, 6 lambda^2 pi M/a per orbit (Newton's
   ring alone gives none), so one dial ties all three. Without the square (separate weights mu on EP and nu on PP) the perihelion factor is
   (2 mu^2 + nu)/3 and the delay's is nu: the square is what makes them agree.
4. Twist: the Sun's spin adds a delay linear in a = J/(Mc), opposite for rays passing on the two sides; for an
   equatorial ray it is about 8 M a / b one way.
Run: python3 echo_predict.py  (writes echo_prediction.json, about a minute)
"""
import datetime
import hashlib
import json
import os

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 40

# constants (SI): the Sun's GM, radius and angular momentum; the astronomical unit
GM_SUN = 1.32712440018e20      # m^3/s^2 (IAU nominal)
C = 2.99792458e8
R_SUN = 6.957e8                # m (IAU nominal)
AU = 1.495978707e11            # m
J_SUN = 1.94e41                # kg m^2/s, helioseismology (Pijpers 1998)
M_SUN_KG = GM_SUN/6.67430e-11
M_GEOM = GM_SUN/C**2           # 1476.6 m


# ------------------------------------------------------------------ 1. the leg time for the dial family (symbolic)
def leg_formula():
    r, r0, M, lam, pr, Ev, b2 = sp.symbols("r r0 M lambda p_r E_v B2", positive=True)
    H2 = -Ev**2 + pr**2 + b2/r**2 - (2*M/r)*(Ev + lam*pr)**2        # 2H, null; E = -p_t, L^2 = b^2 E^2
    tdot = (-sp.diff(H2, Ev)/2).subs(Ev, 1)
    rdot = (sp.diff(H2, pr)/2).subs(Ev, 1)
    H2 = H2.subs(Ev, 1)
    pr0 = sp.solve(sp.Eq(rdot.subs(r, r0), 0), pr)[0]
    b2_r0 = sp.solve(H2.subs({r: r0, pr: pr0}), b2)[0]
    out = {"b^2(r0)": str(sp.factor(b2_r0))}
    s = sp.sqrt(r**2 - r0**2)
    want_out = lam**2/s + 2*lam/r + (r + 2*r0)/((r + r0)*s)                # outgoing leg, first order in M
    roots = sp.solve(H2.subs(b2, b2_r0), pr)
    coeffs = []
    for p_ in roots:
        dtdr = sp.simplify((tdot/rdot).subs(pr, p_))
        c1 = sp.series(dtdr, M, 0, 2).removeO().coeff(M, 1)
        coeffs.append(sp.simplify(c1))
    ok_out = any(sp.simplify(c_ - want_out) == 0 for c_ in coeffs)
    ok_in = any(sp.simplify(c_ + (want_out - 4*lam/r)) == 0 for c_ in coeffs)  # incoming: dr < 0, clock term flips
    out["outgoing_first_order_is_1+lam2_log_plus_clock"] = ok_out
    out["incoming_mirror"] = ok_in
    # antiderivatives of the three pieces
    Lg = sp.log((r + s)/r0)
    pts = [{r: sp.Rational(7, 2) + k_, r0: sp.Rational(3, 2) + k_/sp.Integer(3)} for k_ in range(4)]
    num0 = lambda e: all(abs(sp.N(e.subs(p_), 40)) < sp.Float("1e-35") for p_ in pts)   # 40-digit spot checks
    out["d/dr ln((r+s)/r0) = 1/s"] = num0(sp.diff(Lg, r) - 1/s)
    out["d/dr sqrt((r-r0)/(r+r0)) = r0/((r+r0) s)"] = num0(sp.diff(sp.sqrt((r - r0)/(r + r0)), r) - r0/((r + r0)*s))
    # Newton's corpuscle: v^2 = 1 + 2M/r, dt/dr = 1/v_r, same turning point (b^2 = r0^2 + 2 M r0)
    vr = sp.sqrt((r - r0)*(r + r0 + 2*M))/r
    cN = sp.simplify(sp.series(1/vr, M, 0, 2).removeO().coeff(M, 1))
    out["newton_corpuscle_first_order_is_-log_plus_sqrt"] = num0(cN - (-1/s + r0/((r + r0)*s)))
    return out


# ------------------------------------------------------------------ 2. perihelion and bending with the dial
def perihelion_dial():
    u, M, mu, nu, x, v2 = sp.symbols("u M mu nu x v2", positive=True)
    E2, L2 = sp.symbols("E2 L2", positive=True)
    r = 1/u
    A = 1 - 2*M*nu/r
    B = -2*M*mu*sp.sqrt(E2)/r
    Cc = 1 - E2 - 2*M*E2/r + L2/r**2
    Fs = sp.expand((B**2 - A*Cc)*r**4*u**4/L2)                            # (du/dphi)^2
    sol = sp.solve([Fs.subs(u, x), sp.diff(Fs, u).subs(u, x)], [E2, L2], dict=True)[0]
    k = sp.simplify(-sp.diff(Fs, u, 2).subs(u, x).subs(sol)/2)            # (Omega_r/Omega_phi)^2
    ser = sp.expand(sp.series(k.subs(x, v2/M), v2, 0, 2).removeO())
    lam = sp.symbols("lambda", positive=True)
    return {"ratio_first_order": str(ser),
            "general_is_1-2(2mu^2+nu)v^2": sp.simplify(ser - (1 - 2*(2*mu**2 + nu)*v2)) == 0,
            "square_dial_is_1-6lam^2v^2": sp.simplify(ser.subs({mu: lam, nu: lam**2}) - (1 - 6*lam**2*v2)) == 0}


def R_dial(r, b, lam, M=1):
    """Null radial function (times r^4) for k.p = E + lambda P, no spin."""
    return 4*M**2*lam**2*r**2 + (r**2 - 2*M*lam**2*r)*(r**2 + 2*M*r - b**2)


def bending_dial(lam, b, M=1):
    """Deflection for k.p = E + lambda P, no spin: 2 int_0^u0 b du / sqrt(F(u)) - pi, F(u) = u^4 R(1/u)."""
    lam, b, M = mp.mpf(lam), mp.mpf(b), mp.mpf(M)
    # F(u) = 4 M^2 lam^2 u^2 + (1 - 2 M lam^2 u)(1 + 2 M u - b^2 u^2)  (cubic in u)
    c = [2*M*lam**2*b**2, 4*M**2*lam**2 - b**2 - 4*M**2*lam**2, 2*M - 2*M*lam**2, mp.mpf(1)]  # u^3, u^2, u, 1
    F = lambda u: ((c[0]*u + c[1])*u + c[2])*u + c[3]
    u0 = mp.findroot(F, 1/b)
    # F(u) = (u - u0) * (q0 u^2 + q1 u + q2)
    q0 = c[0]; q1 = c[1] + q0*u0; q2 = c[2] + q1*u0
    G = lambda u: -(q0*u**2 + q1*u + q2)                                  # F = (u0 - u) G
    f = lambda t: 4*b*mp.sqrt(u0)/mp.sqrt(G(u0*(1 - t**2)))
    return mp.quad(f, [0, 1]) - mp.pi


def bending_checks():
    rows = []
    for lam in ("0", "0.5", "1"):
        bb = 1e5
        d = bending_dial(lam, bb)
        rows.append({"lambda": float(lam), "b": bb, "deflection_x_b_over_2M": float(d*bb/2),
                     "predicted_1+lam2": 1 + float(lam)**2})
    return rows


# ------------------------------------------------------------------ 3. exact light times at the Sun (no spin)
def exact_leg_excess(lam, r0, r1, M=1):
    """Exact light time from the closest approach r0 out to r1 for the dial family, minus the flat sqrt(r1^2 - r0^2)
    and the clock relabelling 2 lambda M ln(r1/r0)."""
    lam, r0, r1, M = mp.mpf(lam), mp.mpf(r0), mp.mpf(r1), mp.mpf(M)
    b2 = r0**2*(r0 + 2*M - 2*M*lam**2)/(r0 - 2*M*lam**2)
    # R(r) = 4M^2 lam^2 r^2 + (r^2 - 2 M lam^2 r)(r^2 + 2 M r - b2): factor out (r - r0)
    # R = r^4 + (2M - 2M lam^2) r^3 - b2 r^2 + 2 M lam^2 b2 r   (the 4 M^2 lam^2 r^2 terms cancel)
    c = [mp.mpf(1), 2*M - 2*M*lam**2, -b2, 2*M*lam**2*b2, mp.mpf(0)]
    qc = [c[0]]
    for ci in c[1:-1]:
        qc.append(ci + qc[-1]*r0)

    def f(s):
        r = r0 + s**2
        A = 1 - 2*M*lam**2/r
        q = qc[0]*r**3 + qc[1]*r**2 + qc[2]*r + qc[3]                     # R/(r - r0)
        sq = s*mp.sqrt(q)                                                  # sqrt(R)
        p = (2*M*lam/r + sq/r**2)/A                                        # outgoing p_r
        tdot = 1 + (2*M/r)*(1 + lam*p)
        # dt/dr = tdot / (sqrt(R)/r^2); times dr = 2 s ds
        return tdot*r**2*2/mp.sqrt(q)
    T = mp.quad(f, [0, mp.sqrt(r0)*mp.mpf("0.01"), mp.sqrt(r0), mp.sqrt(r1 - r0)])
    return T - mp.sqrt(r1**2 - r0**2) - 2*lam*M*mp.log(r1/r0)             # remove the clock relabelling


def first_order_leg_excess(lam, r0, r1, M=1):
    lam, r0, r1 = mp.mpf(lam), mp.mpf(r0), mp.mpf(r1)
    s = mp.sqrt(r1**2 - r0**2)
    return M*((1 + lam**2)*mp.log((r1 + s)/r0) + mp.sqrt((r1 - r0)/(r1 + r0)))


def sun_exact_checks():
    r0 = mp.mpf(1.6*R_SUN/M_GEOM)
    r1 = mp.mpf(1.016*AU/M_GEOM)
    out = []
    for lam in ("0", "1"):
        ex = exact_leg_excess(lam, r0, r1)
        fo = first_order_leg_excess(lam, r0, r1)
        out.append({"lambda": float(lam), "exact_excess_m": float(ex*M_GEOM), "first_order_excess_m": float(fo*M_GEOM),
                    "difference_m": float((ex - fo)*M_GEOM), "M2_over_r0_m": float(M_GEOM/r0)})
    return out


# ------------------------------------------------------------------ 4. the twist: linear in the Sun's spin
def kerr_equatorial_time(M, a, r0, r1, r2):
    M, a, r0, r1, r2 = (mp.mpf(x) for x in (M, a, r0, r1, r2))

    def R(r, L):
        P = r**2 + a**2 - a*L
        return P**2 - (r**2 - 2*M*r + a**2)*(L - a)**2
    L = mp.findroot(lambda L: R(r0, L), r0 + M)
    K = (L - a)**2
    c = [mp.mpf(1), mp.mpf(0), 2*a**2 - 2*a*L - K, 2*M*K, (a**2 - a*L)**2 - a**2*K]
    qc = [c[0]]
    for ci in c[1:-1]:
        qc.append(ci + qc[-1]*r0)

    def leg(re):
        def f(s):
            r = r0 + s**2
            P = r**2 + a**2 - a*L
            D = r**2 - 2*M*r + a**2
            num = (r**2 + a**2)*P/D + a*(L - a)
            q = qc[0]*r**3 + qc[1]*r**2 + qc[2]*r + qc[3]
            return 2*num/mp.sqrt(q)
        return mp.quad(f, [0, mp.sqrt(r0)*mp.mpf("0.01"), mp.sqrt(r0), mp.sqrt(re - r0)])
    return leg(r1) + leg(r2)


def twist_delay():
    rows = []
    for r0, r1, r2 in ((1000, 1e6, 1e6), (2000, 1e6, 3e6), (5000, 1e7, 1e7)):
        a = mp.mpf("0.1")
        d = (kerr_equatorial_time(1, a, r0, r1, r2) - kerr_equatorial_time(1, -a, r0, r1, r2))/2
        rows.append({"r0": r0, "r1": r1, "r2": r2, "half_difference_over_Ma/r0": float(d/(a/r0))})
    a_sun = J_SUN/(M_SUN_KG*C)                                              # metres
    b = 1.6*R_SUN
    one_way = 8*M_GEOM*a_sun/b/C
    return {"rows": rows, "coefficient": -8, "sun_a_m": a_sun, "sun_a_over_M": a_sun/M_GEOM,
            "sun_one_way_s_at_1.6Rsun": one_way, "sun_round_trip_s_at_1.6Rsun": 2*one_way}


# ------------------------------------------------------------------ 5. the push dial from the shadow alone (round 28)
def shadow_radius(lam):
    lam = mp.mpf(lam)
    # double root of R_dial in r: R = R' = 0, unknowns r, b
    def eqs(r, b):
        return [R_dial(r, b, lam), mp.diff(lambda x: R_dial(x, b, lam), r)]
    r_, b_ = mp.findroot(eqs, (mp.mpf(3)*lam**2 + mp.mpf("0.5"), mp.mpf(5.2)*lam + mp.mpf("0.3")))
    return float(r_), float(b_)


def shadow_dial():
    shd = json.load(open(os.path.join(HERE, "shadow_split.json")))["eht"]
    tg = [shd["GRAVITY 2019 (VLTI)"]["theta_g_uas"], shd["Do et al. 2019 (Keck)"]["theta_g_uas"]]
    d, e = shd["EHT_2022_shadow_diameter_uas"]
    bc_lo, bc_hi = (d - e)/(2*max(tg)), (d + e)/(2*min(tg))
    table = []
    for lam in [0.6 + 0.05*i for i in range(15)]:
        try:
            table.append((lam, shadow_radius(lam)[1]))
        except Exception:
            pass
    # invert b_c(lambda) by bisection
    def lam_for(bc):
        lo, hi = 0.6, 1.3
        for _ in range(60):
            mid = 0.5*(lo + hi)
            if shadow_radius(mid)[1] < bc:
                lo = mid
            else:
                hi = mid
        return 0.5*(lo + hi)
    l_lo, l_hi = lam_for(bc_lo), lam_for(bc_hi)
    return {"shadow_radius_bc_over_M_range": [bc_lo, bc_hi], "bc_lambda_1": shadow_radius(1)[1],
            "bc_table": table, "lambda_range": [l_lo, l_hi], "gamma_range_lambda2": [l_lo**2, l_hi**2]}


# ------------------------------------------------------------------ 6. the prediction
def prediction(tw):
    # illustrative geometry: Earth at 1.016 AU (June), far station 7.4 AU, ray passing 1.6 solar radii from centre
    r1, r2, b = 1.016*AU, 7.4*AU, 1.6*R_SUN
    log_ = mp.log(4*r1*r2/b**2)
    t_unit = GM_SUN/C**3
    one = float(t_unit*log_)                                               # one-way, per unit of (1 + lambda^2)
    return {"gamma_predicted": 1.0, "rule": "gamma = lambda^2, dials lambda = tau = 1 from rounds 7-8",
            "gamma_from_shadow_alone": "lambda^2 within the shadow_dial range (no use of rounds 7-8's setting)",
            "one_dial_relation": "perihelion factor (2 + 2 gamma - beta)/3 = lambda^2 and gamma = lambda^2, so "
                                 "beta + gamma = 2 for any push strength, if the field stays one square",
            "delay_split": {"time_part": 0.5, "push": 0.5, "twist": "linear in the Sun's spin, opposite on the two "
                            "sides, about 8 GM J/(c^4 M b) one way (well below a part in a million of the total)"},
            "illustrative_geometry": {"r_earth_AU": 1.016, "r_far_AU": 7.4, "b_solar_radii": 1.6},
            "round_trip_delay_s": {"time_part": 2*one, "push": 2*one, "total": 4*one,
                                   "twist_round_trip": tw["sun_round_trip_s_at_1.6Rsun"]},
            "doppler_note": "Cassini-style tracking measures d/dt of the round-trip delay; its relativistic part is "
                            "proportional to (1 + gamma) = 1 + lambda^2, so the prediction for it is the same gamma."}


def main():
    res = {"written_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
           "uses_any_measured_delay_or_gamma": False,
           "dial_settings": {"lambda_push": 1, "tau_twist": 1,
                             "set_by": "rounds 7-8 (the river derived from the shift; exact Kerr), used unchanged in "
                                       "rounds 25-28 (river_term.py, equator_split.py, kn_split.py, shadow_split.py)"}}
    res["leg_formula"] = leg_formula()
    print(res["leg_formula"], flush=True)
    res["perihelion_dial"] = perihelion_dial()
    print(res["perihelion_dial"], flush=True)
    res["bending_dial"] = bending_checks()
    print(res["bending_dial"], flush=True)
    res["sun_exact"] = sun_exact_checks()
    print(res["sun_exact"], flush=True)
    res["twist"] = twist_delay()
    print(res["twist"], flush=True)
    res["shadow_dial"] = shadow_dial()
    print(res["shadow_dial"], flush=True)
    res["prediction"] = prediction(res["twist"])
    res["script_sha256"] = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    json.dump(res, open(os.path.join(HERE, "echo_prediction.json"), "w"), indent=1, default=str)
    print(res["prediction"], flush=True)


if __name__ == "__main__":
    main()
