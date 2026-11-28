"""DRAFT. Round 35: is there a reason in the river's picture why field energy must pull like mass?

Round 34 found that what fixes relativity's 41/16 is a consistency rule: the field's own energy pulls, pushes and drags
exactly as mass does. This round asks whether the river's picture gives a reason for that rule, and tests the first
rung of it.

The river's picture says everything rides the river the same way. Taken at its word, that includes the river's own
energy: a body held together by gravity (negative field energy) must fall exactly like a body of the same total
energy held together any other way. That is the strong equivalence principle. Here it is tested at the first level.

1. The weights (NE1). From the river's first-level metric, time part g00 = -1 + 2U - 2 beta U^2 (beta = how strongly
   field energy pulls, rounds 13-15) and push g_ij = (1 + 2 gamma U) delta_ij (gamma = how strongly the push acts,
   round 29), a body of mass m has L/m = v^2/2 + U + v^4/8 + (gamma + 1/2) U v^2 + (1/2 - beta) U^2.
2. A bound pair falling (NE2). Put a pair of bodies (masses ma, mb, separation r, internal speeds ua, ub) far from a mass
   M at distance R, so each body feels U = M/R + (partner)/r. Collect everything that multiplies M/R: the pair's
   falling weight is m + (2 gamma + 1) T + 2 (2 beta - 1) Omega, with T its internal kinetic energy and Omega = -ma mb/r
   its field (binding) energy. Its inertia is its total energy m + T + Omega (the first-level rules are boost
   invariant, a standard fact we don't recheck here). For a pair in a steady orbit, 2T + Omega = 0, and the ratio of
   falling weight to inertia is 1 + eta Omega/m with eta = 4 beta - gamma - 3: Nordtvedt's formula (1968),
   rederived here from the river's two strengths.
3. What it means. With the river's push, gamma = 1, "the river's own energy rides the river too" (eta = 0) forces
   beta = 1, which is exactly the strength rounds 13-15 had to put in by hand for field energy pulling at the first
   level. So at the first level the consistency rule is not an extra rule: it follows from universality.
4. Measured (NE3). Lunar laser ranging watches whether the Earth's gravitational self-energy makes it fall differently
   from the Moon toward the Sun: eta = (4.4 +/- 4.5) x 10^-4 (Will 2014 review). With Cassini's gamma - 1 =
   (2.1 +/- 2.3) x 10^-5 (Bertotti, Iess, Tortora 2003), beta - 1 = (eta + gamma - 1)/4 = (1.2 +/- 1.1) x 10^-4, within 2 sigma of 0.

    python3 own_energy_falls.py     # writes own_energy_falls.json (seconds)
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
beta, gamma = sp.symbols('beta gamma', real=True)


def weights():
    """Test-body Lagrangian from the river's first-level metric, expanded to first level (eps = 1/c^2 bookkeeping)."""
    U, v, e = sp.symbols('U v epsilon', positive=True)
    g00 = -1 + 2*e*U - 2*beta*e**2*U**2
    gij = 1 + 2*gamma*e*U
    L = -sp.sqrt(-(g00 + gij*e*v**2))          # per unit mass, with v^2 counted as one order like U
    Ls = sp.expand(sp.series(L, e, 0, 3).removeO())
    want = -1 + e*(v**2/2 + U) + e**2*(v**4/8 + (gamma + sp.Rational(1, 2))*U*v**2 + (sp.Rational(1, 2) - beta)*U**2)
    return {"lagrangian": str(Ls), "matches": sp.simplify(Ls - sp.expand(want)) == 0,
            "coef_Uv2": str(sp.Rational(1, 2) + gamma), "coef_U2": str(sp.Rational(1, 2) - beta)}


def bound_pair():
    """Falling weight and inertia of a bound pair; returns eta."""
    ma, mb, ua, ub, r, M, R = sp.symbols('m_a m_b u_a u_b r M R', positive=True)
    def body(m, u, partner):
        Uf = M/R + partner/r
        return m*(u**2/2 + Uf + u**4/8 + (gamma + sp.Rational(1, 2))*Uf*u**2 + (sp.Rational(1, 2) - beta)*Uf**2)
    L = sp.expand(body(ma, ua, mb) + body(mb, ub, ma))
    fall = sp.expand(sp.diff(L, M).subs(M, 0)*R)          # what multiplies M/R, at linear order in M
    T = ma*ua**2/2 + mb*ub**2/2
    Om = -ma*mb/r
    m = ma + mb
    fall_expected = m + (2*gamma + 1)*T + 2*(2*beta - 1)*Om
    weight_ok = sp.simplify(fall - sp.expand(fall_expected)) == 0
    # steady orbit: 2T + Omega = 0; first order in Omega/m
    t, W, mm = sp.symbols('T Omega m', real=True)
    x = sp.symbols('x')
    fallW = mm + (2*gamma + 1)*t + 2*(2*beta - 1)*W
    inert = mm + t + W
    ratio = (fallW/inert).subs(t, -W/2).subs(W, x*mm)
    eta = sp.simplify(sp.diff(sp.series(ratio, x, 0, 2).removeO(), x))
    return {"falling_weight": str(fall_expected), "weight_ok": bool(weight_ok),
            "inertia": "m + T + Omega", "eta": str(sp.expand(eta)),
            "eta_is_nordtvedt": sp.simplify(eta - (4*beta - gamma - 3)) == 0, "eta_expr": eta}


def main():
    w = weights()
    bp = bound_pair()
    eta = bp.pop("eta_expr")
    river = eta.subs({beta: 1, gamma: 1})
    forced_beta = sp.solve(sp.Eq(eta.subs(gamma, 1), 0), beta)
    # measurements
    eta_m, eta_s = 4.4e-4, 4.5e-4            # Will 2014, Living Reviews in Relativity 17, 4 (lunar laser ranging)
    gm1_m, gm1_s = 2.1e-5, 2.3e-5            # Bertotti, Iess, Tortora 2003 (Cassini)
    b1_m = (eta_m + gm1_m)/4
    b1_s = (eta_s**2 + gm1_s**2)**0.5/4
    out = {"weights": {k: (bool(v) if k == "matches" else v) for k, v in w.items()},
           "bound_pair": {k: (bool(v) if isinstance(v, (bool, sp.logic.boolalg.BooleanAtom)) else v)
                          for k, v in bp.items()},
           "eta_river_beta1_gamma1": str(river),
           "beta_forced_by_universality_with_gamma1": [str(s) for s in forced_beta],
           "measured": {"eta": [eta_m, eta_s], "gamma_minus_1": [gm1_m, gm1_s],
                        "beta_minus_1": [b1_m, b1_s],
                        "beta_minus_1_within_2sigma_of_0": abs(b1_m) <= 2*b1_s},
           "sources": ["Nordtvedt 1968, Phys. Rev. 169, 1017 and 1014", "Will 2014, Living Rev. Relativ. 17, 4",
                       "Bertotti, Iess, Tortora 2003, Nature 425, 374"]}
    json.dump(out, open(os.path.join(HERE, "own_energy_falls.json"), "w"), indent=2)
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
