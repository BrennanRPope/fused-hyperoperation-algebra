"""DRAFT. Round 34: do the river's own rules fix how strongly the push and twist meet the time part in motion?

Part 27 pinned relativity's pi^2 to one term, kappa (a . grad) F with kappa = 41/16 (equivalently -kappa (D1 + D2)^2 F),
F = (pi^2/4) G^3 m1^2 m2^2 / r12. The question is whether anything the river already has fixes kappa.

1. Symmetry can't. The river's jet rule is about how a body's field looks when it moves (boosts). Under a boost every
   velocity shifts by the same u, and (D1 + D2) F = (v1 - v2) . grad F does not change, because F depends only on
   y1 - y2. So kappa (D1 + D2)^2 F is boost-invariant on its own: boost symmetry ties velocity terms to the static
   terms one level down (rational, round 14), but it leaves kappa free. A shape that is not boost-invariant, such as
   alpha (D1^2 + D2^2) F, does change, which is how boosts fixed round 31's shape and nothing more.
2. The rules the river has give kappa = 0: the four rules give only plain fractions at the third level (round 30), and
   the time part's field energy with the lag gives a total time derivative (round 32).
3. What does fix it is a consistency rule: the field's own energy pulls, pushes and drags exactly as mass does, at
   every order, and stays conserved. That leaves no free number: once the energy density is fixed, conservation fixes
   the energy current, the current must drag frames like a mass current, and so on up the orders. Done all the way,
   this rule is known to give exactly Einstein's equations (Kraichnan 1955, Deser 1970), and so 41/16. Here the first
   link of that chain is checked: for any moving Newtonian field, the field-energy density (grad U)^2/8pi and the
   current S = -(dU/dt) grad U /4pi obey d/dt(density) + div S = 0 wherever there is no mass. The current is fixed by
   the density: the strength with which moving field energy meets the twist is not a new number.

    python3 rules_fix_strength.py     # writes rules_fix_strength.json (seconds)
"""
import json, os
import sympy as sp
import lagged_field_energy as lg

HERE = os.path.dirname(os.path.abspath(__file__))
u = sp.Matrix(sp.symbols('ux uy uz', real=True))


def boosted(expr):
    """Shift both velocities by u (a boost at this order)."""
    sub = {lg.v1[i]: lg.v1[i] + u[i] for i in range(3)}
    sub.update({lg.v2[i]: lg.v2[i] + u[i] for i in range(3)})
    return expr.xreplace(sub)


def boost_checks():
    F, D1, D2 = lg.F, lg.D1, lg.D2
    rel = D1(D1(F)) + D2(D2(F)) + 2*D1(D2(F))          # (D1 + D2)^2 F, relativity's velocity form / -kappa
    lop = D1(D1(F)) + D2(D2(F))                         # a shape that is not boost-invariant
    return {'relativity_form_boost_invariant': bool(sp.simplify(boosted(rel) - rel) == 0),
            'lopsided_form_changes': bool(sp.simplify(boosted(lop) - lop) != 0)}


def conservation_check():
    """Two bodies moving with constant velocities: U = m1/|x - y1(t)| + m2/|x - y2(t)|.
    density e = (grad U)^2/8pi, current S = -(dU/dt) grad U/4pi: de/dt + div S = -(dU/dt) lap U /4pi = 0 off the bodies."""
    t = sp.symbols('t', real=True)
    X = sp.Matrix(sp.symbols('X Y Z', real=True))
    m1, m2 = sp.symbols('m1 m2', positive=True)
    p1 = lg.y1 + lg.v1*t; p2 = lg.y2 + lg.v2*t
    U = m1/sp.sqrt(((X - p1).T*(X - p1))[0]) + m2/sp.sqrt(((X - p2).T*(X - p2))[0])
    grad = sp.Matrix([sp.diff(U, c) for c in X])
    e = (grad.T*grad)[0]/(8*sp.pi)
    S = -sp.diff(U, t)*grad/(4*sp.pi)
    res = sp.diff(e, t) + sum(sp.diff(S[i], X[i]) for i in range(3))
    import random
    rng = random.Random(34)
    syms = list(X) + list(lg.y1) + list(lg.y2) + list(lg.v1) + list(lg.v2) + [m1, m2, t]
    vals = []
    for _ in range(4):
        pt = {s: sp.Rational(rng.randint(-20, 20), rng.randint(3, 9)) for s in syms}
        pt[m1] = sp.Rational(rng.randint(1, 9), 4); pt[m2] = sp.Rational(rng.randint(1, 9), 4)
        vals.append(abs(sp.N(res.subs(pt), 30)))
    return bool(max(vals) < 1e-20)


def main():
    lfe = json.load(open(os.path.join(HERE, 'lagged_field_energy.json')))
    tlp = json.load(open(os.path.join(HERE, 'third_level_prediction.json')))
    out = {
        'boost': boost_checks(),
        'kappa_from_symmetry': 'free',
        'kappa_four_rules_round30': 0,
        'kappa_lagged_time_part_round32': lfe['K_time_with_lag'],
        'kappa_relativity': '41/16',
        'kappa_consistency_rule': '41/16 (Kraichnan 1955, Deser 1970: gives Einstein\'s equations)',
        'field_energy_conserved_off_bodies': bool(conservation_check()),
        'round30_claim': tlp['claim'],
    }
    json.dump(out, open(os.path.join(HERE, 'rules_fix_strength.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
