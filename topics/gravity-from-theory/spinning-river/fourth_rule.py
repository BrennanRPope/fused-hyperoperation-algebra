"""Round 14 (DRAFT): does adding gravity's own field energy as a fourth rule close the 2PN gap?

Part 9 found that the three rules (conformal time-time and space-space factor, the moving source's push, the
multiplied clocks) miss the exact 2PN redshift of two bodies on a circle by a clean polynomial,
    miss_rules = m1 m2 (23 m1^2 - 50 m1 + 36)/3   (at x^3, total mass 1).
Here we add, piece by piece, what general relativity's 2PN harmonic field has and the rules lack: the field's
own energy acting as a source (the Feynman-Deser bootstrap idea) in its 2PN form.

Stages, all fed through the same redshift pipeline as Part 9 (next_order.redshift with a hook):
  1. rules only                                   -> miss_rules
  2. + closed-form field-energy pieces             (fourth_rule_closed.py: Fock's function ln(r1+r2+r12))
  3. + static field-on-field pieces alpha..delta   (fourth_rule_integrals.py: four 3D integrals, Hadamard finite part)
  4. what is left, and the single point-like term that removes it.

Run:  python3 fourth_rule.py            (uses fourth_rule_cache/ when present; else about 6 min per mass ratio)
      python3 fourth_rule.py --fresh    (recomputes the closed-form pieces)
Writes fourth_rule.json.
"""
import sympy as sp, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import next_order as no

RATIOS = ['1/3', '1/2', '1/4', '2/3', '1/5', '2/5']     # the last one is held out of every fit
XKEYS = ('Xa', 'Xb', 'Xc', 'Xd', 'Xe', 'Xf1c', 'Xf2c', 'Xacc')
CACHE = os.path.join(HERE, 'fourth_rule_cache')
al, be, ga, de = sp.symbols('alpha beta gamma delta')
F_VALUES = {al: sp.Rational(-1, 2), be: sp.Integer(1), ga: sp.Integer(0), de: sp.Integer(-2)}  # fourth_rule_integrals.py
m = sp.symbols('m1', positive=True)


def _key(k): return k if isinstance(k, str) else ','.join(map(str, k))
def _unkey(s): return tuple(int(p) if p.isdigit() else p for p in s.split(',')) if ',' in s else s


def closed_pieces(ratio, fresh=False):
    path = os.path.join(CACHE, f'closed_{ratio.replace("/", "_")}.json')
    if os.path.exists(path) and not fresh:
        raw = json.load(open(path))
        return {int(b): {_unkey(k): tuple(eval_srepr(v) for v in pair)
                         for k, pair in d.items()} for b, d in raw.items()}
    import fourth_rule_closed as fc
    out = fc.at_bodies(sp.Rational(ratio))
    os.makedirs(CACHE, exist_ok=True)
    json.dump({b: {_key(k): [sp.srepr(v[0]), sp.srepr(v[1])] for k, v in d.items()} for b, d in out.items()},
              open(path, 'w'), indent=0)
    return out


def eval_srepr(s):
    return eval(s, {k: getattr(sp, k) for k in dir(sp) if not k.startswith('_')})


def make_hook(ex, fextra=None):
    """Add 4 Wc to g_ij, -8 Rc to g_0i, 8 Xcross + 16 V1k V2k (+ 8 f-terms) to g_00, at both bodies."""
    def hook(res):
        new = {}
        for b in (1, 2):
            G, dG, sqmg = res[b]; G = G.copy(); dG = dG.copy(); E = ex[b]
            for w, tgt in ((0, G), (1, dG)):
                add00 = 8*sum(E[k][w] for k in XKEYS) + 16*E['VV'][w]
                if fextra is not None: add00 += fextra[b][w]
                tgt[0, 0] += add00
                for i in range(3):
                    tgt[0, i+1] += -8*E[('Rc', i)][w]; tgt[i+1, 0] += -8*E[('Rc', i)][w]
                    for j in range(i, 3):
                        tgt[i+1, j+1] += 4*E[('Wc', i, j)][w]
                        if i != j: tgt[j+1, i+1] += 4*E[('Wc', i, j)][w]
            sqmg = sqmg + 2*sum(E[('Wc', k, k)][0] for k in range(3))
            new[b] = (G, dG, sqmg)
        return new
    return hook


def f_terms(m1, m2):
    """Static field-on-field part of 8 X-hat at each body: value and x-gradient, alpha..delta left symbolic.
    At body 1:  X-hat = (alpha m1^2 m2 + beta m1 m2^2)/r12^3,  dX/dx = (gamma m1^2 m2 + delta m1 m2^2)/r12^4."""
    e = no.e; R1, R2 = sp.symbols('R1 R2', positive=True); D = R1 + R2
    M1, M2 = e**2*m1, e**2*m2
    f = {1: ((al*M1**2*M2 + be*M1*M2**2)/D**3, (ga*M1**2*M2 + de*M1*M2**2)/D**4),
         2: ((al*M2**2*M1 + be*M2*M1**2)/D**3, -(ga*M2**2*M1 + de*M2*M1**2)/D**4)}
    return {b: (8*f[b][0], 8*f[b][1]) for b in (1, 2)}


def newtonian_map(m1):
    """The x^3 change in z1 from the alpha..delta terms, worked out by hand at Newtonian order (independent check)."""
    m2 = 1 - m1; nu = m1*m2
    return sp.expand(-sp.Rational(4, 3)*m2*nu*(ga + de)
                     - m2**2*(sp.Rational(8, 3)*nu*(ga + de) - 4*m1*(ga*m1 + de*m2))
                     - 4*m1*m2*(al*m1 + be*m2))


def main(fresh=False):
    out = {'ratios': RATIOS, 'f_values': {str(k): str(v) for k, v in F_VALUES.items()}, 'per_ratio': {}}
    t0 = time.time()
    for r in RATIOS:
        m1 = sp.Rational(r); m2 = 1 - m1
        ex = closed_pieces(r, fresh)
        tr = no.truth(m1)[3]
        miss_rules = sp.nsimplify(no.mismatch_formula(m1))
        cs_closed = no.redshift(m1, 6, res_hook=make_hook(ex))
        miss_closed = sp.simplify(cs_closed[3] - tr)
        cs_f = no.redshift(m1, 6, res_hook=make_hook(ex, fextra=f_terms(m1, m2)))
        miss_sym = sp.expand(cs_f[3] - tr)
        lin = sp.expand(miss_sym - miss_closed)
        miss_f = sp.nsimplify(miss_sym.subs(F_VALUES))
        point = {al: F_VALUES[al], be: F_VALUES[be] + sp.Rational(1, 4), ga: F_VALUES[ga], de: F_VALUES[de] - sp.Rational(1, 4)}
        out['per_ratio'][r] = {
            'truth_x3': str(tr), 'miss_rules': str(miss_rules), 'miss_closed': str(miss_closed),
            'miss_with_f_symbolic': str(miss_sym), 'miss_with_f': str(miss_f),
            'miss_with_point_term': str(sp.nsimplify(miss_sym.subs(point))),
            'newtonian_map_matches': bool(sp.expand(lin - newtonian_map(m1)) == 0)}
        print(r, 'rules', miss_rules, '| +closed', miss_closed, '| +f', miss_f, '| +point',
              out['per_ratio'][r]['miss_with_point_term'], '| map ok', out['per_ratio'][r]['newtonian_map_matches'],
              f'({time.time()-t0:.0f}s)', flush=True)
    # fit the closed-stage miss as nu * quadratic on the first three ratios, test on the other three
    c0, c1, c2 = sp.symbols('c0 c1 c2')
    fit_on = RATIOS[:3]
    eqs = [sp.Eq(sp.Rational(r)*(1 - sp.Rational(r))*(c0 + c1*sp.Rational(r) + c2*sp.Rational(r)**2),
                 sp.sympify(out['per_ratio'][r]['miss_closed'])) for r in fit_on]
    sol = sp.solve(eqs, [c0, c1, c2])
    closed_formula = sp.expand(m*(1 - m)*(sol[c0] + sol[c1]*m + sol[c2]*m**2))
    resid_formula = m*(1 - m)**2*(3 - m)/3
    out['closed_formula'] = str(sp.factor(closed_formula))
    out['closed_formula_fit_on'] = fit_on
    out['closed_formula_holds'] = {r: bool(sp.simplify(closed_formula.subs(m, sp.Rational(r)) -
                                                       sp.sympify(out['per_ratio'][r]['miss_closed'])) == 0) for r in RATIOS}
    out['residual_formula'] = 'm1*(1-m1)**2*(3-m1)/3'
    out['residual_holds'] = {r: bool(sp.simplify(resid_formula.subs(m, sp.Rational(r)) -
                                                 sp.sympify(out['per_ratio'][r]['miss_with_f'])) == 0) for r in RATIOS}
    # residual check by algebra: closed formula + Newtonian map at the computed f-values
    out['residual_by_algebra'] = str(sp.factor(sp.expand(closed_formula + newtonian_map(m).subs(F_VALUES))))
    out['point_term'] = 'X-hat at body a gains (1/4) m_a m_b^2 / (r12^2 r_b): beta -> 5/4, delta -> -9/4'
    out['point_term_closes'] = {r: out['per_ratio'][r]['miss_with_point_term'] == '0' for r in RATIOS}
    json.dump(out, open(os.path.join(HERE, 'fourth_rule.json'), 'w'), indent=1)
    print('closed-stage formula', out['closed_formula'], out['closed_formula_holds'])
    print('residual', out['residual_formula'], out['residual_holds'], 'by algebra', out['residual_by_algebra'])
    print('point term closes', out['point_term_closes'])
    return out


if __name__ == '__main__':
    main(fresh='--fresh' in sys.argv)
