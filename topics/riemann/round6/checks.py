"""Round 6 checks: each reads the saved results and states what must hold."""
import json, numpy as np
LZ = json.load(open('lattice_zeros.json')); EP = json.load(open('euler_pairs.json'))
LP = json.load(open('lattice_pairs.json')); RB = json.load(open('robustness.json'))
out = []
def check(name, ok, detail=''):
    out.append((name, bool(ok), detail)); print(('PASS' if ok else 'FAIL'), name, detail)
Z = lambda run, p: np.array(LZ[run][p]['zeros'])
z = Z('square', '2.0')
check('C1 square p = 2 (zeta L(chi_-4)): every zero on Re 1/2, first zeros 6.0209 and 14.1347',
      np.max(np.abs(z[:, 0] - 0.5)) < 1e-6 and abs(z[0, 1] - 6.0209489) < 1e-5 and np.min(np.abs(z[:, 1] - 14.1347251)) < 1e-5)
z = Z('square', 'inf')
check('C2 square p = inf (8 zeta(2s-1)): every zero on Re 3/4, first at height 14.1347/2',
      np.max(np.abs(z[:, 0] - 0.75)) < 1e-6 and abs(z[0, 1] - 14.1347251/2) < 1e-5)
z2, z4 = Z('hexagon', '2.0'), Z('hexagon', '4.0')
check('C3 hexagon face power: p = 4 gives exactly the p = 2 zeros, all on Re 1/2',
      len(z2) == len(z4) and np.max(np.abs(z2 - z4)) < 1e-6 and np.max(np.abs(z4[:, 0] - 0.5)) < 1e-6)
check('C4 hexagon p = inf (6 zeta(2s-1)): every zero on Re 3/4', np.max(np.abs(Z('hexagon', 'inf')[:, 0] - 0.75)) < 1e-6)
mid = [k for k in LZ['square'] if k not in ('2.0', 'inf')]
check('C5 square, every p strictly between 2 and inf: some zero at least 0.1 off Re 1/2',
      all(np.max(np.abs(Z('square', k)[:, 0] - 0.5)) > 0.1 for k in mid), f'{len(mid)} values of p')
midh = [k for k in LZ['hexagon'] if k not in ('2.0', '4.0', 'inf')]
check('C6 hexagon, every other p: some zero off Re 1/2 by more than 0.02',
      all(np.max(np.abs(Z('hexagon', k)[:, 0] - 0.5)) > 0.02 for k in midh), f'{len(midh)} values of p')
check('C7 ratio dial: the hexagon lattice scaled by 1.7 has the same zeros at p = 3',
      np.max(np.abs(Z('hexagon_scaled_1.7', '3.0') - Z('hexagon', '3.0'))) < 1e-6)
fc = RB['fit_change']
check('C8 changing the small-x fit and x0 moves no zero by more than 1e-5 (square 2.1, 2.5; hexagon 3)',
      all(v['max_dist'] < 1e-5 and v['n1'] == v['n2'] for v in fc.values()))
r20 = RB['p2_curvature_exact_check']['disc-20_square']
check('C9 x^2 + 5y^2: engine zeros agree with the exact genus formula zeta L(chi_-20) + L(chi_-4) L(chi_5) to 1e-5',
      max(r['diff'] for r in r20) < 1e-5)
off = [r['genus_formula'] for r in r20 if abs(r['genus_formula'][0] - 0.5) > 1e-4]
check('C10 x^2 + 5y^2 has a zero off the line at 0.93297 + 15.66825i', len(off) == 1 and abs(off[0][0] - 0.93297) < 1e-4 and abs(off[0][1] - 15.66825) < 1e-4)
r23 = RB['p2_curvature_exact_check']['disc-23_voronoi']['rows']
check('C11 skew lattice, p = 2: engine zeros agree with Epstein exact theta inversion to 1e-5, and two have Re s > 1',
      max(r['diff'] for r in r23) < 1e-5 and sum(r['exact'][0] > 1 for r in r23) == 2)
cr = RB['curvature_refit']
check('C12 curvature runs at p = 3, 4 are stable under a refit (moves < 2e-4); p = inf is not and is set aside',
      all(cr[n][p]['vs_first_run']['max_dist'] < 2e-4 for n in cr for p in ('3.0', '4.0')) and
      any(cr[n]['inf']['vs_first_run']['max_dist'] > 1e-2 for n in cr))
check('C13 odds and ratio scaling: Delta -> t^2 Delta (sign kept), point stays on its parabola d = c k^2',
      EP['td_disc_scaling'] == '-t**2*(4*d - k**2)' and EP['td_parabola'] == '0')
check('C14 curve 37a with odds t: the root line runs 1/2 + log t / log p, spread > 0.25 for t = 1.25 and 0.8',
      EP['td_root_lines']['1.25']['max'] - EP['td_root_lines']['1.25']['min'] > 0.25 and EP['td_root_lines']['0.8']['max'] - EP['td_root_lines']['0.8']['min'] > 0.25
      and abs(EP['td_check_p101']['abs_roots'][0] - EP['td_check_p101']['predicted']) < 1e-12)
sd = EP['sd_detector_matches_delta_negative']
check('C15 height on a degree-2 local factor: fires exactly when Delta < 0 (all grid points), and at every good prime of 37a',
      int(sd['agree']) == sd['points'] and EP['sd_37a']['detector_fires'] == EP['sd_37a']['primes'] and EP['sd_hstar_formula'] == '(-w2 + x**2)/(2*x**2)')
kq = EP['kq_identity_s2']
check('C16 k + q identity at s = 2: the direct sum closes on it, gap halving per doubling of N',
      all(1.6 < g['gaps'][0]/g['gaps'][1] < 2.4 and 1.6 < g['gaps'][1]/g['gaps'][2] < 2.4 and g['gaps'][2] > 0 for g in kq.values()))
ok = True
for key, v in EP['kq_leading_power'].items():
    rows = np.array(v['rows']); fit = np.polyfit(1/np.log(1/rows[:, 0]), rows[:, 1], 1); ok &= abs(fit[1] - v['predicted']) < 0.01
check('C17 k + q leading power: log Z / log zeta extrapolates to 1/(1 - q) within 0.01 for two k', ok)
kr = EP['kr']
check('C18 k + r: zeta zeros kept exactly while the wall is left of 1/2 (|Z| near the first zero = |zeta prime H| * 1e-6)',
      all(abs(v['abs_Z_at_rho1_pm_1e-6'][0]/v['abs_zeta_prime_times_1e-6'] - 1) < 1e-3 for v in kr.values() if 'abs_Z_at_rho1_pm_1e-6' in v)
      and kr['k=1.5,rho=1.3']['wall'] > 0.5)
hc = EP['height_collapse_on_kr']; sf = LP['sf_height_a_points']
check('C19 height collapse: the singular points close on the zeros only as s -> 1 (zeta H and two round lattices)',
      abs(hc['0.8']['off_line']) > abs(hc['0.95']['off_line']) > 0.01 and
      all(sf[k]['0.8']['max_off_line'] > sf[k]['0.95']['max_off_line'] > 0.01 and sf[k]['1.0']['max_off_line'] < 1e-9 for k in ('square_p2.0', 'hexagon_p4.0')))
m = EP['multiplicativity']
check('C20 Euler product: representation counts multiplicative at p = 2 (square, hexagon), not at p = 4, 6 on the square or p = 6 on the hexagon',
      m['square_p2 (a^2+b^2)']['failures'] == 0 and m['hex_p2 (a^2+ac+c^2)']['failures'] == 0 and
      all(m[k]['failures'] > 0 for k in ('square_p4 (a^4+b^4)', 'square_p6 (a^6+b^6)', 'hex_p6 ((a+c)^6+a^6+c^6)/2')))
check('C21 design rule: face power is round at even p <= 2n - 2 (square 2; hexagon 2, 4; octagon 2, 4, 6); hexagon p = 4 is the square of p = 2',
      EP['face_power_round_at'] == {'square': [2], 'hexagon': [2, 4], 'octagon': [2, 4, 6]} and EP['hex_p4_is_square_of_p2'] == '0')
lw = EP['light_wick_hexagon']
check('C22 light: Wick-rotated hexagon form is real at every even p, a power of x^2 - y^2 only at p = 2, 4',
      all(lw[p]['real'] for p in lw) and lw['2']['equals_power_of_x2_minus_y2'] and lw['4']['equals_power_of_x2_minus_y2']
      and not lw['6']['equals_power_of_x2_minus_y2'] and not lw['8']['equals_power_of_x2_minus_y2'])
dq = EP['dq_prime_fractal_dimension']
check('C23 order dial on Q reproduces round 4 (D = 1.39943), and each quadratic field has its own wall', abs(dq['Q'] - 1.39943) < 1e-4 and len(set(round(v, 3) for v in dq.values())) == 4)
fl = LP['flux_discrepancy']
pred = {'square_p2.0': 0.5, 'square_p3.0': 2/3, 'square_p4.0': 0.75, 'square_pinf': 1.0, 'hexagon_p4.0': 0.5, 'hexagon_pinf': 1.0}
check('C24 flux: discrepancy exponents within 0.06 of 1/2 (round), 1 - 1/p (square), 1 (polytope)',
      all(abs(fl[k]['rms_exponent'] - v) < 0.06 for k, v in pred.items()), str({k: round(fl[k]['rms_exponent'], 3) for k in pred}))
n = sum(o[1] for o in out); print(f'{n} of {len(out)} checks pass')
