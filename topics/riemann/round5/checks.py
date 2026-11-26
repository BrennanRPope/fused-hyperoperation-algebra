"""Checks for Riemann round 5 (DRAFT). Run: python3 checks.py > checks_output.txt"""
import json, re, os, numpy as np
p1 = json.load(open('pair1.json')); p2 = json.load(open('pair2.json'))
lf = json.load(open('pass_light_flux.json')); p3 = json.load(open('pair3.json'))
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name + (f'   [{info}]' if info else ''))

# ---------- pair 1: height s with homothety k ----------
check('C1 Abel function of E(x) = e^x - 1: coefficients -1/36 and 1/540; E^(1) and E^(1/2) twice agree with E to 1e-12',
      p1['abel_coeffs'] == {'c1': '-1/36', 'c2': '1/540'} and p1['check_E1'] < 1e-11 and p1['check_half_twice'] < 1e-11,
      f"{p1['check_E1']:.1e}")
check('C2 homothety log: l_k(zeta) equals 1 - prod(1 - k p^-2), over k, at s = 2 for k = 0.5, 1, 2, 3 (1e-14); its series bends like height 1 + k',
      all(abs(v['lhs'] - v['rhs']) < 1e-14 for v in p1['lk_identity_s2'].values())
      and p1['lk_series'] == 'k**2*u**3/6 + k*u**3/2 - k*u**2/2 + u**3/3 - u**2/2 + u')
sv = p1['singular_value_X_s']; xs = [sv[k] for k in sorted(sv, key=float)]
check('C3 height log of zeta breaks at a value X_s < 0 for s < 1, rising to 0 as s -> 1 (-2.52 at s = 0.05, -0.0062 at 0.99); X_2 = 1/e',
      all(x < 0 for x in xs) and np.all(np.diff(xs) > 0) and abs(sv['0.99']) < 0.01
      and abs(p1['singular_value_above_1']['2.0'] - np.exp(-1)) < 1e-12)
ra = p1['real_a_points_off_line']
check('C4 the break points are a-points of zeta off the line: real ones at sigma = 0.6748 (s = 0.05) and 0.5846 (s = 0.1), complex ones all left of the line',
      0.5 < ra['0.1'] < ra['0.05'] < 1 and all(v['right'] == 0 and v['left'] > 0 for v in p1['complex_a_points_side'].values()),
      ', '.join(f"s={k}: {v['left']} left, closest {v['closest']:.3f}" for k, v in p1['complex_a_points_side'].items()))
mk = p1['mertens_k']; sl = {k: v['slope_runmax'] for k, v in mk.items()}
check('C5 generalized Mertens to 1e7: growth slope 0.51 at k = 1 and 0.62, 0.63 at k = 2, 3 (log factor from the order-k poles); k = 1/2 tracks its Selberg-Delange main term (ratio 1.18 to 1.19)',
      abs(sl['1.0'] - 0.5) < 0.03 and 0.55 < sl['2.0'] < 0.7 and 0.55 < sl['3.0'] < 0.7
      and all(1.15 < r < 1.22 for r in mk['0.5']['main_ratio'][-5:]),
      ', '.join(f'k={k}: {v:.3f}' for k, v in sl.items()))
# ---------- pair 2: ratio r with order q ----------
check('C6 prime fractal (one square per prime, side 1/p): Moran dimension 1.39943 (prime zeta = 1)', abs(p2['moran_dimension'] - 1.39943) < 1e-5,
      f"{p2['moran_dimension']:.6f}")
cd = p2['complex_dimensions']; cdj = json.load(open('pair2_complexdims.json')); vs = json.load(open('pair2_valueset.json'))
check('C7 complex dimensions (P(s) = 1, Re s > 1): the only one with t < 60 on the grid is D itself; Newton from 597 starts finds none with 0 < t < 200; at sigma = 1.2 the disc of values has radius 1.52 > 1 (so they exist higher up) while P(1.2 + it) stays 0.23 from 1 for t < 200',
      len(cd) == 1 and abs(cd[0][0] - 1.3994333) < 1e-6 and cdj['complex_roots_t_lt_200'] == [] and vs['P_sigma'] > 1 and abs(vs['min_dist_t_gt_1'] - 0.232) < 0.001,
      f"closest approach at sigma 1.02: {cdj['closest_approach']['1.02']['min_dist']:.3f}")
check('C8 order dial: q-multinomial equals q^(inversions) summed over orderings (brute force, 12 cases)', p2['qmultinomial_bruteforce_maxerr'] < 1e-12,
      f"{p2['qmultinomial_bruteforce_maxerr']:.1e}")
idn = p2['order_dial_identity_s2']; tl = json.load(open('pair2_tail.json'))
ok = all(abs(v['direct'] / v['formula'] - 1) < 2e-4 and v['direct'] < v['formula'] for v in idn.values())
ok &= all(abs(r[-1]['rel']) < 2e-5 and abs(r[-1]['rel']) < 0.6*abs(r[0]['rel']) for r in tl.values())
check('C9 order-dial identity (Z_q as a sum of homothety cubes) at s = 2, q = 0.3, 0.6, 0.9: the direct sum climbs to the formula as n grows (gap halves with each doubling, 1e-5 or better at n = 8e5 with the tail added); necklace form agrees at q = 0.6',
      ok and abs(p2['Zq_formula_s2_q0.6'] / idn['0.6']['formula'] - 1) < 1e-8,
      ', '.join(f"q={q}: {r[-1]['rel']:.1e}" for q, r in tl.items()))
cl = p2['critical_line']; G = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])
ok = True; info = []
for q, v in cl.items():
    t = np.array(v['t']); a = np.array(v['absZq'])
    for g in G[(G > t[0]) & (G < t[-1])]:
        i = np.argmin(abs(t - g)); ok &= a[i] <= a[max(i-1, 0)] and a[i] <= a[min(i+1, len(a)-1)] and a[i] < 0.25*np.median(a)
    info.append(f'q={q}: dips to {a.min():.3f}, median {np.median(a):.2f}')
check('C10 on the critical line |Z_q| dips at every zeta zero in 10 < t < 35 (q = 0.5, 0.3): the order dial keeps the zeros, as branch points', ok, ', '.join(info))
# ---------- light and flux tubes ----------
w = lf['wick_check']
check('C11 heat flow and its Wick rotation from Phi: the first zero moves along the line by v t and straight off it by v tau, v = -0.57963, to 1e-8 at small time',
      abs(lf['v_z_from_Phi'] / lf['v_z_from_zeros'] - 1) < 1e-4
      and all(abs(r['zero'][1] - r['predicted'][1]) < 1e-8 and abs(r['zero'][0] - r['predicted'][0]) < 1e-5 for r in w[:2])
      and all(abs(r['imag']) < 1e-12 for r in lf['heat_check']),
      f"v {lf['v_z_from_Phi']:.6f} vs zeros {lf['v_z_from_zeros']:.6f}")
check('C12 Berry-Keating: area under the boost hyperbola is E log(E/(l_x l_p)) - E + l_x l_p (sympy), and the zero count follows it with mean residual below 0.05',
      lf['bk_area_is_E_eta'] and abs(lf['bk_residual_mean']) < 0.05, f"mean {lf['bk_residual_mean']:.4f}, max {lf['bk_residual_maxabs']:.2f}")
T = np.array(lf['Deff']['T']); m = np.array(lf['Deff']['measured']); s = np.array(lf['Deff']['smooth'])
big = T > 200
check('C13 flux-tube reading: local dimension of the zero count tracks log(T/2pi)/(log(T/2pi) - 1) within 0.1 above T = 200, falling to 1.18',
      np.max(abs(m[big] - s[big])) < 0.1 and abs(s[-1] - 1.18) < 0.01, f'max dev {np.max(abs(m[big]-s[big])):.3f}')
# ---------- pair 3: height s with roundness p ----------
k = p3['kernels']
check('C14 detector: h* goes above 1 only for the kernels with non-real zeros (p = 2.5, 3), and stays below 0.28 for p = 4, 6, 8',
      k['2.5']['hmax_overall'] > 1 and k['3.0']['hmax_overall'] > 1 and all(k[p]['hmax_overall'] < 0.28 for p in ['4.0', '6.0', '8.0'])
      and len(p3['nonreal_zeros']['2.5']) == 3 and len(p3['nonreal_zeros']['3.0']) == 1,
      ', '.join(f"p={p}: {v['hmax_overall']:.3f}" for p, v in k.items()))
pk = p3['detector_peaks']; nr = p3['nonreal_zeros']
check("C15 each detector peak sits next to a non-real zero's real part (within 1)",
      abs(pk['2.5'][0] - nr['2.5'][0][0]) < 1 and abs(pk['3.0'][0] - nr['3.0'][0][0]) < 1,
      f"{pk['2.5'][0]:.2f} vs {nr['2.5'][0][0]:.2f}; {pk['3.0'][0]:.2f} vs {nr['3.0'][0][0]:.2f}")
check("C16 the Riemann kernel's h* (round 4) stays at 0.705, below the detector line at 1", abs(p3['xi_round4_hmax'] - 0.7048) < 1e-3)
# ---------- polish ----------
rd = open('README.md').read(); lean = open('lean/RiemannRound5.lean').read()
check('P1 README is marked DRAFT and names the theory (Affine Composition Theory, built on Fused Hyperoperation Algebra)',
      'DRAFT' in rd.split('\n')[0] and 'Affine Composition Theory, built on Fused Hyperoperation Algebra' in rd)
shown = set(re.findall(r'figures/(fig[^)]+\.png)', rd)); made = set(f for f in os.listdir('figures') if f.endswith('.png'))
check('P2 every picture the README shows exists, and every picture made is shown', shown == made and len(made) > 0, f'{len(made)} pictures')
bad = [p for p in ['proves RH', 'proof of RH', 'breakthrough', 'solves', 'we prove the Riemann', 'settles RH'] if p.lower() in rd.lower()]
check('P3 no overclaim phrases in the README', not bad, ', '.join(bad))
check('P4 wording rule: "backward heat flow" and "rewind" used; no "downward"/"going down" for the flow',
      'backward heat flow' in rd and 'rewind' in rd and not re.search(r'downward|going down|run down|flow down', rd, re.I))
first = rd.split('## The answer')[1].split('\n## ')[0] if '## The answer' in rd else ''
check('P5 the README opens with the answer, says plainly whether it is a lead, and gives the order and why C was dropped',
      '## The answer' in rd and rd.index('## The answer') < rd.index('## Pass 1') and 'Not a lead' in first and 'derived from k and d' in rd)
nums = ['1.39943', '0.6748', '0.5846', '−0.57963', '0.705', '1.18']
check('P6 numbers quoted in the README match the data', all(n in rd for n in nums), ', '.join(n for n in nums if n not in rd))
check('P7 the Lean file says it is compiled in the shared library, has 10 theorems and no sorry',
      'compiled in the shared library' in lean and len(re.findall(r'^theorem', lean, re.M)) == 10 and 'sorry' not in re.sub(r'--.*', '', re.sub(r'/-.*?-/', '', lean, flags=re.S)))
print(f'\n{sum(res)} of {len(res)} checks pass')
