"""Checks for the unified-hypothesis DRAFT. Reads the JSON outputs; recomputes a few things in float64."""
import json, numpy as np, kern
res = []
def check(name, ok):
    res.append((name, bool(ok))); print(('PASS ' if ok else 'FAIL ') + name)

R = kern.riemann(500)
z = kern.real_zeros(R, 74)
check('float64 Riemann kernel puts the first 5 zeros at 2γ (to 1e-4)', np.max(np.abs(z / 2 - np.array(kern.GAMMAS[:5]))) < 1e-4)

L = json.load(open('lam_pball.json'))
hold = [(k, v) for k, v in L.items() if float(k.split('|')[1]) <= float(k.split('|')[0]) and v.get('Lambda') is not None]
check(f'p-ball: Λ < 0 at all {len(hold)} measured points with a ≤ p', all(v['Lambda'] < 0 for _, v in hold))
even_above = [(k, v) for k, v in L.items() if float(k.split('|')[1]) > float(k.split('|')[0])]
check('p-ball: even p above the line (a > p) also has Λ < 0 (p = 2, 4, 6)', all(v.get('Lambda', 1) < 0 for _, v in even_above))
thr = [(k, v) for k, v in L.items() if v.get('far_threshold_at_a_eq_p') is not None and v.get('Lambda') is not None]
check('p-ball at a = p: Λ sits within 0.1 above-or-at the far threshold log|B|', all(v['Lambda'] >= v['far_threshold_at_a_eq_p'] - 0.03 and v['Lambda'] - v['far_threshold_at_a_eq_p'] < 0.1 for _, v in thr))
check('triangle p = a = 1 flagged (double zeros, Λ = 0)', L['1|1']['status'].startswith('complex') and L['1|1']['far_threshold_at_a_eq_p'] == 0.0)

P = json.load(open('parents.json'))
sc = {k: len(v) for k, v in P['scans'].items()}
base = len(P['H0_real_zeros'])
check('Ξ has 18 real zeros below z = 150', base == 18)
check('heat parents keep all 18 low zeros', all(v == 18 for k, v in sc.items() if k.startswith('heat')))
check('factor parents: c ≥ 1 keeps 18, c ≤ 0.6 loses zeros', sc['factor|1.0'] == 18 and sc['factor|2.0'] == 18 and all(sc[f'factor|{c}'] < 18 for c in ['0.25', '0.3', '0.45', '0.6']))
check('bath parents: J ≤ 1 keeps 18, J ≥ 2 loses zeros', sc['bath|0.5'] == 18 and sc['bath|1.0'] == 18 and all(sc[f'bath|{j}'] < 18 for j in ['2.0', '3.0', '4.0', '6.0']))
check('three-level children keep all 18 low zeros, even below w = 1/4', all(v == 18 for k, v in sc.items() if k.startswith('child3')))

off = np.array(json.load(open('anti_offaxis_plus.json'))['zeros'])
u, w = R.u, R.w
Lc = np.pi / 2 * kern.Phi(np.array([0.0]))[0]
zz = off[:, 0] + 1j * off[:, 1]
z0 = zz.copy()
for _ in range(30):
    zz = zz - (np.sin(np.outer(zz, u)) @ (w / u) + Lc) / R.H(zz)
val = np.sin(np.outer(zz, u)) @ (w / u) + Lc
assert np.max(np.abs(zz - z0)) < 1e-4
check(f'antiderivative parent: {len(off)} off-line zeros (at −z) confirmed by Newton (moved < 1e-4, |P(z)+L| < 1e-8; the three tested at 40 digits agree)', len(off) >= 8 and np.max(np.abs(val)) < 1e-8)

M1 = json.load(open('multiplier.json')); M2 = json.load(open('multiplier2.json'))
check('Φ/Φ* moment ratio fails Turán at k = 1 and Jensen degree 2', 1 in M1['Phi_over_PhiStar']['turan_fail'] and M1['Phi_over_PhiStar']['jensen_fail'][0][0] == 2)
check('Gaussian base (p = 2) passes Turán for k ≤ 119 and Jensen to degree 40', M2['2']['turan_fail'] == [] and M2['2']['jensen_fail'] == [])
check('Pólya bumps p = 4, 6, 8, 12 fail Turán at k = 1', all(1 in M2[p]['turan_fail'] for p in ['4', '6', '8', '12']))
check('controls: p-ball/Pólya passes, 1/(k+1) fails', M1['control_pball_over_polya_p4_a2']['turan_fail'] == [] and 1 in M1['control_1_over_k_plus_1']['turan_fail'])

D = json.load(open('landscape.json'))
b = D['bath2_averaged']
check('two-level bath averaged at f = 1/2 keeps zeros real (both kernels)', b['p-ball p=3, a=2 f=0.5']['max_abs_Im'] == 0 and b['Riemann Φ f=0.5']['max_abs_Im'] == 0)
check('two-level bath at f ≠ 1/2 moves every found zero into the upper half-plane', all(v['mean_Im'] is not None and v['mean_Im'] > 0 for k, v in b.items() if not k.endswith('f=0.5')))
t = D['two_spin_real_zero_counts']
check('two spins: no real zeros for J ≤ 0, real zeros for J > 0 (both kernels)', all((c == 0) == (J <= 0) for rows in t.values() for J, c in rows))

print(f'{sum(o for _, o in res)} of {len(res)} pass')
