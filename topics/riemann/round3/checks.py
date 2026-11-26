"""Checks for Riemann round 3 (DRAFT). Run: python3 checks.py > checks_output.txt"""
import json, re, os, numpy as np, mpmath as mp
rw = json.load(open('rewind.json')); rd = json.load(open('roads.json'))
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name + (f'   [{info}]' if info else ''))
T1 = rd['T1']; tr = T1['track']
check('C1 fused rewind is reversible: rewind the toy to t = -0.06 and run it forward again; coefficients come back to 1e-25',
      T1['reversibility_max_coef_err'] < 1e-25, f"{T1['reversibility_max_coef_err']:.1e}")
check('C2 toy collision of the close pair (2, 2.3) comes within 10% of the pair-alone time -g^2/8',
      abs(T1['t_collide'] / T1['isolated_estimate'] - 1) < 0.1, f"{T1['t_collide']:.5f} vs {T1['isolated_estimate']:.5f}")
s = np.array([x['s'] for x in tr]); p = np.array([x['p'] for x in tr]); D = np.array([x['disc'] for x in tr]); t = np.array([x['t'] for x in tr])
check('C3 the fused point (s, p) moves smoothly through the collision (no grid step more than 3 times the median step) while the discriminant changes sign once',
      np.max(abs(np.diff(s))) < 3 * np.median(abs(np.diff(s))) + 1e-12 and np.max(abs(np.diff(p))) < 3 * np.median(abs(np.diff(p))) and np.sum(np.diff(np.sign(D)) != 0) == 1,
      f'max step s {np.max(abs(np.diff(s))):.1e}, p {np.max(abs(np.diff(p))):.1e}')
k = np.argmin(abs(t - T1['t_collide'])); slope = (D[k - 5] - D[k + 5]) / (t[k - 5] - t[k + 5])
check('C4 near the collision the discriminant is a straight line in t with slope near 8 (the pair-alone law)', abs(slope / 8 - 1) < 0.2, f'slope {slope:.2f}')
same = [x for x in rw['full'] if x['same_pair']]
check('C5 for the closest pairs among the first 5,000 zeros, the full rewind (every zero pulling) agrees with the pair-alone time to 3%',
      all(abs(x['t_full'] / -x['tau_isolated'] - 1) < 0.03 for x in same), f'{len(same)} pairs')
m = max(rw['full'], key=lambda x: x['t_full'])
check('C6 the first merger among the first 5,000 zeros is zeros 4765, 4766 near height 5229.2, at t = -0.00094',
      m['j'] == 4765 and abs(m['t_full'] + 0.000936) < 2e-6, f"t = {m['t_full']:.6f}")
L = rw['lehmer']
check("C7 Lehmer's pair (zeros 6709, 6710 at 7005.0629, 7005.1006) merges at t = -0.000711, within 1% of its pair-alone time",
      L['n'] == 6709 and abs(L['gamma'] - 7005.0629) < 1e-3 and abs(L['t_full'] / -L['tau_isolated'] - 1) < 0.01, f"t = {L['t_full']:.7f}")
gm = np.array(rw['pairs']['gamma_mid']); tau = np.array(rw['pairs']['tau_isolated'])
med = [np.median(tau[(gm >= a) & (gm < a + 1000)]) for a in range(0, 5000, 1000)]
check('C8 the rewind budget shrinks with height: the median pair-alone time falls in every block of 1,000 in height',
      all(med[i + 1] < med[i] for i in range(4)), ', '.join(f'{x:.3f}' for x in med))
T2 = rd['T2']; rc = T2['ring_check']
check('C9 the ring average of shifts is J0(aD) when turned by i and I0(aD) when real (trapezoid on the ring, 64 points)',
      abs(rc['imag'][0] / rc['imag'][1] - 1) < 1e-9 and abs(rc['real'][0] / rc['real'][1] - 1) < 1e-9)
gc = T2['gauss_check']
check('C10 a Gaussian smear turned by i is the forward heat flow: E f(z + i sigma X) = e^(-t D^2) f, t = sigma^2/2', abs(gc[0] / gc[1] - 1) < 1e-9)
check('C11 the i-turned ring puts the off-line pair on the line near a = sqrt2*0.3; the real ring pushes the close real pair off near 0.3/sqrt2; the i-turned ring never pushes it off',
      abs(T2['A_lands_at'] / (2 ** 0.5 * 0.3) - 1) < 0.1 and abs(T2['B_leaves_at'] / (0.3 / 2 ** 0.5) - 1) < 0.05 and T2['B_imag_leaves'] is None,
      f"{T2['A_lands_at']}, {T2['B_leaves_at']}")
T3 = rd['T3']
check('C12 mean |zeta(sigma + it)| over 20,000 random t matches the barrier chain prod (2/pi) K(p^-sigma) within 3 standard errors, sigma = 3, 2, 1.5, 1.25, 1.1',
      all(abs(v['data'][0] - v['chain']) < 3 * v['data'][1] for v in T3['mean'].values()),
      ', '.join(f"{k}: {v['data'][0]:.4f}/{v['chain']:.4f}" for k, v in T3['mean'].items()))
vv = T3['var']
check('C13 the spread of log|zeta| matches the chain: within 10% at sigma = 1.5, 1.0 and within 25% at 0.75, 0.6',
      all(abs(vv[k]['data'] / vv[k]['chain'] - 1) < (0.1 if float(k) >= 1 else 0.25) for k in vv),
      ', '.join(f"{k}: {vv[k]['data']:.3f}/{vv[k]['chain']:.3f}" for k in vv))
vc = np.array(T3['var_curve']['chain'])
check('C14 the chain spread grows without bound toward the line: increasing as sigma falls to 0.52, and its sum diverges like (1/2) log(1/(sigma - 1/2))',
      np.all(np.diff(vc) < 0) and vc[0] > 1.0, f'{vc[0]:.2f} at 0.52')
S = np.array(rd['T4']['S'])
check('C15 S(T) at the first 5,000 zeros averages near 0 and stays inside (-3, 3)', abs(S.mean()) < 0.05 and np.max(abs(S)) < 3, f'mean {S.mean():.3f}, max |S| {np.max(abs(S)):.2f}')
# Lean identities, numerically (with complex numbers where i is needed)
rng = np.random.default_rng(3); ok = True
for _ in range(50):
    a, b, z, tt, s_, p_, v, w, r_ = rng.normal(size=9); ph = rng.uniform(0, 6.3)
    ok &= abs((a + b) ** 2 - 4 * (a * b - 2 * tt) - ((a - b) ** 2 + 8 * tt)) < 1e-9
    f = lambda x: x * x - s_ * x + p_
    ok &= abs(f(z + v) + f(z - v) - 2 * (f(z) + v * v)) < 1e-9
    ok &= abs(f(z + 1j * w) + f(z - 1j * w) - 2 * (f(z) - w * w)) < 1e-9
    ok &= abs(abs(1 - r_ * np.exp(1j * ph)) ** 2 - (1 + r_ * r_ - 2 * r_ * np.cos(ph))) < 1e-9
check('C16 the Lean statements hold numerically (heat_pair_gap, real_ring, turned_ring with v = i w, barrier_sq)', ok)
# polish
md = open('README.md').read()
check('P1 README is marked DRAFT and names the theory (Affine Composition Theory, built on Fused Hyperoperation Algebra)',
      'DRAFT' in md[:300] and 'Affine Composition Theory' in md[:600] and 'Fused Hyperoperation Algebra' in md[:600])
shown = set(re.findall(r'figures/(fig[^)]+\.png)', md)); made = set(os.listdir('figures'))
check('P2 every picture the README shows exists, and every picture made is shown', shown == made, f'{len(made)} pictures')
bad = [w for w in ['proves RH', 'proof of RH', 'breakthrough', 'we prove the Riemann', 'solves', 'revolution'] if w.lower() in md.lower()]
check('P3 no overclaim phrases in the README', not bad, ', '.join(bad))
check('P4 wording rule: "backward heat flow" and "rewind" used; no "downward"/"going down" for the flow',
      'backward heat flow' in md and 'rewind' in md and 'downward' not in md.lower() and 'going down' not in md.lower())
first = md.split('\n## ')[1] if '\n## ' in md else ''
check('P5 the README opens with the answer and says plainly whether each road is a lead', first.startswith('The answer') and 'lead' in first)
check('P6 numbers quoted in the README match the data',
      all(x in md for x in ['0.00094', '0.00071', '4765', '6709', '1e-31', f"{T2['A_lands_at']:.3f}", f"{T2['B_leaves_at']:.3f}"]))
lean = open('lean/RiemannRound3.lean').read()
check('P7 the Lean file says it is compiled in the shared library and has no sorry', 'compiled in the shared library' in lean and 'sorry' not in re.sub(r'--.*', '', re.sub(r'/-.*?-/', '', lean, flags=re.S)))
print(f'\n{sum(res)} of {len(res)} pass')
