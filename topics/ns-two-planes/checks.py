"""DRAFT. Checks for ns-two-planes, from algebra_results.json, twist_sectors.json and runs/*.json.
Run: python3 checks.py  -> checks_output.txt, checks_results.json (a few seconds; run algebra.py, batches first)."""
import json, os, numpy as np
L = lambda f: json.load(open(f))
res = []
def check(name, ok, note):
    res.append(dict(name=name, ok=bool(ok), note=note))
alg = L('algebra_results.json')
check('C1 exact algebra (algebra.py)', all(v['pass'] for v in alg.values()), f"{sum(v['pass'] for v in alg.values())} of {len(alg)} pass")
runs = lambda nu, m, k='1': [L(f'runs/n64_nu{nu}_{m}_k{k}_s{s}.json')['rows'] for s in range(3)]
flat = [r for nu in ['0.01', '0.005'] for m in ['A', 'B'] for r in runs(nu, m)]
grow = max(max(x['Z'] for x in r)/r[0]['Z'] for r in flat)
check('C2 a flat flow alone never grows its total spin', grow <= 1 + 1e-9 and all(r[-1]['off_frac'] < 1e-9 for r in flat),
      f'largest peak/start over 12 flat runs = {grow:.6f}; stays exactly flat')
two = {(nu, m): runs(nu, m) for nu in ['0.01', '0.005'] for m in ['twisted', 'untwisted']}
ratios = {k: [max(x['Z'] for x in r)/r[0]['Z'] for r in v] for k, v in two.items()}
check('C3 two crossing flat flows grow total spin in all 12 runs', all(x > 1.0 for v in ratios.values() for x in v),
      '; '.join(f"{k[0]} {k[1]}: {', '.join(f'{x:.2f}' for x in v)}" for k, v in ratios.items()))
offs = [r[-1]['off_frac'] for v in two.values() for r in v]
check('C4 most energy leaves the two planes by t = 8', min(offs) > 0.5, f'off-plane share at t=8: {min(offs):.2f} to {max(offs):.2f}')
def peak_to(f, T=4.0):
    r = [x for x in L(f)['rows'] if x['t'] <= T + 1e-9]; return max(x['Z'] for x in r)/r[0]['Z']
p9 = {m: [peak_to(f'runs/n64_nu0.005_{m}_k1_s{s}.json') for s in range(9)] for m in ['twisted', 'untwisted']}
diff = [a - b for a, b in zip(p9['twisted'], p9['untwisted'])]
check('C5 the twist does not tame the crossing (9 seeds, viscosity 0.005, peak by t = 4)',
      min(diff) < 0 < max(diff) and abs(np.mean(diff)) < 0.1*np.mean(p9['untwisted']),
      f'twisted lower in {sum(d < 0 for d in diff)} of 9; differences {", ".join(f"{d:+.2f}" for d in diff)}; '
      f'means {np.mean(p9["twisted"]):.2f} vs {np.mean(p9["untwisted"]):.2f}')
kd = {}
for m in ['twisted', 'untwisted']:
    kd[m] = []
    for k in ['0', '0.1', '0.25', '0.5', '1']:
        r = L(f'runs/n64_nu0.005_{m}_k{k}_s0.json')['rows']; kd[m].append((max(x['Z'] for x in r)/r[0]['Z'], r[-1]['off_frac']))
mono = all(all(b[0] >= a[0] - 1e-3 for a, b in zip(v, v[1:])) for v in kd.values())
check('C6 k setting on the crossing: peak spin rises with k (k = 0, 0.1, 0.25, 0.5, 1)', mono,
      '; '.join(f"{m}: {', '.join(f'{a:.2f}' for a, _ in v)}" for m, v in kd.items()))
k0 = [kd[m][0] for m in kd]
check('C7 k = 0: planes stay nearly flat and spin does not grow', all(z < 1.02 and o < 0.15 for z, o in k0),
      f"peak/start {k0[0][0]:.3f}, {k0[1][0]:.3f}; off-plane share {k0[0][1]:.3f}, {k0[1][1]:.3f} (the leak is the out-of-plane part fed through the shared line)")
a, b = L('runs/n64_nu0.005_twisted_k1_s0.json')['rows'], L('runs/n96_nu0.005_twisted_k1_s0.json')['rows']
za, zb = max(x['Z'] for x in a), max(x['Z'] for x in b)
ta = np.array([x['t'] for x in a]); tb = np.array([x['t'] for x in b])
gap = max(abs(np.interp(t, ta, [x['Z'] for x in a]) - z)/z for t, z in zip(tb, [x['Z'] for x in b]))
check('C8 grid check: 64³ and 96³ agree on total spin within 3% at every time', gap < 0.03,
      f'peak {za:.3f} vs {zb:.3f}; largest gap {100*gap:.1f}%; peak spin {max(x["wmax"] for x in a):.1f} vs {max(x["wmax"] for x in b):.1f}')
ts = L('twist_sectors.json')
sh = ts['random field: energy share by sector (1, i, -1, -i)']
cl = ts['nonlinear term of sector 0: share landing in sectors 0..3']
check('C9 twisted part: about a quarter of a random flow, closed under the nonlinear term', 0.2 < sh[0] < 0.35 and cl[0] > 1 - 1e-6,
      f'share {sh[0]:.3f}; sector -1 self-interaction lands in twisted part: {ts["nonlinear term of sector 2: share landing in sectors 0..3"][0]:.3f}')
if os.path.exists('runs/klein_planes.json'):
    kr = L('runs/klein_planes.json')['rows']
    fff = max(abs(x['P']['FFF']) for x in kr); tot = max(abs(sum(x['P'].values())) for x in kr)
    check('C10 Klein pair: the flat part alone does no stretching; all of it involves the bend', fff < 1e-6*max(tot, 1),
          f'max |P(flat,flat,flat)| = {fff:.1e} vs max total {tot:.2f}; flat share of spin {kr[0]["plane"]:.2f} -> {min(x["plane"] for x in kr):.2f}')
lines = [f"{'PASS' if r['ok'] else 'FAIL'} {r['name']}: {r['note']}" for r in res]
lines.append(f"{sum(r['ok'] for r in res)} of {len(res)} pass")
open('checks_output.txt', 'w').write('\n'.join(lines) + '\n'); json.dump(res, open('checks_results.json', 'w'), indent=1)
print('\n'.join(lines))
