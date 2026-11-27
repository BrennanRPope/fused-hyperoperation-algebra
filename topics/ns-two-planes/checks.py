"""DRAFT. Checks for ns-two-planes, from algebra_results.json, algebra3_results.json, twist_sectors.json and runs/*.json.
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

# ---- Round 2: three perpendicular planes (algebra3.py, threeplanes.py, batch4.sh) ----
if os.path.exists('algebra3_results.json'):
    a3 = L('algebra3_results.json')
    check('C11 three-plane algebra (algebra3.py): exact split, spins add up, all stretching on crossing lines, no cyclic twist on the cube',
          all(v['pass'] for v in a3.values()), f"{sum(v['pass'] for v in a3.values())} of {len(a3)} pass")
if os.path.exists('runs/p3_two-un_s0.json'):
    def lines_of(T):
        T = np.array(T); return np.array([T[1, :, 2].sum()+T[2, :, 1].sum(), T[0, :, 2].sum()+T[2, :, 0].sum(), T[0, :, 1].sum()+T[1, :, 0].sum()])
    P3 = {m: [L(f'runs/p3_{m}_s{s}.json')['rows'] for s in range(3)] for m in ['two-tw', 'two-un', 'three-tw', 'three-un']}
    def stats(r):
        t = np.array([x['t'] for x in r]); Z = np.array([x['Z'] for x in r]); ip = int(np.argmax(Z))
        comp = np.array([np.array(x['T']).sum(axis=(1, 2)) for x in r]); ln = np.array([lines_of(x['T']) for x in r])
        g = np.trapezoid(comp[:ip+1], t[:ip+1], axis=0); G = np.trapezoid(ln[:ip+1], t[:ip+1], axis=0)
        self_ = max(abs(np.array(x['T'])[i, :, i].sum()) for x in r for i in range(3)); tot = max(abs(np.array(x['T']).sum()) for x in r)
        return dict(peak=Z[ip]/Z[0], gain=g/g.sum(), line=G/G.sum(), early=ln[2]/ln[2].sum(), self_=self_/tot,
                    compgrow=np.array(r[ip]['comp_Z'])/np.maximum(np.array(r[0]['comp_Z']), 1e-12), endshare=np.array(r[-1]['comp_Z'])/r[-1]['Z'])
    S = {m: [stats(r) for r in v] for m, v in P3.items()}
    allS = [s for v in S.values() for s in v]
    check('C12 runs: no plane stretches its own spin; all stretching sits on the three crossing lines (12 runs)',
          max(s['self_'] for s in allS) < 1e-6, f"largest own-plane term / total = {max(s['self_'] for s in allS):.1e}")
    two = S['two-tw'] + S['two-un']
    check('C13 two-plane starts: the third plane (y, no spin at start) takes at least ~0.3 of the spin gain and ends with ~1/3 of the spin',
          all(s['gain'][1] > 0.3 and s['endshare'][1] > 0.3 for s in two),
          'gain share of plane y: ' + ', '.join(f"{s['gain'][1]:.2f}" for s in two) + '; spin share at t=4: ' + ', '.join(f"{s['endshare'][1]:.2f}" for s in two))
    check('C14 two-plane starts: stretching starts on the shared line y (>0.7 at t=0.4), then spreads to all three lines',
          all(s['early'][1] > 0.7 and min(s['line']) > 0.14 for s in two),
          'line y early: ' + ', '.join(f"{s['early'][1]:.2f}" for s in two) + '; smallest line share up to the peak: ' + ', '.join(f"{min(s['line']):.2f}" for s in two))
    pk = {m: [s['peak'] for s in v] for m, v in S.items()}
    check('C15 three planes grow more than two at the same energy (every three-plane seed above every two-plane seed)',
          min(pk['three-tw'] + pk['three-un']) > max(pk['two-tw'] + pk['two-un']),
          '; '.join(f"{m}: {', '.join(f'{x:.2f}' for x in v)}" for m, v in pk.items()))
    d3 = [a - b for a, b in zip(pk['three-tw'], pk['three-un'])]
    check('C16 the twist does not lower growth with three planes counted (3 seeds)', np.mean(d3) > -0.05 and sum(d < 0 for d in d3) < 3,
          f"twisted minus untwisted: {', '.join(f'{d:+.2f}' for d in d3)}; means {np.mean(pk['three-tw']):.2f} vs {np.mean(pk['three-un']):.2f}")
    three = S['three-tw'] + S['three-un']
    check('C17 three-plane starts: every plane\'s own spin grows (no plane stays quiet), each plane ~1/3 of the gain',
          all(min(s['compgrow']) > 1.5 and max(abs(s['gain'] - 1/3)) < 0.06 for s in three),
          'smallest per-plane growth at the peak: ' + ', '.join(f"{min(s['compgrow']):.1f}x" for s in three)
          + '; gain shares within ' + f"{max(max(abs(s['gain'] - 1/3)) for s in three):.2f} of 1/3")
lines = [f"{'PASS' if r['ok'] else 'FAIL'} {r['name']}: {r['note']}" for r in res]
lines.append(f"{sum(r['ok'] for r in res)} of {len(res)} pass")
open('checks_output.txt', 'w').write('\n'.join(lines) + '\n'); json.dump(res, open('checks_results.json', 'w'), indent=1)
print('\n'.join(lines))
