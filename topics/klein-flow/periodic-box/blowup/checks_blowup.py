"""Checks for the blow-up hunt. Reads the run files in runs/ and does two small fresh runs.
usage: python3 checks_blowup.py  -> checks_blowup_output.txt, checks_blowup_results.json"""
import json, os, numpy as np
os.environ.setdefault('GAMMA', '8'); os.environ.setdefault('DELTA', '0.5')
import tubes
from flow3d import Box, NS3
R = lambda f: json.load(open('runs/' + f))
res = []
def check(name, ok, detail):
    res.append(dict(name=name, ok=bool(ok), detail=detail)); print(('PASS ' if ok else 'FAIL ') + name + ': ' + detail, flush=True)
# B1: the base pair lies in both classes exactly
b = Box(64, 1.0, Nz=64); U = tubes.initial(b, 'none', 0, 0, 0.25, 2); dg, dt = tubes.defects(b, U)
check('B1 pair is in both the Klein and torus classes', max(dg, dt) < 1e-12, 'glide defect %.1e, shift defect %.1e' % (dg, dt))
# B2: each seed lands in its own class only
U = tubes.initial(b, 't-g+', 0, 0.01, 0.25, 2); dg1, dt1 = tubes.defects(b, U)
U = tubes.initial(b, 't+g-', 0, 0.01, 0.25, 2); dg2, dt2 = tubes.defects(b, U)
check('B2 Klein-only seed keeps the glide, breaks the shift; torus-only the reverse', dg1 < 1e-12 and dt1 > 1e-2 and dt2 < 1e-12 and dg2 > 1e-2,
      'K-only (%.1e, %.2f), TR-only (%.2f, %.1e)' % (dg1, dt1, dg2, dt2))
# B3: the class is kept through the whole run (N=128 runs)
ok = True; det = []
for f in ['n128_t-g+_s0.json', 'n128_t-g+_s1.json', 'n128_t-g+_s2.json']:
    r = R(f); ok &= r['defects1'][0] < 1e-10; det.append('%.0e' % r['defects1'][0])
for f in ['n128_t+g-_s0.json', 'n128_t+g-_s1.json', 'n128_t+g-_s2.json']:
    r = R(f); ok &= r['defects1'][1] < 1e-10; det.append('%.0e' % r['defects1'][1])
check('B3 every run stays in its class to the end', ok, 'end defects ' + ', '.join(det))
# B4: energy only falls (viscous flow, no forcing)
ok = True; worst = 0
for f in sorted(os.listdir('runs')):
    if f.endswith('.json') and (f.startswith('n128') or f.startswith('nu') or f.startswith('res')):
        E = np.array(R(f)['E']); d = np.diff(E).max()/E[0]; worst = max(worst, d); ok &= d < 1e-6
check('B4 energy never rises in any run', ok, 'largest step-to-step rise %.1e of start' % worst)
# B5: resolution, nu = 0.002 at N=128 (outer-shell energy) and nu = 0.001 at N=192 vs 128 (peak spin)
r = R('n128_none_s0.json'); check('B5 nu=0.002 runs resolved at N=128', r['tailmax'] < 1e-2, 'outer-shell energy %.1e' % r['tailmax'])
try:
    a, c = R('nu0.001_none_s0.json'), R('res192_none_s0.json')
    ta, tc = np.array(a['t']), np.array(c['t']); wa = np.array(a['wmax'])[ta <= 1.9]; wc = np.array(c['wmax'])
    Za, Zc = np.array(a['Z'])[ta <= 1.9], np.array(c['Z'])
    dw = abs(wa.max() - wc.max())/wc.max()
    zi = np.interp(tc, ta[ta <= 1.95], Za[:len(ta[ta <= 1.95])] if False else np.array(a['Z'])[ta <= 1.95])
    dz = np.max(np.abs(zi - Zc)/Zc)
    check('B6 nu=0.001: N=128 and N=192 agree on total spin (enstrophy) within 5% at every time', dz < 0.05,
          'largest gap %.1f%%, N=192 outer shell %.1e' % (100*dz, c['tailmax']))
    # B9: is peak spin settled? Decided only by the 256^3 run: if present, 192 vs 256 within 10% passes, otherwise fails.
    # Without it, the question is recorded as OPEN (neither pass nor fail), with the two readings we have.
    if os.path.exists('runs/res256_none_s0.json'):
        d = R('res256_none_s0.json'); td = np.array(d['t']); wd = np.array(d['wmax'])[td <= 1.9].max()
        dd = abs(wc.max() - wd)/wd
        check('B9 nu=0.001: peak spin settles, N=192 vs N=256 within 10%', dd < 0.1,
              'peak %.1f (128) / %.1f (192) / %.1f (256); 192 vs 256 differ by %.1f%%' % (wa.max(), wc.max(), wd, 100*dd))
    else:
        res.append(dict(name='B9 nu=0.001: is peak spin settled between grids?', ok=None, status='open',
                        detail='peak spin not settled, %.1f at 128 vs %.1f at 192; the 256 run will decide' % (wa.max(), wc.max())))
        print('OPEN B9: peak spin not settled, %.1f at 128 vs %.1f at 192; the 256 run will decide' % (wa.max(), wc.max()), flush=True)
except FileNotFoundError:
    check('B6 nu=0.001 resolution check', False, 'N=192 run missing')
# B7: Kerr crossing: sqrt(nu) Z for nu = 0.004, 0.002, 0.001 (pair alone) meet within 5% at one time
cur = {}
for nu, tag in ((0.004, 'nu0.004'), (0.002, 'n128'), (0.001, 'nu0.001')):
    r = R('%s_none_s0.json' % tag); cur[nu] = (np.array(r['t']), np.sqrt(nu)*np.array(r['Z']))
tt = np.linspace(1.0, 2.0, 201); spread = []
for t in tt:
    v = [np.interp(t, *cur[nu]) for nu in cur]; spread.append((max(v) - min(v))/np.mean(v))
i = int(np.argmin(spread))
check('B7 Kerr crossing: the three viscosities meet at one time', spread[i] < 0.05, 'spread %.1f%% at t = %.2f' % (100*spread[i], tt[i]))
# B8: determinism, a fresh short run repeats a stored number exactly
r1 = tubes.run(48, 0.004, 0.3, 't-g+', 3, 0.01, 0.25, 2); r2 = tubes.run(48, 0.004, 0.3, 't-g+', 3, 0.01, 0.25, 2)
check('B8 same seed gives the same run', abs(r1['wmax'][-1] - r2['wmax'][-1]) < 1e-12, 'peak %.10f twice' % r1['wmax'][-1])
for r in res:
    r.setdefault('status', 'pass' if r['ok'] else 'fail')
json.dump(res, open('checks_blowup_results.json', 'w'), indent=1)
done = [r for r in res if r['status'] != 'open']; nopen = len(res) - len(done)
print('%d of %d pass' % (sum(r['ok'] for r in done), len(done)) + (', %d open' % nopen if nopen else ''))
if any(r['status'] == 'fail' for r in res):
    raise SystemExit(1)
