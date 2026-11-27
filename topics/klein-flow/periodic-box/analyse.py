import json, glob, numpy as np
R = sum((json.load(open(f)) for f in sorted(glob.glob('runs/wall_s*.json'))), [])
hs = sorted({r['h'] for r in R}); As = sorted({r['a'] for r in R}); seeds = sorted({r['seed'] for r in R})
def get(s, h, a, key):
    return np.array([r[key] for r in R if r['surf'] == s and r['h'] == h and r['a'] == a])
def wall(s, h, seed, thr=1.01):
    """smallest a with G > thr, interpolated in log a (nan if never)."""
    rs = sorted([r for r in R if r['surf'] == s and r['h'] == h and r['seed'] == seed], key=lambda r: r['a'])
    prev = None
    for r in rs:
        if r['G'] > thr:
            if prev is None: return r['a']
            la0, la1 = np.log(prev['a']), np.log(r['a']); g0, g1 = prev['G'], r['G']
            return float(np.exp(la0 + (thr - g0)/(g1 - g0)*(la1 - la0)))
        prev = r
    return np.nan
if __name__ == '__main__':
    print('runs', len(R), 'heights', hs, 'amps', As, 'seeds', seeds)
    for h in hs:
        print(f'h={h}')
        for a in As:
            gk, gt = get('K', h, a, 'G'), get('TR', h, a, 'G')
            if len(gk) == 0: continue
            print(f'  a={a:6}: G  K {gk.mean():.3f}±{gk.std():.3f}  TR {gt.mean():.3f}±{gt.std():.3f}   K>TR in {int((gk>gt).sum())}/{len(gk)}')
        wk = [wall('K', h, s) for s in seeds]; wt = [wall('TR', h, s) for s in seeds]
        print('  wall a*  K', np.round(wk, 4), ' TR', np.round(wt, 4))
