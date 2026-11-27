"""DRAFT. (1) Catalog tally of Floquet classes and index pairs, with a precision margin.
(2) Isolation test: do orbits started near a fixed point stay near it (invariant tori) or leave?
Writes analysis.json."""
import json, numpy as np, multiprocessing as mp
from collections import Counter
import cb3, section as S

PAIR_SIGNS = {'E-E': (1, 1), 'H-H': (1, 1), 'IH-IH': (1, 1), 'CU': (1, 1), 'E-H': (-1, -1), 'E-IH': (1, -1), 'H-IH': (-1, 1)}

def margin(o):
    m = []
    for tr, ti in o['tau']:
        if abs(ti) > 1e-9: m.append(abs(ti))
        else: m += [abs(2 - tr), abs(2 + tr)]
    m.append(abs(o['disc']) ** 0.5)
    return min(m)

def stay(args):
    """max distance (in units of r) reached within N returns, for a start at distance r in direction u."""
    s, E, k, r, u, N = args
    x = s + r * u
    dmax = 0.0
    for i in range(N):
        try:
            out = S.F(x, E, k)
        except RuntimeError:
            return np.inf
        if out is None: return np.inf
        x = out[0]; d = np.linalg.norm(x - s) / r
        dmax = max(dmax, d)
        if dmax > 50: return dmax
    return dmax

if __name__ == '__main__':
    d = [o for o in json.load(open('floquet.json')) if o['ok']]
    rob = [o for o in d if margin(o) > 0.01 and o['rel_pal'] < 1e-6]
    res = dict(n_all=len(d), n_robust=len(rob), classes_all=Counter(o['klass'] for o in d),
               classes_robust=Counter(o['klass'] for o in rob))
    pp = [o for o in rob if tuple(o['idx'][:2]) == (1, 1)]
    res['plusplus'] = Counter(o['klass'] for o in pp)
    # parity law check on every robust orbit, iterates 1..4
    def law(o, n):
        t = [complex(*x) for x in o['tau']]; ty = [cb3.pair_type(x) for x in t]
        return (-1) ** ty.count('H') * ((-1) ** (n + 1)) ** ty.count('IH')
    res['parity_law_holds'] = all(o['idx'][n - 1] == law(o, n) for o in rob for n in (1, 2, 3, 4))
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=dict))
    # isolation test: figure-eight (E-E), Simo's H orbit (E-H), I.A4 (H-H), II.C24 (CU)
    cur = {o['cls']: o for o in json.load(open('curated.json')) if o['ok']}
    tw = {(o['name'], o['which']): o for o in json.load(open('twins.json'))}
    targets = {'figure-eight I.A1 (E-E)': tw[('I.A1', 'cat')], "its twin, Simo's H orbit (E-H)": tw[('I.A1', 'twin')],
               'I.A4 (H-H)': cur['I.A4'], 'II.C24 (CU)': cur['II.C24']}
    rng = np.random.default_rng(7)
    radii = [1e-4, 1e-3, 1e-2, 3e-2]
    iso = {}
    for name, o in targets.items():
        if 's' in o: s = np.array(o['s']); E, k = o['E'], o['k']
        else:
            E = o['E']; k = o['k']; s0 = np.array([2, 0, o['v'][0], o['v'][0]]); s, _ = S.refine(s0, E, k)
        jobs = []
        for r in radii:
            for j in range(24):
                u = rng.normal(size=4); u /= np.linalg.norm(u)
                jobs.append((s, E, k, r, u, 300))
        with mp.Pool(4) as p: out = p.map(stay, jobs)
        out = np.array(out).reshape(len(radii), 24)
        iso[name] = dict(radii=radii, frac_stay=[float(np.mean(row < 50)) for row in out],
                         median_excursion=[float(np.median(np.minimum(row, 50))) for row in out])
        print(name, iso[name])
    res['isolation'] = iso
    json.dump(res, open('analysis.json', 'w'), indent=1, default=dict)
