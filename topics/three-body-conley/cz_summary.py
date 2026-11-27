"""DRAFT. Summary of the Conley-Zehnder iteration test on the catalog (reads cz.json, krein.json; writes
cz_summary.json): does the iterate data separate stable from unstable, how many iterates it takes, and how the
round-one fixed-point index sits inside it."""
import json, os, numpy as np
from collections import Counter
import cz
HERE = os.path.dirname(os.path.abspath(__file__))

def margin(r):
    m = []
    for tr, ti in r['tau']:
        if abs(ti) > 1e-9: m.append(abs(ti))
        else: m += [abs(2 - tr), abs(2 + tr)]
    return min(m)

def k_needed(mu, r):
    for K in range(1, len(mu) + 1):
        if cz.decode_r(mu[:K], grid=40001) == r: return K
    return None

if __name__ == '__main__':
    d = json.load(open(os.path.join(HERE, 'cz.json')))
    kr = json.load(open(os.path.join(HERE, 'krein.json')))
    rob = [r for r in d if margin(r) > 0.01 and r['rel_pal'] < 1e-6]
    pp = [r for r in rob if r['klass'] in ('E-E', 'H-H', 'IH-IH', 'CU')]
    out = dict(
        n=len(d), parity_ok=sum(r['parity_ok'] for r in d),
        class_vs_r={f'{k}|{v}': c for (k, v), c in Counter((r['klass'], r['r']) for r in d).items()},
        decode_ok=sum(r['r'] == r['r_decoded'] for r in d),
        plusplus_robust={k: Counter(r['r'] for r in pp if r['klass'] == k) for k in ('E-E', 'H-H', 'IH-IH', 'CU')},
        k_needed={k: Counter(k_needed(r['muhat'], r['r']) for r in rob if r['klass'] == k)
                  for k in ('E-E', 'E-H', 'E-IH')},
        d2={k: Counter(int(r['muhat'][1] - 2 * r['muhat'][0]) for r in d if r['klass'] == k)
            for k in ('E-E', 'E-H', 'E-IH', 'H-H', 'H-IH', 'IH-IH', 'CU')},
        krein_fate=kr['fate'], krein_types=Counter(x['krein'] for x in kr['orbits']),
        figure_eight=[x for x in kr['orbits'] if x['cls'] == 'I.A1'],
    )
    out['plusplus_robust'] = {k: dict(v) for k, v in out['plusplus_robust'].items()}
    out['k_needed'] = {k: {str(a): b for a, b in v.items()} for k, v in out['k_needed'].items()}
    out['d2'] = {k: {str(a): b for a, b in v.items()} for k, v in out['d2'].items()}
    json.dump(out, open(os.path.join(HERE, 'cz_summary.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))
