"""Fixed-point indices of the stable orbit and its same-pattern twin in all 65 logistic windows
(periods 2 to 9) from entropy-knots/windows.json. Writes windows_index.json."""
import json
import numpy as np
from orbits import load_windows, orbits_of_exact_period, pattern

rows = []
for w in load_windows():
    n, r = w['n'], w['r']
    orbs = orbits_of_exact_period(r, n)
    info = [(pattern(o, r), m) for o, m in orbs]
    stable = [i for i, (p, m) in enumerate(info) if abs(m) < 1]
    s = stable[0]
    ps, ms = info[s]
    twins = [(p, m) for i, (p, m) in enumerate(info) if i != s and p == ps]
    rows.append({'n': n, 'r': r, 'stable_pattern': list(ps), 'stable_mult': ms,
                 'stable_index': int(np.sign(1 - ms)) if ms != 1 else 0,
                 'twin_mults': [m for _, m in twins],
                 'twin_indices': [int(np.sign(1 - m)) for _, m in twins],
                 'n_orbits': len(orbs)})
json.dump(rows, open('windows_index.json', 'w'), indent=1)
print(len(rows), 'windows written')
