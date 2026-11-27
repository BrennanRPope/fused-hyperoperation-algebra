import json, mpmath as mp, math
from windows import orbit
from braids import orbit_braid, braid_entropy
mp.mp.dps = 40
R = json.load(open('windows.json'))
for row in R:
    xs = [float(x) for x in orbit(mp.mpf(row['r']), row['n'])]
    word, n = orbit_braid(xs)
    row['braid'] = word
    row['h_braid'] = braid_entropy(word, n, iters=400)
json.dump(R, open('windows.json', 'w'), indent=1)
for row in R:
    print(row['n'], round(row['r'], 5), row['perm'], 'hP=%.4f hB=%.4f' % (row['h_pattern'], row['h_braid']), 'len', len(row['braid']))
