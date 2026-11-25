"""B_s(1,1) for s from -1 to 2 (fused operation at height s applied to 1 and 1). Uses constants_hunt machinery at lower precision."""
import json, mpmath as mp, importlib.util, sys
src = open('/mnt/project-files/math-theory/pi/constants_hunt.py').read()
src = src.split('# sanity: E^1')[0].replace('DPS = 130','DPS = 40').replace('NIT = 20000','NIT = 3000').replace('NC = 60','NC = 30')
ns = {}
exec(src, ns)
Es = ns['Es']; mp.mp.dps = 40
rows = []
for i in range(-20, 41):
    s = mp.mpf(i) / 20
    if s == 0:
        v = mp.mpf(2)
    else:
        v = Es(s, 2 * Es(-s, mp.mpf(1)))
    rows.append([float(s), float(v)])
print(rows[20], rows[30], rows[40], rows[60])
json.dump(rows, open('half_curve.json', 'w'))
