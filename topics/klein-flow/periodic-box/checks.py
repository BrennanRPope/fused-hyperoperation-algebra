"""Checks for the periodic-box round. Writes checks_output.txt and checks_results.json."""
import json, time, numpy as np
from flow3d import Box, NS3, symmetrize, glide_defect, glide
from wall import initial, run
out, res = [], {}
def check(name, ok, msg):
    res[name] = dict(ok=bool(ok), msg=msg); out.append(f"{'PASS' if ok else 'FAIL'} {name}: {msg}"); print(out[-1], flush=True)

# P1 inviscid energy is kept by the scheme on both surfaces
for s in ('K', 'TR'):
    b = Box(32, 1.0); ns = NS3(b, 0.0, s); U = initial(b, s, 0, 0.3); E0 = ns.energy(U)
    for _ in range(40): U = ns.step(U, 0.01)
    d = abs(ns.energy(U)/E0 - 1)
    check(f'P1-{s}', d < 1e-8, f'inviscid energy drift over 40 steps = {d:.1e}')
# P2 the glide symmetry is kept by the equations (K run, 3D, full cube)
r = run(32, 0.01, 2.0, 1, 'K', 1.0, 0.3)
check('P2', r['glide_defect'] < 1e-12, f"glide defect after t = 2: {r['glide_defect']:.1e}")
# P3 glide applied twice is the identity on the cover grid (period 2 pi)
f = np.random.default_rng(0).standard_normal((32, 32, 16))
check('P3', np.array_equal(glide(glide(f, 32), 32), f), 'glide(glide(f)) == f on the grid')
# P4 pure 2D flow: enstrophy never rises (the 2D law)
for s in ('K', 'TR'):
    r = run(32, 0.01, 3.0, 2, s, 1.0, 0.0)
    check(f'P4-{s}', np.all(np.diff(r['Z']) <= 1e-12) and r['G'] == 1.0, f'2D enstrophy monotone, gain {r["G"]:.6f}')
# P5 the shared 2D base is identical on both surfaces
bK = Box(32, 1.0); UK = initial(bK, 'K', 0, 0.0, joint=True); UT = initial(bK, 'TR', 0, 0.0, joint=True)
check('P5', np.abs(UK - UT).max() < 1e-12, f'max |base_K - base_TR| = {np.abs(UK - UT).max():.1e}')
# P6 divergence-free to round-off
b = Box(32, 0.5); U = initial(b, 'K', 3, 0.3)
div = np.abs(b.KX*U[0] + b.KY*U[1] + b.KZ*U[2]).max()/np.abs(U).max()
check('P6', div < 1e-12, f'max |k.u| / max |u| = {div:.1e}')
# P7 thin slab: a 3D part as large as the 2D flow dies (thin-domain behaviour)
r = run(48, 0.01, 3.0, 0, 'K', 0.0625, 1.0)
check('P7', r['E3'][-1] < 1e-3*r['E3'][0] and r['G'] == 1.0, f"h = 1/16, a = 1: 3D share {r['E3'][0]:.2f} -> {r['E3'][-1]:.1e}, gain {r['G']:.3f}")
# P8 resolution: the cube's gain and 3D growth at N = 64 vs N = 96
t0 = time.time()
r64 = run(64, 0.01, 2.5, 0, 'K', 1.0, 1.0); r96 = run(96, 0.01, 2.5, 0, 'K', 1.0, 1.0)
dG = abs(r64['G'] - r96['G'])/r96['G']
check('P8', dG < 0.02, f"cube, a = 1, t <= 2.5: gain {r64['G']:.4f} (N=64) vs {r96['G']:.4f} (N=96), diff {100*dG:.2f}% ({time.time()-t0:.0f}s)")
r64 = run(64, 0.01, 5.0, 1, 'K', 1.0, 0.001, True); r96 = run(96, 0.01, 5.0, 1, 'K', 1.0, 0.001, True)
g64 = r64['E3abs'][-1]/r64['E3abs'][0]; g96 = r96['E3abs'][-1]/r96['E3abs'][0]
check('P9', abs(g64 - g96)/g96 < 0.03, f'cube, small seed: 3D growth {g64:.3f} (N=64) vs {g96:.3f} (N=96)')
n = sum(v['ok'] for v in res.values())
out.append(f'{n} of {len(res)} pass')
open('checks_output.txt', 'w').write('\n'.join(out) + '\n'); json.dump(res, open('checks_results.json', 'w'), indent=1)
print(out[-1])
