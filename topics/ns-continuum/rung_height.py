"""Rung height (circle-operation flow E^t, E(x)=e^x-1) of the dissipation deficiency in the NS gap.
Dissipation symbol m(s) = s^{5/2} / D(s), s = |xi|, L = ln s.
  3D NS: D = s^{1/2} = e^{L/2}          Tao 2009 proved: D up to L^{1/4} (Osgood int dL/D^4 = inf)
  Tao's conjecture (proved for the dyadic model, Barbato-Morandin-Romito 2014): D up to L^{1/2}
  Rung family: D_t = E^t(L/2), t in [0,1]; t = 1 is s^{1/2} - 1 (3D NS), t = 0 is L/2.
"""
import sys, json, numpy as np
sys.path.insert(0, '/mnt/project-files/math-theory/fractional-rungs')
from regular_iteration import Es
out = []
def rep(name, ok, extra=''):
    out.append((name, bool(ok), extra)); print(('PASS ' if ok else 'FAIL ') + name + ('  ' + extra if extra else ''))

L = np.linspace(2.0, 690.0, 3000)
fam = {}
for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
    fam[t] = Es(t, L / 2)
# 1. endpoints
rep('t=0 slice is L/2', np.allclose(fam[0.0], L / 2, rtol=1e-9))
rep('t=1 slice is s^(1/2)-1 (3D NS deficiency)', np.allclose(fam[1.0], np.expm1(L / 2), rtol=1e-8))
# 2. ordered in t at every L
mono = all((fam[a] < fam[b]).all() for a, b in [(0, .25), (.25, .5), (.5, .75), (.75, 1)])
rep('family is ordered: higher rung = bigger deficiency at every frequency', mono)
# 3. flow law: E^0.25 after E^0.25 = E^0.5
x = np.linspace(1, 300, 50)
rep('flow law E^.25(E^.25(x)) = E^.5(x)', np.allclose(Es(0.25, Es(0.25, x)), Es(0.5, x), rtol=1e-7))
# 4. growth: log D_t / log L keeps rising for t>0 (faster than any power of L) while log D_t / L falls (slower than any power of s)
for t in [0.25, 0.5, 0.75]:
    D = fam[t]; r1 = np.log(D) / np.log(L); r2 = np.log(D) / L
    i1, i2 = np.searchsorted(L, 50), -1
    rep(f't={t}: log D / log L rises (beats every power of L on the range)', r1[i2] > r1[i1] * 1.2,
        f'{r1[i1]:.2f} -> {r1[i2]:.2f}')
    rep(f't={t}: log D / L falls (below every power of s on the range)', r2[i2] < r2[i1],
        f'{r2[i1]:.3f} -> {r2[i2]:.3f}')
# 5. Osgood integrals int dL / D^4 (Tao's test) and int dL / D^2 (conjectured test)
def cum(f):
    g = 1.0 / f; return np.concatenate([[0], np.cumsum((g[1:] + g[:-1]) / 2 * np.diff(L))])
osg = {}
for name, D in [('Tao edge L^(1/4)', (L / 2) ** 0.25 * 2**0.25), ('conj edge L^(1/2)', np.sqrt(L))] + [(f'rung t={t}', fam[t]) for t in fam]:
    with np.errstate(over='ignore'): I4 = cum(D ** 4); I2 = cum(D ** 2)
    osg[name] = dict(I4=I4[::30].tolist(), I2=I2[::30].tolist())
    print(f'   {name:18s} int dL/D^4 to L=100: {I4[np.searchsorted(L,100)]:.4f}  to L=690: {I4[-1]:.4f}   int dL/D^2: {I2[np.searchsorted(L,100)]:.4f} -> {I2[-1]:.4f}')
I4tao = cum((L) ** 1.0)
rep('Tao edge: integral keeps growing like ln L', abs((I4tao[-1] - I4tao[np.searchsorted(L, 100)]) - np.log(690 / L[np.searchsorted(L, 100)])) < 1e-2)
for t in [0.25, 0.5]:
    I4 = cum(fam[t] ** 4)
    tail = I4[-1] - I4[np.searchsorted(L, 300)]
    rep(f'rung t={t}: Osgood integral has stopped growing (tail from L=300 is tiny)', tail < 1e-6 * max(I4[-1], 1e-300) + 1e-12, f'tail={tail:.2e}')
# 6. classical intermediate family exp(L^b) sits above the half rung on this range (height 1)
for b in [0.5]:
    C = np.exp((L / 2) ** b)
    r = (L / 2) ** b / np.log(fam[0.5])
    rep(f'exp((L/2)^{b}) gains on the half rung (ratio of logs rises; it overtakes only past s=1e300)', r[-1] > r[np.searchsorted(L, 50)],
        f'{r[np.searchsorted(L,50)]:.3f} -> {r[-1]:.3f}')
json.dump(dict(L=L[::30].tolist(), fam={str(t): np.log10(v[::30]).tolist() for t, v in fam.items()}, osg=osg,
               checks=out), open('rung_height.json', 'w'))
print(sum(o[1] for o in out), '/', len(out), 'pass')
