"""Round 9: the string probe.  Fit settings on zeros 1-500, test on zeros 501-5500.

Primary closeness score: rms of (E_n - gamma_n) / mean spacing on the held-out zeros.
Signatures: bias (smooth error), osc (oscillating error), Gram-interval hit rate,
nearest-zero hit rate, spacing KS (model vs zeros), number variance error, pair
correlation error, missing levels (no crossing found)."""
import json, time
import numpy as np
from math import pi, log
from stringlib import (String, theta, gram, mean_spacing, true_S, score, primes_upto)

gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
TR = np.arange(0, 500)
TE = np.arange(500, NZ)
OUT = {}
t0 = time.time()


def run(label=None, **kw):
    st = String(**kw)
    E = st.levels(NZ)
    r = dict(train=score(E, gam, TR), test=score(E, gam, TE))
    if label:
        print(f"{label:38s} train {r['train']['rms']:.4f}  test {r['test']['rms']:.4f}  "
              f"gram {r['test']['gram_hit']:.3f}  near {r['test']['near']:.3f}  ks {r['test']['ks']:.3f}")
    return r, E


# ---------- 0. baselines ----------
print("== baselines ==")
EG = gram(np.arange(1, NZ + 1) - 1.5)               # smooth (Gram midpoint): S = 0
OUT['smooth'] = dict(train=score(EG, gam, TR), test=score(EG, gam, TE))
print("smooth law only (Gram midpoint)       train %.4f test %.4f" % (OUT['smooth']['train']['rms'], OUT['smooth']['test']['rms']))

# ---------- 1. settings: cutoff X and taper ----------
print("== settings: cutoff and taper (dials at truth) ==")
grid = {}
for taper in ['sharp', 'fejer', 'gauss']:
    for X in ([5, 10, 20, 50, 100, 200, 500, 1000, 3000] if taper != 'gauss' else [5, 10, 20, 50, 100, 200, 500, 1000]):
        r, _ = run(f"taper={taper} X={X}", X=float(X), taper=taper)
        grid[f"{taper}|{X}"] = r
OUT['cutoff_grid'] = grid
best = min(grid, key=lambda k: grid[k]['train']['rms'])
print("best on training zeros:", best, grid[best]['train']['rms'], "-> test", grid[best]['test']['rms'])
OUT['best_in_grid'] = best
# The score falls monotonically with the cutoff for every taper, so the cutoff is a
# regulator (more orbits is always better), not a setting with an optimum.  For the
# sweeps we fix one moderate regulator: Gaussian taper, X = 200.
bt, bX = 'gauss', 200.0
OUT['best_setting'] = dict(taper=bt, X=bX)
base = dict(taper=bt, X=bX)

# ---------- 2. continuous settings around the fit ----------
print("== settings: extension angle, amplitude, phase ==")
sweeps = {}
for name, vals in [('alpha', [-pi, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, pi]),
                   ('A', [0.0, 0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]),
                   ('phi', [-pi / 2, -0.3, -0.1, 0.0, 0.1, 0.3, pi / 2, pi])]:
    sweeps[name] = {}
    for v in vals:
        kw = dict(base); kw[name] = v
        r, _ = run(f"{name}={v:.3f}", **kw)
        sweeps[name][f"{v:.4f}"] = r
OUT['setting_sweeps'] = sweeps

# ---------- 3. growth law (classical dynamics) ----------
print("== settings: growth law (classical dynamics) ==")
gl = {}
g500 = gam[499]
for g in ['xp', 'xp_lead', 'oscillator', 'box', 'gammaC', 'nopi']:
    kw = dict(base, growth=g)
    if g == 'oscillator':
        kw['growth_c'] = 500 / g500     # fitted so that N(gamma_500) = 500
    if g == 'box':
        kw['growth_c'] = 500 / np.sqrt(g500)
    r, _ = run(f"growth={g}", **kw)
    gl[g] = r
OUT['growth'] = gl

# ---------- 4. the seven dials ----------
print("== dials ==")
dials = {'k': [-1.0, -0.5, 0.0, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0],
         't': [0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2],
         's': [0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 1.0],
         'd': [-0.25, -0.1, -0.05, 0.0, 0.05, 0.1, 0.25, 0.5, 1.0],
         'q': [0.0, 0.05, 0.1, 0.25, 0.5, 1.0],
         'r': [0.9, 0.97, 0.99, 0.999, 1.0, 1.001, 1.01, 1.03, 1.1]}
truth = dict(k=1.0, t=1.0, s=0.5, d=0.0, q=0.0, r=1.0)
dres = {}
for dname, vals in dials.items():
    dres[dname] = {}
    for v in vals:
        kw = dict(base); kw[dname] = v
        if dname == 'q' and v > 0:
            kw['X'] = min(bX, 1000.0)
        r, _ = run(f"{dname}={v}", **kw)
        dres[dname][f"{v}"] = r
OUT['dials'] = dres


# ROUNDNESS p: at p = 2 the lattice zeta is 4 zeta(s) L(s, chi_-4); add the levels of chi_-4
def chi4(n):
    n = np.asarray(n, int)
    return np.where(n % 2 == 0, 0, np.where(n % 4 == 1, 1, -1))


class UnionP2(String):
    """zeta * L(chi_-4): smooth count of both, orbits of both (chi_-4 weights)."""
    def __init__(self, **kw):
        super().__init__(**kw)
        self.wL = self.w * chi4(np.round(np.exp(self.L)).astype(int))

    def Nbar(self, E):
        E = np.asarray(E, float)
        thL = E / 2 * np.log(4 * E / (2 * pi * np.e)) + pi / 8   # odd character, to O(1/E)
        return theta(E) / pi + 1 + thL / pi + 0.5

    def S(self, E, chunk=400):
        E = np.atleast_1d(np.asarray(E, float))
        out = np.empty_like(E)
        for i in range(0, len(E), chunk):
            e = E[i:i + chunk, None]
            out[i:i + chunk] = -((self.w + self.wL) * np.sin(e * self.L)).sum(1) / pi
        return out


u = UnionP2(**dict(base, growth='union'))
Eu = u.levels(NZ)
OUT['dial_p2'] = dict(train=score(Eu, gam, TR), test=score(Eu, gam, TE),
                      density_ratio=float((u.Nbar(3000.0) - u.Nbar(2999.0))))
print("p=2 (zeta * L(chi_-4)) test rms %.3f" % OUT['dial_p2']['test']['rms'])

# ---------- 5. sensitivities (central differences at the truth, on the fitted settings) ----------
print("== sensitivities ==")
steps = dict(k=0.05, t=0.02, s=0.02, d=0.02, r=0.001, alpha=0.1, A=0.05, phi=0.05, q=0.05)
sens = {}
for nme, h in steps.items():
    vals = {}
    for sgn in (-1, 1):
        kw = dict(base)
        v0 = truth.get(nme, {'alpha': 0.0, 'A': 1.0, 'phi': 0.0}.get(nme))
        if nme == 'q':
            if sgn < 0:
                continue
            kw['X'] = min(bX, 1000.0)
        kw[nme] = v0 + sgn * h
        r, _ = run(None, **kw)
        vals[sgn] = r
    r0 = run(None, **base)[0] if nme != 'q' else run(None, **dict(base, X=min(bX, 1000.0)))[0]
    if nme == 'q':
        gtr = (vals[1]['train']['rms'] - r0['train']['rms']) / h
        gte = (vals[1]['test']['rms'] - r0['test']['rms']) / h
    else:
        gtr = (vals[1]['train']['rms'] - vals[-1]['train']['rms']) / (2 * h)
        gte = (vals[1]['test']['rms'] - vals[-1]['test']['rms']) / (2 * h)
    sens[nme] = dict(train=gtr, test=gte, step=h)
    print(f"d score / d {nme}: train {gtr:+.4f}  test {gte:+.4f}")
OUT['sensitivity'] = sens

# ---------- 6. best next adjustment: gradient descent on training, then test ----------
print("== next adjustment ==")
cont = ['s', 'A', 'alpha', 'phi', 't', 'k']
x = dict(s=0.5, A=1.0, alpha=0.0, phi=0.0, t=1.0, k=1.0)
hstep = dict(s=0.01, A=0.02, alpha=0.05, phi=0.02, t=0.01, k=0.02)


def trscore(x):
    return run(None, **dict(base, **x))[0]


path = []
r0 = trscore(x)
path.append(dict(x=dict(x), train=r0['train']['rms'], test=r0['test']['rms']))
for it in range(6):
    gvec = {}
    for nme in cont:
        xp = dict(x); xp[nme] += hstep[nme]
        xm = dict(x); xm[nme] -= hstep[nme]
        gvec[nme] = (trscore(xp)['train']['rms'] - trscore(xm)['train']['rms']) / (2 * hstep[nme])
    # scale-free step: normalise each coordinate by its step size
    gn = np.sqrt(sum((gvec[n] * hstep[n]) ** 2 for n in cont))
    if gn < 1e-6:
        break
    best_lr, best_r = None, None
    for lr in [4, 2, 1, 0.5]:
        xt = {n: x[n] - lr * gvec[n] * hstep[n] ** 2 / gn for n in cont}
        rt = trscore(xt)
        if best_r is None or rt['train']['rms'] < best_r['train']['rms']:
            best_lr, best_r, best_x = lr, rt, xt
    if best_r['train']['rms'] >= path[-1]['train']:
        break
    x = best_x
    path.append(dict(x=dict(x), grad=gvec, train=best_r['train']['rms'], test=best_r['test']['rms']))
    print(f"step {it}: train {best_r['train']['rms']:.4f} test {best_r['test']['rms']:.4f}  "
          + " ".join(f"{n}={x[n]:.3f}" for n in cont))
OUT['descent'] = path

# ---------- 7. string-theory checks ----------
print("== string theory ==")
st_out = {}
# (a) Freund-Witten: archimedean place supplies the smooth law, p-adic places the orbits.
#     Drop one p-adic place at a time.
r_full = run(None, **base)[0]
drops = {}
for pdrop in [2, 3, 5, 7, 11, 13, 29, 97]:
    s_ = String(**base)
    pp = np.round(np.exp(s_.L)).astype(int)
    mask = np.ones(len(pp), bool)
    for i, n in enumerate(pp):
        m = n
        while m % pdrop == 0:
            m //= pdrop
        if m == 1:
            mask[i] = False
    s_.w = s_.w[mask]; s_.L = s_.L[mask]
    E = s_.levels(NZ)
    r = dict(train=score(E, gam, TR), test=score(E, gam, TE))
    drops[pdrop] = r
    print(f"drop p={pdrop}: test {r['test']['rms']:.4f} (full {r_full['test']['rms']:.4f})")
st_out['drop_place'] = drops
st_out['archimedean'] = {g: OUT['growth'][g]['test']['rms'] for g in ['xp', 'gammaC', 'nopi']}
# (b) primon-gas statistics = homothety dial: Bose k=1, Fermi k=-1, Boltzmann k=0
st_out['statistics'] = {lab: OUT['dials']['k'][v]['test']['rms'] for lab, v in
                        [('Bose (k=1)', '1.0'), ('Boltzmann (k=0)', '0.0'), ('Fermi (k=-1)', '-1.0')]}
# (c) Connes' absorption sign: phase pi flips the orbit sign
st_out['sign_flip'] = OUT['setting_sweeps']['phi'][f"{pi:.4f}"]['test']['rms']
# (d) time reversal (orbit + reversed orbit add coherently): A = 2
st_out['time_reversal'] = OUT['setting_sweeps']['A']['2.0000']['test']['rms']
OUT['string_theory'] = st_out
print(json.dumps(st_out, indent=0, default=str)[:800])

OUT['elapsed'] = time.time() - t0
json.dump(OUT, open('probe.json', 'w'), default=float)
print("elapsed", OUT['elapsed'])
