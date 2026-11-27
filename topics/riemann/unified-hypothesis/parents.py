"""Real-rooted parent hunt for Riemann's kernel Phi.
A parent of K is K_P = K / U with U a Pólya universal factor; if K_P's transform has only real zeros, so does K's.
We scan heat parents (U = e^{lam u^2}), factor parents (U = 1 + u^2/c^2), two-level-bath parents (U = cosh(J u)),
the antiderivative parent, fused (∘_k) Riemann-sum parents, and the averaged-away bath children."""
import mpmath as mp, json, time, sys
from multiprocessing import Pool
import riemann_kernel as rk
mp.mp.dps = 45
ZMAX, DZ = 150.0, 0.05
GAMMA = [14.134725141734693, 21.022039638771555, 25.010857580145688, 30.424876125859513, 32.935061587739189,
         37.586178158825671, 40.918719012147495, 43.327073280914999, 48.005150881167159, 49.773832477672302,
         52.970321477714460, 56.446247697063394, 59.347044002602353, 60.831778524609809, 65.112544048081606,
         67.079810529494173, 69.546401711173979, 72.067157674481907, 75.704690699083933]   # 2*gamma <= 150: 19 zeros below 75

def parent_weight(kind, par):
    if kind == 'heat':
        return lambda u: mp.e ** (-par * u * u)          # parent = e^{-par u^2} Phi (rewound by par)
    if kind == 'factor':
        return lambda u: 1 / (1 + u * u / par ** 2)
    if kind == 'bath':
        return lambda u: 1 / mp.cosh(par * u)
    if kind == 'child3':                                   # three-level symmetric bath averaged away, outer weight w
        w, J = par
        return lambda u: (1 - 2 * w) + 2 * w * mp.cosh(J * u)
    if kind == 'none':
        return None

def eval_real(args):
    kind, par, z = args
    return float(mp.sign(rk.H(z, 0, parent_weight(kind, par)))) if True else None

def real_zeros(kind, par, pool):
    zs = [i * DZ + DZ / 2 for i in range(int(ZMAX / DZ))]
    sg = pool.map(eval_real, [(kind, par, z) for z in zs], chunksize=50)
    out = [0.5 * (zs[i] + zs[i + 1]) for i in range(len(zs) - 1) if sg[i] != sg[i + 1]]
    return out

def antideriv(z):
    xs, ws = rk.nodes()
    z = mp.mpmathify(z)
    return mp.fsum(w * mp.sin(z * u) / u for u, w in zip(xs, ws))

def eval_anti(z):
    return complex(antideriv(z) - LANTI)

def fused(z, k):
    xs, ws = rk.nodes()
    z = mp.mpmathify(z)
    pr = mp.mpf(1)
    for u, w in zip(xs, ws):
        pr *= 1 + k * w * mp.cos(z * u)
    return (pr - 1) / k

def eval_fused(args):
    z, k = args
    return complex(fused(z, k))

def eval_H(args):
    z, kind, par = args
    return complex(rk.H(z, 0, parent_weight(kind, par)))

def winding(vals):
    import numpy as np
    ph = np.unwrap(np.angle(np.array(vals)))
    return int(round((ph[-1] - ph[0]) / (2 * np.pi)))

def rect(x0, x1, Y, h=0.05):
    import numpy as np
    nb = int((x1 - x0) / h) + 1; nv = int(2 * Y / h) + 1
    return list(np.concatenate([np.linspace(x0, x1, nb) - 1j * Y, x1 + 1j * np.linspace(-Y, Y, nv)[1:],
                                np.linspace(x1, x0, nb)[1:] + 1j * Y, x0 + 1j * np.linspace(Y, -Y, nv)[1:]]))

LANTI = None
if __name__ == '__main__':
    t0 = time.time()
    rk.nodes()
    LANTI = mp.pi / 2 * rk.Phi(0)
    out = {}
    with Pool(3, initializer=rk.nodes) as pool:
        base = real_zeros('none', None, pool)
        out['H0_real_zeros'] = base
        print('H0 real zeros up to 150:', len(base), [round(b / 2, 3) for b in base[:5]], f'{time.time()-t0:.0f}s', flush=True)
        scans = [('heat', v) for v in [0.02, 0.05, 0.1, 0.2, 0.4]] + \
                [('factor', v) for v in [0.25, 0.3, float(mp.pi / 8), 0.45, 0.6, 1.0, 2.0]] + \
                [('bath', v) for v in [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]] + \
                [('child3', (w, 1.0)) for w in [0.05, 0.15, 0.25, 0.35]]
        out['scans'] = {}
        for kind, par in scans:
            rz = real_zeros(kind, par, pool)
            out['scans'][f'{kind}|{par}'] = rz
            print(f'{kind:7s} {str(par):22s} real zeros up to 150: {len(rz):3d} (H0: {len(base)})  first {[round(r,2) for r in rz[:3]]}  {time.time()-t0:.0f}s', flush=True)
            json.dump(out, open('parents.json', 'w'), default=str)
        # argument principle: H0 itself on [0.5,148]x[-6,6]
        R = 0.5 * (base[-1] + base[-2])
        vals = pool.map(eval_H, [(z, 'none', None) for z in rect(0.0, R, 6.0)], chunksize=50)
        out['H0_rect'] = dict(R=R, total=winding(vals), real=sum(1 for b in base if b < R))
        print('H0 rectangle', out['H0_rect'], flush=True)
        # antiderivative parent P(z) - L: left half-plane [-R,0] x [-6,6]
        import numpy as np
        av = pool.map(eval_anti, rect(-R, 0.0, 6.0), chunksize=50)
        out['anti_left'] = dict(R=R, total=winding(av))
        xs = np.arange(-R, -0.5, 0.25)
        neg_vals = pool.map(eval_anti, list(xs), chunksize=50)
        out['anti_left']['real_sign_changes'] = int(np.sum(np.diff(np.sign(np.real(neg_vals))) != 0))
        out['anti_left']['value_far_left'] = float(np.real(neg_vals[0])); out['anti_left']['minus_2L'] = float(-2 * LANTI)
        print('antiderivative parent, left half', out['anti_left'], flush=True)
        av2 = pool.map(eval_anti, rect(0.0, R, 6.0), chunksize=50)
        out['anti_right'] = dict(total=winding(av2))
        print('antiderivative parent, right half', out['anti_right'], flush=True)
        # fused parent (prod(1 + k h_i) - 1)/k on [0,R] x [-Y,Y]
        out['fused'] = {}
        for k in [0.5, 1.0, 2.0]:
            for Y in [3.0, 6.0]:
                fv = pool.map(eval_fused, [(z, k) for z in rect(0.0, R, Y)], chunksize=50)
                out['fused'][f'{k}|{Y}'] = winding(fv)
                print(f'fused k={k} Y={Y}: zeros in rectangle {out["fused"][f"{k}|{Y}"]} (H0 real zeros {out["H0_rect"]["real"]})', flush=True)
        json.dump(out, open('parents.json', 'w'), default=str)
    print('elapsed', time.time() - t0)
