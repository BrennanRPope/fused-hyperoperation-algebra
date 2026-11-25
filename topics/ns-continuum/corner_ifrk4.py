"""Corner dial with dissipation, integrated by integrating-factor RK4 (Lawson) so the
dissipation k^(2 alpha) is exact and only the nonlinear rate sets the step.
Same model and outcome test as corner_viscous.py."""
import numpy as np, sys, json, time
sys.path.insert(0, '/mnt/project-files/math-theory/ns-rungs')
from helical_shell import make_rhs, init

def outcome(eps, alpha, N, nu=1e-3, tmax=40.0, theta=1e-2, seed=0, cfl=0.05):
    rhs0, k = make_rhs(N, eps)
    d = nu * k ** (2 * alpha); D = np.concatenate([d, d, d, d])
    y = init(N, seed); t = 0.0; target = N - 3
    lvl_t = theta * k ** (-2 / 3)
    deepest = -1; reached = False
    kk = np.concatenate([k, k])
    while t < tmax:
        u = (y[:2*N] + 1j*y[2*N:]).reshape(2, N)
        rate = np.max(kk * np.abs(u).ravel()) * 6 + 1e-12
        h = min(cfl / rate, tmax - t)
        E1 = np.exp(-D * h / 2); E2 = E1 * E1
        k1 = rhs0(t, y)
        k2 = rhs0(t, E1 * (y + h / 2 * k1))
        k3 = rhs0(t, E1 * y + h / 2 * k2)
        k4 = rhs0(t, E2 * y + h * E1 * k3)
        y = E2 * y + h / 6 * (E2 * k1 + 2 * E1 * (k2 + k3) + k4)
        t += h
        En = np.abs(u[0])**2 + np.abs(u[1])**2
        above = np.where(En >= lvl_t)[0]
        if len(above): deepest = max(deepest, int(above.max()))
        if En[target] >= lvl_t[target]: reached = True; break
        if not np.isfinite(y).all(): break
    return dict(reached=reached, t=t, deep=deepest)

if __name__ == '__main__':
    epss = [float(x) if x != 'inf' else float('inf') for x in sys.argv[1].split(',')]
    alphas = [float(x) for x in sys.argv[2].split(',')]
    TMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 40.0
    Ns = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [24, 32, 40]
    out = []
    for eps in epss:
        for a in alphas:
            row = []
            for N in Ns:
                t0 = time.time(); o = outcome(eps, a, N, tmax=TMAX); o.update(eps=eps, alpha=a, N=N, secs=time.time() - t0)
                out.append(o); row.append(o)
            print(f"eps={eps} alpha={a}: " + " | ".join(f"N={o['N']} {'REACH' if o['reached'] else 'stop'} t={o['t']:.3f} deep={o['deep']} ({o['secs']:.0f}s)" for o in row), flush=True)
    json.dump(out, open(f"corner_ifrk4_{sys.argv[1].replace(',','_')}_{sys.argv[2].replace(',','_')}.json", 'w'))
