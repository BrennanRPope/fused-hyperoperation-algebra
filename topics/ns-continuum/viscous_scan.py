"""Viscous dyadic model along the dimension continuum at true dissipation alpha.
gam = alpha/theta, theta = 1+n/2, lam = 2^theta. Renormalized time, LSODA.
Outcome per (n, alpha, N, nu): 'reached' if the cascade reaches shell N-6 at 1e-3 of the
K41 level before t=Tmax, else 'stopped'. Blow-up signature: reached for every N, with
arrival times converging as N grows."""
import numpy as np, json, sys, time
from scipy.integrate import solve_ivp

def make(N, lam, nu, gam):
    lj = lam ** np.arange(N + 1); diss = nu * lam ** (2 * gam * np.arange(N + 1))
    def g(t, y):
        a = y[:-1]
        ap = np.zeros_like(a)
        ap[1:] += lj[:-1] * a[:-1] ** 2
        ap[:-1] -= lj[:-1] * a[:-1] * a[1:]
        ap -= diss * a
        return np.concatenate([ap, [1.0]])
    def jac(t, y):
        a = y[:-1]; N1 = len(a); J = np.zeros((N1 + 1, N1 + 1)); i = np.arange(N1)
        J[i, i] = -diss; J[i[:-1], i[:-1]] -= lj[:-1] * a[1:]; J[i[:-1], i[1:]] -= lj[:-1] * a[:-1]
        J[i[1:], i[:-1]] += 2 * lj[:-1] * a[:-1]
        return J
    g.jac = jac
    return g, lj

def outcome(n, alpha, N, nu, Tmax=20.0):
    theta = 1 + n / 2; lam = 2 ** theta; gam = alpha / theta
    g, lj = make(N, lam, nu, gam)
    y0 = np.zeros(N + 2); y0[0] = 1.0
    jt = N - 6
    ev1 = lambda s, y: y[jt] * lj[jt] ** (1/3) - 1e-3; ev1.terminal = True; ev1.direction = 1
    sol = solve_ivp(g, (0, Tmax), y0, method='Radau', jac=g.jac, rtol=1e-8, atol=1e-40, events=[ev1])
    reached = len(sol.t_events[0]) > 0
    tfin = sol.y[-1, -1]
    A = sol.y[:-1]
    kmax = int(np.max(np.where((np.abs(A) * lj[:, None] ** (1/3)).max(1) > 1e-3)[0])) if True else 0
    return dict(reached=bool(reached), t=float(tfin), deepest=kmax, status=int(sol.status))

if __name__ == '__main__':
    alpha = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    nus = [1e-2, 1e-4]
    res = []
    for n in [2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]:
        for nu in nus:
            row = []
            for N in [16, 24, 32]:
                t0 = time.time()
                o = outcome(n, alpha, N, nu)
                o.update(n=n, alpha=alpha, N=N, nu=nu, gam=alpha / (1 + n / 2), secs=time.time() - t0)
                res.append(o); row.append(o)
            print(f"n={n} gam={alpha/(1+n/2):.3f} nu={nu:g}: " + " | ".join(
                f"N={o['N']} {'REACH' if o['reached'] else 'stop'} t={o['t']:.5f} deep={o['deepest']} ({o['secs']:.0f}s)" for o in row), flush=True)
    json.dump(res, open(f'viscous_scan_alpha{alpha}.json', 'w'))
