"""The cross-handed corner dial with fractional dissipation.
Helical shell model from ns-rungs/helical_shell.py (golden-ratio shells, 8 sign corners;
weight 1 on the two same-handed corners, eps on the six cross-handed ones), plus
dissipation -nu k_n^(2 alpha) u_n. Outcome: does energy reach shell N-3 at a fixed
fraction of the Kolmogorov level before tmax, for N = 24 and 32?"""
import numpy as np, sys, json, time
sys.path.insert(0, '/mnt/project-files/math-theory/ns-rungs')
from helical_shell import make_rhs, init
from scipy.integrate import solve_ivp

def outcome(eps, alpha, N, nu=1e-3, tmax=60.0, theta=1e-2, seed=0):
    rhs0, k = make_rhs(N, eps)
    d = nu * k ** (2 * alpha)
    D = np.concatenate([d, d, d, d])
    rhs = lambda t, y: rhs0(t, y) - D * y
    y0 = init(N, seed); target = N - 3
    def ev(t, y):
        u = (y[:2*N] + 1j*y[2*N:]).reshape(2, N); return np.sum(abs(u[:, target])**2) - theta * k[target]**(-2/3)
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0, tmax), y0, method='LSODA', rtol=1e-8, atol=1e-14, events=ev)
    U = (sol.y[:2*N] + 1j*sol.y[2*N:]).reshape(2, N, -1); En = np.sum(abs(U)**2, axis=0)
    lvl = (En / (theta * k[:, None]**(-2/3))).max(1)
    deep = int(np.max(np.where(lvl >= 1)[0])) if (lvl >= 1).any() else -1
    return dict(reached=len(sol.t_events[0]) > 0, t=float(sol.t[-1]), deep=deep, E_end=float(En[:, -1].sum()))

if __name__ == '__main__':
    epss = [float(x) if x != 'inf' else float('inf') for x in sys.argv[1].split(',')]
    alphas = [float(x) for x in sys.argv[2].split(',')]
    out = []
    for eps in epss:
        for a in alphas:
            row = []
            for N in [24, 32, 40]:
                t0 = time.time(); o = outcome(eps, a, N); o.update(eps=eps, alpha=a, N=N, secs=time.time() - t0)
                out.append(o); row.append(o)
            print(f"eps={eps} alpha={a}: " + " | ".join(f"N={o['N']} {'REACH' if o['reached'] else 'stop'} t={o['t']:.3f} deep={o['deep']} ({o['secs']:.0f}s)" for o in row), flush=True)
    json.dump(out, open(f"corner_viscous_{sys.argv[1].replace(',','_')}.json", 'w'))
