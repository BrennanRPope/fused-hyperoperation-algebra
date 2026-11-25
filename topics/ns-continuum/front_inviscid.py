"""Inviscid dyadic front: run until energy reaches the top shells (event), record
per-shell peak amplitude and peak time; fit the front exponents."""
import numpy as np, json
from scipy.integrate import solve_ivp
from dyadic import rhs_factory

def front(n, N=30, thresh=1e-3):
    theta = 1 + n / 2; lam = 2 ** theta
    f, _ = rhs_factory(N, lam, 0.0, 0.5)
    a0 = np.zeros(N + 1); a0[0] = 1.0
    # event: shell N-3 reaches thresh times its K41 level lam^{-(N-3)/3}
    jt = N - 3
    ev = lambda t, a: a[jt] - thresh * lam ** (-jt / 3)
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(f, (0, 10), a0, method='DOP853', rtol=1e-12, atol=1e-40, events=ev)
    t, A = sol.t, sol.y
    # first-arrival time: time shell j first reaches 'thresh' of K41 level; amplitude at that point is fixed by def
    # instead use time when a_j first reaches half of its running max before T end, and local max value
    tj, aj = [], []
    for j in range(N - 2):
        x = A[j]; k = np.argmax(x)
        tj.append(t[k]); aj.append(x[k])
    return lam, t[-1], np.array(tj), np.array(aj), sol.status

res = {}
for n in [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]:
    lam, tend, tj, aj, st = front(n)
    # T*: extrapolate geometric convergence of t_j from last shells
    j = np.arange(len(tj))
    d = np.diff(tj)
    ratio = d[-6:-1] / d[-7:-2]
    q = np.median(ratio)
    Tstar = tj[-1] + d[-1] * q / (1 - q)
    # fit in the middle range
    sel = slice(6, len(tj) - 2)
    beta = -np.polyfit(j[sel], np.log(aj[sel]) / np.log(lam), 1)[0]
    tau = -np.polyfit(j[sel], np.log(Tstar - tj[sel]) / np.log(lam), 1)[0]
    res[n] = dict(lam=lam, Tstar=Tstar, q=q, beta=beta, tau=tau, tend=tend, status=st,
                  tj=tj.tolist(), aj=aj.tolist())
    print(f"n={n:3.1f} lam={lam:7.3f} T*~{Tstar:.10f} step-ratio={q:.4f} beta(amp exp)={beta:.4f} "
          f"tau(time exp)={tau:.4f} beta+tau={beta+tau:.4f} gam*=(1-beta)/2={(1-beta)/2:.4f}")
json.dump(res, open('front_inviscid.json', 'w'))
