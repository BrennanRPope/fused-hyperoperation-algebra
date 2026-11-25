"""Inviscid (or viscous) dyadic model in renormalized time ds = M dt, M = |(lam^j a_j)|_2,
so steps stay O(1) as the cascade speeds up. Records the active shell and its amplitude."""
import numpy as np, json, sys
from scipy.integrate import solve_ivp

def make(N, lam, nu, gam):
    lj = lam ** np.arange(N + 1); diss = nu * lam ** (2 * gam * np.arange(N + 1))
    def g(s, y):
        a = y[:-1]
        ap = np.zeros_like(a)
        ap[1:] += lj[:-1] * a[:-1] ** 2
        ap[:-1] -= lj[:-1] * a[:-1] * a[1:]
        ap -= diss * a
        M = np.sqrt(np.sum((lj * a) ** 2) + (diss * 0 + 1e-300).sum())
        return np.concatenate([ap / M, [1.0 / M]])
    return g, lj

def front(n, N=40, nu=0.0, gam=0.5, S=400.0, a00=1.0):
    lam = 2 ** (1 + n / 2)
    g, lj = make(N, lam, nu, gam)
    y0 = np.zeros(N + 2); y0[0] = a00
    jt = N - 6
    ev = lambda s, y: y[jt] * lj[jt] ** (1/3) - 1e-3; ev.terminal = True; ev.direction = 1
    sol = solve_ivp(g, (0, S), y0, method='DOP853', rtol=1e-10, atol=1e-80, events=ev, max_step=0.5)
    A = sol.y[:-1]; t = sol.y[-1]
    act = np.argmax(lj[:, None] * np.abs(A), axis=0)
    return lam, sol, A, t, act

if __name__ == '__main__':
    res = {}
    for n in [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]:
        lam, sol, A, t, act = front(n)
        N = A.shape[0] - 1
        # per shell: amplitude when it is the active shell (max over those times)
        js, amp, tt = [], [], []
        for j in range(3, N - 6):
            m = act == j
            if m.any():
                k = np.where(m)[0]; kk = k[np.argmax(A[j, k])]
                js.append(j); amp.append(A[j, kk]); tt.append(t[kk])
        js = np.array(js); amp = np.array(amp); tt = np.array(tt)
        sel = js >= 8
        beta = -np.polyfit(js[sel], np.log(amp[sel]) / np.log(lam), 1)[0]
        # T* from final t plus remaining ~ 1/M
        res[n] = dict(lam=lam, beta=beta, js=js.tolist(), amp=amp.tolist(), t=tt.tolist(), tend=float(t[-1]))
        print(f"n={n:3.1f} lam={lam:7.3f} t_end={t[-1]:.10f} beta={beta:.4f}  gam*=(1-beta)/2={(1-beta)/2:.4f}"
              f"  -> alpha*(n)={(1-beta)/2*(1+n/2):.4f}", flush=True)
    json.dump(res, open('front_inviscid.json', 'w'))
