"""Round 8, item 1: where the Fisher zeros of the order-dial prime gas go as q moves.

Z_q(beta) = (q;q)_inf sum_k q^k/(q;q)_k G_k,  log G_k = sum_N b_k(N)/N log zeta(N beta)   (round 5/7).
Near a Riemann zero rho write w = zeta(beta).  G_k ~ c_k w^{q^k/(1-q)}, c_k regular, c_1/c_0 -> 1 as q -> 0.
The two leading terms balance when  w + q w^{q} ~ 0, i.e. w^{1-q} = -q:  to first order the Fisher
zeros are the a-points of zeta with a = -q:  zeta(beta) = -q,  beta ~ rho - q/zeta'(rho).
This script tracks the true Fisher zeros of Z_q and compares them with the a-points and the linear guess.
"""
import json, math, numpy as np, mpmath as mp
from sympy import divisors, mobius
mp.mp.dps = 20

def bk_table(q, K, NMAX):
    return [[sum(int(mobius(N//m))*q**(k*m)/(1-q**m) for m in divisors(N)) if N > 0 else 0
             for N in range(NMAX+1)] for k in range(K+1)]
def qp(q, n):
    r = mp.mpf(1)
    for j in range(1, n+1): r *= (1-q**j)
    return r

def logzeta_path(beta, scale=1, start=4.0, nstep=600):
    s0 = mp.mpf(start); sig = mp.re(beta); t = mp.im(beta)
    path = [s0 + (sig - s0)*j/nstep for j in range(nstep+1)]
    vals = [mp.zeta(scale*mp.mpc(s, t)) for s in path]
    ph = np.unwrap([float(mp.arg(v)) for v in vals])
    return mp.log(abs(vals[-1])) + 1j*mp.mpf(ph[-1])

class Local:
    """continued log zeta near rho, horizontal-cut convention (cut runs left from rho)."""
    def __init__(self, rho):
        self.rho = rho; self.zp = mp.zeta(rho, derivative=1)
        d = mp.mpf('0.05')
        Lr = logzeta_path(rho + d)
        self.Lrho = Lr - mp.log(d) - mp.log(mp.zeta(rho+d)/(self.zp*d))
    def __call__(self, beta):
        x = beta - self.rho
        return self.Lrho + mp.log(x) + mp.log(mp.zeta(beta)/(self.zp*x))

def Zq(beta, q, loc, NMAX=30):
    K = int(36/(-math.log(abs(q)))) + 2
    B = bk_table(q, K, NMAX)
    L = {1: loc(beta), 2: logzeta_path(beta, 2, nstep=200)}
    for N in range(3, NMAX+1): L[N] = mp.log(mp.zeta(N*beta))
    tot = 0
    for k in range(K+1):
        e = sum(B[k][N]/N * L[N] for N in range(1, NMAX+1))
        tot += (q**k/qp(q, k))*mp.e**e
    return qp(q, 300)*tot

if __name__ == '__main__':
    # sanity: local log agrees with the path log at a point right of rho
    rho = mp.zetazero(1); loc = Local(rho)
    b = rho + mp.mpc(0.2, 0.05)
    print("local vs path log at rho+0.2+0.05i:", mp.nstr(loc(b), 10), mp.nstr(logzeta_path(b), 10))
    # round-7 value check: q=-0.5 first Fisher zero was 0.905308+14.068964i
    r = mp.findroot(lambda x: Zq(x, -0.5, loc), mp.mpc(0.9, 14.07))
    print("q=-0.5 first Fisher zero (round 7: 0.905308+14.068964i):", mp.nstr(r, 10))
    qs_neg = [-0.02, -0.05, -0.1, -0.2, -0.35, -0.5, -0.65, -0.8, -0.9]
    qs_pos = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8]
    out = {}
    for n in (1, 2, 3):
        rho = mp.zetazero(n); loc = Local(rho); zp = loc.zp
        rows = []
        for qs in (qs_neg, qs_pos):
            guess = None
            for q in qs:
                lin = rho - q/zp
                try:
                    ap = mp.findroot(lambda x: mp.zeta(x) + q, lin)
                except Exception:
                    ap = None
                g = guess if guess is not None else (ap if ap is not None else lin)
                try:
                    fz = mp.findroot(lambda x: Zq(x, q, loc), g)
                    if abs(fz - rho) > 0.8: fz = None
                except Exception:
                    fz = None
                if fz is not None: guess = fz
                row = {'q': q, 'linear': [float(mp.re(lin)), float(mp.im(lin))],
                       'apoint': None if ap is None else [float(mp.re(ap)), float(mp.im(ap))],
                       'fisher': None if fz is None else [float(mp.re(fz)), float(mp.im(fz))]}
                rows.append(row)
                print(f"rho_{n}  q={q:+.2f}  Fisher {row['fisher']}  a-point {row['apoint']}  linear {row['linear']}")
        out[str(n)] = {'rho': [0.5, float(mp.im(rho))], 'zeta_prime': [float(mp.re(zp)), float(mp.im(zp))], 'rows': rows}
    json.dump(out, open('qgas.json', 'w'))
    print('saved')
