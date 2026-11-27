"""RH line = l2 edge (mean-square abscissa) of the level weights; transition = l1 edge. Three models."""
import numpy as np, json, math
from knauf import levels
out = {}
X = np.unique(np.logspace(1, 4.2, 40).astype(int))
# zeta: a_n=1, line 1/2
n = np.arange(1, X.max()+1); S = np.cumsum(1.0/n)
out['zeta'] = dict(X=X.tolist(), S=S[X-1].tolist())
# Knauf: a_q = phi(q) (Euler phi), line 3/2 ; use exact totient
N = X.max(); ph = np.arange(N+1, dtype=float)
for p in range(2, N+1):
    if ph[p] == p: ph[p::p] *= (1-1/p)
q = np.arange(1, N+1); Sk = np.cumsum(ph[1:]**2/q**3.0)
out['knauf'] = dict(X=X.tolist(), S=Sk[X-1].tolist())
# check chain multiplicities equal Euler phi for q <= 2^n-ish (q<=n+1 safe region)
qq, cc = levels(20); m = qq <= 21
out['chain_mult_is_totient'] = bool(np.allclose(cc[m], ph[qq[m].astype(int)]))
# Ramanujan tau via eta^24
M = 6001
eta = np.zeros(M, dtype=object); eta[:] = 0
k = 0
while True:
    done = True
    for kk in (k, -k) if k else (0,):
        e = kk*(3*kk-1)//2
        if e < M: eta[e] += (-1)**abs(kk); done = False
    if done and k > 0: break
    k += 1
def mul(a, b): return np.convolve(a, b)[:M]
e2 = mul(eta, eta); e4 = mul(e2, e2); e8 = mul(e4, e4); e16 = mul(e8, e8); e24 = mul(e16, e8)
tau = np.zeros(M, dtype=object); tau[1:] = e24[:M-1]   # Delta = q * eta^24
assert tau[1] == 1 and tau[2] == -24 and tau[3] == 252 and tau[12] == -370944
Xt = X[X < M]; nt = np.arange(1, M)
St = np.cumsum([float(tau[i])**2/float(i)**12 for i in nt])
out['tau'] = dict(X=Xt.tolist(), S=St[Xt-1].tolist())
for key in ('zeta', 'knauf', 'tau'):
    Xs = np.array(out[key]['X'], float); Ss = np.array(out[key]['S'])
    sl = np.polyfit(np.log(Xs[Xs > 300]), Ss[Xs > 300], 1)[0]
    out[key]['log_slope'] = float(sl); print(key, "slope of S vs ln X:", round(sl, 4))
C = np.prod([1-2/p**2+1/p**3 for p in range(2, 200000) if all(p % d for d in range(2, int(p**0.5)+1))]) if False else None
json.dump(out, open("l2edge.json", "w")); print("chain multiplicity = totient:", out['chain_mult_is_totient'])
