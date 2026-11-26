"""Round 9: what the zeros teach the string, and the inverse-spectral trap.

(a) Landau's formula: sum over zeros of cos(gamma * ell) peaks at ell = log p^m.
(b) Free weights: fit one weight per prime power (lengths fixed at log p^m) to S(t) on the
    training range; compare with the truth p^{-m/2}/(m pi); test the levels.
(c) Free frequencies: fit K sines with frequencies on a fixed grid (lengths NOT given);
    training error falls to zero as K grows, held-out error grows: the trap.
"""
import json
import numpy as np
from math import pi, log
from stringlib import String, true_S, score, primes_upto, gram, theta

gam = np.array(json.load(open('zeros.json'))['gamma'])
NZ = len(gam)
TR = np.arange(0, 500); TE = np.arange(500, NZ)
OUT = {}

# (a) Landau: F(ell) = sum_{gamma<=T} cos(gamma ell) weighted; peaks at log n with -T Lambda(n)/(2 pi sqrt n)
ell = np.linspace(0.05, 4.2, 4000)
for lab, G in [('train', gam[:500]), ('all', gam)]:
    w = np.exp(-(G / G[-1]) ** 2 * 2)          # gentle taper
    Fl = (np.cos(np.outer(ell, G)) * w).sum(1)
    OUT[f'landau_{lab}'] = Fl.tolist()
    # find the 8 deepest minima
    idx = np.argsort(Fl)[:400]
    mins = []
    for i in idx:
        if all(abs(ell[i] - ell[j]) > 0.06 for j in mins):
            mins.append(i)
        if len(mins) == 10:
            break
    found = sorted(float(ell[i]) for i in mins)
    OUT[f'landau_minima_{lab}'] = found
    print(f"Landau sum over {lab} zeros: deepest dips at", [round(v, 3) for v in found])
OUT['ell'] = ell.tolist()
OUT['log_pm'] = sorted(float(log(n)) for n in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61])

# grid of t for fitting S
tt_tr = np.arange(12.0, gam[499] + 0.3, 0.05)
S_tr = true_S(tt_tr, gam)

# (b) free weights on prime powers n <= 100
base = String(X=100.0, taper='sharp')
L = base.L; n_int = np.round(np.exp(L)).astype(int)
Phi = -np.sin(np.outer(tt_tr, L)) / pi
coef, *_ = np.linalg.lstsq(Phi, S_tr, rcond=None)
m_of = np.array([round(log(n) / log(min(p for p in primes_upto(n) if n % p == 0))) for n in n_int])
truthw = n_int ** -0.5 / m_of
OUT['free_weights'] = dict(n=n_int.tolist(), fitted=coef.tolist(), truth=truthw.tolist())
rel = coef / truthw
print("fitted weight / p^{-m/2}/m, first primes:", [(int(n), round(float(r), 2)) for n, r in zip(n_int[:10], rel[:10])])
OUT['free_weight_ratio_median'] = float(np.median(rel))
st = String(X=100.0, taper='sharp'); st.w = coef.copy()
E = st.levels(NZ)
OUT['free_weight_score'] = dict(train=score(E, gam, TR), test=score(E, gam, TE))
st0 = String(X=100.0, taper='sharp'); E0 = st0.levels(NZ)
OUT['truth_weight_score_X100_sharp'] = dict(train=score(E0, gam, TR), test=score(E0, gam, TE))
print("free weights: train %.4f test %.4f | truth weights same orbits: train %.4f test %.4f" % (
    OUT['free_weight_score']['train']['rms'], OUT['free_weight_score']['test']['rms'],
    OUT['truth_weight_score_X100_sharp']['train']['rms'], OUT['truth_weight_score_X100_sharp']['test']['rms']))

# (c) free frequencies: K frequencies on [0.02, 4.6], sin and cos, least squares (tiny ridge)
trap = []
for K in [5, 10, 20, 40, 80, 160, 320, 640, 1280]:
    f = np.linspace(0.02, 4.6, K)
    A = np.hstack([np.sin(np.outer(tt_tr, f)), np.cos(np.outer(tt_tr, f))])
    lam = 1e-8
    c = np.linalg.solve(A.T @ A + lam * np.eye(2 * K), A.T @ S_tr)
    amps = np.concatenate([c[:K], c[K:]])
    freqs = np.concatenate([f, f])
    ph = np.concatenate([np.zeros(K), np.full(K, pi / 2)])
    s_ = String(X=10.0); s_.free = (amps, freqs, ph)
    E = s_.levels(NZ)
    sc = dict(train=score(E, gam, TR), test=score(E, gam, TE))
    fit_err = float(np.sqrt(np.mean((A @ c - S_tr) ** 2)))
    trap.append(dict(K=K, params=2 * K, S_fit_rms=fit_err, train=sc['train'], test=sc['test']))
    print(f"K={K:5d} ({2*K} settings): S fit rms {fit_err:.3f}  levels train {sc['train']['rms']:.4f}  "
          f"test {sc['test']['rms']:.4f}  missing {sc['test']['missing']:.3f}")
OUT['trap'] = trap

# (d) the exact trap: a diagonal operator with the 500 training zeros as eigenvalues
#     fits training perfectly and says nothing about the rest; its best guess beyond is
#     the smooth law, which is the Gram midpoint.
json.dump(OUT, open('learn.json', 'w'))
