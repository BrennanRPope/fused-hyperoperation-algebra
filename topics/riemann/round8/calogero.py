"""Round 8, item 2: the backward heat flow (rewind) on zeros, read as a Calogero-Moser system.

Convention (de Bruijn): dH/dt = -d^2H/dz^2.  t increasing = forward heat (zeros become real);
t decreasing = rewind / backward heat flow.  At a simple zero, P'(z_k) zdot + dP/dt = 0 gives
zdot_k = P''/P' (z_k) = 2 sum_{j!=k} 1/(z_k - z_j)            (first-order rational CM, Sokal)
Claim checked here: zddot_k = -8 sum_{j!=k} 1/(z_k-z_j)^3      (second-order rational CM, g^2 = -4)
and the CM Lax matrix  L_kk = zdot_k,  L_jk = 2/(z_j - z_k)... is nilpotent along the flow,
so z(t) = eigenvalues of Z0 + t L0.
"""
import json, numpy as np, sympy as sp, mpmath as mp
from numpy.polynomial import polynomial as Pl

out = {}
# --- symbolic identity for n = 3, 4
for n in (3, 4):
    z = sp.symbols(f'z1:{n+1}')
    v = [2*sum(1/(z[k]-z[j]) for j in range(n) if j != k) for k in range(n)]
    acc = [sum(-2*(v[k]-v[j])/(z[k]-z[j])**2 for j in range(n) if j != k) for k in range(n)]
    cm = [-8*sum(1/(z[k]-z[j])**3 for j in range(n) if j != k) for k in range(n)]
    ok = all(sp.simplify(acc[k]-cm[k]) == 0 for k in range(n))
    print(f"symbolic: zddot = -8 sum 1/(zk-zj)^3 for n={n}: {ok}")
    out[f'cm_identity_n{n}'] = bool(ok)

def heat_poly(c, t):
    """apply exp(-t D^2) to polynomial coefficients c (ascending)."""
    n = len(c)-1
    res = np.zeros(n+1, dtype=complex)
    term = np.array(c, dtype=complex)
    fact = 1.0
    for m in range(0, n//2+1):
        if m > 0:
            term = Pl.polyder(term, 2)
            fact *= m
        res[:len(term)] += ((-t)**m/fact) * term
    return res

def lax(z):
    n = len(z)
    L = np.zeros((n, n), dtype=complex)
    for k in range(n):
        for j in range(n):
            if j != k:
                L[k, k] += 2/(z[k]-z[j])
                L[j, k] = 2/(z[j]-z[k])
    return L

# --- test on zeros near the first zeta zeros: symmetric polynomial with roots +-gamma_1..gamma_M
M = 12
gam = [float(mp.im(mp.zetazero(n))) for n in range(1, M+1)]
roots0 = np.array(sorted([-g for g in gam] + gam))
c0 = Pl.polyfromroots(roots0)
# scale-free check of Lax nilpotency and the eigenvalue formula
L0 = lax(roots0.astype(complex))
ev = np.linalg.eigvals(L0)
print("Lax eigenvalues at t=0, max |.| =", np.abs(ev).max(), " (||L0|| =", np.abs(L0).max(), ")")
out['lax_eig_max'] = float(np.abs(ev).max())
Lp = np.linalg.matrix_power(L0, len(roots0))
print("||L0^n|| =", np.abs(Lp).max())
out['lax_pow_n'] = float(np.abs(Lp).max())
# zeros of exp(-tD^2)P versus eigenvalues of Z0 + t L0
errs = []
for t in [0.5, -0.5, 2.0, -2.0, -6.0]:
    r_heat = np.sort_complex(Pl.polyroots(heat_poly(c0, t)))
    r_eig = np.sort_complex(np.linalg.eigvals(np.diag(roots0) + t*L0))
    # match
    e = max(min(abs(a-b) for b in r_eig) for a in r_heat)
    errs.append([t, float(e)])
    print(f"t = {t:5}: max |zero - eigenvalue of Z0 + t L0| = {e:.2e}")
out['eig_formula_errors'] = errs

# --- rewind trajectories on a window of zeta zeros (centred, both signs), collision times
M = 30
gam = [float(mp.im(mp.zetazero(n))) for n in range(1, M+1)]
roots0 = np.array(sorted([-g for g in gam] + gam))
L0 = lax(roots0.astype(complex))
ts = np.linspace(0, -12, 1201)
traj = []
for t in ts:
    traj.append(np.linalg.eigvals(np.diag(roots0) + t*L0))
traj = np.array(traj)
# order by continuity
ordered = [np.sort_complex(traj[0])]
for i in range(1, len(ts)):
    prev = ordered[-1]; cur = list(traj[i]); row = []
    for p in prev:
        j = int(np.argmin([abs(p-c) for c in cur])); row.append(cur.pop(j))
    ordered.append(np.array(row))
ordered = np.array(ordered)
im_max = np.abs(ordered.imag).max(axis=1)
first = ts[np.argmax(im_max > 1e-6)] if (im_max > 1e-6).any() else None
gaps = np.diff(gam)
k = int(np.argmin(gaps))
print(f"closest pair among first {M}: gamma_{k+1}, gamma_{k+2}, gap {gaps[k]:.4f}; naive pair estimate d^2/8 = {gaps[k]**2/8:.4f}")
print(f"rewind: first pair leaves the real line at t = {first}")
out.update({'gammas30': gam, 'rewind_t': ts.tolist(), 'traj_re': ordered.real.tolist(), 'traj_im': ordered.imag.tolist(),
            'first_complex_t': None if first is None else float(first), 'closest_gap': float(gaps[k]), 'closest_index': k+1})
json.dump(out, open('calogero.json', 'w'))
print('saved')
