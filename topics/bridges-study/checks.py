# DRAFT. Bridges study (Nov 28, 2026): checks behind the shared-shape claims.
# Each check is a small exact or numeric identity; no new physics is claimed by passing them.
import json, numpy as np, sympy as sp
rng = np.random.default_rng(28)
def ts(p, q):
    p = np.asarray(p, float).ravel(); p = p[p > 1e-300]
    return float(-(p*np.log(p)).sum()) if abs(q-1) < 1e-12 else float((1-(p**q).sum())/(q-1))
def re(p, a):
    p = np.asarray(p, float).ravel(); p = p[p > 1e-300]
    return float(-(p*np.log(p)).sum()) if abs(a-1) < 1e-12 else float(np.log((p**a).sum())/(1-a))
def rand_p(n): x = rng.random(n); return x/x.sum()
def rand_rho(d):
    X = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)); R = X@X.conj().T; return R/np.trace(R).real
res = {}; out = []
def log(name, ok, note): res[name] = bool(ok); out.append(f"{name} {'PASS' if ok else 'FAIL'}  {note}")
qs = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
# B1 Tsallis entropy of independent pieces composes by the law with k = 1 - q (Tsallis 1988).
err = max(abs(ts(np.outer(a, b), q) - (ts(a, q)+ts(b, q)+(1-q)*ts(a, q)*ts(b, q)))
          for q in qs for a, b in [(rand_p(3), rand_p(4)) for _ in range(30)])
log("B1", err < 1e-12, f"Tsallis composes by a+b+k*ab, k=1-q; max err {err:.1e}")
# B2 Renyi entropy adds on independent pieces at every order (the law's unlifted coordinate).
err = max(abs(re(np.outer(a, b), q) - re(a, q) - re(b, q)) for q in qs for a, b in [(rand_p(3), rand_p(5)) for _ in range(30)])
log("B2", err < 1e-12, f"Renyi adds at every order; max err {err:.1e}")
# B3 The law's lift is the replica trace: 1 + (1-n) S_T,n = Tr rho^n = exp((1-n) S_R,n), quantum states.
err = 0.0
for n in [0.5, 2.0, 3.0]:
    for _ in range(20):
        ev = np.clip(np.linalg.eigvalsh(rand_rho(4)), 0, None); tr = (ev**n).sum()
        err = max(err, abs(1+(1-n)*ts(ev, n)-tr), abs(np.exp((1-n)*re(ev, n))-tr))
log("B3", err < 1e-12, f"lift 1+kS_T = Tr rho^n = Z_n/Z_1^n = exp(k S_R); max err {err:.1e}")
# B4 Dong 2016: cosmic-brane tension T_n = (n-1)/(4nG) and deficit 2pi(n-1)/n = 8piG T_n; zero only at n = 1.
n, G, k = sp.symbols('n G k', positive=True)
T = (n-1)/(4*n*G)
ok4 = sp.simplify(8*sp.pi*G*T - 2*sp.pi*(n-1)/n) == 0 and sp.solve(sp.Eq(T, 0), n) == [1]
Tk = sp.simplify(T.subs(n, 1-k))
log("B4", ok4, f"Dong tension zero only at n=1; in the law's dial k=1-n: T = {Tk}")
# B5 Temperature from the i-turn: smooth Euclidean horizon needs period 4pi/f'(r+) -> T = kappa/2pi.
r, M, a = sp.symbols('r M a', positive=True)
f = 1 - 2*M/r; rp = 2*M; kappa = sp.diff(f, r).subs(r, rp)/2
ok5a = sp.simplify(4*sp.pi/sp.diff(f, r).subs(r, rp) - 8*sp.pi*M) == 0 and sp.simplify(kappa/(2*sp.pi) - 1/(8*sp.pi*M)) == 0
rpk = M + sp.sqrt(M**2-a**2); rmk = M - sp.sqrt(M**2-a**2)
Tkerr = (rpk-rmk)/(4*sp.pi*(rpk**2+a**2))
ok5b = sp.simplify(Tkerr.subs(a, 0) - 1/(8*sp.pi*M)) == 0
log("B5", ok5a and ok5b, "Hawking/Unruh T = kappa/2pi from imaginary-time period; Kerr T reduces to 1/(8 pi M) at a=0")
# B6 Jacobson-style: horizon entropy gamma*A^delta gives 1/(4 G_eff) = dS/dA, constant only at delta = 1.
A, gam, dl = sp.symbols('A gamma delta', positive=True)
Geff = 1/(4*sp.diff(gam*A**dl, A))
ok6 = sp.simplify(sp.diff(Geff.subs(dl, 1), A)) == 0 and sp.simplify(sp.diff(Geff, A)) != 0
log("B6", ok6, f"G_eff = {sp.simplify(Geff)}: constant only at delta=1 (area law)")
# B7 Setup for candidate C1: for a correlated pair the law's composition picks up a cross term.
#   k_eff from S_AB = S_A + S_B + k_eff S_A S_B equals 1-q for products, and moves once A,B are correlated.
def keff(P, q):
    sa, sb, sab = ts(P.sum(1), q), ts(P.sum(0), q), ts(P, q)
    return (sab-sa-sb)/(sa*sb)
q = 2.0; prod = np.outer([.5, .5], [.5, .5]); dev = []
for c in [0.0, 0.1, 0.3, 0.6]:
    P = (1-c)*prod + c*np.diag([.5, .5]); dev.append(keff(P, q)-(1-q))
ok7 = abs(dev[0]) < 1e-12 and all(abs(d) > 1e-3 for d in dev[1:]) and all(np.diff(np.abs(dev)) > 0)
log("B7", ok7, "k_eff = 1-q exactly for independent pieces; drifts as correlation grows: " + ", ".join(f"{d:+.3f}" for d in dev))
res["B7_drift"] = dev
# B8 COW: interferometer phase m^2 g lambda A /(2 pi hbar^2) = (m g H)(L/v)/hbar, phase = action/hbar.
m, g, H, L, v, hb = sp.symbols('m g H L v hbar', positive=True)
lam = 2*sp.pi*hb/(m*v)
ok8 = sp.simplify(m**2*g*lam*(H*L)/(2*sp.pi*hb**2) - m*g*H*(L/v)/hb) == 0
log("B8", ok8, "COW phase = (potential difference x time)/hbar; Newtonian gravity inside Schrodinger, no G-hbar mixing")
n_pass = sum(1 for k2, v2 in res.items() if isinstance(v2, bool) and v2)
n_all = sum(1 for v2 in res.values() if isinstance(v2, bool))
out.append(f"{n_pass}/{n_all} passed")
print("\n".join(out))
json.dump(res, open("results.json", "w"), indent=1)
