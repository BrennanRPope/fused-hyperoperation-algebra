"""Part D and E: three kinds of 'non-repeating but ordered' and what each one's
spectrum is built from.
 (1) the hat patch: Bragg peaks (pure point), on a finitely generated set of positions;
 (2) the Fibonacci chain (Dyson's Pisot example, the 1D cousin of the hat): peaks at (a + b*phi)/sqrt5-type points;
 (3) the square-free integers (a sieve by every prime): peaks at fractions a/q, height a product over p | q;
 (4) the zeta zeros: the sum of cos(gamma x) dips at x = log p^m (the explicit formula)."""
import json, numpy as np, sympy as sp
import hatgen as h

out = {}
# ---------- (1) hat patch diffraction ----------
def keep_hats(level):
    L = h.levels(level); patch = h.constructPatch(*L[level])
    KEEP = set(h.NEW_H_CH + h.NEW_P_CH + h.NEW_F_CH + h.NEW_T_CH)
    return [T for i, (T, g) in enumerate(patch.children) if i in KEEP for T, k in g.hats(T)]

def hexcentre_coords(Ts):
    anc = np.array([h.transPt(T, h.HAT[0]) for T in Ts]); o = anc[0]
    xy = np.round(np.array([h.to_hex(p - o) for p in anc])).astype(np.int64)
    X2, Y2 = 2*xy[:, 0] + xy[:, 1], xy[:, 1] - xy[:, 0]
    assert (X2 % 6 == 0).all() and (Y2 % 6 == 0).all()
    return np.stack([X2 // 6, Y2 // 6], 1)

def diffraction(U, M=1024):
    U = U - U.min(0)
    assert U.max() < M
    g = np.zeros((M, M)); np.add.at(g, (U[:, 0], U[:, 1]), 1.0)
    F = np.fft.fft2(g)
    return np.abs(F)**2 / len(U)**2   # Bragg peak heights tend to a constant (pure point)

peaks = {}
for lev, M in ((3, 512), (4, 1024)):
    U = hexcentre_coords(keep_hats(lev))
    I = diffraction(U, M)
    np.save(f'hat_diffraction_L{lev}.npy', I.astype(np.float32)) if lev == 4 else None
    # strongest peaks away from the origin (h,k as fractions of a reciprocal cell)
    flat = np.argsort(I.ravel())[::-1][:40]
    pk = [(float(i // M) / M, float(i % M) / M, float(I.ravel()[i])) for i in flat]
    peaks[lev] = dict(n=len(U), M=M, top=pk)
out['hat_peaks'] = {str(k): v for k, v in peaks.items()}
# pure-point test: a Bragg peak keeps the same height I/N^2 as the patch grows ~7x.
# Find the strongest peaks of the big patch away from the central (patch-shape) peak,
# refine each to sub-grid accuracy by direct sums, then measure both patches there.
U3, U4 = hexcentre_coords(keep_hats(3)), hexcentre_coords(keep_hats(4))
I4 = np.load('hat_diffraction_L4.npy'); M4 = 1024
hh, kk_ = np.meshgrid(np.arange(M4), np.arange(M4), indexing='ij')
dh = np.minimum(hh, M4 - hh); dk = np.minimum(kk_, M4 - kk_)
mask = (dh**2 + dk**2) > 40**2
cand = np.argsort(np.where(mask, I4, 0).ravel())[::-1]
chosen = []
for i in cand:
    a, b = divmod(int(i), M4)
    if all(min(abs(a - c), M4 - abs(a - c))**2 + min(abs(b - d), M4 - abs(b - d))**2 > 30**2 for c, d in chosen):
        chosen.append((a, b))
    if len(chosen) == 8: break
def I_at(U, q):
    ph = np.exp(-2j*np.pi*(U[:, 0]*q[0] + U[:, 1]*q[1]))
    return float(np.abs(ph.sum())**2 / len(U)**2)
pp_rows = []
for a, b in chosen:
    q = np.array([a/M4, b/M4]); step = 1.0/M4
    for _ in range(8):  # local refinement on a 7x7 grid, shrinking
        grid = [(q[0] + i*step/3, q[1] + j*step/3) for i in range(-3, 4) for j in range(-3, 4)]
        q = np.array(max(grid, key=lambda t: I_at(U4, t))); step /= 3
    # which number field do the peak coordinates live in? best fit (a + b sqrt5)/8 vs best rational c/d, d <= 60
    fits = []
    for t in q:
        f5 = min(((a_ + b_*5**0.5)/8 for a_ in range(-16, 17) for b_ in range(-8, 9)), key=lambda v: abs(v - t))
        fr = min((c_/d_ for d_ in range(1, 61) for c_ in range(0, d_ + 1)), key=lambda v: abs(v - t))
        fits.append(dict(err_sqrt5_over_8=abs(f5 - t), err_best_rational_d60=abs(fr - t)))
    pp_rows.append(dict(q=[float(q[0]), float(q[1])], L3=I_at(U3, q), L4=I_at(U4, q), fits=fits))
out['pure_point_test'] = pp_rows
# fraction of total diffracted intensity in the 200 strongest peaks (Bragg weight)
tot = I4.sum(); out['share_in_top200_L4'] = float(np.sort(I4.ravel())[::-1][:200].sum() / tot)

# ---------- (2) Fibonacci chain ----------
phi = (1 + 5**0.5) / 2
n = np.arange(1, 200001)
fib_pts = n + (1/phi) * np.floor(n / phi)          # long/short gaps 1+1/phi and 1 (a cut-and-project set)
kk = np.linspace(0.02, 2.0, 6000)
def amp(pts, k): return np.abs(np.exp(-2j*np.pi*np.outer(k, pts[:20000])).sum(1)) / 20000
fib_I = amp(fib_pts, kk)**2
# ---------- (3) square-free integers ----------
N = 10**6
mu = np.ones(N + 1, dtype=np.int64); isq = np.ones(N + 1, bool)
for p in sp.primerange(2, int(N**0.5) + 1): isq[p*p::p*p] = False
sqf = np.nonzero(isq[1:])[0] + 1
def sqf_amp(a, q):
    return float(np.exp(2j*np.pi*sqf*a/q).sum().real / N)
def sqf_formula(q):
    f = sp.factorint(q)
    if any(e >= 3 for e in f.values()): return 0.0
    val = 6/np.pi**2
    for p in f: val *= -1.0/(p*p - 1)
    return val
checks = []
for (a, q) in [(1, 2), (1, 3), (1, 4), (3, 4), (1, 6), (1, 9), (2, 9), (1, 12), (5, 12), (1, 18), (1, 25), (7, 36), (1, 8), (1, 27), (1, 5), (1, 7)]:
    checks.append(dict(a=a, q=q, measured=sqf_amp(a, q), formula=sqf_formula(q)))
out['squarefree_amplitudes'] = checks
out['squarefree_density'] = dict(measured=len(sqf)/N, formula=6/np.pi**2)
ks = np.linspace(0.0, 0.5, 6001)[1:]
sqf_I = np.abs(np.exp(-2j*np.pi*np.outer(ks, sqf[:20000])).sum(1) / sqf[20000-1])**2   # density-normalised
# ---------- (4) zeta zeros ----------
g = np.load('../primes/zeta_zeros_1000.npy')
xs = np.linspace(0.3, 3.6, 8000)
S = np.cos(np.outer(xs, g)).sum(1)
# local minima, deepest first, matched to log n for prime powers n
mins = [i for i in range(1, len(xs) - 1) if S[i] < S[i-1] and S[i] < S[i+1] and S[i] < -60]
pp = [m for m in range(2, 40) if len(sp.factorint(m)) == 1]
match = []
for i in mins:
    m = min(pp, key=lambda m: abs(np.log(m) - xs[i]))
    match.append(dict(x=float(xs[i]), depth=float(S[i]), nearest_prime_power=int(m), log=float(np.log(m)),
                      err=float(abs(np.log(m) - xs[i]))))
out['zeros_dips'] = match
# log 2, log 3, log 5, log 7 are independent over Q (unique factorisation): a small integer search finds no relation
rel = [(a, b, c, d) for a in range(-6, 7) for b in range(-6, 7) for c in range(-6, 7) for d in range(-6, 7)
       if (a, b, c, d) != (0, 0, 0, 0) and abs(a*np.log(2) + b*np.log(3) + c*np.log(5) + d*np.log(7)) < 1e-9]
out['log_relations_found'] = rel
json.dump(out, open('spectra.json', 'w'), indent=1)
np.savez_compressed('spectra_curves.npz', kk=kk, fib_I=fib_I, ks=ks, sqf_I=sqf_I, xs=xs, S=S)
print('pure point test (refined Bragg peaks, I/N^2 for L3 and L4 patches):'); [print(r) for r in out['pure_point_test']]
print('share in top 200 peaks', out['share_in_top200_L4'])
print('squarefree density', out['squarefree_density'])
for c in checks: print(c)
print('zero dips:'); [print(m) for m in match]
print('log relations', rel)
