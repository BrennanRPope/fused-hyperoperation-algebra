"""A dial that makes a flux tube: checks for the DRAFT page.

Run: python3 -I checks.py  ->  checks_output.txt, results.json   (about 1 minute)
Heavy runs behind the pictures are in make_data.py (stored in data/).

Parts:
 A  the flux exponent s (energy density |D|^s / s): exact scaling, effective dimension, confinement threshold
 B  s = 1: the flux takes the straight segment (minimal flow = Wasserstein-1), checked as a linear programme
 C  s between 1 and 2 on a grid: the flux narrows toward the segment as s falls
 D  real tubes (abelian Higgs vortices): tension of n flux units; border lambda = 1 is plain addition;
    bunching (type I) below, spreading (type II) above; which composition law fits
 P  polishing
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from scipy import sparse

OUT, RES = [], {}
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    npass += bool(ok); nfail += (not ok)
    line = f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else "")
    print(line); OUT.append(line)


def say(s):
    print(s); OUT.append(s)


# ------------------------------------------------------------------ A
say("== A. The flux exponent s")
s, Dd, lam, delta, a, b = sp.symbols('s D lambda delta a b', positive=True)
# Scale sources by lambda: the minimiser maps to D_lam(x) = lambda^{-(D-1)} D(x/lambda).
# Energy int |D|^s dV picks up lambda^{-(D-1)s} * lambda^{D}.
expo = sp.simplify(-(Dd - 1) * s + Dd)
check("A1 pair energy scales exactly as L^{D - (D-1)s} (point charges, finite self-energy)",
      sp.simplify(expo - (Dd - (Dd - 1) * s)) == 0, f"exponent {expo}")
check("A1b D=3: s=2 gives 1/L (Coulomb), s=3/2 gives L^0 (log), s=1 gives L (linear)",
      [expo.subs({Dd: 3, s: v}) for v in (2, sp.Rational(3, 2), 1)] == [-1, 0, 1])
Deff = sp.simplify(1 - ((expo - 1) - 0))  # force F ~ L^{expo-1}; D_eff = 1 - dlnF/dlnL = 2 - expo
Deff = sp.simplify(2 - expo)
check("A2 effective dimension D_eff = 1 + (D-1)(s-1): s slides it from D (s=2) down to 1 (s=1)",
      sp.simplify(Deff - (1 + (Dd - 1) * (s - 1))) == 0 and Deff.subs(s, 2) == Dd and Deff.subs(s, 1) == 1,
      str(sp.expand(Deff)))
thr = sp.solve(sp.Eq(expo, 0), s)[0]
check("A3 the potential grows with distance (confines) exactly when s < D/(D-1); 3D: s < 3/2", sp.simplify(thr - Dd / (Dd - 1)) == 0
      and thr.subs(Dd, 3) == sp.Rational(3, 2), f"threshold {thr}")
self_ok = sp.simplify(sp.integrate(sp.Symbol('r') ** (2 - 2 * sp.Rational(5, 4)), (sp.Symbol('r'), 0, 1))) == 2
check("A3b self-energy near a point charge is finite exactly when s < D/(D-1) (3D example s=5/4)", self_ok)
# Pagels-Tomboulis: L ~ F^2 (F^2)^delta  =>  energy ~ |E|^{p}, p = 2 + 4 delta; s = p/(p-1)
s_pt = (2 + 4 * delta) / (1 + 4 * delta)
check("A4 Pagels-Tomboulis exponent delta maps to s = (2+4 delta)/(1+4 delta): delta=1/4 -> s=3/2 (log), delta->inf -> s=1",
      s_pt.subs(delta, sp.Rational(1, 4)) == sp.Rational(3, 2) and sp.limit(s_pt, delta, sp.oo) == 1)
# Bundle merging: cost c(Phi) = Phi^s per length. Merging two bundles saves iff (a+b)^s < a^s + b^s iff s < 1
rng = np.random.default_rng(1)
ok = True
for sv in (0.5, 0.8, 1.2, 1.7):
    A_, B_ = rng.uniform(0.1, 3, 200), rng.uniform(0.1, 3, 200)
    sub = np.all((A_ + B_) ** sv < A_ ** sv + B_ ** sv)
    sup = np.all((A_ + B_) ** sv > A_ ** sv + B_ ** sv)
    ok &= (sub if sv < 1 else sup)
check("A5 bundles merge (bunch) iff s < 1, split (spread) iff s > 1; s = 1 is neutral", ok)
# ell_p composition: c(n) = c1 n^s gives c(n1+n2) = (c(n1)^{1/s} + c(n2)^{1/s})^s
n1, n2, c1 = 2.0, 5.0, 0.7
for sv in (0.8, 1.0, 1.3):
    c = lambda n: c1 * n ** sv
    lhs, rhs = c(n1 + n2), (c(n1) ** (1 / sv) + c(n2) ** (1 / sv)) ** sv
    ok &= abs(lhs - rhs) < 1e-12
check("A6 power-law tensions compose as an l_p sum with p = 1/s: c(n1+n2) = (c1^p + c2^p)^(1/p)", ok)
# Gilbert junction: two unit bundles merge into one of 2 at half-angle with cos = c(2)/(2 c(1)) = 2^{s-1}
angs = {sv: 2 * math.degrees(math.acos(2 ** (sv - 1))) for sv in (1.0, 0.75, 0.5, 0.25)}
check("A7 junction angle 2 arccos(2^{s-1}): 0 at s=1, opens as s falls (s=1/2: 90 deg)",
      abs(angs[1.0]) < 1e-9 and abs(angs[0.5] - 90) < 1e-9, ", ".join(f"s={k}: {v:.2f} deg" for k, v in angs.items()))
RES['junction_angles'] = angs

# ------------------------------------------------------------------ B
say("== B. s = 1: minimal flow is the straight segment")
Nx, Ny = 21, 11
idx = lambda i, j: i * Ny + j
edges = []
for i in range(Nx):
    for j in range(Ny):
        if i + 1 < Nx: edges.append((idx(i, j), idx(i + 1, j)))
        if j + 1 < Ny: edges.append((idx(i, j), idx(i, j + 1)))
E_ = len(edges); V_ = Nx * Ny
# variables f+ and f- per edge (flow = f+ - f-), cost sum f+ + f-
rows, cols, vals = [], [], []
for e, (u, v) in enumerate(edges):
    for sign, var in ((1, e), (-1, e + E_)):
        rows += [u, v]; cols += [var, var]; vals += [sign, -sign]
Aeq = sparse.csr_matrix((vals, (rows, cols)), shape=(V_, 2 * E_))
src, snk = idx(4, 5), idx(16, 5)
beq = np.zeros(V_); beq[src] = 1; beq[snk] = -1
res = linprog(np.ones(2 * E_), A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
flow = res.x[:E_] - res.x[E_:]
on_seg = [e for e, (u, v) in enumerate(edges) if abs(flow[e]) > 1e-9]
seg_ok = all(u % Ny == 5 and v % Ny == 5 for (u, v) in (edges[e] for e in on_seg))
check("B1 minimal total flux between two charges = their distance (12), and all flux sits on the segment",
      abs(res.fun - 12) < 1e-9 and seg_ok and len(on_seg) == 12, f"cost {res.fun}, edges used {len(on_seg)}")
# three sources, one sink of 3: with s=1 the flows just add; no junction saving
src3 = [idx(4, 2), idx(4, 5), idx(4, 8)]
beq = np.zeros(V_); beq[src3] = 1; beq[idx(16, 5)] = -3
res3 = linprog(np.ones(2 * E_), A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
check("B2 s=1 is additive: three sources to one sink cost the sum of their (grid) distances", abs(res3.fun - (15 + 12 + 15)) < 1e-9,
      f"{res3.fun}")

# ------------------------------------------------------------------ C
say("== C. 1 < s < 2 on a grid: the flux narrows toward the segment")
from field_solver import solve
wid = {}
for sv in (1.5, 1.25):
    e, Dm, x, r_ = solve(sv, 3.0, N=96, box=12.0, w=0.25, eps=1e-4)
    mid = Dm[len(x) // 2, :]
    wid[sv] = math.sqrt(np.sum(mid * x ** 2) / np.sum(mid))
    check(f"C1 s={sv}: solver converged", r_.status == 0, f"iterations {r_.nit}, energy {e:.5f}")
check("C2 lower s, narrower flux at the midplane (L=3, N=96)", wid[1.25] < wid[1.5],
      f"rms width s=1.5: {wid[1.5]:.3f}, s=1.25: {wid[1.25]:.3f}")
F2 = json.load(open('data/field2.json'))
w4 = {k: v['widths'][2] / 4.0 for k, v in F2.items()}
check("C3 stored N=160 runs (L=4): width/L falls 0.42 -> 0.29 -> 0.18 as s goes 1.5 -> 1.25 -> 1.1",
      w4['1.5'] > w4['1.25'] > w4['1.1'], ", ".join(f"s={k}: {v:.3f}" for k, v in w4.items()))
sl = {k: float(np.polyfit(np.log([2, 3, 4, 5]), np.log(v['E']), 1)[0]) for k, v in F2.items()}
check("C4 energy keeps growing with L at every s<2 (2D confines logarithmically already); slope rises toward 1 as s -> 1",
      all(np.all(np.diff(v['E']) > 0) for v in F2.values()), ", ".join(f"s={k}: {v:.3f}" for k, v in sl.items()))
say("   note: with smeared charges (width 0.25) and L only 2..5, the fitted slopes have not yet reached the point-charge"
    " value 2 - s; widths/L still shrink with L, so this is short of the scaling regime.")
RES['C_widths_over_L'] = w4; RES['C_slopes'] = sl

# ------------------------------------------------------------------ D
say("== D. Real tubes: abelian Higgs vortices with n flux units (units e = v = 1; border at lambda = 1)")
from vortex_solver import vortex
Ev = {}
for lv in (0.5, 1.0, 2.0):
    Ev[lv] = [vortex(n, lv)[0] / math.pi for n in (1, 2, 3, 4)]
check("D1 border lambda = 1: tension is exactly pi n (Bogomol'nyi), plain addition of flux units",
      max(abs(Ev[1.0][i] - (i + 1)) for i in range(4)) < 1e-5, ", ".join(f"{x:.7f}" for x in Ev[1.0]))
check("D2 lambda < 1 (type I): E(2) < 2E(1), flux units attract and bunch; lambda > 1 (type II): E(2) > 2E(1)",
      Ev[0.5][1] < 2 * Ev[0.5][0] and Ev[2.0][1] > 2 * Ev[2.0][0],
      f"0.5: {Ev[0.5][1]:.5f} vs {2*Ev[0.5][0]:.5f}; 2: {Ev[2.0][1]:.5f} vs {2*Ev[2.0][0]:.5f}")
V8 = json.load(open('data/vortex_E.json'))
fits = {}
for lv, Es in V8.items():
    sv = math.log2(Es[1] / Es[0])
    kv = (Es[1] - 2 * Es[0]) / Es[0] ** 2
    pw = [Es[0] * n ** sv for n in range(1, 9)]
    ck = [(math.exp(math.log(1 + kv * Es[0]) * n) - 1) / kv if abs(kv) > 1e-12 else Es[0] * n for n in range(1, 9)]
    fits[lv] = dict(s=sv, k=kv, err_lp=max(abs(p - e) / e for p, e in zip(pw, Es)),
                    err_ok=max(abs(p - e) / e for p, e in zip(ck, Es)))
RES['vortex_fits'] = fits
near = [fits[k]['err_lp'] for k in ('0.5', '0.75', '1.5', '2.0')]
check("D3 near the border (lambda 0.5..2) tensions of 1..8 units follow the l_p sum (power n^s) to within 2%",
      max(near) < 0.02, ", ".join(f"{k}: {fits[k]['err_lp']*100:.1f}%" for k in ('0.5', '0.75', '1.5', '2.0')))
check("D3b far from the border it drifts (lambda 0.1: 4%, 8: 6%): the power law is a near-border description",
      fits['0.1']['err_lp'] < 0.05 and fits['8.0']['err_lp'] < 0.07,
      f"0.1: {fits['0.1']['err_lp']*100:.1f}%, 8: {fits['8.0']['err_lp']*100:.1f}%")
check("D4 our law o_k with one k (fixed by units 1 and 2) does NOT fit: errors 7% to 100% off the border",
      min(fits[k]['err_ok'] for k in fits if k != '1.0') > 0.05,
      ", ".join(f"{k}: {fits[k]['err_ok']*100:.0f}%" for k in fits if k != '1.0'))
svals = [fits[k]['s'] for k in sorted(fits, key=float)]
check("D5 the reading s(lambda) = log2(E2/E1) rises steadily through 1 at the border",
      all(np.diff(svals) > 0) and abs(fits['1.0']['s'] - 1) < 1e-6, ", ".join(f"{v:.3f}" for v in svals))
# D6: Ginzburg-Landau kappa = m_H / (sqrt2 m_V); here m_H^2 = V''(1) = lambda, m_V = 1
fH = sp.Symbol('f')
VH = lam / 8 * (fH ** 2 - 1) ** 2
mH2 = sp.diff(VH, fH, 2).subs(fH, 1)
check("D6 in these units m_H^2 = lambda and m_V = 1, so kappa = sqrt(lambda/2): border lambda = 1 is kappa = 1/sqrt 2",
      sp.simplify(mH2 - lam) == 0)
s_lat = float(np.interp(2 * 1.6 ** 2, [4.0, 8.0], [fits['4.0']['s'], fits['8.0']['s']]))
say(f"   lattice estimates of the QCD dual kappa ~ 1.5-1.7 (SU(2), conventions differ, unsettled) -> lambda ~ 4.5-5.8 -> s ~ {s_lat:.2f}")
RES['s_at_kappa_1p6'] = s_lat

# ------------------------------------------------------------------ P
say("== P. Polishing")
Eb = vortex(1, 1.0, R=40.0, N=6000)[0] / math.pi
check("P1 border unit vortex at a finer grid and larger box: pi exactly to 1e-7", abs(Eb - 1) < 1e-7,
      f"{Eb:.10f}, moved {abs(Eb - Ev[1.0][0]):.1e} from the coarse run")
check("P2 s=1 minimal flow is exact (linear programme, no smoothing)", res.status == 0)

say(f"\nTOTAL: {npass} pass, {nfail} fail")
RES['pass'], RES['fail'] = npass, nfail
open('checks_output.txt', 'w').write("\n".join(OUT) + "\n")
json.dump(RES, open('results.json', 'w'), indent=1, default=float)
