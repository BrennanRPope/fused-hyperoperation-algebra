"""Checks for the proof of the q -> 1 slope rule for Bradley's q-MZVs.  (DRAFT)
Run: python3 slope_proof_checks.py   (sympy, mpmath, numpy)

Asked by Brennan, Nov 27 05:38Z: "Yes" to "Should I try to prove the q -> 1 slope rule at every depth?"

Theorem (proved in README section 7). Let s = (s1, ..., sr) be positive integers with s1 >= 2, not starting (2, 1).
As q -> 1-,  zeta_q(s) = zeta(s) + (q - 1) D(s) + o(1 - q), with
  D(s) = (w/2) zeta(s) + sum_i (s_i/2 - 1) zeta(s with s_i lowered by 1) + [s1 = 2] (1/2) zeta(s2, ..., sr),
where a lowered 1 becomes a letter 0 (n^0 = 1) and zeta() = 1.
"""
import json, numpy as np, mpmath as mp, sympy as sp
mp.mp.dps = 30
res = []
def ok(name, c, note=""):
    res.append({"id": name, "pass": bool(c), "note": note}); print(("PASS " if c else "FAIL ") + name)

# P1 the letter in scaling form: with q = e^-t and u = (1 - q)/t,  f_s(n) = u^s n^-s h_s(nt), where
#    h_s(x) = ((x/2)/sinh(x/2))^s e^((1 - s/2)x)
x, t, n, s = sp.symbols("x t n s", positive=True)
q = sp.exp(-t)
P1 = True
for S in range(1, 6):
    for N in range(1, 6):
        f = q**((S-1)*N) / sum(q**j for j in range(N))**S
        hS = ((N*t/2)/sp.sinh(N*t/2))**S * sp.exp((1 - sp.Rational(S, 2))*N*t)
        rhs = ((1 - q)/t)**S * sp.Integer(N)**(-S) * hS
        P1 &= abs(float((f - rhs).subs(t, sp.Rational(3, 10)))) < 1e-14 and abs(float((f - rhs).subs(t, sp.Rational(1, 7)))) < 1e-14
ok("P1 Bradley's letter is u^s n^-s h_s(nt) with h_s(x) = ((x/2)/sinh(x/2))^s e^((1 - s/2)x)  (s, n = 1..5)", P1)

# P2 h_s(x) = 1 + (1 - s/2) x + O(x^2), and h_2(x) = 1 - x^2/12 + ...; u = 1 - t/2 + O(t^2)
P2 = True
for S in range(1, 7):
    hS = ((x/2)/sp.sinh(x/2))**S * sp.exp((1 - sp.Rational(S, 2))*x)
    ser = sp.series(hS, x, 0, 2).removeO()
    P2 &= sp.simplify(ser - (1 + (1 - sp.Rational(S, 2))*x)) == 0
h2 = ((x/2)/sp.sinh(x/2))**2 * sp.exp(0*x)
P2 &= sp.simplify(sp.series(h2, x, 0, 3).removeO() - (1 - x**2/12)) == 0
P2 &= sp.simplify(sp.series((1 - sp.exp(-t))/t, t, 0, 2).removeO() - (1 - t/2)) == 0
ok("P2 h_s(x) = 1 + (1 - s/2)x + O(x^2) (s = 1..6); h_2 = 1 - x^2/12 + ...; u = 1 - t/2 + O(t^2)", P2)

# P3 the bounds used: |h_s(x) - 1 - (1 - s/2)x| <= C min(x, x^2) on a fine grid, h_s bounded for s >= 2, h_1(x) <= 1 + x
X = np.concatenate([np.linspace(1e-4, 1, 2000), np.linspace(1, 200, 20000)])
def h(S, X): return ((X/2)/np.sinh(X/2))**S * np.exp((1 - S/2)*X) if S > 1 else X/(1 - np.exp(-X))
P3 = True; consts = {}
for S in range(1, 7):
    r = np.abs(h(S, X) - 1 - (1 - S/2)*X) / np.minimum(X, X**2)
    consts[S] = float(r.max()); P3 &= np.isfinite(r).all() and r.max() < 2
    if S >= 2: P3 &= h(S, X).max() <= 1 + 1e-12
P3 &= (h(1, X) <= 1 + X + 1e-12).all()
ok("P3 |h_s(x) - 1 - (1 - s/2)x| <= C min(x, x^2) with C < 2 (s = 1..6); 0 < h_s <= 1 for s >= 2; h_1(x) <= 1 + x", P3,
   "largest ratio per s: " + ", ".join(f"{k}: {v:.3f}" for k, v in consts.items()))

# P4 the residue as an area: g(x) = (h_2(x) - 1)/x^2 = 1/(4 sinh^2(x/2)) - 1/x^2 is increasing from -1/12 to 0,
#    and its integral over (0, oo) is exactly -1/2 (antiderivative 1/x - coth(x/2)/2)
g = lambda X: 1/(4*np.sinh(X/2)**2) - 1/X**2
G = np.linspace(1e-2, 60, 200000)
P4a = np.all(np.diff(g(G)) > 0) and abs(g(np.array([1e-2]))[0] + 1/12) < 1e-5
I = mp.quad(lambda X: 1/(4*mp.sinh(X/2)**2) - 1/X**2, [0, 1, 10, mp.inf])
anti = lambda X: 1/X - sp.coth(X/2)/2
P4b = abs(float(I) + 0.5) < 1e-20 and (sp.limit(anti(x), x, sp.oo) == -sp.Rational(1, 2) and sp.series(anti(x), x, 0, 2).removeO() == -x/12
       and sp.simplify(sp.diff(anti(x), x) - (1/(4*sp.sinh(x/2)**2) - 1/x**2)) == 0)
ok("P4a g(x) = (h_2(x) - 1)/x^2 is increasing, from -1/12 at 0 to 0 at infinity (so its Riemann sums converge)", P4a)
ok("P4b the integral of g over (0, infinity) is exactly -1/2: the 1/2 in the slope is this area", P4b, f"quad = {mp.nstr(I, 25)}")
# Riemann sums t * sum g(n t) -> -1/2
rs = [t_ * np.sum(g(t_*np.arange(1, int(80/t_)))) - 1/(t_*int(80/t_)) for t_ in (1e-2, 1e-3, 1e-4)]   # plus the tail of -1/x^2 beyond x = 80
ok("P4c Riemann sums t * sum_n g(nt) (tail beyond 80 added) -> -1/2: " + ", ".join(f"{v:.6f}" for v in rs), abs(rs[-1] + 0.5) < 1e-3)

# P5 the convergence condition: for s1 >= 2, "s1 + ... + s_i >= i + 2 for all i >= 2" holds exactly when s does not start (2, 1)
import itertools
P5 = True
for r in range(1, 5):
    for w in itertools.product(range(1, 5), repeat=r):
        if w[0] < 2: continue
        H = all(sum(w[:i]) >= i + 2 for i in range(2, r + 1))
        P5 &= H == (not (r >= 2 and w[0] == 2 and w[1] == 1))
ok("P5 the proof's hypothesis is the same as 'does not start with 2, 1' (all words of depth <= 4, letters <= 4)", P5)

# P6 the theorem against the q-series, depth 2 and 3 (values from known evaluations; slopes measured at q -> 1-)
z = {k: float(mp.zeta(k)) for k in range(2, 8)}
Zv = {(2,): z[2], (3,): z[3], (4,): z[4], (2, 1): z[3], (2, 2): 0.75*z[4], (3, 1): z[4]/4, (2, 1, 1): z[4],
      (4, 1): 2*z[5] - z[2]*z[3], (3, 2): 3*z[2]*z[3] - 11*z[5]/2, (2, 3): 4.5*z[5] - 2*z[2]*z[3],
      (3, 1, 1): 2*z[5] - z[2]*z[3], (2, 2, 1): 3*z[2]*z[3] - 11*z[5]/2, (2, 2, 2): float(mp.pi**6/5040), (): 1.0}
Zv[(3, 0)] = z[2] - z[3]; Zv[(4, 0)] = z[3] - z[4]; Zv[(2, 0)] = None
Zv[(3, 1, 0)] = Zv[(3, 0)] - Zv[(3, 1)]
Zv[(3, 0, 1)] = Zv[(2, 1)] - Zv[(3, 0)] - Zv[(3, 1)]
Zv[(2, 2, 0)] = Zv[(2, 1)] - Zv[(2, 2)]
def D(w):
    tot = sum(w)/2 * Zv[w]
    for i, si in enumerate(w):
        if si == 2: continue                                 # coefficient s/2 - 1 = 0 (a front 2 gets its residue below)
        low = w[:i] + (si - 1,) + w[i+1:]
        if si - 1 == 0 and i == len(w) - 1 and len(w) == 1: low = ()
        tot += (si/2 - 1) * Zv[low]
    if w[0] == 2: tot += 0.5 * Zv[w[1:]]
    return tot
def letters(S, q, Nmax):
    nn = np.arange(1, Nmax + 1, dtype=np.float64); qn = np.exp(nn*np.log(q))
    return (1 - q)**S * qn**(S - 1) / (1 - qn)**S
def brad(w, q):
    Nmax = int(90/(1 - q)) + 10; acc = np.ones(Nmax)
    for S in reversed(w):
        L = letters(S, q, Nmax) * acc; acc = np.concatenate([[0.0], np.cumsum(L)[:-1]])
    return float(np.sum(L))
def slope(w):
    hs = [4e-3, 2e-3, 1e-3, 5e-4]; sl = [(Zv[w] - brad(w, 1 - hh))/hh for hh in hs]
    return 2*sl[-1] - sl[-2]
rows = []; P6 = True
for w in [(2,), (3,), (4,), (2, 2), (3, 1), (2, 3), (3, 2), (4, 1), (2, 2, 1), (3, 1, 1), (2, 2, 2)]:
    m, p = slope(w), D(w); rows.append((w, m, p)); P6 &= abs(m - p) < 2e-5
ok("P6 slope rule vs measured slope on 11 words of depth 1 to 3 (largest gap %.1e)" % max(abs(m - p) for _, m, p in rows), P6,
   "; ".join(f"{w}: measured {m:.6f}, rule {p:.6f}" for w, m, p in rows))

# P7 control: dropping the front-2 residue breaks the depth-3 words that start with 2
P7 = all(abs(m - (p - 0.5*Zv[w[1:]])) > 0.3 for w, m, p in rows if w[0] == 2)
ok("P7 negative control: without the front-2 term every word starting with 2 misses by at least 0.3", P7)

json.dump({"checks": res, "slopes": [{"word": list(w), "measured": m, "rule": p} for w, m, p in rows]},
          open("slope_proof_checks.json", "w"), indent=1)
print(f"{sum(r['pass'] for r in res)}/{len(res)} pass")
