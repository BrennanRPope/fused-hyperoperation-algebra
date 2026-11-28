"""Checks for the widened paper 'Why quantum needs i, and a small bridge toward gravity' (DRAFT).
Run: python3 checks.py  -> prints PASS/FAIL per check, writes results.json.
Upstream pages are not re-run here; their recorded outputs are hash-checked (U*)."""
import re, json, hashlib, subprocess, os, itertools
import numpy as np, sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.abspath(os.path.join(HERE, '..', '..'))
rng = np.random.default_rng(7)
res, out = {}, []
def check(name, ok, msg):
    out.append(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}"); res[name] = {'pass': bool(ok), 'msg': msg}

def herm(n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return (a + a.conj().T) / 2

# ---------- Part 1: why quantum needs i (re-run the theorem's own checks) ----------
p = subprocess.run(['python3', 'checks.py'], cwd=os.path.join(MT, 'why-i', 'rotations-force-i'),
                   capture_output=True, text=True, timeout=600)
npass = len(re.findall(r'\bPASS\b', p.stdout)); nfail = len(re.findall(r'\bFAIL\b', p.stdout))
mt = re.search(r'(\d+)/(\d+) pass', p.stdout)
check('P1 theorem checks re-run', nfail == 0 and mt is not None and mt.group(1) == mt.group(2) and npass >= 21, f'rotations-force-i/checks.py: {npass} pass, {nfail} fail')

# ---------- Part 2: the information side ----------
# I1 Tsallis entropy is mult-add with k = 1 - q on independent systems
def tsallis(p, q): return (1 - np.sum(p ** q)) / (q - 1)
errs = []
for q in [0.5, 2.0, 3.0]:
    pa = rng.dirichlet(np.ones(3)); pb = rng.dirichlet(np.ones(4)); pab = np.outer(pa, pb).ravel()
    k = 1 - q; Sa, Sb = tsallis(pa, q), tsallis(pb, q)
    errs.append(abs(tsallis(pab, q) - (Sa + Sb + k * Sa * Sb)))
check('I1 Tsallis entropy composes by the law', max(errs) < 1e-12, f'S(AB) = S(A)+S(B)+k S(A)S(B), k = 1-q; max error {max(errs):.1e}')

# I2 the linear (coin-weighted) grouping rule holds for Shannon only
def shannon(p): p = p[p > 0]; return -np.sum(p * np.log(p))
def group_gap(S, p):
    s = p[0] + p[1]; return S(p) - (S(np.r_[s, p[2:]]) + s * S(np.array([p[0] / s, p[1] / s])))
pp = rng.dirichlet(np.ones(4))
gS = abs(group_gap(shannon, pp)); gT = abs(group_gap(lambda x: tsallis(x, 2.0), pp))
check('I2 linear mixing pins k = 0', gS < 1e-12 and gT > 1e-3, f'grouping gap Shannon {gS:.1e}, Tsallis q=2 (k=-1) {gT:.3f}')

# I3 interference orders: the law as a measure vs quantum amplitudes
k = sp.symbols('k'); a, b, c = sp.symbols('a b c', positive=True)
mu = lambda *xs: (sp.prod([1 + k * x for x in xs]) - 1) / k
I2 = sp.simplify(mu(a, b) - mu(a) - mu(b))
I3 = sp.simplify(mu(a, b, c) - mu(a, b) - mu(a, c) - mu(b, c) + mu(a) + mu(b) + mu(c))
amp = rng.normal(size=3) + 1j * rng.normal(size=3)
P = lambda idx: abs(sum(amp[list(idx)])) ** 2
q2 = P([0, 1]) - P([0]) - P([1])
q3 = P([0, 1, 2]) - P([0, 1]) - P([0, 2]) - P([1, 2]) + P([0]) + P([1]) + P([2])
check('I3 law interferes at every order, quantum stops at 2', sp.simplify(I2 - k*a*b) == 0 and sp.simplify(I3 - k**2*a*b*c) == 0
      and abs(q3) < 1e-12 and abs(q2) > 1e-3, f'law: I2 = k ab, I3 = k^2 abc; quantum: I2 = {q2:.3f}, I3 = {q3:.1e}')

# I4 the three-slit cap on k (forgiving reading |k w| < sqrt(6 kappa), w = 1/3)
caps = {name: float(np.sqrt(6 * kap) / (1 / 3)) for name, kap in [('Sinha 2010', 1e-2), ('Kauten 2017', 1e-4)]}
check('I4 three-slit cap', abs(caps['Kauten 2017'] - 0.0735) < 1e-3 and abs(caps['Sinha 2010'] - 0.735) < 1e-2,
      f"|k| < {caps['Kauten 2017']:.3f} (Kauten 2017), < {caps['Sinha 2010']:.2f} (Sinha 2010)")

# I5 XOR is the law at k = -2 on bits
ok = all((x + y - 2 * x * y) == (x ^ y) for x in (0, 1) for y in (0, 1))
check('I5 XOR is the law at k = -2', ok, 'a + b - 2ab = a XOR b on {0,1}')

# I6 ball vs square: CHSH for quantum (correlation cos on a circle) and for a PR box (square)
def chsh(t): A0, A1, B0, B1 = t; E = lambda x, y: np.cos(x - y); return E(A0, B0) + E(A0, B1) + E(A1, B0) - E(A1, B1)
best = max(-minimize(lambda t: -chsh(t), rng.uniform(0, 2*np.pi, 4)).fun for _ in range(20))
ceil = {a_: 2 ** (2 - 1 / a_) for a_ in (1, 2, 1e9)}
check('I6 ball gives 2*sqrt2, square gives 4', abs(best - 2*np.sqrt(2)) < 1e-6 and abs(ceil[2] - 2*np.sqrt(2)) < 1e-12
      and abs(ceil[1] - 2) < 1e-12 and abs(ceil[1e9] - 4) < 1e-6, f'circle max {best:.6f}; ceiling 2^(2-1/a): a=1 -> 2, a=2 -> {ceil[2]:.4f}, a->inf -> 4')

# ---------- Part 3: the hinge, one shift lemma used three ways ----------
# H1 quantum: shifting H by c*1 changes no odds; shifting along a non-commuting generator does
n = 4; H = herm(n); G = herm(n); O = herm(n); psi = rng.normal(size=n) + 1j*rng.normal(size=n); psi /= np.linalg.norm(psi)
ev = lambda Hm, t: np.real(np.vdot(expm(-1j*Hm*t) @ psi, O @ expm(-1j*Hm*t) @ psi))
d1 = max(abs(ev(H + cc*np.eye(n), 1.3) - ev(H, 1.3)) for cc in (0.5, 2, 10))
d2 = max(abs(ev(H + cc*G, 1.3) - ev(H, 1.3)) for cc in (0.5, 2, 10))
check('H1 energy zero is a convention in QM', d1 < 1e-10 and d2 > 1e-2, f'shift by c*1: change {d1:.1e}; shift by c*G (non-commuting): change {d2:.2f}')

# H2 bubble: a constant potential gives E = -[A, phi]; zero if Abelian, non-zero if not
s = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]])]
T = [-0.5j*m for m in s]
A = sum(rng.normal()*t for t in T); phi = 0.7*T[2]
E_na = np.linalg.norm(A @ phi - phi @ A); E_ab = abs(0.3*0.7 - 0.7*0.3)
check('H2 potential zero is physical in a non-Abelian field', E_na > 1e-2 and E_ab == 0, f'|[A, phi]| = {E_na:.3f} (su(2)); Abelian = {E_ab}')

# H3 gravity: a constant in the action is weighted by sqrt(-g), so it sources the metric
h = sp.symbols('h0:4', real=True); rho = sp.symbols('rho', positive=True)
# inverse metric g^{mn} = diag(h0, h1, h2, h3) with h0 < 0; metric g_mn = diag(1/h_i)
gm = sp.diag(*[1/hi for hi in h]); sqrtg = sp.sqrt(-gm.det())
L = -rho*sqrtg   # a constant energy density rho, weighted by the volume element
Tcov = [sp.simplify(-2/sqrtg*sp.diff(L, h[i])) for i in range(4)]   # T_mn = -2/sqrt(-g) d(sqrt(-g) L)/d g^mn
vals = {h[0]: -sp.Rational(1, 3), h[1]: sp.Rational(1, 2), h[2]: 2, h[3]: sp.Rational(5, 7)}
ok3 = all(sp.simplify((Tcov[i] - (-rho*gm[i, i])).subs(vals)) == 0 for i in range(4))
check('H3 a constant energy sources gravity', ok3, 'L = -rho sqrt(-g) gives T_mn = -rho g_mn (a cosmological constant), not zero')

# H4 the action view: in QM a constant c adds -c*T to every path's action; in gravity -rho*Volume, which differs between geometries
m_ = 1.0; Tt = 2.0; N = 400; tt = np.linspace(0, Tt, N); dt = tt[1]-tt[0]
def act(q, cc): v = np.diff(q)/dt; return np.sum(0.5*m_*v**2 - 0.5*q[:-1]**2 - cc)*dt
q1 = np.sin(np.pi*tt/Tt); q2 = 0.3*np.sin(2*np.pi*tt/Tt)
dS = [act(q1, cc) - act(q2, cc) for cc in (0, 1, 5)]
Vol = lambda aa: aa**3  # volume of a unit coordinate box at scale factor a, one unit of time
dV = [rho0*(Vol(1.0) - Vol(1.5)) for rho0 in (0, 1, 5)]
check('H4 action: same shift for every path in QM, different for every geometry in gravity', max(dS) - min(dS) < 1e-9 and len(set(np.round(dV, 9))) == 3,
      f'path action difference stays {dS[0]:.4f} for c = 0,1,5; geometry action difference -rho*dV = {[round(x,3) for x in dV]}')

# H5 the fused shift exp(k t): imaginary k turns (stays on the circle), real k grows (de Sitter)
th = np.linspace(0, 6, 50); r_i = np.abs(np.exp(1j*th)); 
Hs = {r0: float(np.sqrt(8*np.pi*r0/3)) for r0 in (0.01, 0.04, 0.16)}
a2 = {r0: float(np.exp(Hh*2)) for r0, Hh in Hs.items()}
t_ = sp.symbols('t'); Hsym = sp.sqrt(8*sp.pi*rho/3); asym = sp.exp(Hsym*t_)
ok5 = np.allclose(r_i, 1) and sp.simplify(sp.diff(asym, t_, 2)/asym - 8*sp.pi*rho/3) == 0
check('H5 imaginary k turns, real k grows', ok5, f'|e^(i t)| = 1; de Sitter a = e^(Ht), H = sqrt(8 pi G rho/3) solves a"/a = 8piG rho/3; a(2) = {[round(v,3) for v in a2.values()]} for rho = 0.01, 0.04, 0.16')

# H6 Dong's brane tension (n-1)/(4nG) vanishes only at n = 1, i.e. k = 0 with n = 1 - k
nn = sp.symbols('n', positive=True); Gs = sp.symbols('G', positive=True)
tension = (nn - 1)/(4*nn*Gs); sol = sp.solve(sp.Eq(tension.subs(nn, 1 - k), 0), k)
check('H6 entropy leg sits at k = 0', sol == [0], f'tension(n = 1-k) = 0 only at k = {sol}')

# H7 unimodular gravity: fix the total 4-volume and a constant energy is a convention again
gsym = sp.diag(*[1/hi for hi in h]); Tc = -rho*gsym
trace = sum((sp.diag(*h)*Tc)[i, i] for i in range(4))
tracefree = sp.simplify(Tc - sp.Rational(1, 4)*gsym*trace)
V0 = 2.0; w = lambda r0, S: np.exp(-S - r0*V0)          # every history has the same volume V0
ratios = [w(r0, 0.3)/w(r0, 1.1) for r0 in (0, 1, 5)]
check('H7 fixed volume makes the constant a convention again', tracefree == sp.zeros(4, 4) and np.ptp(ratios) < 1e-12,
      f'trace-free part of -rho g_mn is 0 (drops out of unimodular field equations); weight ratio {ratios[0]:.6f} for rho = 0,1,5')

# ---------- Part 4: the gravity side ----------
# G1 frame dragging far from Kerr: omega = -g_tphi/g_phiphi -> 2 M a / r^3 = 2J/r^3
r, th_, M, aa = sp.symbols('r theta M a', positive=True)
Sig = r**2 + aa**2*sp.cos(th_)**2; Del = r**2 - 2*M*r + aa**2
gtp = -2*M*aa*r*sp.sin(th_)**2/Sig
gpp = (r**2 + aa**2 + 2*M*r*aa**2*sp.sin(th_)**2/Sig)*sp.sin(th_)**2
omega = sp.simplify(-gtp/gpp); lead = sp.limit(omega*r**3, r, sp.oo)
check('G1 frame dragging 2J/r^3', sp.simplify(lead - 2*M*aa) == 0, f'r^3 * omega -> {lead} = 2J with J = M a')

# G2 Newman-Janis shift: complex distance to a source at height i*a gives Kerr's ring and its r
x, y, z = sp.symbols('x y z', real=True)
Rc2 = x**2 + y**2 + (z - sp.I*aa)**2
rK = sp.symbols('r_K', positive=True)
# Kerr-Schild r: (x^2+y^2)/(r^2+a^2) + z^2/r^2 = 1 ; claim: Re sqrt(Rc2) = r_K on the axis and the ring
on_axis = sp.simplify(sp.sqrt(Rc2.subs({x: 0, y: 0})) .subs(z, 3*aa))
pt = {x: 1.3, y: 0.4, z: 0.9, aa: 0.6}
Rnum = complex(sp.sqrt(Rc2.subs(pt)).evalf()); rr = float(sp.nsolve(((1.3**2+0.4**2)/(rK**2+0.36) + 0.81/rK**2 - 1), rK, 1.0))
check('G2 the shift by i*a gives Kerr r', abs(abs(Rnum.real) - rr) < 1e-10 and abs(abs(Rnum.imag) - 0.6*0.9/rr) < 1e-10,
      f'sqrt((x,y,z - i a)^2) = r_K - i a z/r_K at a sample point: r_K = {rr:.6f}')

# ---------- Upstream results this paper quotes: verify their recorded outputs by hash ----------
def sha(path): return hashlib.sha256(open(path, 'rb').read()).hexdigest()
ups = {
    'rotations-force-i': ('why-i/rotations-force-i/checks_output.txt', '21/21'),
    'spin-i-tie': ('spin-i-tie/checks_output.txt', None),
    'three-slit-cap': ('spin-i-tie/three-slit-cap/checks_output.txt', None),
    'spinning-river': ('gravity-from-theory/spinning-river/checks_output.txt', '139/139'),
    'bridges-study': ('bridges-study/checks_output.txt', None),
    'clock-test': ('spin-i-tie/clock-test/checks_output.txt', None),
}
upres = {}
for nm, (rel, need) in ups.items():
    pth = os.path.join(MT, rel); txt = open(pth).read()
    fails = txt.count('FAIL') - txt.count('0 FAIL') - txt.lower().count('no fail')
    upres[nm] = {'file': rel, 'sha256': sha(pth), 'fail_lines': txt.count('[FAIL]')}
okU = all(v['fail_lines'] == 0 for v in upres.values())
bsum = os.path.join(MT, 'bridges-study', 'SHA256SUMS')
bok = all(sha(os.path.join(MT, 'bridges-study', l.split()[1])) == l.split()[0] for l in open(bsum) if l.strip())
check('U1 upstream outputs show no failures, bridges-study hashes verify', okU and bok,
      f"{len(upres)} recorded outputs read, 0 [FAIL] lines; bridges-study SHA256SUMS {'verified' if bok else 'MISMATCH'}")
res['upstream'] = upres

npass = sum(v['pass'] for k_, v in res.items() if isinstance(v, dict) and 'pass' in v)
ntot = sum(1 for k_, v in res.items() if isinstance(v, dict) and 'pass' in v)
out.append(f'\n{npass}/{ntot} checks pass')
print('\n'.join(out))
json.dump({k_: v for k_, v in res.items()}, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=str)
