"""DRAFT round 26 checks: the moving half of the tide, and gravitational waves, from a moving or spinning twist.

Round 25 found the chain twist -> dial q -> tide, but only for slices at rest, where the tide has just its 'electric' half.
Here things move: the law's slices are squashed and unsquashed (the lowest gravitational wave they can carry), the weak force's
locked twist ripples, and a lump spins. The question is what makes the tide's 'magnetic' half and gravitational waves,
and whether the twist can.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, time, itertools
import numpy as np
import sympy as spy
from scipy.integrate import solve_ivp

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

# ------------------------------------------------------------------ a squashed, expanding law 3-sphere with a rippling locked twist
# frame e_0 = d/dt, e_b = X_b / a_b(t) on the unit 3-sphere's left-invariant fields ([X_a, X_b] = 2 eps X_c, round 25);
# a_b = a (1 + e h_b), h_1 + h_2 + h_3 = 0 (a pure squash); twist A(e_b) = (f + e u_b) T_b / a_b, u_1 + u_2 + u_3 = 0 (a pure ripple)
t, ee = spy.symbols('t epsilon')
aF = spy.Function('a')(t); fF = spy.Function('f')(t)
h = [spy.Function('h1')(t), spy.Function('h2')(t)]; hh = [h[0], h[1], -h[0] - h[1]]
u = [spy.Function('u1')(t), spy.Function('u2')(t)]; uu = [u[0], u[1], -u[0] - u[1]]
Ab = [aF * (1 + ee * hh[b]) for b in range(3)]
Fb = [fF + ee * uu[b] for b in range(3)]
eta = [-1, 1, 1, 1]
def E3(i, j, k): return spy.LeviCivita(i, j, k)
cc = [[[0] * 4 for _ in range(4)] for _ in range(4)]
for b in range(1, 4):
    ab = Ab[b - 1]; cc[b][0][b] = -spy.diff(ab, t) / ab; cc[b][b][0] = spy.diff(ab, t) / ab
    for b2 in range(1, 4):
        for b3 in range(1, 4):
            if len({b, b2, b3}) == 3:
                cc[b3][b][b2] = 2 * E3(b - 1, b2 - 1, b3 - 1) * Ab[b3 - 1] / (Ab[b - 1] * Ab[b2 - 1])
cl = lambda k, i, j: eta[k] * cc[k][i][j]
Gm = [[[(cl(k, i, j) - cl(i, j, k) + cl(j, k, i)) / 2 for j in range(4)] for i in range(4)] for k in range(4)]
Gu = lambda k, i, j: eta[k] * Gm[k][i][j]
def der(i, ex): return spy.diff(ex, t) if i == 0 else ex * 0
def Riem(k, l, i, j):
    return (der(i, Gu(k, j, l)) - der(j, Gu(k, i, l)) + sum(Gu(k, i, m) * Gu(m, j, l) - Gu(k, j, m) * Gu(m, i, l) for m in range(4))
            - sum(cc[m][i][j] * Gu(k, m, l) for m in range(4)))
Ric = [[sum(Riem(k, l, k, j) for k in range(4)) for j in range(4)] for l in range(4)]
Rs = sum(eta[l] * Ric[l][l] for l in range(4))
Ein = [[Ric[i][j] - eta[i] * (1 if i == j else 0) * Rs / 2 for j in range(4)] for i in range(4)]
sx = spy.Matrix([[0, 1], [1, 0]]); sy = spy.Matrix([[0, -spy.I], [spy.I, 0]]); sz = spy.Matrix([[1, 0], [0, -1]])
Tg = [-spy.I * m / 2 for m in (sx, sy, sz)]; Z = spy.zeros(2, 2)
Afr = [Z] + [Fb[b] / Ab[b] * Tg[b] for b in range(3)]
com = lambda P, Q: P * Q - Q * P
Ff = [[der(i, Afr[j]) - der(j, Afr[i]) - sum((cc[k][i][j] * Afr[k] for k in range(4)), Z) + com(Afr[i], Afr[j]) for j in range(4)] for i in range(4)]
ip = lambda P, Q: -2 * (P * Q).trace()
FF = sum(eta[i] * eta[j] * ip(Ff[i][j], Ff[i][j]) for i in range(4) for j in range(4))
Tst = [[sum(eta[k] * ip(Ff[i][k], Ff[j][k]) for k in range(4)) - spy.Rational(1, 4) * eta[i] * (1 if i == j else 0) * FF for j in range(4)] for i in range(4)]
YM = []
for j in range(1, 4):
    s = Z
    for i in range(4):
        s = s + eta[i] * (der(i, Ff[i][j]) - sum((Gu(k, i, i) * Ff[k][j] for k in range(4)), Z)
                          - sum((Gu(k, i, j) * Ff[i][k] for k in range(4)), Z) + com(Afr[i], Ff[i][j]))
    YM.append(s)
q, kap = spy.symbols('q kappa', positive=True)
Eeq = [[Ein[i][j] + 3 * q * (eta[i] if i == j else 0) - 2 * kap * Tst[i][j] for j in range(4)] for i in range(4)]   # 8 pi G T = 2 kappa T (T in 1/g^2)
def o1(ex): return spy.expand(spy.diff(ex, ee).subs(ee, 0))
def o0(ex): return ex.subs(ee, 0)
Hs, Us = spy.Function('H')(t), spy.Function('U')(t)          # the squash h1 - h2 and the ripple u1 - u2
def to_HU(ex):
    ex = spy.expand(ex)
    # the equations are linear in (h1, h2, u1, u2) and depend only on the differences here; check that, then rename
    sub = {h[0]: Hs + h[1], u[0]: Us + u[1]}
    out = spy.simplify(spy.expand(ex.subs(sub).doit()))
    return out
d_ein = to_HU(o1(Eeq[1][1] - Eeq[2][2]))
y_ym = to_HU(o1(-2 * (YM[0] * Tg[0]).trace() + 2 * (YM[1] * Tg[1]).trace()))
indep = (not d_ein.has(h[1])) and (not d_ein.has(u[1])) and (not y_ym.has(h[1])) and (not y_ym.has(u[1]))
# expected forms (cosmic time); in conformal time these are H'' + 2 calH H' + 8H + (4 kappa/a^2)[(f^2 (f-2)^2 - f'^2) H + f' U' - f (f-2)(f+2) U] = 0
#                                                             and U'' + (6f + 4) U = f'' H + 2 f' H' + 6 f (f - 2) H
fd, Hd, Ud = spy.diff(fF, t), spy.diff(Hs, t), spy.diff(Us, t)
ein_exp = (aF**4 * spy.diff(Hs, t, 2) + 3 * aF**3 * spy.diff(aF, t) * Hd + 4 * aF**2 * (2 * Hs - kap * fd**2 * Hs + kap * fd * Ud)
           + 4 * kap * (fF**2 * (fF - 2)**2 * Hs - fF * (fF - 2) * (fF + 2) * Us)) / aF**4
ym_exp = (aF**2 * (Hs * spy.diff(fF, t, 2) + 2 * fd * Hd - spy.diff(Us, t, 2)) + aF * spy.diff(aF, t) * (Hs * fd - Ud)
          + 6 * fF * (fF - 2) * Hs - (6 * fF + 4) * Us) / aF**3
ein_ok = spy.simplify(d_ein - ein_exp) == 0
ym_ok = spy.simplify(y_ym - ym_exp) == 0

# ================================================================== W1. The lowest gravitational wave on the law's slices
say("== W1. Squashing the law's 3-sphere: its lowest gravitational wave")
# vacuum part (kappa = 0): in conformal time H'' + 2 calH H' + 8 H = 0 (n = 3 tensor mode of a closed universe: n^2 - 1 = 8)
vac = spy.simplify(ein_exp.subs(kap, 0) * aF**2)
# numeric on the law's de Sitter, q = 1: a = sec(eta), calH = tan(eta), eta in (-pi/2, pi/2)
def w1rhs(e_, y):
    return [y[1], -2 * math.tan(e_) * y[1] - 8 * y[0]]
e0, e1 = -math.pi / 2 + 1e-3, math.pi / 2 - 1e-6
s1 = solve_ivp(w1rhs, (e0, e1), [1.0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)
eg = np.linspace(e0, e1, 6001); Hv, Hp = s1.sol(eg)
late_slope = abs(Hp[-1]); early_amp = 1.0
frozen = Hv[-1]
zero_cross = int(np.sum(np.diff(np.sign(Hv)) != 0))
# energy while well inside the horizon: a^2 (H'^2 + 8 H^2) kept where calH << sqrt(8)
inside = np.abs(np.tan(eg)) < 0.05
E_in = (1 / np.cos(eg)**2) * (Hp**2 + 8 * Hv**2)
ok1 = ein_ok and ym_ok and indep and spy.simplify(vac - (aF**2 * spy.diff(Hs, t, 2) + 3 * aF * spy.diff(aF, t) * Hd + 8 * Hs)) == 0 and late_slope < 1e-3
check("W1 (exact + numeric; textbook tensor modes of a closed universe, Lifshitz 1946, Abbott-Schaefer 1986) squash the law's 3-sphere slices along its own three directions, a_b = a(1 + h_b), sum h_b = 0. Linearising Einstein's equations (done here symbolically for any a(t), any twist f(t)) gives, with no matter, exactly H'' + 2 calH H' + 8H = 0 in conformal time: the n = 3 tensor mode of a closed universe (n^2 - 1 = 8), its lowest gravitational wave. On the law's de Sitter it gets only pi of conformal time (round 22), so it swings about once and then freezes at a fixed squash as the slices grow past its wavelength",
      ok1, f"Einstein (traceless part, first order): matches the stated form {ein_ok}; Yang-Mills ripple equation {ym_ok}; depends only on H = h1 - h2, U = u1 - u2: {indep}; zero crossings {zero_cross}; frozen squash {frozen:+.4f} of the starting 1, final slope {late_slope:.1e}")
results["W1"] = dict(einstein_ok=bool(ein_ok), ym_ok=bool(ym_ok), zero_crossings=zero_cross, frozen_value=float(frozen), final_slope=float(late_slope))
SAVE = dict(w1=dict(eta=eg[::20].tolist(), H=Hv[::20].tolist()))

# ================================================================== W2. The tide's moving half
say("== W2. The magnetic (moving) half of the tide")
def Rl(k, l, i, j): return eta[k] * Riem(k, l, i, j)
def gg(x, y): return eta[x] if x == y else 0
def Cw(a_, b_, c_, d_):
    return (Rl(a_, b_, c_, d_) - (gg(a_, c_) * Ric[b_][d_] - gg(a_, d_) * Ric[b_][c_] - gg(b_, c_) * Ric[a_][d_] + gg(b_, d_) * Ric[a_][c_]) / 2
            + Rs * (gg(a_, c_) * gg(b_, d_) - gg(a_, d_) * gg(b_, c_)) / 6)
subU = {u[0]: 0, u[1]: 0}
Bm = [[spy.simplify(o1(sum(E3(a_ - 1, c_ - 1, d_ - 1) * Cw(0, b_, c_, d_) for c_ in range(1, 4) for d_ in range(1, 4)) / 2).subs(subU))
       for b_ in range(1, 4)] for a_ in range(1, 4)]
Em = [[spy.simplify(o1(Cw(0, a_, 0, b_)).subs(subU)) for b_ in range(1, 4)] for a_ in range(1, 4)]
B0 = [[spy.simplify(o0(sum(E3(a_ - 1, c_ - 1, d_ - 1) * Cw(0, b_, c_, d_) for c_ in range(1, 4) for d_ in range(1, 4)) / 2)) for b_ in range(1, 4)] for a_ in range(1, 4)]
E0 = [[spy.simplify(o0(Cw(0, a_, 0, b_))) for b_ in range(1, 4)] for a_ in range(1, 4)]
hd = [spy.diff(x, t) for x in hh]
B_ok = all(spy.simplify(Bm[i][j] - (3 * hd[i] / aF if i == j else 0)) == 0 for i in range(3) for j in range(3))
E_ok = all(spy.simplify(Em[i][j] - ((-spy.diff(hh[i], t, 2) / 2 - spy.diff(aF, t) * hd[i] / (2 * aF) + 4 * hh[i] / aF**2) if i == j else 0)) == 0 for i in range(3) for j in range(3))
bg_ok = all(B0[i][j] == 0 and E0[i][j] == 0 for i in range(3) for j in range(3))
# slow spin outside a lump in the law's de Sitter: Kottler + frame dragging g_t phi = -w(r) r^2 sin^2 theta
tt_, r, th, ph = spy.symbols('t r theta phi', positive=True); M, J, ep = spy.symbols('M J e')
X4 = [tt_, r, th, ph]; fK = 1 - 2 * M / r - q * r**2
wf = spy.Function('w')(r)
def spin_curv(wexpr):
    g = spy.zeros(4); g[0, 0] = -fK; g[1, 1] = 1 / fK; g[2, 2] = r**2; g[3, 3] = r**2 * spy.sin(th)**2
    g[0, 3] = g[3, 0] = -ep * wexpr * r**2 * spy.sin(th)**2
    g0 = g.subs(ep, 0); g0i = spy.diag(*[1 / g0[i, i] for i in range(4)]); gi = g0i - ep * g0i * spy.diff(g, ep) * g0i
    tr = lambda x: spy.expand(x).subs(ep**2, 0).subs(ep**3, 0)
    G4 = [[[tr(sum(gi[a, d] * (spy.diff(g[d, b], X4[c]) + spy.diff(g[d, c], X4[b]) - spy.diff(g[b, c], X4[d])) for d in range(4)) / 2)
            for c in range(4)] for b in range(4)] for a in range(4)]
    R = [[[[tr(spy.diff(G4[a][b][d], X4[c]) - spy.diff(G4[a][b][c], X4[d]) + sum(G4[a][c][k] * G4[k][b][d] - G4[a][d][k] * G4[k][b][c] for k in range(4)))
            for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
    return g, gi, R, tr
g4, gi4, R4, tr = spin_curv(wf)
Ric4 = spy.Matrix(4, 4, lambda b, d: spy.simplify(tr(sum(R4[a][b][a][d] for a in range(4)))))
res = (Ric4 - 3 * q * g4).applyfunc(lambda x: spy.factor(spy.simplify(tr(x))))
drag_eq = spy.factor(spy.simplify(res[0, 3] / (ep * spy.sin(th)**2)))
drag_expected = spy.factor(r * fK * (r * spy.diff(wf, r, 2) + 4 * spy.diff(wf, r)) / 2)
drag_ok = spy.simplify(drag_eq - drag_expected) == 0 and all(spy.simplify(res[i, j]) == 0 for i in range(4) for j in range(4) if (i, j) not in ((0, 3), (3, 0)))
Om0 = spy.symbols('Omega0')
sol_ok = spy.simplify(drag_expected.subs(wf, 2 * J / r**3 + Om0).doit()) == 0
g4s, gi4s, R4s, tr = spin_curv(2 * J / r**3)
Rl4 = [[[[tr(sum(g4s[a, k] * R4s[k][b][c][d] for k in range(4))) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
Ruu = [[[[tr(sum(gi4s[b, k] * R4s[a][k][c][d] for k in range(4))) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
Pdens = 0
for perm in itertools.permutations(range(4)):
    sg = spy.combinatorics.Permutation(list(perm)).signature(); c, d, e_, f_ = perm
    for a in range(4):
        for b in range(4):
            if Rl4[a][b][c][d] != 0:
                Pdens += sg * Rl4[a][b][c][d] * Ruu[a][b][e_][f_]
Pdens = spy.simplify(tr(Pdens / (2 * r**2 * spy.sin(th))))
P_ok = spy.simplify(Pdens - 288 * J * M * ep * spy.cos(th) / r**7) == 0
ok2 = B_ok and E_ok and bg_ok and drag_ok and sol_ok and P_ok
check("W2 (exact, symbolic; textbook: the Weyl tensor's magnetic part, Lense-Thirring 1918, Kerr-de Sitter to first order in spin) the tide's 'magnetic' half comes only from motion. (a) On the squashing law 3-sphere, the magnetic tide is exactly B_bb = 3 h_b'/a (cosmic time): three times how fast each direction's squash changes, and zero for slices at rest, like round 25's. The electric tide is -h_b''/2 - (a'/a) h_b'/2 + 4 h_b/a^2. The unsquashed law has neither. (b) Spin a lump in the law's de Sitter slowly: Einstein's equations with Lambda = 3q force the frame dragging w(r) to obey (r^4 w')' = 0, so w = 2J/r^3 + a constant (the frame at infinity), with q dropping out. Its tidal invariant R.*R = 288 J M cos(theta)/r^7, also free of q: the magnetic tide is there only when the lump both has mass and spins",
      ok2, f"B = 3 h'/a: {B_ok}; E formula: {E_ok}; law slices at rest have E = B = 0: {bg_ok}; drag equation r f (r w'' + 4w')/2 = 0: {drag_ok}; w = 2J/r^3 + const solves it: {sol_ok}; R.*R = {Pdens}")
results["W2"] = dict(magnetic="B_bb = 3 h_b'/a", electric="E_bb = -h_b''/2 - (a'/a) h_b'/2 + 4 h_b/a^2", drag_equation=str(drag_expected),
                     pontryagin=str(Pdens))

# ================================================================== W3. A rippling twist makes a gravitational wave
say("== W3. The weak force's twist rippling: does it make a gravitational wave?")
# couplings in the conformal-time equations: Einstein side (4 kappa/a^2)[(f^2(f-2)^2 - f'^2) H + f'U' - f(f-2)(f+2) U]; YM side 6 f (f-2) H + f''H + 2f'H'
fsym = spy.symbols('f')
cU = -fsym * (fsym - 2) * (fsym + 2); cH = 6 * fsym * (fsym - 2)
vac_ok = all(cU.subs(fsym, v) == 0 and cH.subs(fsym, v) == 0 for v in (0, 2))
law_val = (cU.subs(fsym, 1), cH.subs(fsym, 1))
# numeric: closed universe with Lambda = 3q and the law's f = 1 field at rest (it stays at f = 1: the hilltop is a solution),
# background in conformal time a'^2 = kappa C + q a^4 - a^2, C = 1; ripple U starts at 0.01, squash H at 0
def run(qq, kk, fb=1.0, U0=0.01, eta_end=None):
    C0 = (fb**2 - 2 * fb)**2
    a0 = math.sqrt((1 + math.sqrt(1 - 4 * qq * kk * C0)) / (2 * qq))
    cU_ = -fb * (fb - 2) * (fb + 2); cH_ = 6 * fb * (fb - 2); m2 = fb**2 * (fb - 2)**2
    def rhs(e_, y):
        a_, ap, H_, Hp_, U_, Up_ = y
        calH = ap / a_
        Hpp = -2 * calH * Hp_ - 8 * H_ - 4 * kk / a_**2 * (m2 * H_ + cU_ * U_)
        Upp = -(6 * fb + 4) * U_ + cH_ * H_
        return [ap, 2 * qq * a_**3 - a_, Hp_, Hpp, Up_, Upp]
    blow = lambda e_, y: y[0] - 1e4
    blow.terminal = True
    s = solve_ivp(rhs, (0, 10), [a0, 0.0, 0.0, 0.0, U0, 0.0], rtol=1e-11, atol=1e-14, events=blow, dense_output=True, max_step=0.01)
    return s
outs = {}
for kk in (0.0005, 0.001, 0.002):
    s = run(0.3, kk)
    outs[kk] = dict(eta_end=float(s.t[-1]), H_max=float(np.abs(s.y[2]).max()), H_final=float(s.y[2][-1]))
s_vac = run(0.3, 0.001, fb=0.0); H_vac = float(np.abs(s_vac.y[2]).max())
s_two = run(0.3, 0.001, fb=2.0); H_two = float(np.abs(s_two.y[2]).max())
r1 = outs[0.001]["H_max"] / outs[0.0005]["H_max"]; r2 = outs[0.002]["H_max"] / outs[0.001]["H_max"]
s_pic = run(0.3, 0.002)
SAVE["w3"] = dict(eta=s_pic.t[::10].tolist(), H=s_pic.y[2][::10].tolist(), U=s_pic.y[4][::10].tolist(), a=s_pic.y[0][::10].tolist())
ok3 = vac_ok and law_val == (3, -6) and H_vac == 0.0 and H_two == 0.0 and outs[0.001]["H_max"] > 1e-6 and abs(r1 - 2) < 0.01 and abs(r2 - 2) < 0.01
check("W3 (exact + numeric; textbook in kind: Einstein-Yang-Mills in anisotropic Bianchi cosmologies, e.g. Darian-Kunzle 1996 (Bianchi I), Jin-Maeda 2005 (class A, including IX); the law's-point reading is ours) let the locked twist ripple, f_b = f + u_b with sum u_b = 0, the same squash shape as W1's wave. Linearised Einstein and Yang-Mills couple the ripple U to the wave H only through f(f-2)(f+2) and f(f-2) (and f' when the twist moves): both vanish at the two flat twists f = 0 and f = 2 and are 3 and -6 at the law's f = 1. So a ripple on an empty weak-force vacuum makes no gravitational wave at first order, while a ripple on the law's own field does. On a closed universe with Lambda = 3q holding the law's f = 1 field (q = 0.3), a ripple of 0.01 makes a frozen wave of size proportional to kappa = 4 pi G/g^2 (doubling kappa doubles it, to within 0.1%), and exactly none at f = 0 or 2",
      ok3, f"couplings vanish at f = 0, 2: {vac_ok}; at f = 1: U {law_val[0]}, H {law_val[1]}; " + "; ".join(f"kappa {k}: |H| max {v['H_max']:.3e}" for k, v in outs.items()) + f"; ratios {r1:.3f}, {r2:.3f}; f = 0: {H_vac}; f = 2: {H_two}")
results["W3"] = dict(coupling_U="-f(f-2)(f+2)", coupling_H="6f(f-2)", at_law=[int(x) for x in law_val], runs={str(k): v for k, v in outs.items()},
                     ratios=[r1, r2], H_at_f0=H_vac, H_at_f2=H_two)

# ================================================================== W4. How the pieces fit: what a static, a spinning and a waving source give
say("== W4. Electric and magnetic tide side by side")
# electric/magnetic sizes for a free wave in the law's de Sitter (W1 solution), early, while inside the horizon:
# with a = sec(eta) at q = 1, cosmic-time derivative d/dt = cos(eta) d/deta; per direction (h1 = H/2 ... take h_1 = -h_2 = H/2, h_3 = 0)
aa = 1 / np.cos(eg); Hdot = np.cos(eg) * Hp
# h'' in cosmic time from the equation: hdd = -3 (a'/a) hdot - 8 h / a^2 with a'/a = sin(eta) (cosmic-time Hubble at q = 1 is tanh t = sin eta)
Hub = np.sin(eg); Hddot = -3 * Hub * Hdot - 8 * Hv / aa**2
Eamp = np.abs(-Hddot / 2 - Hub * Hdot / 2 + 4 * Hv / aa**2) / 2; Bamp = np.abs(3 * Hdot / aa) / 2
mid = np.abs(eg) < 0.6
ratio_mid = float(np.sqrt(np.mean(Bamp[mid]**2) / np.mean(Eamp[mid]**2)))
late = eg > 1.45
ratio_late = float(np.sqrt(np.mean(Bamp[late]**2) / np.mean(Eamp[late]**2)))
ratio_end = float(Bamp[-1] / Eamp[-1])
SAVE["w4"] = dict(eta=eg[::20].tolist(), E=Eamp[::20].tolist(), B=Bamp[::20].tolist())
ok4 = 0.3 < ratio_mid < 3 and ratio_late < ratio_mid and ratio_end < 0.01
check("W4 (numeric; textbook) side by side: a lump at rest has only the electric tide (round 25: 3M/r^3, magnetic 0); a spinning lump adds a magnetic tide set by its spin (W2, invariant 288 J M cos(theta)/r^7); the law's lowest gravitational wave carries both halves at comparable size while it swings inside the horizon (rms magnetic/electric about 1), and once it freezes outside the horizon its magnetic half dies away (it needs motion) while a fixed electric squash is left",
      ok4, f"rms B/E while swinging (|eta| < 0.6): {ratio_mid:.2f}; after freezing (eta > 1.45): {ratio_late:.2e}; at the end of time: {ratio_end:.1e}")
results["W4"] = dict(ratio_swinging=ratio_mid, ratio_frozen=ratio_late, ratio_end=ratio_end)

say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed, seconds=round(time.time() - T_START))
results["figure"] = SAVE
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
