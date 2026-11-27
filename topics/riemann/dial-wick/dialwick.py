"""Dial-Wick: turn one of our dials so the Riemann zero question becomes a shock (Burgers) question,
look at it there, and turn back.  DRAFT, Nov 27 2026.

Conventions (as on the earlier riemann pages):
  H_t(z) = int e^{t u^2} Phi(u) cos(z u) du,  dH/dt = -H_zz,  zeros at z = 2*gamma.
  RH  <=>  H_0 has only real zeros  <=>  Lambda <= 0.

Three turns are tested.
  1. Quarter turn (d -> -d in the plane law, i.e. x -> i x, tan <-> tanh):  G_t(w) = H_t(i w).
  2. Half turn in the metric only (the d value of the gradient metric dx^2 + d dy^2, called kappa here).
  3. Viscosity phase (complex zoom by e^{i theta}):  nu(theta) = -e^{-2 i theta}.
"""
import json, math
import numpy as np
from numpy.polynomial import polynomial as Pn
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

out = {}

# ---------------------------------------------------------------- kernel and the two frames
def Phi(u):
    s = np.zeros_like(u)
    for n in range(1, 10):
        s += (2 * np.pi**2 * n**4 * np.exp(9 * u) - 3 * np.pi * n**2 * np.exp(5 * u)) * np.exp(-np.pi * n**2 * np.exp(4 * u))
    return s

U = 1.4
xg, wg = np.polynomial.legendre.leggauss(1200)
uq = (xg + 1) * U / 2; wq = wg * U / 2
PHI = Phi(uq)

def H(t, z):            # original frame: backward heat in t (dH/dt = -H_zz)
    z = np.atleast_1d(np.asarray(z, dtype=complex))
    return 2 * (np.cos(np.outer(z, uq)) * (np.exp(t * uq**2) * PHI * wq)).sum(1)

def G(t, w, deriv=0):   # rotated frame G_t(w) = H_t(i w) = 2 int e^{t u^2} Phi cosh(w u) du
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    k = [np.cosh, np.sinh][deriv % 2](np.outer(w, uq)) * uq**deriv
    return 2 * (k * (np.exp(t * uq**2) * PHI * wq)).sum(1)

# 1a. the quarter turn is exact: G_t(w) = H_t(iw), and G obeys FORWARD heat dG/dt = +G_ww
ws = np.array([0.3, 2.0 + 1.0j, 7.5, 11.0 - 4.0j])
rot_err = float(np.max(np.abs(G(0.1, ws) - H(0.1, 1j * ws)) / np.abs(G(0.1, ws))))
h = 1e-4; e = 1e-2
dGdt = (G(0.1 + h, ws) - G(0.1 - h, ws)) / (2 * h)
Gww = (G(0.1, ws + e) - 2 * G(0.1, ws) + G(0.1, ws - e)) / e**2
out["quarter_turn"] = {"G_equals_H_of_iw_relerr": rot_err,
                       "forward_heat_relerr": float(np.max(np.abs(dGdt - Gww) / np.abs(Gww)))}

# 1b. in the rotated frame phi = G_t(x) > 0 on the real line for EVERY t (integrand positive since Phi > 0)
xs = np.linspace(-60, 60, 601)
tgrid = np.linspace(-1.0, 1.0, 21)
minG = min(float(np.min(G(t, xs).real)) for t in tgrid)
out["physical_every_t"] = {"min_G_over_grid": minG, "Phi_min_on_quad": float(PHI.min()),
                           "t_range": [-1, 1], "x_range": [-60, 60]}

# 1c. Burgers field u = -2 G_x / G solves u_t + u u_x = u_xx (viscosity +1, forward): checked by differences
def ub(t, x):
    return (-2 * G(t, x, 1) / G(t, x)).real
x0 = np.array([0.5, 3.0, 10.0, 25.0]); ht = 1e-4; hx = 1e-3
ut = (ub(ht, x0) - ub(-ht, x0)) / (2 * ht)
ux = (ub(0, x0 + hx) - ub(0, x0 - hx)) / (2 * hx)
uxx = (ub(0, x0 + hx) - 2 * ub(0, x0) + ub(0, x0 - hx)) / hx**2
res = ut + ub(0, x0) * ux - uxx
out["burgers_residual"] = {"max_abs": float(np.max(np.abs(res))), "scale": float(np.max(np.abs(uxx)))}

# 1d. poles of u in the rotated frame = zeros of G = i * (zeros of H): on the imaginary axis at 2*gamma
ys = np.linspace(1, 70, 3000)
Hy = H(0.0, ys).real
roots = [brentq(lambda y: H(0.0, y)[0].real, ys[i], ys[i + 1]) for i in range(len(ys) - 1) if Hy[i] * Hy[i + 1] < 0]
gam_known = np.array([14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588])
out["poles_on_centerline"] = {"heights": roots[:5], "two_gamma": (2 * gam_known).tolist(),
                              "max_err": float(np.max(np.abs(np.array(roots[:5]) - 2 * gam_known)))}

# 1e. magnetisation reading: G(h) = int e^{h u} mu(du) with mu symmetric, so G'/G = mean spin m(h), u_Burgers = -2 m(h)
hs = np.linspace(0, 40, 401)
mag = {}
for t in [-0.3, 0.0, 0.3]:
    m = (G(t, hs, 1) / G(t, hs)).real
    mag[str(t)] = m.tolist()
mono = all(np.all(np.diff(np.array(v)) > 0) for v in mag.values())
out["magnetisation"] = {"increasing_in_h_all_t": bool(mono), "m_at_h40": {k: v[-1] for k, v in mag.items()},
                        "spin_range_max_u": U}

# complex-plane picture of |G_0| (rotated frame) and |H_0| near the first zeros
WX, WY = np.meshgrid(np.linspace(-12, 12, 161), np.linspace(-70, 70, 321))
Wg = (WX + 1j * WY).ravel()
logG = np.log10(np.abs(G(0.0, Wg))).reshape(WX.shape)

# ---------------------------------------------------------------- 2. exact toy: H = c + e^t cos z  (Lambda = ln c)
c = 1.5; Lam = math.log(c)
# also a two-mode check that it solves dH/dt = -H_zz is immediate; record the zero / pole tracks
ttoy = np.linspace(-0.4, 1.0, 141)
def toy_zero(t):      # zero of c + e^t cos z nearest z = pi (upper/right branch)
    r = c * math.exp(-t)
    return math.pi + (1j * math.acosh(r) if r > 1 else math.acos(r))
zt = np.array([toy_zero(t) for t in ttoy])
wt = -1j * zt                                   # rotated frame: pole at w = -i z
out["toy"] = {"c": c, "Lambda": Lam,
              "off_line_before": bool(np.all(np.abs(zt[ttoy < Lam - 1e-9].imag) > 0)),
              "real_after": bool(np.all(np.abs(zt[ttoy >= Lam + 1e-9].imag) < 1e-12)),
              "pole_column_x_at_t_minus0p4": float(abs(wt[0].real)),
              "columns_merge_at_x0": bool(abs(toy_zero(Lam).imag) < 1e-6)}
# spin reading: spin value 0 with weight c, +-1 with weight e^t/2 each -> Z(h) = c + e^t cosh h;
# Lee-Yang (zeros only at imaginary h) iff t >= ln c.  Numeric check of LY failure/success:
def toy_LY(t):
    r = c * math.exp(-t)          # cosh h = -r :  h = i pi +- acosh(r) (off axis) if r > 1 else i(pi +- acos r)
    return r <= 1
out["toy"]["lee_yang_at"] = {str(t): toy_LY(t) for t in [0.0, 0.3, Lam, 0.6]}
np.savez("toy_tracks.npz", t=ttoy, z=zt, w=wt)

# ---------------------------------------------------------------- 1f. pole ODE in the rotated frame (polynomial, exact)
def heat_poly(cf, s):     # e^{s D^2} on a polynomial
    outc = np.zeros_like(cf); term = cf.copy(); k = 0
    while True:
        outc[:len(term)] += term * s**k / math.factorial(k)
        if len(term) <= 2: break
        term = Pn.polyder(term, 2); k += 1
    return outc
w0 = np.array([0.0 + 3j, 0.0 - 3j, 0.0 + 5.5j, 0.0 - 5.5j, 0.6 + 1.2j, 0.6 - 1.2j, -0.6 + 1.2j, -0.6 - 1.2j])
cf = Pn.polyfromroots(w0).real
def vel(w):
    D = w[:, None] - w[None, :]; np.fill_diagonal(D, 1.0); R = 1 / D; np.fill_diagonal(R, 0)
    return -2 * R.sum(1)
s1 = 0.05
r_exact = Pn.polyroots(heat_poly(cf, s1))
sol = solve_ivp(lambda t, v: np.concatenate([vel(v[:8] + 1j * v[8:]).real, vel(v[:8] + 1j * v[8:]).imag]),
                [0, s1], np.concatenate([w0.real, w0.imag]), rtol=1e-12, atol=1e-14, method="DOP853")
r_ode = sol.y[:8, -1] + 1j * sol.y[8:, -1]
err = max(min(abs(r - q) for q in r_exact) for r in r_ode)
out["rotated_pole_ode"] = {"max_err": float(err), "s": s1,
                           "rule": "dw/dt = -2 sum 1/(w_j - w_k)  (viscous Burgers pole law, forward)"}
# merger of the two side columns (they start at x = +-0.6) under forward heat in the rotated frame
tt = np.linspace(0, 0.6, 601)
colx = []
for t in tt:
    r = Pn.polyroots(heat_poly(cf, t))
    colx.append(float(np.max(np.abs(r.real))))
colx = np.array(colx)
tm = float(tt[np.argmax(colx < 1e-6)]) if np.any(colx < 1e-6) else None
out["rotated_merger"] = {"merge_time": tm}
np.savez("merger.npz", t=tt, colx=colx)

# ---------------------------------------------------------------- 3. half turn in the metric: kappa = d of the gradient metric
# law: v_j = 2 sum (dx + i kappa dy)/r^2 = (1-kappa)/2 * Hol + (1+kappa)/2 * conj(Hol),  Hol = 2 sum 1/(z_j - z_k)
rng = np.random.default_rng(3)
zz = rng.normal(0, 2, 6) + 1j * rng.normal(0, 2, 6)
def v_kappa_direct(z, kap):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, 1.0)
    V = 2 * (D.real + 1j * kap * D.imag) / np.abs(D)**2; np.fill_diagonal(V, 0); return V.sum(1)
def v_kappa_split(z, kap):
    D = z[:, None] - z[None, :]; np.fill_diagonal(D, 1.0); R = 1 / D; np.fill_diagonal(R, 0)
    Hol = 2 * R.sum(1); return (1 - kap) / 2 * Hol + (1 + kap) / 2 * np.conj(Hol)
out["kappa_split_identity_err"] = float(max(np.max(np.abs(v_kappa_direct(zz, k) - v_kappa_split(zz, k)))
                                            for k in [-1, -0.3, 0, 0.7, 1]))
out["kappa_minus1_is_zeta_law"] = float(np.max(np.abs(v_kappa_direct(zz, -1) - v_kappa_split(zz, -1) )))

# linear growth of a vertical ripple on a crystal of spacing 1 (periodic ring of N, image sums exact)
N = 400
n = np.arange(1, N)
Dinv2 = (np.pi / N)**2 / np.sin(np.pi * n / N)**2
qs = 2 * np.pi * np.arange(1, N // 2) / N
rate_ring = np.array([2 * np.sum((1 - np.cos(q * n)) * Dinv2) for q in qs])   # times kappa
rate_formula = 2 * (np.pi * qs - qs**2 / 2)
out["kappa_dispersion"] = {"max_relerr": float(np.max(np.abs(rate_ring - rate_formula) / rate_formula)),
                           "law": "vertical ripple grows at kappa * 2(pi q - q^2/2); horizontal ripple decays at 2(pi q - q^2/2) for every kappa"}
np.savez("dispersion.npz", q=qs, rate=rate_formula)

# an off-line pair inside a crystal (ring of L = 40, spacing 1), pair at x = 0.5, height 0.3, for kappa from -1 to +1
L = 40
xr = np.arange(L, dtype=float); xr = xr[xr != 20]          # remove one bead, add the pair: count stays L+1
def ring_vel(z, kap):
    D = z[:, None] - z[None, :]
    C = (np.pi / L) / np.tan(np.pi * np.where(np.eye(len(z), dtype=bool), 1.0, D) / L)
    np.fill_diagonal(C, 0)
    Hol = 2 * C.sum(1)
    return (1 - kap) / 2 * Hol + (1 + kap) / 2 * np.conj(Hol)
z0r = np.concatenate([xr.astype(complex), [20.0 + 0.3j, 20.0 - 0.3j]])
pair_tracks = {}
for kap in [-1.0, -0.5, 0.0, 0.5, 1.0]:
    def rhs(t, v): zc = v[:len(v)//2] + 1j * v[len(v)//2:]; vv = ring_vel(zc, kap); return np.concatenate([vv.real, vv.imag])
    def hit(t, v): return v[len(v)//2 + len(xr)] - 1e-3
    hit.terminal = True
    sol = solve_ivp(rhs, [0, 0.12], np.concatenate([z0r.real, z0r.imag]), rtol=1e-10, atol=1e-12,
                    max_step=1e-3, events=hit, method="DOP853")
    yy = sol.y[len(z0r) + len(xr)]
    pair_tracks[str(kap)] = {"t": sol.t.tolist(), "y": yy.tolist()}
out["kappa_pair"] = {k: {"y_end": v["y"][-1], "t_end": v["t"][-1]} for k, v in pair_tracks.items()}
out["kappa_pair"]["capture_time_formula_kappa_-1"] = "lone pair: y^2/2 = 0.045; the crowd shortens it, as t = ln(1+2Sy^2)/(4S) <= y^2/2 on the spin-orbit page"
json.dump(pair_tracks, open("kappa_pairs.json", "w"))

# ---------------------------------------------------------------- 4. viscosity phase: rotate the plane by theta
# G_theta(w) = H(e^{i theta} w)  =>  dG/dt = -e^{-2 i theta} G_ww,  nu(theta) = -e^{-2 i theta}
nu = {str(th): complex(-np.exp(-2j * th)) for th in [0, np.pi / 4, np.pi / 2]}
out["viscosity_phase"] = {k: [v.real, v.imag] for k, v in nu.items()}
th = np.pi / 4; e2 = 1e-2
wsx = np.array([1.0 + 0.5j, 6.0])
Gth = lambda t, w: H(t, np.exp(1j * th) * w)
lhs = (Gth(0.1 + 1e-4, wsx) - Gth(0.1 - 1e-4, wsx)) / 2e-4
rhs_ = -np.exp(-2j * th) * (Gth(0.1, wsx + e2) - 2 * Gth(0.1, wsx) + Gth(0.1, wsx - e2)) / e2**2
out["viscosity_phase_check_relerr"] = float(np.max(np.abs(lhs - rhs_) / np.abs(rhs_)))

# ---------------------------------------------------------------- save
np.savez("frames.npz", WX=WX, WY=WY, logG=logG, hs=hs, **{f"m_{k}": np.array(v) for k, v in mag.items()},
         xs=xs, u0=ub(0.0, xs), um=ub(-0.3, xs), up=ub(0.3, xs), roots=np.array(roots))
json.dump(out, open("dialwick.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
