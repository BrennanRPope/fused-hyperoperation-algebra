"""DRAFT. Depth checks for the wall-shapes study.

Claims tested (proofs are in README.md):
  D1  Rayleigh n=2:  R = (1+k a1)(1+k a2) - (1+k c) is real stable iff 1+kc >= 0
      (Branden's multiaffine criterion: beta*gamma - alpha*delta = k^2 (1+kc)).
  D2  Line test: for real x and positive v, t -> R(x + t v) has only real roots
      (n=2, W>=0), and fails for n=2, W<0 and for n>=3, W != 0.
  D3  Explicit upper-half-plane witnesses for n>=3 (exact, Gaussian integers).
  D4  Sector law: inputs with arg in (0, alpha), depth n.  W>0 is missed iff n*alpha <= 2 pi,
      W<0 is missed iff n*alpha <= pi.
  D5  Lee-Yang two-spin dictionary: Cayley z=(x-i)/(x+i) turns 1+a(z1+z2)+z1 z2 into
      (1+a) x1 x2 - (1-a), i.e. x1 x2 = W with W = (1-a)/(1+a) = tanh J (a = e^{-2J}).
  D6  Pure n-body parity coupling: Z = prod(1+z_j) + T prod(1-z_j), T = tanh J_n.
      Cayley image: prod x_j = -T i^n.  Diagonal zeros all on |z|=1 only for n=2, T>0.
  D7  n=4, J4 -> infinity (T -> 1): diagonal zeros at |z| = sqrt2 +- 1.
  D8  N=4 sites, pair couplings J2 on all pairs + one 4-body J4: diagonal Lee-Yang map,
      plus multivariate line test at sample points.
Uses only numpy and fractions (our own tools).  Writes depth.json.
"""
import itertools, json, cmath, math
from fractions import Fraction as F
import numpy as np

rng = np.random.default_rng(20261127)
out = {}

# ---------- D1 Rayleigh / Branden ----------
def branden_n2(k, c):
    alpha, beta, gamma, delta = -k * c, k, k, k * k   # R = alpha + beta a1 + gamma a2 + delta a1 a2
    return beta * gamma - alpha * delta
d1 = []
for k in [F(1), F(2), F(-3, 2), F(1, 3)]:
    for c in [F(-2), F(-1), F(-1, 2), F(0), F(1), F(5, 2)]:
        lhs = branden_n2(k, c); rhs = k * k * (1 + k * c)
        d1.append(dict(k=str(k), c=str(c), branden=str(lhs), k2_1pkc=str(rhs), equal=(lhs == rhs)))
out["D1_rayleigh_identity_all_equal"] = all(r["equal"] for r in d1)

# ---------- D2 line test ----------
def R_coeffs_on_line(k, W, x, v):
    """coefficients (numpy poly1d) of t -> prod(1+k(x_i+t v_i)) - W."""
    p = np.poly1d([1.0])
    for xi, vi in zip(x, v):
        p = p * np.poly1d([k * vi, 1 + k * xi])
    return p - np.poly1d([W])
def line_test(n, k, W, trials=4000):
    worst = 0.0
    for _ in range(trials):
        x = rng.normal(0, 2, n); v = rng.exponential(1, n) + 1e-3
        r = np.roots(R_coeffs_on_line(k, W, x, v).coeffs)
        worst = max(worst, float(np.max(np.abs(r.imag))) if len(r) else 0.0)
    return worst
d2 = {}
for n in [2, 3, 4]:
    for W in [2.0, 0.5, 0.0, -0.5]:
        d2[f"n={n},W={W}"] = line_test(n, 1.0, W)
out["D2_line_test_max_imag_part"] = d2
out["D2_pass"] = (d2["n=2,W=2.0"] < 1e-6 and d2["n=2,W=0.5"] < 1e-6 and d2["n=2,W=0.0"] < 1e-6
                  and d2["n=2,W=-0.5"] > 1e-3 and all(d2[f"n={n},W={W}"] > 1e-3
                  for n in [3, 4] for W in [2.0, 0.5, -0.5]) and d2["n=3,W=0.0"] < 1e-6)

# ---------- D3 Gaussian-integer witnesses (u_i = 1 + k a_i, all Im u_i > 0) ----------
def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def find_witness(n, sign, box=4):
    pts = [(a, b) for a in range(-box, box + 1) for b in range(1, box + 1)]
    for combo in itertools.combinations_with_replacement(pts, n):
        p = (1, 0)
        for u in combo: p = gmul(p, u)
        if p[1] == 0 and p[0] != 0 and (p[0] > 0) == (sign > 0):
            return dict(u=[f"{a}{'+' if b >= 0 else '-'}{abs(b)}i" for a, b in combo], product=p[0])
    return None
d3 = {f"n={n},W{'>' if s > 0 else '<'}0": find_witness(n, s) for n in [3, 4] for s in [1, -1]}
out["D3_witnesses"] = d3
out["D3_pass"] = all(v is not None for v in d3.values())

# ---------- D4 sector law ----------
def sector_hits(n, alpha, target_arg):
    """exact: products of n points with arg in (0, alpha) reach every arg in (0, n alpha), any modulus."""
    t = target_arg
    while t <= 0: t += 2 * math.pi
    hits = False
    while t < n * alpha - 1e-12:
        hits = True; break
    if hits:   # construct the witness and check it numerically
        u = cmath.rect(1.0, t / n)
        assert 0 < t / n < alpha
        prod = u ** n
        assert abs(cmath.phase(prod) - cmath.phase(cmath.rect(1, target_arg))) < 1e-9
    return hits
d4 = []
for n in [1, 2, 3, 4, 6]:
    for alpha_over_pi in [0.25, 0.4, 0.5, 2 / 3, 0.75, 1.0]:
        a = alpha_over_pi * math.pi
        pos = sector_hits(n, a, 0.0)          # W > 0 sits at arg 2 pi
        neg = sector_hits(n, a, math.pi)      # W < 0 sits at arg pi
        d4.append(dict(n=n, alpha_over_pi=round(alpha_over_pi, 4), Wpos_free=not pos, Wneg_free=not neg,
                       rule_pos=n * alpha_over_pi <= 2 + 1e-12, rule_neg=n * alpha_over_pi <= 1 + 1e-12))
out["D4_sector_rows"] = d4
out["D4_pass"] = all(r["Wpos_free"] == r["rule_pos"] and r["Wneg_free"] == r["rule_neg"] for r in d4)

# ---------- D5 Lee-Yang two-spin Cayley dictionary ----------
def cay(x): return (x - 1j) / (x + 1j)
d5 = []
for J in [-1.0, -0.3, 0.0, 0.4, 1.5]:
    a = math.exp(-2 * J); worst = 0.0
    for _ in range(200):
        x1, x2 = rng.normal(size=2) + 1j * rng.normal(size=2)
        z1, z2 = cay(x1), cay(x2)
        P = (1 + a * (z1 + z2) + z1 * z2) * (x1 + 1j) * (x2 + 1j)
        Q = 2 * ((1 + a) * x1 * x2 - (1 - a))
        worst = max(worst, abs(P - Q))
    W = (1 - a) / (1 + a)
    d5.append(dict(J=J, W=W, tanhJ=math.tanh(J), identity_err=worst))
out["D5_rows"] = d5
out["D5_pass"] = all(r["identity_err"] < 1e-9 and abs(r["W"] - r["tanhJ"]) < 1e-12 for r in d5)

# ---------- partition-function helper ----------
def Zcoef(N, terms):
    """diagonal: coefficients of Z(z) = sum_sigma exp(sum_A J_A prod sigma_A) z^{#minus}."""
    coef = np.zeros(N + 1)
    for s in itertools.product([1, -1], repeat=N):
        E = sum(J * np.prod([s[i] for i in A]) for J, A in terms)
        coef[s.count(-1)] += math.exp(E)
    return coef
def diag_moduli(N, terms):
    return np.sort(np.abs(np.roots(Zcoef(N, terms)[::-1])))

# ---------- D6 pure n-body parity coupling ----------
d6 = []
for n in [2, 4, 6]:
    for J in [0.5, -0.5, 1.2]:
        T = math.tanh(J)
        num = diag_moduli(n, [(J, tuple(range(n)))])
        # closed form: x^n = -T i^n, z = (x-i)/(x+i)
        xs = np.roots([1] + [0] * (n - 1) + [T * (1j ** n)])
        closed = np.sort(np.abs((xs - 1j) / (xs + 1j)))
        all_on_circle = bool(np.all(np.abs(num - 1) < 1e-9))
        d6.append(dict(n=n, J=J, T=T, moduli=[round(float(m), 6) for m in num],
                       closed_form_match=bool(np.max(np.abs(num - closed)) < 1e-8),
                       all_on_unit_circle=all_on_circle))
out["D6_rows"] = d6
out["D6_pass"] = all(r["closed_form_match"] for r in d6) and all(
    r["all_on_unit_circle"] == (r["n"] == 2 and r["J"] > 0) for r in d6)

# ---------- D7 J4 -> infinity ----------
lim = diag_moduli(4, [(12.0, (0, 1, 2, 3))])
out["D7_moduli_J4_12"] = [float(m) for m in lim]
out["D7_target"] = [math.sqrt(2) - 1, math.sqrt(2) + 1]
out["D7_pass"] = bool(abs(lim[0] - (math.sqrt(2) - 1)) < 1e-6 and abs(lim[-1] - (math.sqrt(2) + 1)) < 1e-6)

# ---------- D8 pair + 4-body map (N=4) ----------
pairs = list(itertools.combinations(range(4), 2))
def terms_of(J2, J4): return [(J2, p) for p in pairs] + [(J4, (0, 1, 2, 3))]
J2s = np.round(np.linspace(-0.6, 1.2, 37), 4); J4s = np.round(np.linspace(-1.2, 1.2, 49), 4)
grid = []
for J4 in J4s:
    row = []
    for J2 in J2s:
        m = diag_moduli(4, terms_of(float(J2), float(J4)))
        row.append(float(np.max(np.abs(np.log(m)))))   # 0 means every diagonal zero on |z|=1
    grid.append(row)
out["D8_grid"] = dict(J2=J2s.tolist(), J4=J4s.tolist(), max_abs_log_modulus=grid)

def multivariate_line_test(N, terms, trials=1500):
    """Cayley-transformed Z: real stable iff every line x + t v (x real, v>0) is real-rooted.
       Z(z) multilinear; with z_j = (x_j - i)/(x_j + i), multiply by prod (x_j + i)."""
    configs = list(itertools.product([1, -1], repeat=N))
    w = [math.exp(sum(J * np.prod([s[i] for i in A]) for J, A in terms)) for s in configs]
    worst = 0.0
    for trial in range(trials):
        if trial == 0:   # the diagonal line itself
            x = np.zeros(N); v = np.ones(N)
        else:            # half the lines near the diagonal, half spread out
            x = rng.normal(0, 0.3 if trial % 2 else 2.0, N); v = rng.exponential(1, N) + 1e-3
        poly = np.poly1d([0.0])
        for s, wt in zip(configs, w):
            p = np.poly1d([wt])
            for j in range(N):   # spin up -> factor (x+i), spin down -> z (x+i) = (x - i)
                p = p * np.poly1d([v[j], x[j] + (1j if s[j] == 1 else -1j)])
            poly = poly + p
        r = np.roots(poly.coeffs)
        worst = max(worst, float(np.max(np.abs(r.imag))))
    return worst
d8 = {}
for J2, J4 in [(0.5, 0.0), (0.0, 0.5), (0.0, -0.5), (0.5, 0.3), (0.5, -0.3), (0.3, 0.6), (0.2, -0.6),
               (0.45, -0.6), (0.45, -0.7), (0.15, 1.2), (0.14, 6.0), (0.13, 6.0)]:
    d8[f"J2={J2},J4={J4}"] = dict(line_test_max_imag=multivariate_line_test(4, terms_of(J2, J4)),
                                  diag_max_abs_log=float(np.max(np.abs(np.log(diag_moduli(4, terms_of(J2, J4)))))))
out["D8_points"] = d8
out["D8_pass"] = (d8["J2=0.5,J4=0.0"]["line_test_max_imag"] < 1e-6 and
                  d8["J2=0.0,J4=0.5"]["line_test_max_imag"] > 1e-3 and
                  d8["J2=0.45,J4=-0.7"]["line_test_max_imag"] > 1e-3)

# ---------- D9 large-J4 threshold ----------
# J4 -> +inf keeps only even #minus: Z ~ e^{6J2}(1+z^4) + 6 e^{-2J2} z^2, on |z|=1 iff 6 e^{-8 J2} <= 2.
thr = math.log(3) / 8
d9 = {}
for J2 in [thr - 0.01, thr + 0.01]:
    m = diag_moduli(4, terms_of(J2, 14.0))
    d9[f"J2={J2:.4f}"] = float(np.max(np.abs(np.log(m))))
out["D9_threshold_ln3_over_8"] = thr
out["D9_rows"] = d9
out["D9_pass"] = d9[f"J2={thr-0.01:.4f}"] > 1e-4 and d9[f"J2={thr+0.01:.4f}"] < 1e-7

# ---------- D10 many-body only: fails hot, holds cold (symmetric, exact by Grace-Walsh-Szego) ----------
from math import comb
def e_d(N, kminus, dd):
    return sum(comb(N - kminus, dd - j) * comb(kminus, j) * (-1) ** j for j in range(dd + 1))
def sym_max_abs_log(N, J2, J4):
    coef = [comb(N, kk) * math.exp(J2 * e_d(N, kk, 2) + J4 * e_d(N, kk, 4)) for kk in range(N + 1)]
    return float(np.max(np.abs(np.log(np.abs(np.roots(coef[::-1]))))))
d10 = {}
for N in [5, 6, 8, 10]:
    d10[N] = {str(J): sym_max_abs_log(N, 0.0, J) < 1e-7 for J in [0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0]}
out["D10_all_4_subsets_no_pairs_LY"] = d10
out["D10_pass"] = all(not d10[N]["0.001"] and d10[N]["2.0"] for N in d10)

out["all_pass"] = all(out[k] for k in out if k.endswith("_pass") or k == "D1_rayleigh_identity_all_equal")
json.dump(out, open("depth.json", "w"), indent=1)
for k, v in out.items():
    if k.endswith("pass") or k.startswith("D1_") or k in ("D2_line_test_max_imag_part", "D3_witnesses", "D7_moduli_J4_12", "D8_points", "D9_rows", "D10_all_4_subsets_no_pairs_LY"):
        print(k, v)
print("ALL PASS" if out["all_pass"] else "SOME FAIL")
