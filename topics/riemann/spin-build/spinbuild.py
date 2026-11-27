"""Spin-build: try to build zeta's Lee-Yang spin out of simple Lee-Yang spins joined by ferromagnetic links.  DRAFT, Nov 27 2026.

Lee-Yang (Newman's sense): a symmetric spin distribution mu is Lee-Yang when Z(h) = int e^{h u} mu(du) vanishes only at
imaginary h, i.e. its Fourier transform int cos(z u) mu(du) has only real zeros.  zeta's spin: mu = Phi(u) du, RH <=> Lee-Yang.
The Lee-Yang theorem (Lee-Yang 1952; Newman 1974) needs: (1) symmetric Lee-Yang single spins, (2) ferromagnetic couplings,
(3) the field coupled linearly to the total, and (4) a measure that exists (finite total weight, or a limit of finite ones).
"""
import json, math
import numpy as np
import mpmath as mp
from scipy.optimize import brentq

mp.mp.dps = 30
out = {}

# ---------------------------------------------------------------- zeros of H_0 (z = 2 gamma)
NZ = 40
gam = [float(mp.zetazero(k).imag) for k in range(1, NZ + 1)]
Z = np.array([2 * g for g in gam])
out["first_zeros_z"] = Z[:6].tolist()

# ---------------------------------------------------------------- A. independent simple spins: the 3x rule
# Ising spin +-a: cos(a z), zeros (2m+1) pi / (2a)  -> if x is a zero then 3x, 5x, ... are zeros.
# Uniform spin on [-a, a]: sin(a z)/(a z), zeros m pi / a -> x, 2x, 3x, ...
# Gaussian: no zeros.  A product of such factors has a zero set that is a union of these progressions,
# so EVERY zero x needs 3x (Ising factor) or 2x (uniform factor) to be a zero as well.
def nearest(x):
    i = np.argmin(np.abs(Z - x)); return float(Z[i]), float(abs(Z[i] - x))
rule = []
for x in Z[:8]:
    n2, d2 = nearest(2 * x); n3, d3 = nearest(3 * x)
    gap = 2 * math.pi / math.log(x / (4 * math.pi)) * 2      # mean spacing of z = 2 gamma near 3x (approx)
    rule.append({"x": float(x), "2x": float(2 * x), "nearest_to_2x": n2, "miss_2x": d2, "3x": float(3 * x),
                 "nearest_to_3x": n3, "miss_3x": d3})
out["three_x_rule"] = rule
out["three_x_rule_fails_at_first_zero"] = bool(rule[0]["miss_3x"] > 0.5 and rule[0]["miss_2x"] > 0.5)

# ---------------------------------------------------------------- B. a ferromagnetic link escapes the rule
# two Ising spins, amplitudes a > b, coupling J: Z(z) = e^J cos((a+b) z) + e^{-J} cos((a-b) z)
a, b = 1.0, 0.37
def zeros_two_spin(J, zmax=40.0):
    f = lambda z: math.exp(J) * math.cos((a + b) * z) + math.exp(-J) * math.cos((a - b) * z)
    zs = np.linspace(0, zmax, 40001); fv = np.array([f(z) for z in zs])
    real_roots = [brentq(f, zs[i], zs[i + 1]) for i in range(len(zs) - 1) if fv[i] * fv[i + 1] < 0]
    # count all zeros in |Re z| < zmax via argument principle proxy: compare to expected count (a+b) zmax / pi
    expected = (a + b) * zmax / math.pi
    return real_roots, expected
two = {}
for J in [-0.4, -0.1, 0.0, 0.1, 0.5]:
    rr, ex = zeros_two_spin(J)
    two[str(J)] = {"real_zeros": len(rr), "expected_all": ex}
out["two_spin"] = two
out["two_spin_note"] = "all zeros real (count = (a+b) zmax/pi) iff J >= 0; J < 0 loses real zeros in complex pairs"
rr, _ = zeros_two_spin(0.5)
out["two_spin_J0.5_first_zeros"] = rr[:4]
out["two_spin_J0.5_breaks_3x"] = bool(min(abs(np.array(rr) - 3 * rr[0])) > 0.05)

# ---------------------------------------------------------------- C. zeta's spin is a MIXTURE over n, of a non-positive bump
def phi1(u):
    return (2 * np.pi**2 * np.exp(9 * u) - 3 * np.pi * np.exp(5 * u)) * np.exp(-np.pi * np.exp(4 * u))
def Phi(u):
    s = 0
    for n in range(1, 60):
        s = s + (2 * np.pi**2 * n**4 * np.exp(9 * u) - 3 * np.pi * n**2 * np.exp(5 * u)) * np.exp(-np.pi * n**2 * np.exp(4 * u))
    return s
uu = np.linspace(-1.2, 1.2, 961)
mix = sum(n**-0.5 * phi1(uu + 0.5 * math.log(n)) for n in range(1, 60))
out["mixture_identity_maxerr"] = float(np.max(np.abs(mix - Phi(uu))))
out["Phi_even_err"] = float(np.max(np.abs(Phi(uu) - Phi(-uu))))
out["phi1_min"] = float(phi1(uu).min()); out["phi1_negative_below_u"] = float(uu[np.argmax(phi1(uu) > 0)])
out["Phi_min"] = float(Phi(uu).min())
np.savez("mixture.npz", u=uu, Phi=Phi(uu), phi1=phi1(uu),
         **{f"term{n}": n**-0.5 * phi1(uu + 0.5 * math.log(n)) for n in [1, 2, 3, 4, 5]})

# ---------------------------------------------------------------- D. the symmetric primon gas: exact, symmetric, free, field-linear
# prime spin X_p = k ln p, k in Z, weight p^{-sigma |k|}  (symmetric, discrete Laplace).  Field h couples to X = sum X_p
# = ln(a/b) with a, b coprime; total weight (ab)^{-sigma}.  Z(h) = sum_{(a,b)=1} a^{-(sigma-h)} b^{-(sigma+h)}
#                                                              = zeta(sigma-h) zeta(sigma+h) / zeta(2 sigma)   (|Re h| < sigma - 1)
x_, y_ = 2.5, 3.1
lhs = mp.mpf(0)
for A in range(1, 1500):
    for B in range(1, 1500 // A + 400):
        if math.gcd(A, B) == 1:
            lhs += mp.mpf(A)**(-x_) * mp.mpf(B)**(-y_)
rhs = mp.zeta(x_) * mp.zeta(y_) / mp.zeta(x_ + y_)
out["coprime_identity_relerr"] = float(abs(lhs - rhs) / rhs)
# each prime spin: Z_p(h) = (1 - q^2) / ((1 - q p^h)(1 - q p^-h)), q = p^-sigma: no zeros, poles on Re h = +-sigma
# -> Lee-Yang vacuously.  Completed version: zeros of xi(1/2 + h) xi(1/2 - h) = xi(1/2+h)^2 at h = rho - 1/2:
#    imaginary iff Re rho = 1/2.  So RH <=> the symmetric primon gas at sigma = 1/2 (completed) is Lee-Yang.
# But its total weight sum (ab)^{-sigma} is finite only for sigma > 1: the measure does not exist at sigma = 1/2.
tot = {}
for sg in [1.5, 1.2, 1.05]:
    tot[str(sg)] = float(mp.zeta(sg)**2 / mp.zeta(2 * sg))
out["total_weight_zeta2_over_zeta2s"] = tot
out["total_weight_diverges_at_sigma_1"] = "zeta(sigma)^2/zeta(2 sigma) has a double pole at sigma = 1"

# finite-volume versions (primes <= P) at sigma = 1/2 along the imaginary field h = i t: zero-free for every P
primes = [p for p in range(2, 20000) if all(p % q for q in range(2, int(p**0.5) + 1))]
ts = np.linspace(10, 50, 1601)
def ZP(P, t, sg=0.5):
    lz = np.zeros_like(ts)
    for p in primes:
        if p > P: break
        q = p**-sg
        num = 1 - q * q
        lz += np.log(num) - np.log(np.abs(1 - q * np.exp(1j * t * np.log(p)))) - np.log(np.abs(1 - q * np.exp(-1j * t * np.log(p))))
    return lz       # log Z_P(it) (real, positive Z)
logZ = {P: ZP(P, ts) for P in [10, 100, 1000, 10000]}
zeta2 = np.array([float(mp.log(abs(mp.zeta(mp.mpf(0.5) + 1j * t))**2)) for t in ts])
out["finite_volume"] = {str(P): {"min_logZ": float(v.min()), "corr_with_log_zeta2": float(np.corrcoef(v, zeta2)[0, 1])}
                        for P, v in logZ.items()}
out["finite_volume_note"] = "Z_P(it) = prod |1 - p^{-1/2 - it}|^{-2} times constants: positive (no zeros) for every P"
np.savez("primon.npz", t=ts, zeta2=zeta2, **{f"logZ_{P}": v for P, v in logZ.items()}, gam=np.array(gam))

# ---------------------------------------------------------------- E. the four-condition table
out["conditions"] = {
    "newman_spin_Phi": {"symmetric": "yes", "ferromagnetic": "unknown (= RH)", "field_linear": "yes", "measure_exists": "yes"},
    "symmetric_primon_gas_sigma_half": {"symmetric": "yes", "ferromagnetic": "yes (free, J = 0)", "field_linear": "yes",
                                        "measure_exists": "NO (only for sigma > 1; sigma = 1 is its phase transition)"},
    "knauf_chain": {"symmetric": "yes (Ising spins)", "ferromagnetic": "yes (proved, Knauf 1993)",
                    "field_linear": "NO (s is inverse temperature, not field)", "measure_exists": "yes at finite length"},
    "independent_ising_uniform": {"symmetric": "yes", "ferromagnetic": "yes (J = 0)", "field_linear": "yes",
                                  "measure_exists": "yes", "fails": "3x / 2x rule at the first zero"}}
json.dump(out, open("spinbuild.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
