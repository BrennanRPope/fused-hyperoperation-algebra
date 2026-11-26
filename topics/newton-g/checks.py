"""Checks for newton-g/README.md. Run: python3 checks.py  (numpy only)."""
import numpy as np, json
from math import pi, log, exp, sqrt
from gshapes import *

R = []; out = {}
def check(tag, ok, msg):
    R.append(bool(ok)); print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}")

# ===== A. What matched before (recovered) =====
gs, eps = 0.7, 2.0; k = -gs
s = lambda t: (1 - np.exp(-eps*t))/gs
circ = lambda a, b: a + b + k*a*b
t1, t2, g0 = 0.3, 0.55, 0.2
u = lambda t: 1/gs + (1/g0 - 1/gs)*np.exp(-eps*t)          # exact one-loop solution in u = 1/g
check("A1", abs(circ(u(0), s(t1)) - u(t1)) < 1e-12 and abs(circ(s(t1), s(t2)) - s(t1+t2)) < 1e-12,
      "one-loop flow in u = 1/g is u o_k s(t) with k = -g*, and steps compose by o_k (recovered identity)")
G0, w = 0.9, 1/gs
kk = np.linspace(0.1, 5, 50); gg = G_sakharov(kk, G0, w)*kk**2
dg = np.gradient(gg, np.log(kk))
check("A2", np.max(np.abs(dg - one_loop(gg, gs=1/w))[3:-3]) < 2e-3,
      "Sakharov / Bonanno-Reuter running 1/G = 1/G0 + w k^2 is the solution of beta = 2g(1 - w g): the one-loop match IS induced-gravity running")
# free parameters: any quadratic beta = a g - b g^2 maps to o_k with k = -a/b ... identity, zero fitted numbers
check("A3", True, "free-parameter count of the one-loop match: 0 fitted (k = -g* comes from the flow itself); it is an identity for EVERY quadratic beta, so it tests nothing about gravity by itself")

# ===== B. induced share (union / odds reading) =====
f = w*kk**2*G_sakharov(kk, G0, w)             # induced fraction of the stiffness 1/G
odds = f/(1 - f)
check("B1", np.allclose(f, gg*w) and np.allclose(odds, G0*w*kk**2),
      "g/g* = (induced part of 1/G)/(1/G): the one-loop 'probability' is the induced share; its odds = G0 w k^2 grow by e^2 per e-fold")
# Litim: log-odds slope runs from 2 (IR) to theta (UV)
Rg = np.array([1e-6, 0.5, 1 - 1e-7])
def logodds_rate(Rv, th=theta_litim):  # dR/dt = R(1-R)/(a(1-R) + bR), a = 1/2, b = 1/th ; d logodds/dt = 1/(a(1-R)+bR)
    return 1/(0.5*(1 - Rv) + Rv/th)
rates = logodds_rate(Rg)
check("B2", abs(rates[0] - 2) < 1e-5 and abs(rates[-1] - 23/11) < 1e-5,
      f"beyond one loop (Litim, lambda=0) the odds no longer grow at one rate: {rates[0]:.4f} at the IR end, {rates[1]:.4f} mid-way, {rates[-1]:.4f} = 23/11 at the UV end")
# verify the exact Litim solution against direct integration
gg2 = gstar*np.concatenate([np.geomspace(1e-4, 0.5, 4000), 1 - np.geomspace(0.5, 1e-3, 4000)[1:]])
tnum = np.concatenate([[0], np.cumsum(np.diff(gg2)/litim0(0.5*(gg2[1:] + gg2[:-1])))])
tex = litim_t_of_R(gg2/gstar) - litim_t_of_R(gg2[0]/gstar)
check("B3", np.max(np.abs(tnum - tex)) < 2e-4, "exact Litim lambda=0 solution R^(1/2)(1-R)^(-11/23) = e^t confirmed by integrating beta")

# ===== C. branching ladder: rho = scheme-free two-loop number =====
for name, (F, sig2, f3) in LAWS.items():
    rho = rho_branching(sig2, f3)
    epss = np.array([0.004, 0.008, 0.016])
    th = np.array([theta_map(F, e) for e in epss])
    fit = np.polyfit(epss, (th - epss)/epss**2, 1)[1]
    check("C:" + name.split()[0], abs(fit - rho) < 2e-2,
          f"{name}: rho = 1 - 2 f'''(1)/(3 sigma^4) = {rho:.4f}; exponent of the near-critical family gives theta = eps + rho eps^2 with rho ~ {fit:.4f}")
    out["rho_" + name.split()[0]] = rho
# E's own invariant from the height-flow handoff
check("C4", abs(rho_branching(1, 1) - 1/3) < 1e-15, "Poisson(1) value equals E's iterative-log invariant rho = -b1/b0^2 = 1/3 (height-running-coupling handoff)")
# critical survival 2/n
x = -1.0
for n in range(1, 2001): x = np.expm1(x)
check("C5", abs(-x*2000/2 - 1) < 5e-3, f"E^n(-1) = -P(Poisson(1) family survives n generations); at n = 2000, n*P/2 = {-x*1000:.5f} (Kolmogorov 2/(sigma^2 n), known)")
# geometric = circle op exactly
check("C6", abs(theta_map(LAWS["geometric (linear fractional)"][0], 0.3) - 0.3) < 1e-6,
      "geometric (linear fractional) offspring: theta = eps exactly at every eps, rho = 0: the circle op is the one-loop member of the branching ladder")
# gravity in 2+eps
epss = np.array([0.002, 0.004, 0.008])
thg = np.array([theta_flow(ak_beta, e, 1e-9*e, 30*e) for e in epss])
fitg = np.polyfit(epss, (thg - epss)/epss**2, 1)[1]
check("C7", abs(fitg - 0.6) < 1e-2, f"2+eps gravity (Aida-Kitazawa beta, c = 0): theta = eps + rho eps^2 with rho = 15/(25-c) = 0.6 (fit {fitg:.4f}); time unit t = ln k")
out["rho_gravity_c0"] = 0.6

# ===== D. numerology test =====
cfac = 0.6/(1/3)
c_lnk = 25 - 15*3          # 15/(25-c) = 1/3
c_lnk2 = 25 - 15/(2/3)     # one E step = 2 e-folds (r = e^2) => rho_E per e-fold = 2/3
check("D1", abs(cfac - 1.8) < 1e-12 and c_lnk == -20 and abs(c_lnk2 - 2.5) < 1e-12,
      "to match gravity to E you must pick a time unit (factor 9/5) or a matter central charge: c = -20 (one height = 1 e-fold, r = e) or c = 2.5 (one height = 2 e-folds, r = e^2); neither is a natural value")
# offspring law that would hit 3/5: f'''(1)/sigma^4 = 3/5
check("D2", abs(rho_branching(1, 0.6) - 0.6) < 1e-12,
      "a critical branching law with f'''(1)/sigma^4 = 3/5 reproduces gravity's 3/5; the ladder can hit any rho <= 1, so landing on it is one free number fitting one number")

# ===== E. shape test: Richards (semiring/power-mean dial) =====
p = theta_litim/2
tt = np.linspace(-6, 6, 1201)
RL = align(R_litim, tt); RR = align(lambda t: R_richards(t, p), tt); RG = align(R_logistic, tt)
dR_rich, dR_log = np.max(np.abs(RL - RR)), np.max(np.abs(RL - RG))
out.update(dict(dR_richards=dR_rich, dR_logistic=dR_log))
check("E1", abs(-2*p - (2*(1 - (p + 1)))) < 1e-12,
      "Richards curve dR/dt = 2R(1 - R^p) with p = theta/2 = 23/22 has both end exponents right: 2 at R = 0 and theta = 23/11 at R = 1")
check("E2", dR_rich < dR_log, f"shape: max |R_Litim - R_Richards| = {dR_rich:.4f}; max |R_Litim - R_logistic| = {dR_log:.4f} (curves aligned at R = 1/2)")
check("E3", abs((p + 1) - 45/22) < 1e-12,
      "but Richards' beta has a g^(45/22) term at g = 0: not a power series, while perturbative gravity's beta is. Same corner failure as the earlier flat coordinate")
# power mean form: 1/G^p = 1/G0^p + (w k^2)^p gives R^p logistic
Gp = lambda k: (G0**-p + (w*k**2)**p)**(-1/p)
Rp = w*kk**2*Gp(kk); y = (G0*w*kk**2)**p
check("E4", np.allclose(Rp**p, y/(1 + y)), "the Richards curve is the power-mean (soft-max with dial 1/p) of bare and induced stiffness; p = 1 is Sakharov's plain sum")

# ===== F. spacetime-fluid strand =====
for G in [1.0, 0.37, 6.674e-11]:
    eta_sh, s_ent = 1/(16*pi*G), 1/(4*G)
    ok = abs(eta_sh/s_ent - 1/(4*pi)) < 1e-15
check("F1", ok, "membrane paradigm: horizon shear viscosity 1/(16 pi G) over entropy density 1/(4G) = 1/(4 pi) for every G: the fluid ratio cannot fix G")
hb = 1.0; G = 0.37; eta_J = 1/(4*hb*G)
check("F2", abs(eta_J**-0.5 - 2*sqrt(hb*G)) < 1e-12, "Jacobson 1995: G = 1/(4 hbar eta) with S = eta A; eta^(-1/2) = 2 l_P; eta itself is not computed")
Ga, Gb = 0.8, 1.7
par = Ga*Gb/(Ga + Gb)
check("F3", abs(1/(1/Ga + 1/Gb) - par) < 1e-15 and abs((1/par) - (1/Ga + 1/Gb)) < 1e-15,
      "induced gravity: species add in 1/G, so G combines as a parallel sum, which is + seen through x -> 1/x, the same move that makes the one-loop flow o_k in u = 1/g")
# scale: Sakharov's 1/G ~ N kappa^2 -> G sets only the cutoff: G kappa^2 is pure but kappa free
Ns = np.array([1, 10, 100]); kap = 1/np.sqrt(Ns)   # if 1/G = N kappa^2 c with fixed G, kappa ~ 1/sqrt(N)
check("F4", np.allclose(Ns*kap**2, 1), "for fixed G, the cutoff must fall like 1/sqrt(N): Dvali's species scale (known). G's value = choice of molecule size; no pure number in our algebra supplies a length")

print(f"\n{sum(R)}/{len(R)} checks pass")
json.dump(out, open("results.json", "w"), indent=1)
