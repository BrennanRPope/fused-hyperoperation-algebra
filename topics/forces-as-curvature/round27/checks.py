"""DRAFT round 27: hunting for a gravity prediction from the law. Each check classifies one candidate:
does the law force the number or is it chosen, is it known or measured, and what data decide it.
Data values are typed in below with their sources (see README); nothing is fetched at run time."""
import json, math, io, sys, time
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}; FIG = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ---------------- data (typed in; sources in README) ----------------
H0 = 67.4                      # km/s/Mpc, Planck 2018 base LCDM
Om, OL, Orad = 0.315, 0.685, 9.1e-5
OK_DESI = (0.0023, 0.0011)     # DESI DR2 BAO + CMB (o-LCDM), positive = open
OK_PLBAO = (0.001, 0.002)      # Planck 2018 + BAO
c = 299792.458                 # km/s
Mpc = 3.0857e22                # m
hbar_c = 3.1615e-26            # J m
g2 = 0.42                      # SU(2) coupling squared at m_Z (g ~ 0.65)
G = 6.674e-11; c_si = 2.998e8
rho_crit = 3 * (H0 * 1e3 / Mpc)**2 / (8 * math.pi * G) * c_si**2   # J/m^3
R_H = c / H0 * Mpc             # Hubble radius in m

# ================= G1. naive identification: today's slices are the law's S^3 at curvature q = Lambda/3
say("== G1. If q = Lambda/3 also set today's spatial curvature")
Ok_naive = -OL                       # Omega_k = -k c^2/(a0 H0)^2 with 1/a0^2 = Lambda/3 -> -Omega_Lambda
s_desi = (OK_DESI[0] - Ok_naive) / OK_DESI[1]; s_pl = (OK_PLBAO[0] - Ok_naive) / OK_PLBAO[1]
check("G1 (exact + data) the law ties its dial to both Lambda = 3q and the slices' curvature q. If that held today, Omega_k = -Lambda/(3 H0^2) = -Omega_Lambda = -0.685. Measured curvature is a thousand times smaller: the identification is ruled out, so it can hold at most at a moment of rest",
      s_desi > 100 and s_pl > 100, f"law (naive) Omega_k = {Ok_naive:.3f}; DESI DR2+CMB {OK_DESI[0]} +/- {OK_DESI[1]} -> {s_desi:.0f} sigma away; Planck+BAO {OK_PLBAO[0]} +/- {OK_PLBAO[1]} -> {s_pl:.0f} sigma")
results["G1"] = dict(Ok_naive=Ok_naive, sigma_DESI=s_desi, sigma_PlanckBAO=s_pl)

# ================= G2. did our universe ever have the law's moment of rest?
say("== G2. A moment of rest in our past?")
def E2(x, Ok):     # x = 1/a = 1 + z, Lambda fixed by closure
    OLx = 1 - Om - Orad - Ok
    return Orad * x**4 + Om * x**3 + Ok * x**2 + OLx
def min_E2(Ok, xmax=50.0):
    xs = np.linspace(1, xmax, 20001); return float(E2(xs, Ok).min())
Ok_thr = brentq(lambda k: min_E2(k), -3.0, -0.01)      # most open-to-closed value where E^2 first touches 0 for some z >= 0
xs = np.linspace(1, 50, 20001); x_rest = float(xs[np.argmin(E2(xs, Ok_thr))])
no_rest_desi = min_E2(OK_DESI[0] - 5 * OK_DESI[1]) > 0 and min_E2(OK_PLBAO[0] - 5 * OK_PLBAO[1]) > 0
check("G2 (numeric, textbook Friedmann) rounds 16, 23 and 25 tie the law's q to the slices' curvature only at a moment of rest (H = 0). With Planck's matter and Lambda, a moment of rest anywhere in our past (z >= 0) needs Omega_k below the threshold found here; the measured curvature is far above it even 5 sigma toward closed, so our universe never had one: its past is a big bang, not the law's bounce",
      no_rest_desi and Ok_thr < -0.5, f"threshold Omega_k = {Ok_thr:.3f} (rest at 1+z = {x_rest:.2f}); 5-sigma closed edges {OK_DESI[0]-5*OK_DESI[1]:.4f}, {OK_PLBAO[0]-5*OK_PLBAO[1]:.4f} give min H^2/H0^2 = {min_E2(OK_DESI[0]-5*OK_DESI[1]):.3f}, {min_E2(OK_PLBAO[0]-5*OK_PLBAO[1]):.3f}")
results["G2"] = dict(Ok_threshold=Ok_thr, rest_at_1pz=x_rest)
zz = np.linspace(0, 6, 400)
FIG["g2"] = dict(z=zz.tolist(), E2={str(k): E2(1 + zz, k).tolist() for k in (OK_DESI[0], -0.3, round(Ok_thr, 4), -1.2)})

# ================= G3. the weak version: closed slices of any size (sign only)
say("== G3. Sign only: are the slices closed?")
p_closed_desi = float(norm.cdf((0 - OK_DESI[0]) / OK_DESI[1])); p_closed_pl = float(norm.cdf((0 - OK_PLBAO[0]) / OK_PLBAO[1]))
check("G3 (data) if the cosmic rest slices must be the law's 3-spheres (q > 0, as Lambda > 0 says), of whatever size, then Omega_k < 0: closed. This is the one sign the law leans toward. It is chosen, not forced: de Sitter can be sliced flat or open too, and nothing in the law makes matter's rest frame its slicing. DESI DR2 + CMB lean open at about 2 sigma; Planck + BAO are consistent with either",
      0 < p_closed_desi < 0.05 and p_closed_pl > 0.2, f"P(Omega_k < 0): DESI DR2+CMB {p_closed_desi:.3f} ({OK_DESI[0]/OK_DESI[1]:.1f} sigma open); Planck+BAO {p_closed_pl:.2f}")
results["G3"] = dict(P_closed_DESI=p_closed_desi, P_closed_PlanckBAO=p_closed_pl, sigma_open_DESI=OK_DESI[0] / OK_DESI[1])

# ================= G4. the locked weak-force field as cosmic radiation
say("== G4. The law's f = 1 twist filling space: how much energy?")
Ok_closed_bound = OK_PLBAO[0] - 2 * OK_PLBAO[1]                  # most closed allowed at 2 sigma (Planck+BAO)
a0_min = R_H / math.sqrt(abs(Ok_closed_bound))                   # curvature radius today
rho_YM = 3 * 1.0 / (2 * g2) * hbar_c / a0_min**4                  # 8 pi G rho/3 = kappa C/a^4, kappa = 4 pi G/g^2, C = 1 at f = 1
Omega_YM = rho_YM / rho_crit
check("G4 (exact formula from round 25 + data) the law's own f = 1 field locked to a closed universe of curvature radius a0 carries rho = 3/(2 g^2) hbar c/a0^4 (round 25: 8 pi G rho/3 = kappa C/a^4, C = 1). With a0 at least 18 Hubble radii (Planck + BAO, 2 sigma), that is a fixed number but about 125 orders of magnitude below the critical density; and the real weak force is screened by the Higgs beyond 1e-18 m, so such a field could not span space anyway",
      Omega_YM < 1e-100, f"a0 >= {a0_min/R_H:.1f} Hubble radii = {a0_min:.2e} m; rho <= {rho_YM:.1e} J/m^3; Omega <= {Omega_YM:.1e} (critical density {rho_crit:.2e} J/m^3)")
results["G4"] = dict(a0_min_hubble=a0_min / R_H, rho_YM=rho_YM, Omega_YM=Omega_YM)

# ================= G5. the n = 3 wave and its (3, -6) couplings
say("== G5. The law's lowest gravitational wave today")
lam_min = 2 * math.pi * a0_min / math.sqrt(8)                    # physical wavelength of the n = 3 mode
check("G5 (exact + data) round 26's lowest wave has wavelength 2 pi a0/sqrt(8). With a0 >= 18 Hubble radii it is longer than about 40 Hubble radii, far beyond what any telescope or detector can see; its (3, -6) coupling to the twist needs a cosmic weak-force field that the Higgs forbids. Fixed by the law, but untestable",
      lam_min / R_H > 10, f"wavelength >= {lam_min/R_H:.0f} Hubble radii (observable universe radius about 3.3)")
results["G5"] = dict(lambda_min_hubble=lam_min / R_H)

# ================= G6. the ledger
say("== G6. The ledger")
ledger = [
 dict(c="Lambda = 3q and slices' curvature q today", forced="yes, if taken literally", known="measured", verdict="ruled out (G1)", test="done: Omega_k"),
 dict(c="Moment of rest in our past (law's bounce)", forced="chosen", known="measured", verdict="ruled out (G2)", test="done: Omega_k, CMB"),
 dict(c="Closed slices, any size (Omega_k < 0)", forced="chosen", known="open question", verdict="open; data lean open at 2 sigma (G3)", test="DESI final, Euclid, CMB-S4"),
 dict(c="Size of Lambda (q)", forced="no: q is free", known="measured", verdict="not a prediction", test="none"),
 dict(c="Sign of Lambda > 0", forced="yes inside the law (q > 0 closes)", known="measured 1998", verdict="agrees; textbook-equivalent", test="done"),
 dict(c="Energy of the locked f = 1 twist", forced="yes, given the field", known="no", verdict="fixed but ~1e-125 of critical (G4)", test="none possible"),
 dict(c="n = 3 wave, couplings (3, -6)", forced="yes, given the field", known="no", verdict="fixed but beyond the horizon (G5)", test="none possible"),
 dict(c="Frame dragging 2(1+gamma)", forced="gamma = 1 put in", known="measured (GP-B)", verdict="agrees; GR", test="done"),
 dict(c="Fused expansion laws without Lambda", forced="chosen", known="measured", verdict="ruled out (dark-sector page)", test="done: SNe + DESI"),
 dict(c="Mult-add gravity for galaxies", forced="chosen", known="measured", verdict="ruled out: 1 to 5% of the gap (dark-sector page)", test="done: SPARC"),
]
n_open = sum(1 for r in ledger if r["verdict"].startswith("open"))
check("G6 (bookkeeping) every candidate gets: forced or chosen, known or measured, and what decides it. Exactly one is open and decidable (closed slices). The only one that a literal reading of the law forces and that differs from textbook physics (G1) is ruled out; the rest are chosen, textbook, or fixed but untestable",
      n_open == 1 and len(ledger) == 10, f"{len(ledger)} candidates; open and testable: {n_open}")
results["ledger"] = ledger

say(f"\n{passed} passed, {failed} failed  ({time.time() - T0:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed)
results["data"] = dict(H0=H0, Om=Om, OL=OL, Orad=Orad, OK_DESI=OK_DESI, OK_PlanckBAO=OK_PLBAO, g2=g2)
results["figure"] = FIG
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
