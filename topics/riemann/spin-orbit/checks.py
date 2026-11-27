"""Checks for spin-orbit: read spinorbit.json and ../known-tools/rewind_snaps.npz."""
import json, numpy as np
R = json.load(open("spinorbit.json")); res = []
def ck(n, ok, d=""):
    res.append(ok); print(("PASS " if ok else "FAIL ") + n + ("  " + d if d else ""))
L = R["lean"]
ck("every prime 2..31 leans toward 180 deg at the zeros (mean cos < 0)", all(v["mean_cos"] < 0 for v in L.values()))
ck("the lean matches Landau's formula within 0.001 for every prime 2..31",
   all(abs(v["mean_cos"] - v["predicted"]) < 1e-3 for v in L.values()))
ck("the lean is largest near p = 7 (log p / sqrt p peaks at e^2)",
   min(L, key=lambda p: L[p]["mean_cos"]) in ("5", "7", "11"))
ck("no pair of primes locks: |mean cos(gamma log p/q)| < 0.0005 for all six pairs",
   all(abs(v["mean_cos"]) < 5e-4 for v in R["pairs"].values()))
ck("random heights show no lean for p = 2 (|mean| < 3 standard errors)",
   abs(R["control_random_heights_mean_cos2"]) < 3 * L["2"]["se"])
w = R["lean2_windows"]
ck("the p = 2 lean fades with height (first window more negative than last)", w[0]["mean"] < w[-1]["mean"])
ck("the p = 2 lean in every window is within 2 standard errors of -(log 2/sqrt 2)/log(gamma/2pi)",
   all(abs(v["mean"] - v["pred"]) < 2 * v["se"] for v in w))
for c in R["capture"]:
    ck(f"inserted pair y0 = {c['y0']} at gamma {c['x0']/2:.0f}: captured, formula within 10%, inside de Bruijn's bound",
       abs(c["t_capture_nbody"] / c["t_formula"] - 1) < 0.10 and c["t_capture_nbody"] < c["debruijn_bound"],
       f"{c['t_capture_nbody']:.3f} vs {c['t_formula']:.3f}, bound {c['debruijn_bound']:.3f}")
K = np.load("../known-tools/rewind_snaps.npz"); x = K["S"] * K["d0"]**2 / 2
ck("both rewind families are populated (each between 30% and 70% of pairs)", 0.3 < (x < 1).mean() < 0.7, f"{(x<1).mean():.3f}")

# ---- squeeze (LAMBDA.md) ----
Q = json.load(open("squeeze.json")); B = Q["backward"]
tc_ = Q["tcap_at_3e12_beta1"]
ck("S1 forward: capture of a beta = 1 zero at 3e12 is positive and below Platt-Trudgian's 0.2 (heuristic)",
   0 < tc_ < 0.2 and all(0 < v < 0.2 for v in Q["tcap_at_3e12_beta1_S_range"]), f"{tc_:.4f}")
tau = np.array(B["taus"]); fr = np.array(B["broken_fraction"]); pr = np.array(B["gue_small_gap_prediction"])
m = (fr > 0.003) & (fr < 0.1)
slope = np.polyfit(np.log(tau[m]), np.log(fr[m]), 1)[0]
ck("S2 backward: broken share grows like tau^1.5 (GUE small-gap law), slope within 0.3", abs(slope - 1.5) < 0.3, f"slope {slope:.2f}")
ratio = np.median(fr[m] / pr[m])
ck("S3 backward: broken share is about half the GUE prediction (ratio 0.3 to 0.8)", 0.3 < ratio < 0.8, f"ratio {ratio:.2f}")
a, b = Q["broken_fraction_tau1e-3_at_1e3_1e22"]
ck("S4 squeeze: at fixed tau = 0.001 the broken share grows with height", b > 100 * a, f"{a:.1e} -> {b:.3f}")
print(f"\n{sum(res)} of {len(res)} checks pass")
