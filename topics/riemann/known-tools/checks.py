"""Checks for known-tools: read rewind.json, euler_dial.json, parents_lambda.json."""
import json
R = json.load(open("rewind.json")); E = json.load(open("euler_dial.json")); P = json.load(open("parents_lambda.json"))
res = []
def ck(n, ok, d=""):
    res.append(ok); print(("PASS " if ok else "FAIL ") + n + ("  " + d if d else ""))
nb = R["nbody_closest"]
ck("N-body rewind of the closest pair matches the pair formula within 2%",
   abs(nb["t_collision"] / nb["formula"] - 1) < 0.02, f"{nb['t_collision']:.3e} vs {nb['formula']:.3e}")
tf = R["low_formula"]["pairs"][3]["t_c"]; td = R["low_direct"]["t_first_nonreal"]
ck("low pair: formula collision lies in the direct rewind's step (td, td + 0.25]", td < tf <= td + 0.25, f"{tf:.2f}, {td}")
b = R["pairs"]["top10"][0]
ck("closest pair gives a lower bound between -1e-3 and 0", -1e-3 < b["t_c"] < 0, f"{b['t_c']:.2e} at gamma {b['gamma'][0]:.2f}")
F = R["forward"]; ks = ["0.0", "0.01", "0.05", "0.2"]
ck("forward flow: gap variance falls", all(F[a]["var"] > F[c]["var"] for a, c in zip(ks, ks[1:])))
ck("forward flow: smallest gap grows", all(F[a]["min"] < F[c]["min"] for a, c in zip(ks, ks[1:])))
ck("forward flow: no gap below half the mean at t = 0.2", F["0.2"]["frac_below_half"] == 0)
ck("M(1e8) = 1928 (known value)", E["mertens"]["M_at_1e8"] == 1928)
ck("|M(x)|/sqrt x < 1 for 1e3 <= x <= 1e8", E["mertens"]["max_abs_M_over_sqrt_x_from_1e3"] < 1)
ck("squarefree error / x^(1/4) stays below 1.5 from 1e3 to 1e8", E["squarefree"]["max_abs_err_over_x_quarter_above_1e3"] < 1.5)
L = E["second_layer"]
ck("pi(x) < li(x) on the whole sample 1e4..1e8", not L["A_ever_positive_1e4_1e8"])
ck("second layer: log-mean moves from below -1 to within 0.3 of 0 when 1/2 li(sqrt x) is added back",
   L["mean_A_logweighted_1e4_1e8"] < -1 and abs(L["mean_B_logweighted_1e4_1e8"]) < 0.3)
Js = sorted(P, key=float); ts = [P[j] for j in Js]
ck("bath parents: local Lambda never falls as J grows (Lambda_Phi <= Lambda_P)",
   all(t is not None for t in ts) and all(a <= c for a, c in zip(ts, ts[1:])), str(ts))
ck("J = 0 reproduces Phi's low collision within one scan step of the formula", abs(P["0.0"] - tf) <= 0.25)
print(f"\n{sum(res)} of {len(res)} checks pass")
