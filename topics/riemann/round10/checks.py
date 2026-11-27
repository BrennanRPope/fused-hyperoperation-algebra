"""Round 10 checks: read r10.json and jets.json (made by r10.py and jets.py) and test each claim in the README."""
import json, math
r = json.load(open("r10.json")); j = json.load(open("jets.json"))
G = 0.5772156649015329
res = []
def ck(name, ok, detail=""):
    res.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))

for p in ("11", "101"):
    d = r["lpf_share"][p]; rel = abs(d["at_1e8"] / d["limit"] - 1)
    ck(f"lpf share of {p} at 1e8 within 1% of (1/p)prod(1-1/q)", rel < 0.01, f"rel {rel:.4f}")
ck("lpf share of 1009 still above its limit at 1e8 (slow approach)",
   r["lpf_share"]["1009"]["at_1e8"] > r["lpf_share"]["1009"]["limit"])
ck("Mertens: open fraction * log p at 1e5 within 0.5% of e^-gamma",
   abs(r["mertens_at_1e5"] / math.exp(-G) - 1) < 0.005, f"{r['mertens_at_1e5']:.5f}")
ck("wheel overcount at 1e5 within 0.5% of 2e^-gamma", abs(r["wheel_overcount_at_1e5"] / (2*math.exp(-G)) - 1) < 0.005,
   f"{r['wheel_overcount_at_1e5']:.5f} vs {2*math.exp(-G):.5f}")
ck("record gaps above 1e3 stay below (log p)^2", r["record_gap_max_ratio_above_1e3"] < 1,
   f"{r['record_gap_max_ratio_above_1e3']:.3f}")
for p in ("11", "101"):
    v = r["landau_spike_last500"][p]
    ck(f"Landau spike at p = {p} within 2% of 1", abs(v - 1) < 0.02, f"{v:.4f}")
vs = [r["variance_share"][k] for k in ("2", "3", "5", "7", "11", "31", "101", "1009", "10007")]
ck("variance share of log|zeta| grows with the primes kept", all(a < b for a, b in zip(vs, vs[1:])), f"{vs[0]:.3f} -> {vs[-1]:.3f}")
ck("prime 2 alone carries 10-25% of the variance near t = 5000", 0.10 < vs[0] < 0.25)
ck("small-gap fraction: zeros < wheel < primes",
   r["zero_gap_small_frac"] < r["wheel_gap_small_frac"] < r["prime_gap_small_frac"],
   f"{r['zero_gap_small_frac']:.3f} < {r['wheel_gap_small_frac']:.3f} < {r['prime_gap_small_frac']:.3f}")
ck("Chebyshev bias: non-squares ahead for >99.9% of primes to 1e8, mod 4 and mod 3",
   min(r["chebyshev_ahead_fraction"].values()) > 0.999)
h = r["hyper"]
for k in ("1.0", "-0.5", "2.0"):
    d = h[f"rung2_k{k}"]
    ck(f"rung-2 hyper integral -> (e^(k int f)-1)/k at k = {k}", d["err_at_max_n"] / abs(d["target"]) < 1e-4,
       f"rel {d['err_at_max_n']/abs(d['target']):.1e}")
ck("rung-3 (exponent) hyper integral -> exp(int log log(1+f))",
   abs(h["rung3_power_at_max_n"] - h["rung3_power_target_log"]) < 1e-9)
for d in h["primes"]:
    ck(f"rung-2 fold over primes at s = {d['s']} equals zeta(s)/zeta(2s) - 1",
       abs(d["rung2_fold"] - d["zeta(s)/zeta(2s)-1"]) < 1e-6, f"{d['rung2_fold']:.9f}")
    ck(f"rung-3 fold over primes at s = {d['s']} runs to 0 (log hugely negative)", d["rung3_log_of_fold"] < -1e5)
A = j["A"]
ck("matrix rung-2 product -> ordered exponential, error ~ 1/n",
   A["err_vs_ordered_exp"][-1] < 1e-3 and A["err_vs_ordered_exp"][0] / A["err_vs_ordered_exp"][-1] > 100)
ck("det of the matrix product -> exp(int tr S)", A["err_det_vs_exp_int_trace"][-1] < 1e-3)
ck("Q and R jet layers total ~ 1/n (vanish in the limit)",
   A["sum_of_Q_R_parts"][-1] < 0.01 and A["sum_of_Q_R_parts"][0] / A["sum_of_Q_R_parts"][-1] > 100)
B = j["B"]["2.0"]
ck("sum of prime layers at k = 1 equals zeta(2)/zeta(4)", abs(B["sum_k1"] - B["target_k1_zeta2_over_zeta4"]) < 1e-6)
ck("alternating sum at k = -1 equals 1/zeta(2)", abs(B["sum_km1"] - B["target_km1_one_over_zeta2"]) < 1e-6)
ck("layer 2 matches the direct sum over p < q", abs(B["e"][2] - j["B"]["e2_direct_primes_below_2000"]) < 1e-4)
for s in ("2.0", "3.0"):
    e = j["B"][s]["e"]
    ok = all(e[m]**2 >= (1 + 1/m) * e[m-1] * e[m+1] for m in range(1, 10))
    ck(f"Newton's inequalities e_m^2 >= (1+1/m) e_(m-1) e_(m+1) at s = {s}, m = 1..9", ok)
n_ok = sum(ok for _, ok, _ in res)
print(f"\n{n_ok} of {len(res)} checks pass")
