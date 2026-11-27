import json, math
d = json.load(open("spinbuild.json"))
res = []
def chk(name, ok, val):
    res.append(ok); print(f"{'PASS' if ok else 'FAIL'}  {name}: {val}")
r = d["three_x_rule"][0]
chk("S1 first zero x1 = 28.27: neither 2x1 nor 3x1 is a zero (independent Ising/uniform spins ruled out)",
    r["miss_2x"] > 0.5 and r["miss_3x"] > 0.5, (round(r["miss_2x"], 3), round(r["miss_3x"], 3)))
t = d["two_spin"]; full = t["0.5"]["real_zeros"]
chk("S2 two Ising spins: all zeros real for J >= 0, real zeros lost for J < 0",
    t["0.0"]["real_zeros"] == t["0.1"]["real_zeros"] == full >= math.floor(t["0.5"]["expected_all"]) and t["-0.1"]["real_zeros"] < full and t["-0.4"]["real_zeros"] < full,
    {k: v["real_zeros"] for k, v in t.items()})
chk("S3 with a ferromagnetic link the zeros break the 3x rule", d["two_spin_J0.5_breaks_3x"], d["two_spin_J0.5_first_zeros"][:2])
chk("S4 Phi(u) = sum_n n^{-1/2} phi1(u + ln(n)/2) (a mixture of shifts)", d["mixture_identity_maxerr"] < 1e-12, d["mixture_identity_maxerr"])
chk("S5 Phi is even", d["Phi_even_err"] < 1e-12, d["Phi_even_err"])
chk("S6 the single bump phi1 goes negative; only the total Phi is >= 0", d["phi1_min"] < -0.1 and d["Phi_min"] > -1e-12, (d["phi1_min"], d["Phi_min"]))
chk("S7 sum over coprime a,b of a^-x b^-y = zeta(x) zeta(y)/zeta(x+y) (x=2.5, y=3.1; truncated sum)", d["coprime_identity_relerr"] < 1e-4, d["coprime_identity_relerr"])
w = d["total_weight_zeta2_over_zeta2s"]
chk("S8 the symmetric primon gas's total weight blows up as sigma -> 1", w["1.5"] < w["1.2"] < w["1.05"] and w["1.05"] > 200, w)
fv = d["finite_volume"]
chk("S9 every finite gas (primes <= P, sigma = 1/2) is zero-free on the imaginary field axis and tracks log|zeta|^2",
    all(math.isfinite(v["min_logZ"]) for v in fv.values()) and fv["10000"]["corr_with_log_zeta2"] > 0.85,
    {k: (round(v["min_logZ"], 2), round(v["corr_with_log_zeta2"], 3)) for k, v in fv.items()})
print(f"\n{sum(res)}/{len(res)} checks pass")
