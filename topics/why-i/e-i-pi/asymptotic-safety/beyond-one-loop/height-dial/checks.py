"""Height-dial test: does putting e at a different rung height s fix the 4.5% gap? Run: python3 checks.py"""
from hd import *
out = []; ok_all = True
def check(tag, cond, msg):
    global ok_all; ok_all &= bool(cond); out.append(f"[{'PASS' if cond else 'FAIL'}] {tag}: {msg}"); print(out[-1])

# H1  the machinery: e at height s (ours, from e-every-height)
vals = {s: e_s(s) for s in (0, 0.5, 1, 1.5, 2)}
check("H1", abs(vals[0] - 2) < 1e-9 and abs(vals[1] - e) < 1e-9, "e_s: " + ", ".join(f"e_{s} = {v:.4f}" for s, v in vals.items()))

# H2  the fit: which height turns 23/11 into the circle operation's own rate
L_lit = 23/22; s_star = s_for(L_lit); e_star = e_s(s_star)
check("H2", abs(log(e_star) - L_lit) < 1e-9 and abs(e_star - e**(23/22)) < 1e-9,
      f"theta = (d-2) ln e_s hits 23/11 at s* = {s_star:.4f}, e_s* = {e_star:.4f} = e^(23/22)")

# H3  e_s is continuous and increasing, so every ratio above ln 2 has its own s: one match is free
grid = np.linspace(0, 2.5, 26); es = np.array([e_s(s) for s in grid])
check("H3", np.all(np.diff(es) > 0), f"e_s increases on [0, 2.5]: {es[0]:.3f} -> {es[-1]:.3f}; any theta/(d-2) in ({log(2):.3f}, {log(es[-1]):.3f}) is hit by exactly one s")

# H4  the second numbers, each with its own height (all published; ratios theta/(d-2))
v = newton(eh_beta(), [0.7, 0.19]); th = np.linalg.eigvals(-jac(eh_beta(), v))
nums = [("Litim, lambda = 0 (23/11)", 23/22),
        ("Bonanno-Reuter's own beta", 1.18),
        ("two loops, 2+eps at eps = 2", ak_ratio(2.0)),
        ("Litim, full EH, real part", th[0].real/2),
        ("large N (exponents 2 and 4)", 1.0)]
rows = []
for name, L in nums:
    s = s_for(L); pred = 2*L_lit; actual = 2*L
    rows.append((name, L, s, e_s(s), actual, pred, pred/actual - 1))
for r in rows: print(f"     {r[0]:32s} needs ln e_s = {r[1]:.4f} -> s = {r[2]:.3f}, e_s = {r[3]:.4f}; theta = {r[4]:.3f}, s* predicts {r[5]:.3f} ({100*r[6]:+.0f}%)")
seconds = rows[1:]
check("H4", all(abs(r[6]) > 0.04 for r in seconds),
      "s* fixed by 23/11 misses every second published number: " + ", ".join(f"{100*r[6]:+.0f}%" for r in seconds))

# H5  the full EH exponent is complex; a real height gives a real number
check("H5", abs(th[0].imag) > 3, f"full EH (Litim): theta = {th[0].real:.3f} +- {abs(th[0].imag):.3f} i; theta = (d-2) ln e_s is real for real s")

# H6  dimension: the needed height runs with d and returns to 1 at d = 2
ds = [2.05, 2.5, 3, 3.5, 4, 5, 6]
sl = [s_for(litim0_ratio(dd)) for dd in ds]
sa = [s_for(ak_ratio(dd - 2)) for dd in ds[:5]]
print("     d:            " + "  ".join(f"{dd:5.2f}" for dd in ds))
print("     s (Litim l=0):" + "  ".join(f"{x:5.3f}" for x in sl))
print("     s (two loops):" + "  ".join(f"{x:5.3f}" for x in sa))
check("H6", abs(litim0_ratio(4) - 23/22) < 1e-9 and abs(sl[0] - 1) < 0.01 and abs(sa[0] - 1) < 0.07 and max(sl) - min(sl) > 0.1 and max(sa) - min(sa) > 0.5,
      f"needed s goes to 1 as d -> 2 (theta = eps exactly at one loop), is {sl[4]:.3f} (Litim) or {sa[4]:.3f} (two loops) at d = 4, and moves with d: not a constant")

# H7  the low-energy end: G constant needs V_s(a)/a -> a constant as a -> infinity. Only s = 1 does that.
A = [10.0, 100.0, 1e3, 1e4]
tab = {s: [V(s, a)/a for a in A] for s in (0.5, 1, 1.1, 1.5, 2)}
for s, r in tab.items(): print(f"     s = {s:3}: V_s(a)/a at a = 1e1..1e4: " + "  ".join(f"{x:7.4f}" for x in r))
check("H7", abs(tab[1][-1] - 1) < 1e-3 and tab[0.5][-1] < 0.5*tab[0.5][0] and tab[2][-1] > 2*tab[2][0] and tab[1.1][-1] > tab[1.1][1],
      "at a = 1/g -> infinity the rate settles only for s = 1; s < 1 drifts to 0, s > 1 grows without end (log-type), so G would keep running at low energy")

# H8  the high-energy end: the wall rate is 1 for every height
wr = {s: wall_rate(s) for s in (0.5, 1, 1.1, 1.5, 2)}
check("H8", all(abs(x - 1) < 1e-3 for x in wr.values()), "rate at the finite wall E^s(-inf) = 1 for every s: " + ", ".join(f"s={s}: {x:.4f}" for s, x in wr.items()))

# H9  so the ratio (wall rate)/(low-energy rate) is 1 at s = 1 and 0 or infinity otherwise; 23/22 is never reached
drift = {s: abs(V(s, 1e5)/1e5/(V(s, 1e4)/1e4) - 1) for s in (0.5, 1, 1.1, 1.5, 2)}
check("H9", drift[1] < 1e-3 and all(drift[s] > 0.02 for s in (0.5, 1.1, 1.5, 2)),
      "low-energy rate still drifting between a = 1e4 and 1e5: " + ", ".join(f"s={s}: {100*x:.1f}%" for s, x in drift.items()) +
      "; so (wall rate)/(low-energy rate) is a number only at s = 1, where it is 1. 23/22 needs a second, unrelated dial")

# H10 rescaling (1/k) B_s(ka, kb) moves the wall to w_s/k and leaves both rates alone
k = 0.37; s = 1.5; a = 50.0; hb = 1e-6
Bk = lambda a, b: float(Es(s, Es(-s, np.array([k*a])) + Es(-s, np.array([k*b])))[0])/k
Vk_num = (Bk(a, hb) - Bk(a, -hb))/(2*hb)          # generator of (1/k) B_s(ka, kb), read off the law itself
check("H10", abs(Vk_num - V(s, k*a)) < 1e-4*abs(Vk_num),
      f"generator of (1/k)B_s(ka,kb) = V_s(ka) ({Vk_num:.5f} vs {V(s, k*a):.5f}): the k dial multiplies both end rates by k and moves the wall to w_s/k, so their ratio is untouched")

# H11 the exact second dial that does fit: lift u^p (1 + k u) with p = theta/(d-2) - 1 = 1/22 at Litim lambda = 0
p = L_lit - 1
check("H11", abs(p - 1/22) < 1e-12, "the dial that fits is p = theta/(d-2) - 1 = 1/22: it is theta restated, not a prediction")

# H12 large N is exact at s = 1 and only there
check("H12", abs(s_for(1.0) - 1) < 1e-9, f"large N needs ln e_s = 1 exactly, s = {s_for(1.0):.6f}; s* = {s_star:.4f} would give {2*L_lit:.4f} and {4*L_lit:.4f} in place of 2 and 4")

print("ALL PASS" if ok_all else "SOME FAILED")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + ("\nALL PASS\n" if ok_all else "\nSOME FAILED\n"))
