"""Finance identities as the circle operation, checked on Shiller's S&P data (1871-2022).
a (+) b = a + b + ab = (1+a)(1+b) - 1;  a (-) b = (a - b)/(1 + b) is its inverse."""
import numpy as np, itertools, load
madd = lambda a, b: a + b + a * b
msub = lambda a, b: (a - b) / (1 + b)
out = []
def say(s=""): print(s); out.append(s)

A = load.shiller_annual()
say("== 1. Fisher: real rate = nominal (-) inflation ==")
real_exact = msub(A.r_nom, A.infl)
real_add = A.r_nom - A.infl
gap = (real_add - real_exact)
say(f" years {A.index[0]}-{A.index[-1]}; additive Fisher error (pct points): mean |gap| {100*gap.abs().mean():.3f}, max {100*gap.abs().max():.2f} in {gap.abs().idxmax()}")
for y in [1918, 1947, 1974, 1980, 2022]:
    if y in A.index: say(f"  {y}: yield {100*A.r_nom[y]:.2f}%, inflation {100*A.infl[y]:.2f}%: exact real {100*real_exact[y]:.2f}%, additive {100*real_add[y]:.2f}%")
say(f" check: real (+) inflation = nominal: {np.allclose(madd(real_exact, A.infl), A.r_nom)}")
for pi in [0.5, 1.0, 2.0]:
    say(f"  at {100*pi:.0f}% inflation and a 10% real return, nominal = {100*madd(0.10,pi):.0f}% (adding says {100*(0.10+pi):.0f}%)")

say("\n== 2. Thirty years of real S&P returns as a 30-cube (1993-2022) ==")
R = A.R_real.loc[1993:2022].values
# layer k = elementary symmetric polynomial e_k of the 30 returns (sum over corners with k switches on)
e = np.zeros(len(R) + 1); e[0] = 1
for r in R: e[1:] = e[1:] + r * e[:-1]
total = np.prod(1 + R) - 1
say(f" growth multiple {np.prod(1+R):.3f}x; (+) of all 30 returns = {total:.4f}; sum of layers 1..30 = {e[1:].sum():.4f}")
say(" layer totals k=1..8: " + ", ".join(f"{e[k]:.3f}" for k in range(1, 9)))
say(f" layer 1 (adding the returns) = {e[1]:.3f}, i.e. {100*e[1]/total:.0f}% of the real growth; layers 2+ (return on return) = {total-e[1]:.3f}")

say("\n== 3. Volatility drag is the square cut out of the corner ==")
m, s = 0.07, 0.20
say(f" up {100*(m+s):.0f}% then down {100*(s-m):.0f}%: (1+m+s)(1+m-s) = (1+m)^2 - s^2 = {(1+m+s)*(1+m-s):.4f} = {(1+m)**2:.4f} - {s*s:.4f}")
Rall = A.R_real.dropna().values
am, gm = Rall.mean(), np.exp(np.log1p(Rall).mean()) - 1
say(f" S&P real 1871-2022: arithmetic mean {100*am:.2f}%, geometric (the circle operation) mean {100*gm:.2f}%, drag {100*(am-gm):.2f} pts; variance/2 = {100*Rall.var()/2:.2f} pts; variance/(2(1+m)) = {100*Rall.var()/(2*(1+am)):.2f}")

say("\n== 4. Perpetuities and the pole ==")
for i in [0.02, 0.05]:
    say(f" perpetuity at {100*i:.0f}%: price 1/i = {1/i:.1f}; duration (1+i)/i = {(1+i)/i:.1f} = 1/d with d = -( (-)i ) = {i/(1+i):.5f}")
r, g, D0 = 0.07, 0.04, 1.0
P_gordon = D0 * (1 + g) / (r - g)
P_perp = D0 / msub(r, g)
N = 4000; P_sum = sum(D0 * ((1 + g) / (1 + r)) ** k for k in range(1, N + 1))
say(f" Gordon: D0(1+g)/(r-g) = {P_gordon:.4f}; perpetuity at r (-) g = {msub(r,g):.6f}: D0/(r(-)g) = {P_perp:.4f}; summed {N} years = {P_sum:.4f}")
say(" so a growing stream is a flat perpetuity at the circle-operation difference of the rates; it blows up at g = r (a pole of order 1)")
m_ = load.shiller_monthly()
for y in ["1921-01", "1950-01", "1982-01", "2000-01", "2009-03", "2022-12"]:
    row = m_.loc[y].iloc[0]
    dy = row.D / row.P
    say(f"  {y}: dividend yield {100*dy:.2f}% = implied r (-) g; duration 1 + P/D = {1+1/dy:.0f} years")
open("identities_output.txt", "w").write("\n".join(out) + "\n")
