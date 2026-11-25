"""Checks for the closing-rate proof (fused-linear rule).
For each rung k, at b = c_k (ordinary base, fused constant c_k - 1):
  - orbit o_j of s(x) = b^(x-1) from o_0 = 2, with m = k - 3
  - (L) says: if the tower settles then o_{m+1} < log_b 3.5   (must hold at c_k)
  - (U) says: if o_{m-2} <= log_b 3 (plus a tail condition) the tower settles,
    so at c_k we must see o_{m-2} > log_b 3 or the tail condition failing.
  - x2(b): repelling fixed point of s, multiplier q_b = x2 ln b -> ln 4
Also fits theta_k = (2 - c_k)(ln 4)^k."""
import json, math
from scipy.optimize import brentq
R = json.load(open('rate_scan.json'))
L4 = math.log(4)
out = []
def orbit(b, n):
    o = [2.0]
    for _ in range(n):
        o.append(b ** (o[-1] - 1))
    return o
out.append(" k   c_k                 o_{m+1}  log_b3.5  o_{m-2}  log_b3   x2(b)-2     q_b    ")
for r in R:
    k = r['k']; b = float(r['c']); m = k - 3
    if m < 3: continue
    o = orbit(b, m + 2)
    lb35 = math.log(3.5) / math.log(b); lb3 = math.log(3) / math.log(b)
    f = lambda x: b ** (x - 1) - x
    x2 = brentq(f, 2.0 + 1e-15, 2.6)
    q = x2 * math.log(b)
    eps = 2 - b
    tail = 0.02 * 0.7 ** (m - 29) if m > 29 else float('nan')
    out.append(f"{k:2d}  {b:.15f}  {o[m+1]:.4f}   {lb35:.4f}   {o[m-2]:.4f}   {lb3:.4f}  {x2-2:.3e}  {q:.6f}")
out.append("")
out.append("theta_k = (2 - c_k)(ln 4)^k, and successive differences:")
th = [(r['k'], (2 - float(r['c'])) * L4 ** r['k']) for r in R]
for i, (k, t) in enumerate(th):
    d = t - th[i - 1][1] if i else float('nan')
    out.append(f"  k={k:2d} theta={t:.6f}  diff={d:+.6f}")
# geometric extrapolation from the last three (only after the kiss switch at k=19)
t1, t2, t3 = th[-3][1], th[-2][1], th[-1][1]
rr = (t3 - t2) / (t2 - t1)
lim = t3 + (t3 - t2) * rr / (1 - rr)
out.append(f"difference ratio at the end = {rr:.4f}; geometric extrapolation of theta = {lim:.4f} (a guess, not a proof)")
s = "\n".join(out); print(s); open('rate_checks_output.txt', 'w').write(s + "\n")
