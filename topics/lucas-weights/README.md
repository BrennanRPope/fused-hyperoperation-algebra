# DRAFT: Lucas corner weights U_m(k, d)

> **DRAFT, not finished.** Pushed Nov 25 so it isn't lost. The work stopped when Brennan merged the corner-weights work into one stream, led by the size/pair-weights task. No figures or Lean yet. [HANDOFF.md](HANDOFF.md) holds everything found so far, including a Lean plan.

**Early verdict:** for the two-dial plane law (x + y + k·xy)/(1 − d·xy), a chain of n terms weights each rank-m cube corner by the Lucas sequence Uₘ(k, d) on top and −d·U_(m−1) on the bottom. The weights are known (Lucas 1878; Dickson and Chebyshev polynomials), and the corner formula is a two-line corollary of Cayley–Hamilton for a 2×2 matrix with trace k and determinant d. So it is probably folklore: useful as a tool, not a research lead on its own. Named points on the plane give plain addition, the circle operation, inclusion–exclusion, tangent and tanh addition, Fibonacci, Pell, Jacobsthal and Mersenne weights.

**Known vs ours:** the general (k, d) formula written out was not found in print in five searches on Nov 25, so it is "ours as stated, probably folklore". One recent paper (Belbachir, arXiv 2610.24667) is unread and could be a prior statement.

Run `python3 checks.py` (16 symbolic checks pass, sympy; output in `checks_output.txt`).
