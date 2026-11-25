# Dial survey: a quick pass on the five dials (Nov 25, 2026)

Page: https://claude.ai/artifact/LoM7iRCttrFnAi9VGz2kVr (source index.html)

- checks.py -> checks_output.txt: 40 of 40 pass (sympy, mpmath). results.json holds Z_k(2), the C_k curve, direct counts.
- lean/DialSurvey.lean: 18 theorems, core Lean 4.37, no sorry. Merged and compiled in the shared library as
  Topics/DialSurvey (clean build: 111 modules, 2,801 theorems).

## Findings
1. Only two independent numeric dials: k (q = k+1) and the per-element dial t. Prime dial = per-element dial on primes
   with t_p = p^-s. Twist = q x area (twisted boosts), fixed by k and directions. Choice = a setting (partial order), not a number.
2. k: only sign survives rescaling (units_change_weight). Also Tsallis q-algebra with k = 1 - q_T. No lead by itself.
3. Per-element t is odds: union of independent random subsets has odds t ✳ s; complement t -> 1/t. Uniform t = 2^a - 1
   is the fractional exponent through the mult-add ruler (exponent_ruler).
4. New: k acting on per-element dials = how many copies allowed. Per-element factor (1 - kz)^(-1/k), coefficient
   prod_{j<e}(1+jk)/e!. k=-1 subsets, k=-1/m at most m copies (proved), k=0 Poisson, k=1 multisets.
5. On primes: Z_k(s) = prod_p (1 - k p^-s)^(-1/k) = prod_N zeta(Ns)^(M(k,N)/k), M = aperiodic necklaces (cyclotomic
   identity). k=-1 fermionic primon gas zeta(s)/zeta(2s), k=1 bosonic zeta(s). Densities C_k = prod (1-1/p)(1-k/p)^(-1/k):
   C_-1 = 6/pi^2, C_0 = e^(M - gamma) = 0.7293, C_1 = 1. For k > 2 the critical point moves from s=1 to s=log2 k.
   Only k = +-1 use finitely many zeta factors; every other real k has a wall at Re s = 0. Follow-up (wall/README.md): this is a special case of Dahlquist 1952, not new.

## Verdicts
- k alone: not a lead. Twist: geometric anchor, no new prediction (the fusion inside it is an identity of GR).
- Per-element: lead only paired with k. Primes: same dial as per-element, test bed.
- Choice: not a dial; the fractional power set in ZF stays the open lead.
- k on elements / primes: followed up (wall/). The wall at Re s = 0 for k != +-1 is known (Dahlquist 1952), and the
  k = -1/m points meet the Lee-Yang values only at k = -1 (0 is a limit). Both leads close; nothing for the catalog.

Sources: Julia 1990; Spector 1990; Metropolis-Rota 1983; Landau-Walfisz 1920; Shirai-Takahashi 2003; Tsallis 1988;
Thomas 1926; Wigner 1939; Tarski 1924; Howard-Rubin 1998.
