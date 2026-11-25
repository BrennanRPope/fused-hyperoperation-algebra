# The wall at Re s = 0, and the Lee–Yang comparison (Nov 25, 2026)

Follow-up to the dial survey. Z_k(s) = ∏_p (1 − k·p^−s)^(−1/k) = ∏_N ζ(Ns)^(M(k,N)/k),
with M(k,N) the aperiodic necklace counts.

## Results

1. **Wall.** For every real k other than ±1, Re s = 0 is a natural boundary of Z_k.
   k = ±1 are the only values with no wall (ζ(s)/ζ(2s) and ζ(s)); k = 0 is Landau–Walfisz (1920).
2. **Known, not new.** As three secondary sources state it, Dahlquist (Ark. Mat. 1, 1952, 533–554)
   covers any Euler factor h holomorphic at 0 with h(0) = 1 and only isolated singularities in the
   closed unit disc: Re s = 0 is a natural boundary unless h is a finite product of cyclotomic factors.
   Our factor (1 − kz)^(−1/k) fits for every real k, and is cyclotomic only at k = ±1. So the wall is a
   special case of Dahlquist. Caveat: the original paper is behind a bot check; the statement here is
   from Bhowmik and Schlage-Puchta (arXiv 1001.1869 and 1001.1891) and arXiv 2406.18190 (Theorem 9.2).
   Estermann (1928) needs polynomial factors with integer coefficients, so it covers k = −1 only.
3. **Our own route** (short, standard style, shown for transparency, not claimed as new):
   near s0 = ρ/N, Z_k ~ (s − s0)^C with C = (1/k) Σ m(N's0) M(k,N'). That sum is a polynomial in k
   of degree max N' with nonzero leading coefficient, so it vanishes only on a countable set of
   algebraic k. The points ρ/N crowd every point of the imaginary axis. Hence: unconditional wall for
   every k outside a countable algebraic set (all transcendental k), and for every k ≠ 0, ±1 under RH
   (each ρ/N then has a single contributor, and M(k,p) = (k^p − k)/p ≠ 0 at prime p).
   Rational k can cancel: at k = −1/2, M1 + M2 + M3 = 0 (proved in Lean). Dahlquist's argument must
   handle this with extra input; we did not reproduce that step.
4. **Lee–Yang.** The copies points k = −1/m and the Lee–Yang (Beraha) points k = 1 + 2cos(2π/n) meet
   only at k = −1 (m = 1, n = 2). k = 0 is a Lee–Yang point (n = 3) and the limit of the copies
   points; k = 1 (n = 4) is the other no-wall point. No Lee–Yang point lies in (−1, 0), and for m ≥ 2
   −1/m is not an algebraic integer while every Lee–Yang point is. Lead closed.

## Verdict

Both leads close. The wall is known (Dahlquist 1952) and the Lee–Yang match is only at the integer
points. No line for the "helpful to mathematicians" catalog. What is worth keeping is the picture:
the strength of each singularity is a necklace count (rotations of an N-gon), and k = ±1 are exactly
the dial settings where those counts switch off.

## Files

- `wall_checks.py`, `wall_checks_output.txt`, `wall_results.json`: 15 of 15 checks pass.
- `make_figs.py` → `figs/walls_by_k.svg`, `figs/points_pile_up.svg`, `figs/strength_by_k.svg`
  (the `.inline.svg` copies use the page's colour classes).
- `DialWall.lean` (library module `FusedAlgebra.Topics.DialWall`): 15 theorems, core Lean 4.37,
  no sorry. Library build: 116 modules, 2,930 theorems, standard axioms only.
