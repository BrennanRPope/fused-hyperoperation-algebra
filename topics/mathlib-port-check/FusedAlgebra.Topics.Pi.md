# FusedAlgebra.Topics.Pi: review

Theorems (45): SAME 21, GENERAL 10, EQUIV 1, MATHLIB 8, DROPPED-OK 5. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked:
- `G` became `GaussianInt` = `Zsqrtd (-1)`. Multiplication, norm (re² + im²) and conjugation (`star`) agree with G's formulas.
- `gma`/`gpow`/`G.pow` became library `✳` / `mpow` / `^`, with the same recursions.
- `toG x = ⟨den, num⟩`, `tl`, `q4` and `firstCross`: identical.
- `gt`/`greg` became `gregTerm`/`greg` over any field. `greg` is now a `Finset.sum`, with the same recursion (`greg_succ`).
- `machinLo`/`machinHi`: identical ℚ formulas.

The Machin intervals (checked as asked):
- `machinLo N = 16·greg(1/5)(2N) − 4·greg(1/239)(2N+1)` and `machinHi N = 16·greg(1/5)(2N+1) − 4·greg(1/239)(2N)` are unchanged.
- `machinLo_le_machinHi` (was `machin_lo_le_hi`) and `machin_nested` are unchanged.
- `machin_50_digits` uses the same N = 22. I compared both 50-digit bounds (…937510 and …937511, over 10^50) character for character. Still `decide +kernel`.
- Neither the original nor the port proves that the intervals contain π. The port says so plainly. The new `hasSum_gregTerm` (Gregory's series = arctan for ‖x‖ < 1) is only a link to mathlib.

Gaussian identities, rechecked by hand:
- (5+i)⁴(239−i) = 114244(1+i)
- (2+i)(3+i) = 5+5i
- (2+i)²(7−i) = 25+25i
- (3+i)²(7+i) = 50+50i
- `cross_10/100/1000`: same fuel and same answers (32, 315, 3142).

Generalisations:
- `circle_real_points` and `multAdd_self_zero`: ℚ → any ring without zero divisors.
- `parity_twice`: ℚ → any ring.
- The Gregory lemmas: ℚ → any (linearly ordered) field.
- None of these adds a hypothesis that ℚ fails.

Proof-route changes: the `G` ring-law lemmas (`gma_*`, `gpow_lift`, `norm_mul`, `norm_conj`, `mul_assoc'`, `div_le_div_of`, `neg_one_pow_*`) are replaced by library or mathlib lemmas. These were bookkeeping. Your Gregory nesting proofs are kept, in a shorter form.

New port theorems with no original: 2 (`greg_succ`, `hasSum_gregTerm`).
