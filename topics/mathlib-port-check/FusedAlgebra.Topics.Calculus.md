# FusedAlgebra.Topics.Calculus: review

Theorems (29): SAME 18, DROPPED-OK 4, GENERAL 3, MATHLIB 3, EQUIV 1. No WEAKER / DIFFERENT / MISSING / VACUOUS.

Definitions checked: `diff f k = f (k+1) - f k` is mathlib `Δ_[1] f k` (fwdDiff h f x = f (x + h) - f x: forward, step 1). `dif (k+1) f = Δ (dif k f)` is `(Δ_[1])^[k+1] f`. `mint` as a Finset.prod in PreQuasiregular has the same recursion (mint_zero, mint_succ). liftProd, sumTo and isum are range sums/products. Cube.layer is Nat.choose.

Medium-confidence items:
- `dif_shift`, marked MATHLIB as `fwdDiff_iter_comp_add`. I believe mathlib states Δ_[h]^[n] (fun x => f (x + m)) y = Δ_[h]^[n] f (y + m), which covers it with m = 1. The port never uses it by name. A Lean `#check fwdDiff_iter_comp_add` would settle it.
- `dif_one`, DROPPED-OK to `Function.iterate_succ_apply'`. There is no port statement. Δ^[k] f 1 = Δ^[k] f 0 + Δ^[k+1] f 0 follows from iterate_succ_apply' plus unfolding fwdDiff and rearranging, so it is covered, but only after a small step.
- `rate_product` (EQUIV, high): the hypotheses F' = …, G' = … are substituted into the conclusion.

Route changes worth knowing:
- `newton`: the user's hand proof (dif_shift, dif_one, pascal_sum, induction on n) is replaced by one line from mathlib's `shift_eq_sum_fwdDiff_iter`. `pascal_sum` is kept as a statement but is no longer used.
- `ftc_add`: the user's induction is replaced by mathlib `Finset.sum_range_sub`.
- `mint_lift`/`one_add_mint`, `mint_add`, `mint_multAdd`: the user's inductions are replaced by monoid product lemmas (map_prod, prod_range_add, prod_mul_distrib).
- `mpow_mul` now comes from Basic via monoid `pow_mul`.

New port theorems with no original: 1 (mint_zero).
