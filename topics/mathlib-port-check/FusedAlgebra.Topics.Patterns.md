# FusedAlgebra.Topics.Patterns

Theorems (32): SAME 22, EQUIV 4, GENERAL 3, DROPPED-OK 3.

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.
- `step_trace`, `step_det`: EQUIV. They were bare ring identities (the trace and the determinant written out by hand) and are now `Matrix.trace` / `Matrix.det` of the new `stepMatrix = !![1 + k*a, a; -(q*a), 1]`, which really is `step`'s matrix (`stepMatrix_mulVec`). They say more than before, nothing less.
- The three `_pair` helpers are conjunctions of the main lemma at n and n+1; dropping them loses nothing.
- `SimpleLaws.fib` (ℤ) and `nfib` (ℕ) replaced by `Nat.fib`: same values.

Route changes:
- `lucas_add`, `lucas_fib`, `one_add_mfib` (was `mfib_lift`): paired induction through the helper lemmas replaced by `Nat.twoStepInduction`.
- `multAdd_iter`: was derived from the Lucas formula (`iter_lucas` + `multadd_den_one`); now a corollary of the new `iter_multAdd`, so it no longer goes through the Lucas route.
- `step_trace`/`step_det` use mathlib's `trace_fin_two`/`det_fin_two`.
- `mirror_num` cites `MultAdd.wAdd_neg_neg` (his own library lemma).

New port theorems with no original: 7 (stepMatrix_mulVec, iter_zero, iter_succ, iter_multAdd, mfib_zero, mfib_one, mfib_succ_succ). The porter's claim is true: `iter_multAdd` proves what the original `multadd_iter` docstring promised but didn't state.
