# FusedAlgebra.Topics.JetFields: review

Theorems (46): SAME 33, MATHLIB 7, EQUIV 2, GENERAL 2, DROPPED-OK 2. No WEAKER / DIFFERENT / MISSING / VACUOUS. No low-confidence items.

Definitions checked: the WDual and Jet2 products are the same formulas. The weight k now sits in the type (WDual k), which changes no statement. The M3/M2 records map to Matrix with the right indices (vort, symm2 = X + Xᵀ, fsq = trace (X * Xᵀ), quad = w ⬝ᵥ X *ᵥ w, which parses as intended). reStep is the same. dDet uses updateRow and gives the same three determinants. dE2 was redefined as trace X * trace Y - trace (X*Y); expanding it by hand, the diagonal terms cancel and the off-diagonal terms are exactly the original six pairs. d1 and d2 are unchanged, and d1_eq_fwdDiff / d2_eq_fwdDiff are rfl links to fwdDiff 1 (forward, step 1).

New `WDual.dualNumberEquiv : WDual 0 ≃+* DualNumber R` (inl re + inr du) is a correct identification: the ε-part of the product is a*d + b*c + 0*b*d against a*d + b*c. No original theorem depends on it.

Table remarks: `jma` and `jma_lift` are marked "mathlib", but the replacement is the port's own Basic lemma `MultAdd.one_add_multAdd` (cosmetic). `det_one_add_smul` and `det2_step` really are covered by `LinearAlgebraMA.det_one_add_smul_fin_three` / `_fin_two`; for det2_step this is up to unfolding det_fin_two.

Route changes worth knowing:
- `sym_e2_nonpos`: the hand proof (sum of squares via sym_e2_sumsq and `sq_nonneg_int`, then omega) is replaced by `nlinarith`; `sq_nonneg_int` is now mathlib's `mul_self_nonneg`.
- `tr_sq_traceless`, `tr_cube_traceless`, `e2_vorticity_strain`, `sym_tr_cube`, `re_dQ` are now corollaries of new general identities (Newton identities, no trace hypothesis), still proved by expansion.
- `wmul_comm`/`wmul_assoc` are now fields of a CommRing instance; `twoPt_mul` is `map_mul` of a RingHom.
- `cayley_hamilton` is still proved by expansion, not from mathlib's aeval_self_charpoly.

New port theorems with no original: 36 (31 of them are simp projection lemmas for WDual/Jet2).
