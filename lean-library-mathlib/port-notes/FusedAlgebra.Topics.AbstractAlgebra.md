# FusedAlgebra.Topics.AbstractAlgebra: port remarks

- The hand-written records `CommMonoidLaws`, `CommGroupLaws`, `OpIso` are gone; the statements
  that used them are now instances / equivalences in mathlib's hierarchy:
  - `(R, ✳, 0)`: `PreQuasiregular R` with `PreQuasiregular.instCommMonoid` (Basic);
  - the shift: `MultAdd.oneAddEquiv : PreQuasiregular R ≃* R` (Basic);
  - `(R, wAdd k, 0)`: the new type copy `Weighted k` (a one-field structure, like
    `PreQuasiregular`) with `Weighted.instCommMonoid`; `Weighted.scaleHom`, `Weighted.scaleEquiv`,
    `Weighted.commGroupOfPowEqZero`;
  - `(R, fadd, ✳)`: the new type copy `Fused R` with `Fused.instCommRing`, obtained by
    transporting the ring structure of `R` along the shift (`Equiv.commRing`), and
    `Fused.shiftRingEquiv : Fused R ≃+* R`;
  - the Tate quotient: `tateCon q : Con (PreQuasiregular R)` and its quotient monoid.
  `Weighted` is general enough that it could live next to `wAdd` in `FusedAlgebra.Basic`
  (`Weighted 1` is `PreQuasiregular`); it was put here because Basic is not mine to edit.
- `Weighted.commGroupOfPowEqZero` is a reducible definition rather than a global instance,
  because it needs the proof `k ^ m = 0`. Use `letI := Weighted.commGroupOfPowEqZero k m hk`.
- "Has a ✳-inverse" was `∃ c, (1 + a) * c = 1`; it is now `IsUnit (1 + a)`, and
  `multAdd_unit_iff_isQuasiregular` (new) says this is mathlib's `IsQuasiregular a`.
  `wAdd_unit_iff` was an implication despite its name and is renamed
  `isUnit_one_add_mul_of_wAdd_eq_zero`.
- The `Fin n` examples are stated on `ZMod n` (same underlying type for these `n`, and the
  mathlib-standard ring); all are still closed by `decide`. Names `finN_*` became `zmodN_*`.
- New: `exists_slope_wall_iff`, the general rule behind the five mod-p examples (walls of slope
  addition exist mod a prime p iff p % 4 ≠ 3), a two-line corollary of mathlib's
  `ZMod.exists_sq_eq_neg_one_iff`. The original only stated this rule in a comment.
- New small lemmas: `fadd_def`, `Weighted.val_one`, `Weighted.val_mul`, `Weighted.val_scaleHom`,
  `Fused.val_add/val_mul/val_zero/val_one/val_neg`, `Fused.shiftRingEquiv_apply`, `tateCon_apply`.
- The integers acting on `(R, ✳)`: `mpow_multAdd` is `mul_pow` in `PreQuasiregular R`; the other
  three power lemmas of the original are already in Basic (`mpow_mul`, `mpow_one`) or one line
  (`zero_mpow`).
