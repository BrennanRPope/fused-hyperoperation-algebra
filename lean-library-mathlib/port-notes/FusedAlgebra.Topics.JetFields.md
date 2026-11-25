# FusedAlgebra.Topics.JetFields: port remarks

- The port imports `FusedAlgebra.Topics.LinearAlgebra` (the original imported only `Base`) to
  share `LinearAlgebraMA.e2` and `det_one_add_smul_fin_three`; the original had its own copies of
  the 3×3 record, `e2` and the layer expansion.
- `WDual` now takes the weight as a type parameter: `WDual k` is the ring `R[ε]/(ε² - k ε)` with
  a hand-written `CommRing` instance, so `wmul k x y` is `x * y`. It is the quotient
  `AdjoinRoot (X ^ 2 - C k * X)`; that identification was not built. At `k = 0` the exact match
  with mathlib is given: `WDual.dualNumberEquiv : WDual (0 : R) ≃+* DualNumber R`. `twoPt` is a
  bundled `RingHom`, and `WDual.twoPtEquiv` (new) is the ring isomorphism with `K × K` for
  `k ≠ 0` in a field.
- `Jet2 R` is `R[ε]/(ε³)` with a hand-written `CommRing` instance; `jmul` is `*` and `jma` is the
  library's `✳`, so `jma_lift` is `MultAdd.one_add_multAdd`. It was not identified with a
  quotient of `Polynomial R` or with `PowerSeries.trunc 3`: no cheap exact match.
- Both `CommRing` instances use `nsmulRec` / `zsmulRec` and the default casts and powers; no
  `Algebra R` instance was added.
- `cayley_hamilton` is proved by expansion (`ring`), not as a corollary of
  `Matrix.aeval_self_charpoly`; that would need the 3×3 characteristic polynomial as a
  polynomial identity. `DialSweep.eval_charpoly_fin_three` has its evaluated form.
- `dE2` is defined as `trace X * trace Y - trace (X * Y)` instead of the written-out sum (same
  polynomial). `dDet` uses `Matrix.updateRow`; it also equals `trace (adjugate X * Y)` (Jacobi),
  which is not proved here.
- New general identities, of which the traceless originals are corollaries:
  `trace_mul_self_fin_three`, `trace_mul_mul_self_fin_three`, `e2_vorticity_strain_add`.
- `d1`, `d2` are kept as two-variable forward differences, with `d1_eq_fwdDiff`,
  `d2_eq_fwdDiff` linking them to mathlib's `fwdDiff`.
