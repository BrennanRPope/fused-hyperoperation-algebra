# Port notes: FusedAlgebra.Topics.Relativity

- The structure `B2` (four rationals) and `B2.mul` are replaced by `Matrix (Fin 2) (Fin 2) R` and
  matrix multiplication. `boost u = !![1, u; u, 1]`.
- `kN`, `kD` were a mutual recursion; they are now the entries `(boost u ^ n) 0 1` and
  `(boost u ^ n) 0 0`, and the old recursion equations are the theorems `kN_zero`, `kN_succ`,
  `kD_zero`, `kD_succ`. `kicks` is now a `Function.iterate` (`kicks_zero`, `kicks_succ`).
- New statements, not in the original: `boost_mul_boost`
  (`boost u * boost v = (1 + u * v) • boost (vel u v)`), `det_boost` (`= 1 - u ^ 2`), and the
  equation lemmas above.
- Identities hold in any commutative ring or field. Those that divide by `2` or `4` (`dop_back`,
  `dop_wall`, `rung_natural_q`, `rung_parabola`, `rung_raw_walls`) take `[NeZero (2 : K)]`; they
  are false in characteristic 2.
- Hypotheses the proofs never needed were removed (recorded in the table): `velC_units` and
  `landscape_walls` (`c ≠ 0`), `kick_step` (`kD n u + u * kN n u ≠ 0`).
- Velocity addition is not packaged as a group: `vel` is a group law only on the open interval
  `(-1, 1)` of an ordered field, and building that subtype was not cheap. `dop_vel` and
  `ThermoFused.lens_vel` are the two homomorphism statements.
