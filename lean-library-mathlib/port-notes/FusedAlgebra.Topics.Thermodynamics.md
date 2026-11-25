# Port notes: FusedAlgebra.Topics.Thermodynamics

- Identities are stated in any commutative ring `R` or field `K`, inequalities in any linearly
  ordered commutative ring. `pauling_ice` needs characteristic zero and `lens_involution` needs
  `2 ≠ 0` (`[NeZero (2 : K)]`); both are false in characteristic 2.
- `outer` is now `xs.flatMap fun x ↦ ys.map (x * ·)`; `outer_nil`, `outer_cons` (new) recover the
  original defining equations.
- The multiplicative map `f` of `outer_map`, `power_sum_multiplies`, `qlog_qsum` is now any
  `MulHomClass` map (was `f : ℚ → ℚ` with `∀ u v, f (u * v) = f u * f v`). The "logarithm" `L` of
  `free_energy_adds` and `renyi_adds` keeps its hypothesis form, because its target is only an
  additive type.
- Hypotheses that the proofs never needed were removed (each recorded in the table):
  `carnot_series` (`h ≠ 0`), `reverse_is_inverse` (`1 - eta ≠ 0`), `shields` (`1 - r ≠ 0`),
  `wick_num` (`i * i = -1`; the numerator identity holds for every `i`).
- `fermion_lift` and `fermion_cube` restate `MultAdd.one_add_multAddList` and
  `MultAdd.one_add_multAddList_eq_sum_sublists`; kept as one-liners for the physical reading.
- `lens_vel` says the lens is a homomorphism from velocity addition to multiplication. It is not
  packaged as a `MonoidHom` because `vel` is only a group law away from its poles
  (`1 + u * v = 0`), so there is no total monoid structure to hang it on.
