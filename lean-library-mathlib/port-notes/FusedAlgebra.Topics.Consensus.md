# FusedAlgebra.Topics.Consensus: port remarks

- `mpow_lift` (which shadowed the `Base` name) is exactly `MultAdd.one_add_mul_wpow` and
  `mpow_lift1` is `MultAdd.wpow_one` with `MultAdd.one_add_mpow`; both dropped.
- `dbl` is now `(fun x ↦ x ✳ x)^[n] a` (`Function.iterate`); `dbl_lift` → `one_add_dbl`, and
  `dbl_eq_mpow : dbl n a = mpow a (2 ^ n)` is new.
- `depth_useless` is stated for every weight `k` (`wpow k 0 n = 0`), not just `k = 1`. It is a
  general `wpow` lemma that could live in `Basic` as `wpow_zero_left`.
- `additive_split` is `map_nsmul` for the additive map built from `f`.
- `ConsensusFused.maj3` and `QuantumNoise.maj3` are the same polynomial written in a different
  order; each module keeps its own, as in the original.
- The purely numerical examples stay over `ℚ` and are proved by `norm_num`.
