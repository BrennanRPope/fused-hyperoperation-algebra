# FusedAlgebra.Topics.QuantumNoise: port remarks

- `Pauli` is now `Pauli R` for any type `R` (the original is the case `R = ℚ`), with `@[ext]`.
- `Pauli.comp` is multiplication in the group algebra of the Klein four-group and `lx`, `ly`,
  `lz` are its non-trivial characters. The port does not build that bridge to `MonoidAlgebra`;
  it would be the natural next step (then `comp_lx` etc. become "characters are ring maps").
- `f` → `vote`, `f_fixed` → `vote_sub_self`, `f_helps` → `vote_helps` (one-letter name).
- The four Shor-code brackets (`gx_below`, …) stay over `ℚ` and are proved by `norm_num` instead
  of `decide +kernel`.
- `amplitude_damping` and `bit_flip` are the same statements as `InfoFused.survival` and
  `InfoFused.bsc_polarization`; both modules keep their own copy because neither imports the
  other in the original.
