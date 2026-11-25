# FusedAlgebra.Topics.QuantumNoise

Theorems: 21. GENERAL 17, SAME 4. No findings.

- Added hypotheses, all automatic in the original's ℚ:
  - `maj3_polar`, `vote_sub_self` (was `f_fixed`) and `macwilliams_rep3` need `[NeZero (2 : K)]`. This is genuinely needed: each fails in characteristic 2.
  - `undetected_two_ways` needs `[CharZero K]`. Only characteristic ≠ 2 is needed, so this is stronger than necessary, but it is automatic in ℚ.
  - `depolarizing` needs only a Field. The identity holds for any constant c in place of 4/3, so characteristic 3 is harmless.
- The Shor brackets keep the same rationals (806/1000, 807/1000, 850/1000, 851/1000), with `gx`/`gz` built from `vote` (= the old `f`).
- `Pauli` is now `Pauli R`, with identical `comp`, `lx`, `ly`, `lz` formulas.

Route changes: none meaningful (`norm_num` instead of `decide +kernel` for the brackets).

New port theorems with no original: 0.
