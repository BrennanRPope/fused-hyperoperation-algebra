# FusedAlgebra.Topics.Fermat.FreyPieces

- The original imported `Topics.Topology.CubeProofs` only for `IsPrime`, `exists_prime_factor`
  and `MEq`. In the port these are mathlib's `Nat.Prime`, `Nat.exists_prime_and_dvd` and
  `Int.ModEq`, so the port does not import `CubeProofs` at all (only `Topics.Elliptic`).
- `FLT n` is replaced by mathlib's `FermatLastTheoremFor n` (which also assumes `z ≠ 0`);
  `fermatLastTheoremFor_iff` shows the two shapes agree. "For all `n ≥ 3`" is mathlib's
  `FermatLastTheorem`.
- The hand-rolled `Weier` structure and its invariants are mathlib's `WeierstrassCurve`, `b₂`,
  `b₄`, `b₆`, `b₈`, `c₄`, `Δ`. The Frey curve itself is `EllipticLead.freyCurve` (in
  `Topics/Elliptic.lean`).
- `normalize_core` and `frey_normalize` were strengthened so that the near-duplicate proofs in
  `FermatRoute` (`normalize_core'`, `frey_normalize_cop`) are not needed twice. `cop_left` from
  `FermatRoute` lives here as `gcd_eq_one_of_add_pow_eq`, because `frey_normalize` needs it.
- New, with no original: `fermatLastTheoremFor_iff`, `frey_model_equation`,
  `variableChange_freyCurve` (the minimal model as mathlib's `VariableChange` with `u = 2`,
  `s = 1`, over a field with `2 ≠ 0`), `freyModel_map_zmod`, `gcd_eq_one_of_add_pow_eq`.
- Nothing deep is proved: modularity, Mazur, Ribet and the genus formula are not in the file.
  "Minimal model" is used as a name only; minimality of the model is not proved (it was not in
  the original either).
