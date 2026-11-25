# FusedAlgebra.Topics.Fermat.FreyPieces — review

Theorems (46): EQUIV 16, MATHLIB 12, GENERAL 11, SAME 4, DROPPED-OK 2 (pow_mod_two, one_add_div_mul), WEAKER 1 (pow_dvd_bound).
Definitions checked: FLT, Weier and its b2/b4/b6/b8/c4/disc, freyModel, a2of, a4of, nodeX, chordDet, freyRHS, mu2, nu2, nu3, nuInf.

## Findings
- WEAKER (minor), high: pow_dvd_bound.
  - Original: `(hq : 2 ≤ q) (hm : ¬ q ∣ m) (h : q ^ a ∣ q ^ b * m) : a ≤ b`, true for any q ≥ 2.
  - Replacement: Nat.Prime.pow_dvd_iff_le_factorization `(pp : Prime p) (hn : n ≠ 0) : p ^ k ∣ n ↔ k ≤ n.factorization p`. This needs a prime and is a different shape; the porter's note admits it.
  - The original only ever used it for primes (inside exact_exponent_pth_power), which is ported, so nothing downstream is lost.
- DROPPED-OK at medium:
  - pow_mod_two survives only as a local `have` inside coprime_parity. The named replacement does not state it; mathlib's Nat.even_pow' / Nat.pow_mod does.
  - one_add_div_mul is inlined as field_simp. It is a trivial helper.

## Conventions checked (as asked)
- Curves: Weier → mathlib WeierstrassCurve. I checked mathlib's b₂ = a₁² + 4a₂, b₄ = 2a₄ + a₁a₃, b₆ = a₃² + 4a₆, b₈ = a₁²a₆ + 4a₂a₆ − a₁a₃a₄ + a₂a₃² − a₄², c₄ = b₂² − 24b₄ and Δ = −b₂²b₈ − 8b₄³ − 27b₆² + 9b₂b₄b₆. They are identical to the original's, with the same signs and scale.
  - freyModel is ⟨1, a₂, 0, a₄, 0⟩ in both.
  - The Frey curve is EllipticLead.freyCurve = ⟨0, B − A, 0, −AB, 0⟩, which is y² = x(x − A)(x + B).
  - frey_c4 (c₄ = A² + AB + B²) and frey_disc (256Δ = (ABC)²) recomputed: they hold in any commutative ring with the same hypotheses 4a₂ = B − A − 1 and 16a₄ = −AB.
- FLT n vs FermatLastTheoremFor n: same exponent. Mathlib adds c ≠ 0, and fermatLastTheoremFor_iff proves them equal (n = 0 included).
- IsPrime → Nat.Prime: equivalent definitions. MEq q x y (q ∣ x − y) → Int.ModEq (q ∣ y − x): the sign doesn't matter.
- Level-2 counts: mu2/nu2/nu3/nuInf changed from List filters to Finset cards over ZMod 2 and a φ-sum over divisors. Values recomputed: 3, 1, 0, 2, the same as the original.

## Stays an assumption (in original and port alike)
- genus_X0_2: the genus formula 12g = 12 + μ − 3ν₂ − 4ν₃ − 6ν∞ is a hypothesis. The theorem is arithmetic on it, exactly as in the original.
- Modularity, Mazur and Ribet appear nowhere. "Minimal model" is a name only; minimality is not proved (it wasn't in the original either).

## Route changes
- exact_exponent_pth_power: proved from Nat.factorization_pow instead of the original's factor_out / not_dvd_pow / pow_dvd_bound.
- cop_left (from FermatRoute) becomes gcd_eq_one_of_add_pow_eq: IsCoprime algebra instead of a prime-factor argument.
- four_or_odd_prime and the divisor split in flt_reduction come from mathlib (Nat.four_dvd_or_exists_odd_prime_and_dvd_of_two_lt, FermatLastTheoremFor.mono).
- Euclid's lemma (euclid, euclid_int) and similar helpers are mathlib's.

New port theorems with no original: 4 (fermatLastTheoremFor_iff, frey_model_equation, variableChange_freyCurve, freyModel_map_zmod). gcd_eq_one_of_add_pow_eq is the port of FermatRoute.cop_left.
