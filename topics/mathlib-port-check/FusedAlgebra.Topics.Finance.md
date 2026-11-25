# FusedAlgebra.Topics.Finance

Theorems (16): GENERAL 10, EQUIV 2, MATHLIB 4. Definitions checked: gsum (EQUIV; the Finset sum has the same indices, q^1 to q^N), mall, growth (MATHLIB, they match).

No WEAKER / DIFFERENT / MISSING / VACUOUS / low-confidence items.

- `am_hm`. Original: `x * u = 1 → y * v = 1 → ∀ c, c * ((u+v)/2) = 1 → c ≤ (x+y)/2`. Port: `((u+v)/2)⁻¹ ≤ (u⁻¹+v⁻¹)/2`. These are equivalent because u, v > 0 force x = u⁻¹, y = v⁻¹ and c = ((u+v)/2)⁻¹.
- `average_factors`: the same reciprocal-hypothesis change, and equivalent for the same reason.
- Route changes:
  - `am_hm` was a hand-built proof (the key identity and order cases); it is now `field_simp; ring; positivity` on the difference `(u-v)^2/(2uv(u+v))`.
  - His helpers `sq_nonneg'`, `pow_nonneg'`, `one_le_pow`, `pow_le_pow'` are replaced by mathlib lemmas.
  - `mall_growth` is now a corollary of `one_add_multAddList`.
- New port theorems with no original: 0.
