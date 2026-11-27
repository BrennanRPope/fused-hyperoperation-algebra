# Mathlib candidates (DRAFT)

Files written in Mathlib's own style, meant to be offered to Mathlib once the theory is further
along. Nothing here has been sent anywhere: no Mathlib pull request is planned yet (Brennan,
Nov 25, 2026). CI builds this library on every push to main, with the same Lean and Mathlib pins as
the rest of `lean-library-mathlib/`, and `Check.lean` fails if any theorem here uses `sorry` or an
axiom beyond Lean's three standard ones.

## The Tate curve (Nov 27, 2026)

Asked for by Brennan after the FLT tracking work ([`topics/flt-tracking/`](../../topics/flt-tracking/)):
the Tate curve was the one lead there that Mathlib had no definitions for.

| File | What is proved |
|---|---|
| `TateCurve.lean` | `curve R`: the Tate curve `Y² + XY = X³ + a₄X + a₆` over `R⟦q⟧` for any commutative ring `R`, with `a₄ = −5s₃`, `a₆ = −(5s₃ + 7s₅)/12`, `s_k = Σ σ_k(n) qⁿ` |
| | `twelve_dvd`: `12 ∣ 5d³ + 7d⁵`, so `a₆` has integer coefficients and the curve exists over every ring (`twelve_mul_a₆`) |
| | `curve_c₄`, `curve_c₆`: `c₄ = E₄ = 1 + 240s₃` and `c₆ = −E₆ = −(1 − 504s₅)` |
| | `curve_c_relation`: `1728 Δ = E₄³ − E₆²` |
| | `curve_Δ`, `X_dvd_Δ`, `not_X_sq_dvd_Δ`: `Δ = q · (unit)`, so `Δ` has order exactly 1 in `q` |
| | `curve_map_constantCoeff`: at `q = 0` the curve is the node `Y² + XY = X³`, whose two tangent slopes `0` and `−1` are defined over `R` |
| | `hasSplitMultiplicativeReduction`: over a field `k`, the Tate curve over `k⸨q⸩` has split multiplicative reduction with respect to `k⟦q⟧` (and the model is minimal) |
| `TateCurve/ModularForms.lean` | `qExpansion_E₄`, `qExpansion_E₆`: Mathlib's normalised Eisenstein series `ModularForm.E₄`, `ModularForm.E₆` have q-expansions `TateCurve.E₄ ℂ`, `TateCurve.E₆ ℂ` |
| | `curve_c₄_eq_qExpansion`, `curve_c₆_eq_qExpansion`: so the Tate curve's `c₄`, `c₆` over `ℂ⟦q⟧` are those q-expansions |
| | `qExpansion_discriminant`: the Tate curve's discriminant over `ℂ⟦q⟧` is the q-expansion of Mathlib's modular discriminant `Δ` |

41 theorems, all resting only on `propext`, `Classical.choice` and `Quot.sound` (checked locally with `Check.lean`'s `#check_mathlib_candidates`).

### Assumptions

- `R` is any commutative ring for the definitions and the formulas for `c₄`, `c₆`, `Δ`;
  `not_X_sq_dvd_Δ` also needs `R` nontrivial.
- The reduction theorem is over a field `k`, with discrete valuation ring `k⟦q⟧` and fraction
  field `k⸨q⸩` (Laurent series). It does not cover the Tate curve over a general complete
  non-archimedean field at a point `|q| < 1`; that is the Lean FLT project's version.
- `TateCurve.isMinimal_of_valuation_c₄_eq_one` is copied from the Lean FLT project (Kevin Buzzard,
  William Coram, Claude; Apache 2.0), with credit in the file, until it reaches Mathlib.

### Known vs new

All of the mathematics is classical: Tate (1959, published 1995), Silverman, *Advanced Topics in
the Arithmetic of Elliptic Curves*, chapter V; the Eisenstein series formulas are Silverman I.7.
Nothing here is new mathematics, and none of it uses the circle operation.

In Lean: Mathlib has no Tate curve. The Lean FLT project (checked at commit 0198058, Nov 18, 2026)
has the same `a₄`, `a₆` series over `ℤ⟦q⟧`, the Weierstrass identity behind the uniformisation
(proved), and the curve over a non-archimedean field (its uniformisation still `sorry`). As far as
we found, it does not prove `c₄ = E₄`, the order of `Δ` in `q`, the reduction type of the formal
curve, or the link to Mathlib's modular forms; those are what this file adds. The two should be
merged before either goes to Mathlib.

### Not yet here

- The product formula `Δ = q ∏ (1 − qⁿ)²⁴` for the formal curve. Mathlib proves the analytic
  version (`ModularForm.discriminant_eq_q_prod`); turning it into an identity of formal power
  series needs the q-expansion of an infinite product, which Mathlib does not yet have.
- The `j`-invariant `j = 1/q + 744 + 196884q + ⋯` as a Laurent series.
- Tate's uniformisation `K^×/q^ℤ ≅ E_q(K)`.
- Where the circle operation meets this: at `q = 0` the curve is the node, whose smooth points
  form the multiplicative group with chord law the circle operation ∘₁ (`FLTLeads/SplitNode.lean`
  in `topics/flt-tracking/lean/`). The Tate curve deforms that group into `K^×/q^ℤ`. Stating this
  needs the uniformisation, so it stays a lead.

## Build

From `lean-library-mathlib/`: `lake exe cache get`, `lake build`, `lake env lean Check.lean`.
