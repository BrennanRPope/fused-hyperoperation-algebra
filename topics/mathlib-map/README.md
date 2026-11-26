# Mathlib map: what Mathlib already has that is ours

**Verdict.** Mathlib already holds our law at k = 1: it is Mathlib's *circle operation*
`x ∘ y = x + y + xy` (`PreQuasiregular`, `IsQuasiregular`, `quasispectrum`), and about 1,550 theorems
in 53 files are built directly on it (Unitization, the Neumann series, spectrum theory, the whole
non-unital continuous functional calculus). The links here show three new things on top of that:
every dial k ≠ 0 is the same monoid; the walls an element meets as the dial turns are exactly its
quasispectrum, flipped by λ ↦ −1/λ; and the Jacobson radical is exactly the elements that never meet
a wall at any dial.

Page: see `PAGES-INDEX.md` ("Mathlib map").

## Files

All eight modules compile with **no `sorry`** and only the standard axioms
(`propext`, `Classical.choice`, `Quot.sound`); `Check.lean` prints them.

| File | What it links | Main theorems |
|---|---|---|
| `MathlibMap/CircleOperation.lean` | core `CircleOp.circK` ↔ Mathlib `PreQuasiregular` | `circ_eq_quasiregular_mul`, `circK_one_eq_quasiregular_mul`, `circK_scale_quasiregular_mul`, `Dial R k` monoid, `Dial.equivCircle : Dial R k ≃* PreQuasiregular R` (k a unit), `wayBack_iff_isQuasiregular` |
| `MathlibMap/Walls.lean` | walls of the dial ↔ `quasispectrum` | `dialAdd`, `HasWayBack`, `hasWayBack_iff_isQuasiregular`, **`not_hasWayBack_iff_mem_quasispectrum : ¬HasWayBack k a ↔ -k⁻¹ ∈ quasispectrum 𝕜 a`**, `mem_quasispectrum_iff_wall`, `not_hasWayBack_iff_mem_spectrum` (unital case) |
| `MathlibMap/Neumann.lean` | Mathlib's Neumann series ↔ "no wall near k = 0" | `hasWayBack_of_norm_mul_lt_one : ‖k‖‖a‖ < 1 → HasWayBack k a` (complete normed non-unital algebra) |
| `MathlibMap/Jacobson.lean` | `Ring.jacobson` ↔ "never meets a wall" | `wayBack_iff_isUnit`, **`mem_jacobson_iff_wayBack : x ∈ Ring.jacobson R ↔ ∀ k, ∃ b, circK k x b = 0`**, `wayBack_of_isNilpotent`, `wayBack_every_dial_field` |
| `MathlibMap/Lucas.lean` | core `Patterns.lucas` at Q = 1 ↔ `Polynomial.Chebyshev.S` | `lucas_one_eq_S` |
| `MathlibMap/LeeYangBand.lean` | **statement only**: the Lee–Yang p-ball rule | `LeeYang.fieldFn`, `LeeYang.Holds`, `LeeYang.Rule`, `LeeYang.OpenBand` (propositions, not proved) |
| `MathlibMap/LeeYangStep.lean` | **statement only** (DRAFT, Nov 26): the step identity's consequences, from shape-dial-survey section N3 | `LeeYang.StepDown`, `LeeYang.TopStrip`, `LeeYang.TopStripReduces` (propositions, not proved) |
| `MathlibMap/LeeYangP3.lean` | **statement only** (DRAFT, Nov 26): Lee–Yang at roundness 3, from shape-dial-survey section "N3, p = 3" (computer-assisted on paper) | `LeeYang.RuleAtThree`, `LeeYang.StripAtThree` (propositions, not proved) |

19 theorems plus the `Dial` monoid and its equivalence. Names follow the library's 1.0 names
(`CircleOp.circ`, `CircleOp.circK`). When the folder was added to the repo (Nov 26), its references to the
pre-1.0 names were ported, and its first three theorems were renamed from `multAdd_eq_circle`,
`wAdd_one_eq_circle` and `wAdd_scale_circle`.

`MathlibMap/LeeYangBand.lean` is a **statement only**: it defines the Lee–Yang p-ball rule (`Rule`)
and its open band (`OpenBand`) as propositions, using Mathlib's interval integral. It compiles, proves
nothing about them, and has no `sorry`. The proof work stays with the shape-dial thread.

## Build

```
lake exe cache get
lake build
lake env lean Check.lean
```

The lakefile pins Mathlib to tag `v4.37.0-rc1`, the same as `../../lean-library-mathlib`, and points
`FusedAlgebra` at `../../lean-library`. `lake-manifest.json` pins the same Mathlib commit as
`../../lean-library-mathlib`. In the session that wrote it, the Mathlib cache server was out
of reach, so the needed Mathlib modules (interval integrals included) were built from source at Mathlib master commit `348837d`
(Nov 25, 2026, toolchain `v4.37.0-rc1`). The folder is now built on the pinned tag by the repo's Mathlib CI
(`.github/workflows/mathlib-build.yml`, steps "Mathlib map"), which also fails on any `sorry`.

## How the links read

- **k = 1 is the circle operation.** `circK 1 a b = a + b + ab`, Mathlib's `x ∘ y = y + x + xy`
  (commutative case). Scaling by k (`circK_scale`) carries any unit dial to k = 1, so every unit
  dial gives the same monoid (`Dial.equivCircle`).
- **Walls are the quasispectrum.** An element `a` has a way back at dial k (a two-sided
  ∘ₖ-inverse) unless `−1/k` lies in `quasispectrum 𝕜 a`. Turning the dial, the walls sit at
  k = −1/λ for the λ in the quasispectrum. For a matrix that is det(1 + kA) = 0.
- **Neumann.** Mathlib's `IsQuasiregular.of_norm_lt_one` (the geometric series) gives: no wall
  while ‖k‖‖a‖ < 1.
- **Jacobson.** `Ring.mem_jacobson_iff` (x ∈ J ↔ 1 + xy is a unit for all y) becomes: x is in the
  Jacobson radical exactly when it has a way back at every dial. In ℤ/12 these are {0, 6}.
- **Lucas.** Our Lucas sequence at Q = 1 is the Vieta–Fibonacci polynomial `Chebyshev.S`, so
  S_n(2x) = U_n(x) ties our pattern to Chebyshev's second kind.

## Figures

`figures.py` makes `figures/f2..f5`: what is built on the circle operation, walls along the dial,
the Jacobson grid for ℤ/12, and the open-problems map.
