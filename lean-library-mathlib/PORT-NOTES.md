# PORT-NOTES: the fused mult-add library in mathlib form

Folder `lean-library-mathlib/`, on main since Nov 25 (it was written on branch `mathlib-port`); main's CI builds it. This is a restatement of the core-Lean
library `lean-library/` (62 modules, 1,750 theorems, as of main `7a51396`, Nov 25, 2026) on
Lean 4.37.0-rc1 with mathlib v4.37.0-rc1. The core-Lean library is untouched.

## Status

- All 62 modules are ported. `lake build` completes with no errors and no warnings, with
  mathlib's standard linter set switched on in `lakefile.toml`.
- `lake exe runLinter FusedAlgebra` (mathlib's declaration linter) passes.
- `lake env lean Check.lean` reports 1,970 theorems, all using only `propext`,
  `Classical.choice` and `Quot.sound`, and no `sorry`.
- Every one of the 1,750 original theorem names (and every original definition) has a row in a
  name table under `port-notes/`. `python3 tools/port_notes.py` checks that coverage against
  `../lean-library/theorems.tsv`; `lake env lean tools/CheckNames.lean` checks that all 1,720
  "ported" rows name a declaration that really exists in the port.
- Nothing is marked "not ported".
- Also built here since Nov 25: the proof-complexity study (`ProofComplexity/`, 17 proofs of six
  known results, copied unchanged from `topics/proof-complexity/lean/` on main). It is a separate
  library in the same package, with its own options matching where it was first compiled, and
  `Check.lean` checks it the same way. It is not part of the port and has no name table.

These checks alone do not establish that each ported statement says exactly what the original
said. Where a statement was restated, generalised, strengthened or given an extra hypothesis,
the table row says so. A separate line-by-line review against the originals has since been
done; it found four helper lemmas that had lost content, now restored (see "Review fixes" at
the end).

## How to build

```
cd lean-library-mathlib
lake exe cache get        # once, downloads mathlib's prebuilt files
lake build
lake exe runLinter FusedAlgebra
lake env lean Check.lean
```

## The core: mult-add on mathlib's structures

`a ✳ b = a + b + a * b` is the circle operation of ring theory. Mathlib already has it as the
monoid structure on `PreQuasiregular R` (in `Mathlib/Algebra/Algebra/Spectrum/Quasispectrum`),
where `*` means mult-add and `1` means `0`. `FusedAlgebra/Basic.lean` builds on that:

| Idea | In the port |
|---|---|
| mult-add on the ring itself | `MultAdd.multAdd`, notation `a ✳ b`, scoped in `MultAdd` (use `open MultAdd`) |
| `(R, ✳, 0)` is a commutative monoid | `PreQuasiregular.instCommMonoid` (new instance; mathlib only had `Monoid`) |
| the lift `x ↦ 1 + x` carries `✳` to `*` | `MultAdd.oneAddHom : PreQuasiregular R →* R` |
| in a ring the lift is an isomorphism | `MultAdd.oneAddEquiv : PreQuasiregular R ≃* R` |
| the `✳`-invertible elements are a group | `MultAdd.unitsEquiv : (PreQuasiregular R)ˣ ≃* Rˣ` |
| "has a way back" | mathlib's `IsQuasiregular`; in a field it is `1 + u ≠ 0` |
| mult-add power, list fold | `mpow`, `multAddList`: powers and list products in `PreQuasiregular R` |
| the weighted family `a + b + k * a * b` | `wAdd`, `wpow` (plain definitions, commutative semirings) |
| the fused ladder | `MultAdd.rung n a b = hyperoperation n (a + 1) (b + 1) - 1` on mathlib's `hyperoperation` |
| cube layers, factorial | `Nat.choose`, `Nat.factorial` |

The `✳` notation is scoped, so it clashes with nothing. It has the same precedence as `+` (as in
`topics/tesseract/TesseractCoords.lean`), so write `1 + (a ✳ b)` with parentheses.

Lemma names follow mathlib: `multAdd_lift` became `one_add_multAdd`, `mpow_lift` became
`one_add_mpow`, `foldMA` became `multAddList`, and so on. The complete old-to-new tables for
the three core modules are `port-notes/FusedAlgebra.Base.tsv`, `…Rungs.tsv`, `…Cube.tsv`.

## Name tables

One table per original module: `port-notes/<Module>.tsv`, columns
`original / status / new_name / note`. Status is one of

- **ported**: restated in the port (same name, renamed, or generalised);
- **mathlib**: replaced by an existing mathlib declaration, named in `new_name`;
- **dropped**: redundant; `new_name` says what covers it.

Most modules also have `port-notes/<Module>.md` with the porter's remarks and judgment calls.

| Module | Original theorems | Rows | Ported | Mathlib | Dropped | Not ported |
|---|---:|---:|---:|---:|---:|---:|
| `Base` | 46 | 57 | 42 | 14 | 1 | 0 |
| `Cube` | 7 | 9 | 0 | 9 | 0 | 0 |
| `Foundations.FusedPeano` | 65 | 89 | 52 | 29 | 8 | 0 |
| `Foundations.Independence` | 7 | 27 | 11 | 13 | 3 | 0 |
| `Foundations.MultAddRing` | 51 | 78 | 41 | 34 | 3 | 0 |
| `Rungs` | 36 | 42 | 40 | 2 | 0 | 0 |
| `Topics.AbstractAlgebra` | 62 | 74 | 50 | 19 | 5 | 0 |
| `Topics.Calculus` | 29 | 36 | 24 | 8 | 4 | 0 |
| `Topics.Chemistry.ChemConservation` | 10 | 11 | 11 | 0 | 0 | 0 |
| `Topics.Chemistry.ChemPeriodic` | 12 | 16 | 16 | 0 | 0 | 0 |
| `Topics.Chemistry.ChemQuantum` | 8 | 10 | 9 | 1 | 0 | 0 |
| `Topics.ClosingRate` | 7 | 8 | 7 | 1 | 0 | 0 |
| `Topics.Consensus` | 27 | 29 | 27 | 0 | 2 | 0 |
| `Topics.DarkSector` | 17 | 17 | 14 | 0 | 3 | 0 |
| `Topics.DataCenter` | 16 | 17 | 15 | 2 | 0 | 0 |
| `Topics.DialSweep` | 28 | 34 | 29 | 5 | 0 | 0 |
| `Topics.E` | 37 | 44 | 25 | 19 | 0 | 0 |
| `Topics.EarlyScience` | 17 | 24 | 21 | 3 | 0 | 0 |
| `Topics.Elliptic` | 22 | 28 | 28 | 0 | 0 | 0 |
| `Topics.EllipticS3` | 17 | 19 | 17 | 1 | 1 | 0 |
| `Topics.EntropyNoise` | 41 | 49 | 42 | 7 | 0 | 0 |
| `Topics.Fermat.Fermat3` | 108 | 127 | 40 | 81 | 6 | 0 |
| `Topics.Fermat.Fermat4` | 9 | 9 | 5 | 3 | 1 | 0 |
| `Topics.Fermat.FermatRoute` | 9 | 10 | 5 | 4 | 1 | 0 |
| `Topics.Fermat.FreyPieces` | 46 | 64 | 41 | 20 | 3 | 0 |
| `Topics.Finance` | 16 | 19 | 13 | 6 | 0 | 0 |
| `Topics.FractionalRungs` | 8 | 9 | 8 | 0 | 1 | 0 |
| `Topics.FusedConstants` | 20 | 22 | 19 | 2 | 1 | 0 |
| `Topics.GeometryTrig` | 24 | 26 | 26 | 0 | 0 | 0 |
| `Topics.GravityProbeB` | 16 | 17 | 17 | 0 | 0 | 0 |
| `Topics.Infinities` | 12 | 19 | 11 | 5 | 3 | 0 |
| `Topics.InformationTheory` | 52 | 56 | 52 | 3 | 1 | 0 |
| `Topics.JetFields` | 46 | 78 | 51 | 25 | 2 | 0 |
| `Topics.Landscape` | 27 | 29 | 29 | 0 | 0 | 0 |
| `Topics.LeeYang` | 39 | 48 | 37 | 11 | 0 | 0 |
| `Topics.LeeYangCluster` | 10 | 10 | 9 | 1 | 0 | 0 |
| `Topics.LeeYangSpill` | 10 | 10 | 9 | 1 | 0 | 0 |
| `Topics.LinearAlgebra` | 51 | 85 | 46 | 38 | 1 | 0 |
| `Topics.NSContinuum` | 18 | 24 | 21 | 2 | 1 | 0 |
| `Topics.NavierStokes` | 6 | 11 | 10 | 1 | 0 | 0 |
| `Topics.NumberPlane` | 26 | 29 | 28 | 1 | 0 | 0 |
| `Topics.OsgoodInverse` | 53 | 69 | 49 | 18 | 2 | 0 |
| `Topics.Patterns` | 32 | 37 | 33 | 1 | 3 | 0 |
| `Topics.Pi` | 45 | 64 | 41 | 16 | 7 | 0 |
| `Topics.QAcross.Knots` | 19 | 50 | 49 | 1 | 0 | 0 |
| `Topics.QAcross.Partitions` | 13 | 30 | 29 | 1 | 0 | 0 |
| `Topics.QAcross.QCounting` | 13 | 21 | 20 | 1 | 0 | 0 |
| `Topics.QAcross.QDial` | 16 | 17 | 17 | 0 | 0 | 0 |
| `Topics.QuantumNoise` | 21 | 30 | 30 | 0 | 0 | 0 |
| `Topics.Relativity` | 30 | 38 | 36 | 2 | 0 | 0 |
| `Topics.Rulecraft` | 10 | 10 | 9 | 1 | 0 | 0 |
| `Topics.RungShape` | 17 | 19 | 16 | 3 | 0 | 0 |
| `Topics.ShapeDirections` | 36 | 44 | 42 | 1 | 1 | 0 |
| `Topics.SolitonsRungs` | 26 | 29 | 28 | 1 | 0 | 0 |
| `Topics.Statistics` | 98 | 137 | 126 | 8 | 3 | 0 |
| `Topics.Thermodynamics` | 68 | 73 | 70 | 2 | 1 | 0 |
| `Topics.Topology.CubeProofs` | 23 | 29 | 7 | 18 | 4 | 0 |
| `Topics.Topology.HopfInvariantOne` | 17 | 25 | 13 | 10 | 2 | 0 |
| `Topics.Topology.KummerWilson` | 30 | 35 | 20 | 14 | 1 | 0 |
| `Topics.Toric.ProjectivePlane` | 16 | 24 | 19 | 5 | 0 | 0 |
| `Topics.Toric.ToricCharts` | 16 | 24 | 24 | 0 | 0 | 0 |
| `Topics.TwistedBoosts` | 36 | 58 | 44 | 13 | 1 | 0 |
| **Total** | 1750 | 2283 | 1720 | 487 | 76 | 0 |

"Rows" is larger than "Original theorems" because the tables also list definitions, structures
and instances. The port has more theorems than "Ported" because of new bridge lemmas and
generalisations.

## What was replaced by mathlib

The larger replacements (details in each module's table):

- **Base**: `pw` and its laws are `^` with `pow_add`, `pow_mul`, `mul_pow`, `one_pow`;
  `prodL`, `sumL`, `groups` are `List.prod`, `List.sum`, `List.sublists'`.
- **Cube**: `layer` is `Nat.choose`, `fact` is `Nat.factorial`; the whole module's seven
  theorems are mathlib lemmas. The port's `Cube.lean` keeps only the mult-add statements
  (`one_add_mpow_eq_sum_choose`, `card_corners`).
- **Rungs**: `tet`, `pent` and the five named rungs sit on mathlib's `hyperoperation`.
- **Foundations**: a mult-add ring now yields real mathlib `CommRing` / `Field` instances, with
  `MARing R ≃ CommRing R`; a fused-Peano model gets `N ≃+* ℕ`.
- **Topology (CubeProofs, KummerWilson, HopfInvariantOne)**: `IsPrime` is `Nat.Prime`, `MEq` is
  `Int.ModEq`, `bsum` is `Finset.sum`. Lucas, Fermat's little theorem, the freshman's dream,
  Wilson and Kummer are corollaries of mathlib (`Choose.lucas_theorem`, `ZMod.wilsons_lemma`,
  `add_pow_char`, mathlib's Kummer). The cube-route proof of Wilson is no longer the proof; its
  ingredients are kept and the route is described in a docstring.
- **Fermat3 / Fermat4**: the headline theorems are corollaries of mathlib's
  `fermatLastTheoremThree` and `fermatLastTheoremFour`. The hand-written Euler descent
  (about 190 lines) is dropped. The Eisenstein integers are `QuadraticAlgebra ℤ (-1) (-1)`,
  with a new `EuclideanDomain` instance. Euler's cube lemma `key_cube`, which is not in
  mathlib, is kept with its proof.
- **FreyPieces / FermatRoute / Elliptic**: curves are mathlib `WeierstrassCurve`s; `FLT n` is
  mathlib's `FermatLastTheoremFor n`.
- **LinearAlgebra, DialSweep, LeeYang, Relativity, TwistedBoosts, EarlyScience**: hand-rolled
  2×2 and 3×3 records are `Matrix (Fin n) (Fin n) R`, with `Matrix.det`, `Matrix.trace`,
  `Matrix.charpoly`. TwistedBoosts moved from pairs of rationals to `ℂ`.
- **Pi, NumberPlane, GeometryTrig**: Gaussian integers are `GaussianInt`; the number plane is
  bridged to `QuadraticAlgebra` and to `ℂ`.
- **Calculus, JetFields**: finite differences are mathlib's `fwdDiff`; the weight-0 dual number
  type is identified with `DualNumber R`.
- **Statistics, E, Infinities, Chemistry, QDial**: hand-written list recursions are
  `Finset.sum` / `List.sum`; indicators are `IsIdempotentElem`; `fall` is `Nat.descFactorial`.
- **AbstractAlgebra**: the hand-written "laws" records are mathlib instances (`CommMonoid`,
  `CommRing`, `MulEquiv`, `RingEquiv`, `Con`).

Most statements over `ℚ` or `ℤ` were generalised to commutative rings, fields or linearly
ordered fields where the proof was unchanged. Numeric examples stay over `ℚ`, `ℤ` or `ℕ`.

## What did not fit, and other things to know

Honesty first: nothing in the port proves more deep mathematics than the original did, and
where mathlib supplies a classical theorem the docstring and table say so.

- **Fermat.** `FermatRoute` still proves only "no Frey curve ⇒ Fermat's Last Theorem"; the
  modularity step remains an explicit hypothesis. The new converse `noFreyCurve_of_flt` shows
  that hypothesis is equivalent to FLT itself, so the reduction is bookkeeping. `genus_X0_2`
  still takes the genus formula as a hypothesis, and "minimal model" is still only a name.
- **Partitions and q-series stay finite checks.** The pentagonal, triple-product, four-squares
  and Ramanujan statements are kernel computations on truncated coefficient lists, as in the
  original. `partList` is not connected to mathlib's `Nat.Partition`. This is the biggest gap
  to real mathlib form.
- **Knot and q-polynomials stay coefficient lists**, because mathlib's `Polynomial` is not
  computable and the invariants are computed by `decide +kernel`. A bridge
  `QKnots.toPolynomial : Poly → ℤ[X]` is provided. Mathlib has no Gaussian binomials, so `gauss`
  stays a recursive definition.
- **Statistics stays on finite weighted lists**, not mathlib's `PMF` or measure theory: the
  theorems are ring identities with no positivity or normalisation assumed.
- **Toric.** No schemes are built. `P1 K` is `OnePoint K` with an exact equivalence to
  `ℙ K (Fin 2 → K)`; for the plane there is a map to `ℙ K (Fin 3 → K)` but it is not proved
  bijective.
- **Velocity addition is not packaged as a group** (it is a group law only on `(-1, 1)`).
- **Pi.** The Machin intervals are not proved to contain `Real.pi`.
- **Statements that changed shape** (all flagged in the tables): the four independence
  theorems are now 4-part statements over mathlib `AddCommGroup` models instead of 8-part
  conjunctions; several hypotheses written as `x * u = 1` are now `u⁻¹`; `x * x` is `x ^ 2`;
  `% 2 = 1` is `Odd`; `Int.gcd a b = 1` is `IsCoprime a b` in Fermat3.
- **Hypotheses added by generalising** (automatic in the original's type):
  `FusedRate.grow`, `DarkSector.add_le_multAdd`, `DarkSector.sum_le_multAddList` need
  nonnegativity; some field statements need `[NeZero (2 : K)]` or `[CharZero K]`;
  `ChemQuantum.cisd_not_extensive` and `SimpleLaws.constant_product` need `NoZeroDivisors`.
- **Hypotheses removed** (statements got stronger) are listed per row, for example in
  Thermodynamics, Relativity, Statistics (`bayes_odds`, `bp_layers`) and
  `ShapeDirections.periods_outnumber_curves` (now all `g ≥ 2`).
- **`decide +kernel`** is still used for the large numeric checks in Pi, Knots, Partitions and
  QCounting, as in the original. It is kernel evaluation and adds no axioms.
- **Import graph changed.** `EntropyNoise` now imports `Finance`; `JetFields` imports
  `LinearAlgebra`; `Independence` imports `MultAddRing`; `NSContinuum`, `OsgoodInverse` and
  `FreyPieces` no longer import `CubeProofs` / `Pi` because what they used is now mathlib.
- **Root-namespace names.** `fermat4`, `no_x4_y4_z2` and `fourth_power_ne_mpow_four` keep their
  original, non-mathlib-style names in the root namespace.
- **No copyright headers.** The repo has no LICENSE file, so `linter.style.header` is switched
  off in `lakefile.toml` rather than inventing a licence line.
- **Module system.** Files use Lean's `module` / `public import` / `@[expose] public section`
  form, as mathlib does; without it every file warns when importing `Mathlib`.
- **Whole-`Mathlib` import.** `Basic.lean` imports all of `Mathlib` for simplicity, not the
  minimal set mathlib itself would require.
- **Namespaces are the originals'** (`FinanceMultAdd`, `SimpleLaws`, `M2`, …), kept so names map
  one to one. They are not a single mathlib-style hierarchy.

## Candidates to move into `Basic.lean`

Porters left these in topic files because the core was frozen during the port:

- a noncommutative weighted mult-add (`LinearAlgebraMA.wmadd`) and the one-field type
  `AbstractAlgebra.Weighted k` with its `CommMonoid` instance (`Weighted 1` is
  `PreQuasiregular`);
- `wpow k 0 n = 0`, `wAdd_neg_div`, `multAddList_add_one`, `add_le_multAdd`,
  `multAddList_nonneg`, `m2_multAdd` / `m3_multAdd` / `m4_multAdd`, `wAdd_axioms`;
- three small `rung` lemmas for `Rungs.lean`: `rung_zero`, `rung_one_zero_right`,
  `rung_add_two_zero_right`;
- the `EuclideanDomain` and `IsDomain` instances for the Eisenstein integers, if another module
  needs them.

## How the port was done

`Basic`, `Rungs`, `Cube` and `Finance` were ported first as the reference. The other 58
modules were ported in parallel in 13 dependency-ordered groups following `PORTING-GUIDE.md`,
one commit per module, each gated on compiling with no output. The final checks listed under
Status were then run on the whole library.

## Review fixes (Nov 25, 2026)

A line-by-line check of the port at `95bf1f4` against the core-Lean originals found four helper
lemmas that had lost content. All four are restored, each saying at least what the original
said:

| Original | Problem found | Now |
|---|---|---|
| `FusedPeano.add_right_cancel'` | missing; the table pointed to mathlib's `add_right_cancel`, which did not apply because no cancellation instance existed on the models | `MAFound.FusedPeano.add_right_cancel`, proved for every model, plus `add_left_cancel` and an `IsCancelAdd N` instance, so mathlib's lemma now applies as well |
| `Fermat3.measure_lt` | missing (bookkeeping for the exponent-3 descent) | `Fermat3.measure_lt`, original statement. The descent itself is still not ported |
| `FreyPieces.pow_dvd_bound` | weaker: needed `q` prime | `FreyPieces.pow_dvd_bound` with the original hypothesis `2 ≤ q` |
| `OsgoodInverse.const_block` | weaker: only the block starting at `2 ^ j` | `MultAdd.OsgoodInverse.const_block` for any start `s`; `blk_const` is now its special case |

The same review asked for `#check` on seven mathlib names cited in the name tables. All seven
exist in the pinned mathlib (v4.37.0-rc1), so no table entry needed changing:
`Nat.add_one_mul_choose_eq`, `Nat.multiplicity_mul`, `Int.ModEq.sum`, `Nat.two_pow_succ`,
`IsIdempotentElem.iff_eq_zero_or_one`, `neg_one_pow_eq_ite`, `fwdDiff_iter_comp_add`.
One caveat on `Nat.multiplicity_mul`: it takes `[Fact p.Prime]` and `a, b ≠ 0`, as the
KummerWilson table row already says.

While re-checking names, the automatic name check (`tools/CheckNames.lean`) turned up five more
table entries, not on the review's list, whose mathlib name did not resolve. Four were real
errors: mathlib v4.37 no longer has `Monoid.npow`, `Monoid.Pow`, `DivInvMonoid.zpow` or
`AddMonoid.nsmul` (the operations now live in the parent classes). They are corrected to
`NPow.npow`, `Monoid.toNPow`, `ZPow.zpow` and `NSMul.nsmul` in the FusedPeano, MultAddRing,
OsgoodInverse, Pi, Base, Fermat3 and LeeYang tables. The fifth, `mul_zero`, exists but under its full name
`MulZeroClass.mul_zero`, which the table now uses. After these fixes every single-identifier
mathlib name in the tables resolves.
