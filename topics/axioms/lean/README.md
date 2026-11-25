# Mult-add as an axiom: the foundation files

Core Lean 4.37.0-rc1, no mathlib, no `sorry`. `./build.sh` builds the shared
library (if needed) and then these files; `Axioms/Check.lean` confirms every
theorem uses only Lean's three standard axioms. Last clean build: Nov 25, 2026,
157 theorems written here plus 133 re-derived by `#rederive`.

## Files

| File | What it does | Imports the library? |
|---|---|---|
| `Axioms/MultAddRing.lean` | MA-ring axioms (abelian group + M1–M4), proof that every MA-ring is a commutative ring with a·b := a⊕b − a − b and 1 := −w, round trip both ways, MA-field (F1–F2) is a field | no |
| `Axioms/Independence.lean` | a counter-model for each of M1–M4 | no |
| `Axioms/FusedPeano.lean` | Peano with rung 1 and rung 2 as primitives; ⊕ laws by induction; + and × defined; categoricity (every model is a copy of Nat); recursion on any model; the fused ladder from one schema; F n a b + 1 = H n (a+1) (b+1) for every n; derived commutative semiring; + and × definable from S and ⊕ | no |
| `Axioms/Rederive.lean` | `#rederive`: re-derives every generic library theorem for an arbitrary MA-ring/MA-field (and the semiring ones for any fused-Peano model); restates 23 rung laws on any fused-Peano model | yes |
| `Axioms/Check.lean` | standard-axioms check, writes `build/axioms-theorems.tsv` | yes |

## Count against the library (658 theorems on Nov 25; 523 when this task began)

- 99 generic theorems (any commutative semiring, ring or field): all 99 re-derived from the MA axioms, each kernel-checked as `MAFound.Rederived.<name>`. The 34 semiring ones are also re-derived from the fused Peano axioms.
- 23 rung laws from `FusedAlgebra.Rungs` restated and proved on any model of the fused Peano axioms (`MAFound.Rederived.Rungs`).
- Total re-derived one by one: 122 of 658 (119 of the original 523).
- The other 141 Nat theorems hold in every model by the categoricity theorem, but are not restated one at a time.
- 395 theorems about Int, Rat (and 3 others) need the integers and fractions, built from the naturals by the usual constructions; the axioms here do not replace those constructions.

## For the library thread

`MultAddRing`, `Independence` and `FusedPeano` are standalone and could be folded in as `FusedAlgebra/Foundations/`. They declare global instances `[MARing R] → Add R, Mul R, …` and `[FusedPeano N] → Add N, …`; harmless unless a type carries both its own ring and one of these classes. `Rederive` imports the whole library, so it belongs after it, like `Check`.
