# FLT leads, checked in Lean with mathlib

Built Nov 27, 2026 against Lean and mathlib `v4.37.0-rc1` (the same pins as `lean-library-mathlib/`).
`lake build` succeeds with no `sorry`; `lake env lean Axioms.lean` shows every main theorem resting
only on Lean's three standard axioms (`propext`, `Classical.choice`, `Quot.sound`).

This is a separate Lake project, not part of `lean-library-mathlib/`, because it copies code from
the Lean FLT projects and keeps their credit lines together. Both sources are Apache 2.0.

## What is proved

| File | Result | Source |
|---|---|---|
| `FLTLeads/FreyMinimal.lean` | `isMinimal_of_valuation_c₄_eq_one`: an integral model with `v(c₄) = 0` is minimal | copied from upstream FLT (Buzzard, Coram, Claude) |
| | `frey_hasGood_or_hasMultiplicative`: over **every** DVR (every prime, 2 included) the Frey model `Y² + XY = X³ + a₂X² + a₄X` has good or multiplicative reduction | ours |
| | `frey_isMinimal`: so the model is minimal at every prime. This turns "minimal model" in `FreyPieces` from a name into a theorem | ours |
| | `frey_hasMultiplicative_iff`: away from 2, multiplicative exactly when `ABC` lies in the maximal ideal | ours |
| `FLTLeads/SplitNode.lean` | `frey_hasSplit_iff`: where `a₄` reduces to 0 (primes dividing `AB`), mathlib's split multiplicative reduction holds exactly when the node's slope polynomial `T² + T − a₂` has a root mod the prime, which is when the port's circle-operation chord law (`node_chord_multAdd`) applies | ours |
| | `frey_hasSplit_iff_isSquare`: away from 2, split exactly when `B − A` is a square mod the prime | ours |
| `FLTLeads/Torsion.lean` | `group_theory_lemma`: `d`-torsion of size `d^r` for all `d ∣ n` gives `n`-torsion `≅ (ZMod n)^r` | statement upstream FLT (still `sorry` there); proof copied from the futarchy-fi fork (krandder) |
| | `galoisRepresentation`: the Galois action on points is a `DistribMulAction` | statement upstream FLT (four laws still `sorry` there); proofs from the fork |

## Assumptions

- The Frey theorems take `A + B = C`, `gcd A B = 1`, `4a₂ = B − A − 1` and `16a₄ = −AB` as
  hypotheses. The normalisation `A ≡ 3 mod 4`, `16 ∣ B` that supplies such `a₂, a₄` is
  `frey_minimal_model` in the port; it is not repeated here.
- `R` is any discrete valuation ring with fraction field `K`; the model is pushed from ℤ to `R`
  by the canonical map. Nothing is assumed about the residue characteristic except where a
  theorem says `2` is a unit.
- `freyModel`, `frey_c4`, `frey_disc`, `frey_c4_not_dvd` and `frey_semistable` are repeated from
  `lean-library-mathlib/FusedAlgebra/Topics/Fermat/FreyPieces.lean` so this folder builds on
  mathlib alone.

## Known vs new

All of this is known mathematics. Minimality from a unit `c₄` is Silverman, *The Arithmetic of
Elliptic Curves*, VII.1.1 and the remark after it; semistability and the split criterion for the
Frey curve are standard (Silverman VII.5; Serre 1987; Diamond and Shurman ch. 9). The torsion
lemma is classical group theory. What is ours is the Lean: the Frey reduction results stated with
mathlib's reduction classes, the bridge from mathlib's split condition to the port's node law, and
the check that the fork's two proofs still build on a newer mathlib than either FLT project uses.

## Not done

- **Tate curve.** Now started, as a Mathlib candidate in
  [`lean-library-mathlib/MathlibCandidates/`](../../../lean-library-mathlib/MathlibCandidates/)
  rather than here, so CI builds it. It has the formal curve, its invariants and its split
  multiplicative reduction, but not the uniformisation, so the "circle-operation law as a
  q-deformation" lead stays a lead.
- **Upstream contribution.** `group_theory_lemma` and the Galois laws are still `sorry` in
  Buzzard's repository. Sending them upstream would mean a public PR from Brennan's account; not
  done, and the fork may get there first.

## Build

```
lake build
lake env lean Axioms.lean
```

The mathlib download cache (`lake exe cache get`) was unreachable from the machine that built
this, so mathlib was compiled from source (about 5,200 files, roughly an hour on 4 cores). On a
machine with the cache, `lake exe cache get` first makes it a few minutes.
