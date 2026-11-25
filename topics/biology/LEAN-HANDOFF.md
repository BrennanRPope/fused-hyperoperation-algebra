# Biology Lean files for the library

> History note (Nov 25, 2026): all 11 biology files are now in the shared library under `lean-library/FusedAlgebra/Topics/Biology/`, renamed for 1.0 (`circ`, `circK`, notation `◯`). The topic copies in each chapter's `lean/` folder are the library versions. The table below is the hand-off as it was written.

All 11 files are Lean 4.37 core only, with no mathlib and no sorry. Each was compiled against
/tmp/lib (FusedAlgebra.Base, plus LeeYang for molecular biology) with tools/check.sh, and only
standard axioms appear. Each file has its own namespace and uses Brennan's ∘₁ notation
(local infixl).

Intended library path: FusedAlgebra/Topics/Biology/<File>.lean

| # | File | Theorems | Imports |
|---|------|----------|---------|
| 01 | 01-molecular/lean/MolecularBiology.lean | 28 | Base, Topics.LeeYang |
| 02 | 02-metabolism/lean/Metabolism.lean | 19 | Base |
| 03 | 03-cell/lean/CellBiology.lean | 13 | Base |
| 04 | 04-genetics/lean/Genetics.lean | 16 | Base |
| 05 | 05-development/lean/Development.lean | 8 | Base |
| 06 | 06-physiology/lean/Physiology.lean | 11 | Base |
| 07 | 07-neuroscience/lean/NeuroScience.lean | 17 | Base |
| 08 | 08-evolution/lean/Evolution.lean | 18 | Base |
| 09 | 09-population/lean/Population.lean | 12 | Base |
| 10 | 10-ecology/lean/Ecology.lean | 10 | Base |
| 11 | 11-behaviour/lean/Behaviour.lean | 12 | Base |

Total: 164 theorems. Each folder's lean/theorems.tsv lists the theorems with their axioms.

Highlights:
- NeuroScience.twoport_lens: exact lens law for two synapses on any passive network.
- Evolution.path_consistent_iff and consistent_forces_dial: linear global epistasis is
  order-consistent only on the dial.
- Population.bh_compose: two Beverton–Holt seasons are one.
- Genetics.haldane_unlinked and lens_unlinked: the wall at r = ½.
