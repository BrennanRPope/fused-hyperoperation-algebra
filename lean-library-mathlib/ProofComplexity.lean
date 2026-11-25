import ProofComplexity.Framework
import ProofComplexity.FrameworkPlane
import ProofComplexity.Operation
import ProofComplexity.Tangent
import ProofComplexity.Standard
import ProofComplexity.Control

/-!
# Proof complexity

The 17 proofs from `topics/proof-complexity/` on main: six known results, each proved the standard
way, through the circle operation, and as a control with the same plan and the operation taken out.
Built here so the Mathlib build check covers them. `Measure.lean` and `Tree.lean` compute the
yardsticks and write `yard.jsonl` and `tree.tsv` when compiled, so they are kept beside these files
but left out of this build; `topics/proof-complexity/README.md` on main says how to rerun them.
-/
