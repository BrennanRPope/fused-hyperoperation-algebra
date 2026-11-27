import FusedAlgebra
import ProofComplexity
import MathlibCandidates
import Lean

/-!
# Check

Every theorem of the mathlib port uses only Lean's three standard axioms (`propext`,
`Classical.choice`, `Quot.sound`) and no `sorry`. Run with `lake env lean Check.lean` after
`lake build`. It prints the number of theorems and fails if any uses another axiom.
-/

open Lean Elab Command

/-- Check every theorem declared in a module under the given prefix. -/
def checkLibrary (pfx : Name) : CommandElabM Unit := do
  let env ← getEnv
  let ok : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut count := 0
  let mut bad : Array Name := #[]
  for (name, info) in env.constants.toList do
    let some idx := env.getModuleIdxFor? name | continue
    let mod := env.header.moduleNames[idx.toNat]!
    unless pfx.isPrefixOf mod do continue
    unless info matches .thmInfo _ do continue
    if name.isInternal then continue
    count := count + 1
    let axioms ← liftCoreM (collectAxioms name)
    unless axioms.all (ok.contains ·) do bad := bad.push name
  if bad.isEmpty then
    logInfo m!"{pfx}: {count} theorems, all using only the standard axioms, no sorry."
  else
    throwError m!"{pfx}: theorems using sorry or other axioms: {bad}"

/-- Check the port's own library. -/
elab "#check_fused_library" : command => checkLibrary `FusedAlgebra

/-- Check the proof-complexity study. -/
elab "#check_proof_complexity" : command => checkLibrary `ProofComplexity

/-- Check the Mathlib candidates (the Tate curve). -/
elab "#check_mathlib_candidates" : command => checkLibrary `MathlibCandidates

#check_fused_library
#check_proof_complexity
#check_mathlib_candidates
