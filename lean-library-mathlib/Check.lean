import FusedAlgebra
import Lean

/-!
# Check

Every theorem of the mathlib port uses only Lean's three standard axioms (`propext`,
`Classical.choice`, `Quot.sound`) and no `sorry`. Run with `lake env lean Check.lean` after
`lake build`. It prints the number of theorems and fails if any uses another axiom.
-/

open Lean Elab Command

/-- Check every theorem declared in a `FusedAlgebra.*` module. -/
elab "#check_fused_library" : command => do
  let env ← getEnv
  let ok : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut count := 0
  let mut bad : Array Name := #[]
  for (name, info) in env.constants.toList do
    let some idx := env.getModuleIdxFor? name | continue
    let mod := env.header.moduleNames[idx.toNat]!
    unless (`FusedAlgebra).isPrefixOf mod do continue
    unless info matches .thmInfo _ do continue
    if name.isInternal then continue
    count := count + 1
    let axioms ← liftCoreM (collectAxioms name)
    unless axioms.all (ok.contains ·) do bad := bad.push name
  if bad.isEmpty then
    logInfo m!"{count} theorems, all using only the standard axioms, no sorry."
  else
    throwError m!"theorems using sorry or other axioms: {bad}"

#check_fused_library
