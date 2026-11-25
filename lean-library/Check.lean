/-
  Check: every theorem in the library uses only Lean's three standard axioms
  (propext, Classical.choice, Quot.sound) and no `sorry`. Writes the list of
  theorems to build/theorems.tsv (module, name, axioms).
-/
import FusedAlgebra
import Lean

open Lean Elab Command

elab "#check_fused_library" : command => do
  let env ← getEnv
  let ok : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut rows : Array String := #[]
  let mut bad : Array Name := #[]
  for h : i in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[i]
    unless (`FusedAlgebra).isPrefixOf mod do continue
    for n in env.header.moduleData[i]!.constNames do
      let some ci := env.find? n | continue
      unless ci matches ConstantInfo.thmInfo _ do continue
      if n.isInternalDetail then continue
      -- skip lemmas Lean generates itself (equation, match and constructor lemmas):
      -- only declarations written in a source file have a source range
      let some _ ← liftCoreM (findDeclarationRanges? n) | continue
      let axs ← liftCoreM (collectAxioms n)
      if axs.any (fun a => !ok.contains a) then bad := bad.push n
      rows := rows.push s!"{mod}\t{n}\t{", ".intercalate (axs.toList.map toString)}"
  let sorted := rows.qsort (· < ·)
  IO.FS.writeFile "build/theorems.tsv" ("\n".intercalate sorted.toList ++ "\n")
  if bad.isEmpty then
    logInfo m!"{sorted.size} theorems, all using only the standard axioms, no sorry."
  else
    throwError m!"theorems using sorry or other axioms: {bad}"

#check_fused_library
