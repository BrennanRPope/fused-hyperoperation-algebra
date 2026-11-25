/-
  Axioms.Check: every theorem in the axioms modules uses only Lean's three
  standard axioms (propext, Classical.choice, Quot.sound) and no `sorry`,
  including the theorems `#rederive` added. Writes build/axioms-theorems.tsv.
-/
import Axioms.Rederive
import Axioms.Independence
import Lean

open Lean Elab Command

elab "#check_axioms_modules" : command => do
  let env ← getEnv
  let ok : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut rows : Array String := #[]
  let mut bad : Array Name := #[]
  let mut nOwn := 0
  let mut nRed := 0
  for h : i in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[i]
    unless (`Axioms).isPrefixOf mod do continue
    for n in env.header.moduleData[i]!.constNames do
      let some ci := env.find? n | continue
      unless ci matches ConstantInfo.thmInfo _ do continue
      if n.isInternalDetail then continue
      let generated := (`CircFound.Rederived).isPrefixOf n && !(`CircFound.Rederived.Rungs).isPrefixOf n
        || (`CircFound.RederivedPeano).isPrefixOf n
      unless generated do
        let some _ ← liftCoreM (findDeclarationRanges? n) | continue
      let axs ← liftCoreM (collectAxioms n)
      if axs.any (fun a => !ok.contains a) then bad := bad.push n
      if generated then nRed := nRed + 1 else nOwn := nOwn + 1
      rows := rows.push s!"{mod}\t{n}\t{", ".intercalate (axs.toList.map toString)}"
  let sorted := rows.qsort (· < ·)
  IO.FS.writeFile "build/axioms-theorems.tsv" ("\n".intercalate sorted.toList ++ "\n")
  if bad.isEmpty then
    logInfo m!"{nOwn} theorems written in the axioms files and {nRed} re-derived by #rederive, \
      all using only the standard axioms, no sorry."
  else
    throwError m!"theorems using sorry or other axioms: {bad}"

#check_axioms_modules

#check @CircFound.Rederived.CircleOp.mpow_lift
#check @CircFound.Rederived.CircleOp.circ_circInv
#check @CircFound.RederivedPeano.CircleOp.foldCirc_lift
