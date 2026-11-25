/-
  Extract: for every theorem in the shared library, write
    * the engine's picture of its statement (the Reader, run on the elaborated type),
    * the constants its compiled proof term uses (aux lemmas Lean generated are opened),
  one JSON object per line, to build/extract.jsonl. Read-only on the library.
-/
import FusedPictures
open Lean Elab Command Meta FusedPictures

partial def usedConsts (env : Environment) (n : Name) (seen : NameSet) (depth : Nat) : NameSet := Id.run do
  let some ci := env.find? n | return seen
  let some v := ci.value? (allowOpaque := true) | return seen
  let mut s := seen
  for c in v.getUsedConstants do
    if s.contains c then continue
    s := s.insert c
    -- open lemmas Lean made for this proof (auxiliary, match, proof_ terms)
    if depth < 3 && (c.isInternalDetail || c.isInternal) then
      s := usedConsts env c s (depth + 1)
  return s

elab "#extract_fused" : command => do
  let env ← getEnv
  let h ← IO.FS.Handle.mk "build/extract.jsonl" .write
  let mut count := 0
  for hi : i in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[i]
    unless (`FusedAlgebra).isPrefixOf mod do continue
    for n in env.header.moduleData[i]!.constNames do
      let some ci := env.find? n | continue
      unless ci matches ConstantInfo.thmInfo _ do continue
      if n.isInternalDetail then continue
      let some _ ← liftCoreM (findDeclarationRanges? n) | continue
      let shape ← liftTermElabM do
        try
          let sh ← readStatement ci.type
          pure (some sh)
        catch _ => pure none
      let stmt ← liftTermElabM do
        try pure (toString (← ppExpr ci.type)) catch _ => pure ""
      let cs := (usedConsts env n {} 0).toList.map toString
      let j := Json.mkObj [("name", toString n), ("module", toString mod),
        ("type", stmt),
        ("shape", match shape with | some s => s.toJson | none => Json.null),
        ("maxSides", match shape with | some s => toJson s.maxSides | none => Json.null),
        ("marks", match shape with | some s => toJson s.marks.eraseDups | none => Json.null),
        ("consts", toJson cs),
        ("typeConsts", toJson (ci.type.getUsedConstants.toList.map toString))]
      h.putStrLn j.compress
      count := count + 1
  logInfo m!"{count} theorems extracted"

#extract_fused
