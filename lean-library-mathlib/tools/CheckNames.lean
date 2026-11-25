import FusedAlgebra
import Lean

/-!
# CheckNames

Reads every table in `port-notes/*.tsv` and checks that each `ported` row names a declaration that
exists in the port, and reports `mathlib` rows whose `new_name` is a single identifier that does
not resolve. Run from `lean-library-mathlib/` with `lake env lean tools/CheckNames.lean`.
-/

open Lean Elab Command

/-- Check the names in the port-notes tables against the environment. -/
elab "#check_port_names" : command => do
  let env ← getEnv
  let dir : System.FilePath := "port-notes"
  let mut ported := 0
  let mut badPorted : Array String := #[]
  let mut mathlibOk := 0
  let mut badMathlib : Array String := #[]
  for entry in (← dir.readDir) do
    unless entry.fileName.endsWith ".tsv" do continue
    let lines := (← IO.FS.readFile entry.path).splitOn "\n"
    for line in lines.drop 1 do
      let parts := line.splitOn "\t"
      if parts.length < 3 then continue
      let status := parts[1]!
      let new := parts[2]!.trimAscii.toString
      let isIdent := !new.isEmpty && new.all fun c => !c.isWhitespace && c != ',' && c != '(' && c != ';'
      if status == "ported" then
        ported := ported + 1
        unless isIdent && env.contains new.toName do
          badPorted := badPorted.push s!"{entry.fileName}: {parts[0]!} -> {new}"
      else if status == "mathlib" && isIdent then
        if env.contains new.toName then mathlibOk := mathlibOk + 1
        else badMathlib := badMathlib.push s!"{entry.fileName}: {parts[0]!} -> {new}"
  logInfo m!"{ported} ported rows, {badPorted.size} unresolved; {mathlibOk} mathlib names resolved, {badMathlib.size} single-identifier mathlib names unresolved"
  for b in badPorted do logInfo m!"PORTED? {b}"
  for b in badMathlib do logInfo m!"MATHLIB? {b}"

#check_port_names
