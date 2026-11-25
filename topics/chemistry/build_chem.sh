#!/bin/bash
# Compile the chemistry Lean files against the shared library (read-only) and check axioms.
# Usage: ./build_chem.sh   (needs lean 4.37.0-rc1 on PATH; library built in ../lean-library/build)
set -e
cd "$(dirname "$0")"
LIB=../lean-library/build
OUT=${CHEM_BUILD:-/tmp/chembuild}
mkdir -p "$OUT/Chem"
export LEAN_PATH="$LIB:$OUT"
mods=()
for f in */lean/*.lean; do
  m=$(basename "$f" .lean)
  echo "== $m"
  lean -o "$OUT/Chem/$m.olean" "$f"
  mods+=("$m")
done
{
  for m in "${mods[@]}"; do echo "import Chem.$m"; done
  cat <<'L'
import Lean
open Lean Elab Command
elab "#check_chem" : command => do
  let env ← getEnv
  let ok : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut rows : Array String := #[]
  let mut bad : Array Name := #[]
  for h : i in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[i]
    unless (`Chem).isPrefixOf mod do continue
    for n in env.header.moduleData[i]!.constNames do
      let some ci := env.find? n | continue
      unless ci matches ConstantInfo.thmInfo _ do continue
      if n.isInternalDetail then continue
      let some _ ← liftCoreM (findDeclarationRanges? n) | continue
      let axs ← liftCoreM (collectAxioms n)
      if axs.any (fun a => !ok.contains a) then bad := bad.push n
      rows := rows.push s!"{mod}\t{n}\t{", ".intercalate (axs.toList.map toString)}"
  let sorted := rows.qsort (· < ·)
  IO.FS.writeFile "chem-theorems.tsv" ("\n".intercalate sorted.toList ++ "\n")
  if bad.isEmpty then logInfo m!"{sorted.size} chemistry theorems, standard axioms only, no sorry."
  else throwError m!"bad: {bad}"
#check_chem
L
} > "$OUT/CheckChem.lean"
lean "$OUT/CheckChem.lean"
cut -f1 chem-theorems.tsv | sort | uniq -c
