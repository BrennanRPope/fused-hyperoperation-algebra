#!/bin/bash
# Compile the mult-add axioms and re-derive the library from them.
# Needs `lean` 4.37.0-rc1 on PATH (core Lean only, no mathlib).
set -e
cd "$(dirname "$0")"
LIB="$(cd ../../lean-library && pwd)"
if [ ! -f "$LIB/build/Check.olean" ] && [ ! -f "$LIB/build/FusedAlgebra.olean" ]; then
  echo "== building the shared library first"; bash "$LIB/build.sh"
fi
export LEAN_PATH="$LIB/build:$PWD/build"
mkdir -p build/Axioms
for m in MultAddRing Independence FusedPeano Rederive; do
  echo "== Axioms.$m"
  lean -o "build/Axioms/$m.olean" "Axioms/$m.lean"
done
echo "== Axioms.Check"
lean Axioms/Check.lean
echo "All axiom modules compiled."
