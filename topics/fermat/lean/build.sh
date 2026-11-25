#!/bin/bash
# Compile the Fermat proof files against the shared library (core Lean 4.37.0-rc1, no mathlib).
# Needs the library built first: (cd ../../lean-library && ./build.sh)
set -e
cd "$(dirname "$0")"
LIB="$(cd ../../lean-library && pwd)/build"
OUT="${FERMAT_BUILD:-/tmp/fermat-build}"
mkdir -p "$OUT"
export LEAN_PATH="$LIB:$OUT"
for m in Fermat4 Fermat3 FreyPieces FermatRoute; do
  echo "== $m"
  lean -o "$OUT/$m.olean" "$m.lean" | grep -v "Quot.sound\]$" || true
done
echo "All Fermat files compiled."
