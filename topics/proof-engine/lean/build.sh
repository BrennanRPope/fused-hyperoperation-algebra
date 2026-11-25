#!/bin/bash
# Build the fused proof engine against the shared library (core Lean, no lake, no mathlib).
# Usage: bash build.sh   (expects `lean` on PATH and the library built by ../../lean-library/build.sh)
set -e
cd "$(dirname "$0")"
mkdir -p build
export LEAN_PATH="$PWD/../../lean-library/build:$PWD/build"
lean -o build/FusedPictures.olean FusedPictures.lean
lean Demo.lean > build/demo_trace.txt          # the eight proofs; their trace feeds film.py
lean Extract.lean                              # writes build/extract.jsonl for the gallery
echo "built; now run from proof-engine/: python3 -I film.py"
