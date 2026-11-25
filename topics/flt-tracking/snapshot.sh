#!/usr/bin/env bash
# Snapshot the public Lean FLT project (upstream) and the futarchy-fi fork.
# Usage: snapshot.sh [workdir]   -> prints one JSON line per repo; append to snapshots.jsonl
set -eu
W="${1:-/tmp/flt-track}"; mkdir -p "$W"
for spec in "upstream https://github.com/ImperialCollegeLondon/FLT.git" "fork https://github.com/futarchy-fi/FLT.git"; do
  set -- $spec; name=$1; url=$2; d="$W/$name"
  rm -rf "$d"; git clone -q --depth 1 "$url" "$d"
  head=$(git -C "$d" log -1 --format=%h); date=$(git -C "$d" log -1 --format=%cs)
  tc=$(cat "$d/lean-toolchain"); ml=$(grep -A3 'name = "mathlib"' "$d/lakefile.toml" | grep rev | cut -d'"' -f2)
  axioms=$(grep -A2 "depends on axioms" "$d/FermatsLastTheorem.lean" | tr -d '\n' | sed 's/.*axioms:[[:space:]]*//; s/-\/.*//; s/[[:space:]]\+/ /g')
  files=$(find "$d/FLT" -name '*.lean' | wc -l)
  # count `sorry` tokens outside line comments
  sorries=$(find "$d/FLT" -name '*.lean' -print0 | xargs -0 sed 's/--.*//' | grep -ow "sorry" | wc -l || true)
  sfiles=$(find "$d/FLT" -name '*.lean' -print0 | xargs -0 grep -lw "sorry" | wc -l || true)
  assum=$(ls "$d/FLT/Assumptions" 2>/dev/null | grep -c '\.lean$' || true)
  printf '{"repo":"%s","head":"%s","date":"%s","toolchain":"%s","mathlib":"%s","flt_axioms":"%s","lean_files":%s,"sorry_tokens":%s,"files_with_sorry":%s,"assumption_files":%s}\n' \
    "$name" "$head" "$date" "$tc" "$ml" "$axioms" "$files" "$sorries" "$sfiles" "$assum"
done
