#!/bin/bash
# usage: check.sh <subject-dir> <LeanFileStem>   (needs /tmp/lib built with Base + LeeYang)
set -e
export PATH=/opt/lean-4.37.0-rc1-linux/bin:$PATH LEAN_PATH=/tmp/lib/build
d=$1; n=$2; mod=FusedAlgebra.Topics.Biology.$n
cp "$d/lean/$n.lean" /tmp/lib/FusedAlgebra/Topics/Biology/
mkdir -p /tmp/lib/build/FusedAlgebra/Topics/Biology
(cd /tmp/lib && lean -o build/FusedAlgebra/Topics/Biology/$n.olean FusedAlgebra/Topics/Biology/$n.lean)
sed -e "s/MODULE/$mod/g" -e "s#OUT#$d/lean/theorems.tsv#" "$(dirname "$0")/CheckBio.lean.tmpl" > /tmp/CheckBio.lean
(cd /tmp && lean CheckBio.lean)
