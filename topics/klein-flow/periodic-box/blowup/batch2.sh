#!/bin/bash
cd "$(dirname "$0")"
export GAMMA=8 DELTA=0.5 EVERY=0.1
for nu in 0.001 0.004; do for c in none t-g+; do
  python3 tubes.py 128 $nu 2.6 $c 0 runs/nu${nu}_${c}_s0.json 0.01 0.25 2 >> runs/nu.log 2>&1 &
done; done; wait
for s in 1 2; do for c in t-g+ t+g-; do
  python3 tubes.py 128 0.002 2.6 $c $s runs/n128_${c}_s${s}.json 0.01 0.25 2 >> runs/n128.log 2>&1 &
done; done; wait
# resolution check (appended): nu=0.001 pair alone at N=192
WORKERS=4 python3 tubes.py 192 0.001 1.9 none 0 runs/res192_none_s0.json 0.01 0.25 2 >> runs/nu.log 2>&1
