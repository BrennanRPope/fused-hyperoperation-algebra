#!/bin/bash
# usage: batch.sh N nu T eps seed tag cls...
N=$1; NU=$2; T=$3; EPS=$4; SEED=$5; TAG=$6; shift 6
for c in "$@"; do
  GAMMA=8 DELTA=0.5 python3 tubes.py $N $NU $T $c $SEED runs/${TAG}_${c}_s${SEED}.json $EPS 0.25 2 >> runs/${TAG}.log 2>&1 &
done
wait
