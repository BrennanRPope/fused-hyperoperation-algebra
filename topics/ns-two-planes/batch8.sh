#!/bin/bash
# DRAFT. Blow-up hunt runs (hunt.py). Outputs runs/hunt/h_{ic}_n{N}_nu{nu}_th{theta}.json
run() { f=runs/hunt/h_$1_n$2_nu$3_th$4.json; [ -f $f ] || python3 hunt.py $1 $2 $3 $4 $5 $f; }
for N in 64 128; do for ic in tg rings klein; do run $ic $N 0 1 6; done; done            # Euler
for N in 64 128; do for ic in tg rings; do for th in 0 0.5 0.8 0.9 0.95 0.98; do run $ic $N 0 $th 6; done; done; done   # push-back dial
for N in 64 128; do run tg $N 0.005 1 10; run rings $N 0.005 1 4; run klein $N 0.005 1 6; done   # Navier-Stokes
for ic in tg rings klein; do run $ic 256 0 1 6; done
run tg 256 0.005 1 10; run rings 256 0.005 1 4; run klein 256 0.005 1 6
for th in 0.9 0.98; do run tg 256 0 $th 6; run rings 256 0 $th 6; done
