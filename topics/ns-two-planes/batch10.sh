#!/bin/bash
# DRAFT note: this batch was stopped after "tg 256 0.95"; the runs listed after it (tg 256 0.98, rings 256, klein 256 at viscosity 0.005) were not done.
# DRAFT. Blow-up hunt with the faster solver (hunt2.py). Outputs runs/hunt2/h_{ic}_n{N}_nu{nu}_th{theta}.json
run() { f=runs/hunt2/h_$1_n$2_nu$3_th$4.json; [ -f $f ] || DIAG=0.025 python3 hunt2.py $1 $2 $3 $4 $5 $f; }
for N in 64 128; do for ic in tg rings; do for th in 0 0.5 0.8 0.9 0.95 0.98 1; do run $ic $N 0 $th 6; done; done; run klein $N 0 1 6; done
for N in 64 128; do run tg $N 0.005 1 10; run rings $N 0.005 1 4; run klein $N 0.005 1 6; done
for th in 0 0.5 0.8 0.9 1; do run tg 256 0 $th 6; done
run klein 256 0 1 6
for th in 0.95 0.98; do run tg 256 0 $th 6; done
run rings 256 0 1 6; run rings 256 0 0.9 6
run klein 256 0.005 1 6
