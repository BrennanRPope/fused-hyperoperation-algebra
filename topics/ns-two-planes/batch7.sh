# DRAFT batch 7 (round 4): fastest spin growth in 3D vs 4D at the same viscosity, on growing grids. `bash batch7.sh`
run() { [ -f runs/feedd/g_d$1_n$2_nu$3_s$4.json ] || python3 feedd.py $1 $2 $3 $4 runs/feedd/g_d$1_n$2_nu$3_s$4.json 1 ${5:-500}; }
for nu in 0.1 0.05; do for N in 16 24 32 48; do run 3 $N $nu 0; done; done
for nu in 0.1 0.05; do for N in 16 24; do run 4 $N $nu 0; done; done
run 4 16 0.1 1; run 4 16 0.05 1
SAVE=1 python3 feedd.py 4 16 0.05 0 runs/feedd/pic_d4_n16_nu0.05_s0.json 1 500
run 4 32 0.05 0 150      # capped climb: its value is a lower bound on the fastest rate at 32^4
