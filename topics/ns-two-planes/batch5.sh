# DRAFT batch 5 (round 3): fastest feeding at fixed total spin Z = 1. Run with `bash batch5.sh` from this folder.
# 32^3 keeps modes |k_i| <= 8, 64^3 keeps <= 16, 96^3 keeps <= 24. The 32^3 runs at viscosity 0.025 and 0.0125 turned out
# under-resolved (the 64^3 optimum sits at higher wavenumber), so the scaling uses 64^3 with a 96^3 check.
run() { [ -f runs/feed/f_n$1_nu$2_$3_s$4.json ] || python3 feed.py $1 $2 $3 $4 runs/feed/f_n$1_nu$2_$3_s$4.json 1 600; }
for nu in 0.1 0.07 0.05 0.035 0.025 0.0125; do for c in any two one; do for s in 0 1 2; do run 32 $nu $c $s; done; done; done
for nu in 0.07 0.05 0.035; do for c in any two; do run 64 $nu $c 0; done; done
for c in any two; do for s in 0 1; do run 64 0.025 $c $s; done; done
run 64 0.025 one 0
for c in any two; do run 96 0.025 $c 0; done
