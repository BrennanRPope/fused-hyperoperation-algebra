# DRAFT batch 4 (round 2, three planes). Run with `bash batch4.sh` from this folder.
for s in 0 1 2; do for m in two-tw two-un three-tw three-un; do python3 threeplanes.py 64 0.005 4 $m $s runs/p3_${m}_s${s}.json; done; done
