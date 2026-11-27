# DRAFT batch: run with `bash batch1.sh` from this folder
for nu in 0.01 0.005; do for s in 0 1 2; do
 for m in twisted untwisted A B; do python3 twoplanes.py 64 $nu 8 1 $m $s runs/n64_nu${nu}_${m}_k1_s${s}.json; done
done; done
for k in 0 0.1 0.25 0.5; do for m in twisted untwisted; do python3 twoplanes.py 64 0.005 8 $k $m 0 runs/n64_nu0.005_${m}_k${k}_s0.json; done; done
python3 twoplanes.py 128 0.005 8 1 twisted 0 runs/n128_nu0.005_twisted_k1_s0.json
python3 twoplanes.py 128 0.005 8 1 untwisted 0 runs/n128_nu0.005_untwisted_k1_s0.json
