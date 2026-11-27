# DRAFT batch 3: more seeds at viscosity 0.005 for the twisted vs untwisted comparison. Run with `bash batch3.sh`
for s in 3 4 5 6 7 8; do for m in twisted untwisted; do python3 twoplanes.py 64 0.005 4 1 $m $s runs/n64_nu0.005_${m}_k1_s${s}.json; done; done
