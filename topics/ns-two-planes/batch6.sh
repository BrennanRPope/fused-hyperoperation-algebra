# DRAFT batch 6 (round 3): two planes sharing the box diagonal; saved fields for the pictures. Run with `bash batch6.sh`.
for nu in 0.05 0.035; do python3 feed.py 64 $nu twod 0 runs/feed/f_n64_nu${nu}_twod_s0.json 1 600; done
for c in twod any two; do SAVE=1 python3 feed.py 64 0.025 $c 0 runs/feed/f_n64_nu0.025_${c}_s0.json 1 600; done
python3 feed.py 96 0.025 twod 0 runs/feed/f_n96_nu0.025_twod_s0.json 1 600
