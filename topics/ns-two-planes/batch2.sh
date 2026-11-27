# DRAFT batch 2: grid check (same flow on 64 and 96), snapshots, Klein pair reading. Run with `bash batch2.sh`
python3 twoplanes.py 64 0.005 0.3 1 twisted 0 runs/regen_check64.json
python3 twoplanes.py 96 0.005 8 1 twisted 0 runs/n96_nu0.005_twisted_k1_s0.json
SNAP=0,0.5,1,2 python3 twoplanes.py 64 0.005 2.1 1 twisted 0 runs/snap_twisted.json
GAMMA=8 DELTA=0.5 SIG=0.3 WORKERS=4 python3 klein_planes.py 96 0.002 3 runs/klein_planes.json > klein_planes.log
