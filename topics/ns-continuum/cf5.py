import time, sys
from corner_ifrk4 import outcome
for nu in [0.0, 1e-5]:
    t0 = time.time(); o = outcome(0.05, 0.2, 24, nu=nu, tmax=400.0)
    print(f"eps=0.05 alpha=0.2 nu={nu}: N=24 {'REACH' if o['reached'] else 'stop'} t={o['t']:.3f} deep={o['deep']} ({time.time()-t0:.0f}s)", flush=True)
