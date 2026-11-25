"""Test the t* ~ eps^-2 trend: run eps=0.05 to t=200 (prediction from eps=0.1..0.5: t* roughly 60-120)."""
import json
from multiprocessing import Pool
from shell_scan import job
if __name__=="__main__":
    with Pool(3) as p: res=p.map(job,[(0.05,30,s,200.0) for s in range(3)])
    json.dump(res,open("scan_eps005_long.json","w"),indent=1,default=str); print(res)
