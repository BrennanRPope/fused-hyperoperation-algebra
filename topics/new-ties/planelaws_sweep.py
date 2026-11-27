import numpy as np, json
from fractions import Fraction as F
import planelaws as P
s1,s2,s3=F(1,5),F(3,20),F(1,4)
out=[]
for q in [F(3,2),F(2),F(3),F(5),F(8),F(12),F(20),F(40),F(100)]:
    r=P.run([(F(1),q,[s1,s2]),(F(-1),-q,[s3])],n=10,name=f"q={q}")
    out.append(dict(q=float(q),box=r["box"],Dn=(r["Dn_wall1"]+r["Dn_wall2"])/2,D_mult=r["D_mult"],D_sup=r["D_sup"]))
json.dump(out,open("planelaws_sweep.json","w"),indent=1)
