import math, time
def Hf(k,b,x):
    if k==3: return b**x
    if x<=0: return 1+x
    return Hf(k-1,b,Hf(k,b,x-1))
for k in range(4,16):
    t=time.time(); v=Hf(k,1.9,3.5); print(k, v, time.time()-t)
