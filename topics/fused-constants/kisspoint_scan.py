# Where does the kiss happen for high rungs (linear extension)?  float scan, k up to 60
import sys
sys.setrecursionlimit(100000)
def Hf(k,b,x):
    if x>200: return float('inf')
    if k==3: return b**x
    if x<=0: return 1+x
    return Hf(k-1,b,Hf(k,b,x-1))
xs=[1+i/1000 for i in range(1,4000)]
for k in (8,12,16,20,25,30,40,50,60):
    lo,hi=1.5,2.0
    for _ in range(50):
        mid=(lo+hi)/2
        if min(Hf(k,mid,x)-x for x in xs)<=0: lo=mid
        else: hi=mid
    g=[(Hf(k,lo,x)-x,x) for x in xs]
    xm=min(g)[1]
    # second-best local minimum in a different unit interval
    other=min((v,x) for v,x in g if int(x)!=int(xm))
    gap2=min(Hf(k,2.0,x)-x for x in xs)
    arg2=min(xs,key=lambda x:Hf(k,2.0,x)-x)
    print(f"k={k:2d} c_k={lo:.12f} 2-c_k={2-lo:.3e} kiss x={xm:.3f}  next-closest {other[0]:.2e} at x={other[1]:.3f} | base 2: min gap {gap2:.4f} at x={arg2:.3f}")
