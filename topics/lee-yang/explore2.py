from dialgas import *
import mpmath as mp, time, numpy as np, json, pickle
C=pickle.load(open('torus_counts.pkl','rb'))
ks=[8,3+2*2**0.5,4,3,2,1.5,1.0001,1,0.999,0.5,0.25,0.1,0.01,0,-0.25,-0.5,-1,-2,-3]
out={}
for L in (4,5,6):
    for k in ks:
        r=z_zeros(C[L],k,dps=80 if L==6 else 60)
        zs=[complex(x) for x in r]
        mods=[abs(x) for x in zs]
        args=[abs(np.angle(x)) for x in zs]
        out[(L,k)]=zs
        print(f"L={L} k={k:.4f} n={len(zs)} |z| in [{min(mods):.4g},{max(mods):.4g}]  k^-2={k**-2 if k else float('inf'):.4g}  min|arg|={min(args):.4f}  nearest +axis z={zs[int(np.argmin(args))]:.4f}", flush=True)
pickle.dump(out,open('torus_zeros.pkl','wb'))
