import pickle, numpy as np, mpmath as mp, json
from up_step_grid import F
mp.mp.dps=20
res=pickle.load(open('up_step_grid.pkl','rb'))+pickle.load(open('up_step_grid_b.pkl','rb'))
out={}
for k in (1,2):
    R=[r for r in res if r[2]==k]
    xs=sorted(set(r[0] for r in R)); ys=sorted(set(r[1] for r in R))
    A=np.full((len(xs),len(ys)),np.nan,dtype=complex)
    for x,y,kk,v in R: A[xs.index(x),ys.index(y)]=v
    L=np.log(np.abs(A)+1e-300)
    cands=[]
    for i in range(1,len(xs)-1):
        for j in range(0,len(ys)-1):
            nb=[L[i+a,j+b] for a in (-1,0,1) for b in (-1,0,1) if 0<=j+b<len(ys) and (a,b)!=(0,0)]
            if L[i,j]<min(nb): cands.append((xs[i],ys[j]))
    zs=[]
    for x,y in cands:
        try:
            z=mp.findroot(lambda z:F(z,k),mp.mpc(x,y)); 
            if all(abs(complex(z)-w)>1e-6 for w in zs): zs.append(complex(z))
        except Exception as e: print('fail',k,x,y)
    zs.sort(key=lambda z:z.real)
    print(k,[ (round(z.real,3),round(z.imag,3)) for z in zs])
    out[k]=dict(zeros=[[z.real,z.imag] for z in zs], grid_x=xs, grid_y=ys, logabs=L.real.tolist(), real_axis_min=float(np.abs(A[:,0]).min()))
pickle.dump(out,open('up_zeros.pkl','wb'))
import json
json.dump({str(k): v for k, v in out.items()}, open('up_step.json', 'w'))
