from numseries import *
import json,sys
mp.mp.dps=70
N=int(sys.argv[1]); out={}
for mv in ['0','0.2','0.5','0.8']:
    G=slope(mp.mpf(mv),N); V=julia(G,N)
    out[mv]=[(n,float(mp.log10(abs(V[n])))) for n in range(5,N,2) if V[n]!=0]
    print(mv,[ (n,round(l,2)) for n,l in out[mv][::10]],flush=True)
json.dump(out,open('julia_growth.json','w'))
