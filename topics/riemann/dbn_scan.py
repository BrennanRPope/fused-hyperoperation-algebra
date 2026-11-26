import mpmath as mp, numpy as np, time, pickle, sys
mp.mp.dps=45
def Phi(u):
    s=mp.mpf(0)
    for n in range(1,8):
        s+=(2*mp.pi**2*n**4*mp.e**(9*u)-3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u))
    return s
U=mp.mpf('1.3')
pts=mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec)
us=[(x+1)*U/2 for x,w in pts]; ws=[w*U/2 for x,w in pts]
Ph=[Phi(u) for u in us]
XS=[mp.mpf(i)/4 for i in range(0,601)]   # x in [0,150]
t0=time.time()
C=[[mp.cos(x*u) for u in us] for x in XS]
print('cos',time.time()-t0,flush=True)
def Hrow(t):
    A=[w*p*mp.e**(t*u*u) for u,w,p in zip(us,ws,Ph)]
    return [mp.fdot(A,c) for c in C]
TS=[float(a) for a in sys.argv[1:]] if len(sys.argv)>1 else [0.5,0.2,0,-0.5,-1,-2,-4]
out={}
for t in TS:
    row=Hrow(mp.mpf(t))
    s=[mp.sign(v) for v in row]
    ch=[float(XS[i]) for i in range(len(s)-1) if s[i]*s[i+1]<0]
    out[t]=dict(changes=ch, row=[float(mp.log10(abs(v)+mp.mpf(10)**-80))*float(mp.sign(v)) for v in row])
    print(t,len(ch),[round(c,2) for c in ch[:12]],time.time()-t0,flush=True)
pickle.dump(out,open('dbn_fine.pkl' if len(sys.argv)>2 else 'dbn_coarse.pkl','wb'))
