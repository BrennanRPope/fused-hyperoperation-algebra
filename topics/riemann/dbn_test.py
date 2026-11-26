import mpmath as mp, numpy as np, time
mp.mp.dps=45
def Phi(u):
    s=mp.mpf(0)
    for n in range(1,8):
        s+=(2*mp.pi**2*n**4*mp.e**(9*u)-3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u))
    return s
U=mp.mpf('1.3'); NN=360
xg,wg=mp.gauss_legendre if False else (None,None)
# Gauss-Legendre nodes in mpmath precision
import mpmath
nodes=mp.calculus.quadrature.GaussLegendre(mp.mp)
pts=nodes.calc_nodes(9, mp.mp.prec)   # degree 9 -> 3*2^8 = 768 points on [-1,1]
us=[(x+1)*U/2 for x,w in pts]; ws=[w*U/2 for x,w in pts]
print('nodes',len(us))
t0=time.time()
Ph=[Phi(u) for u in us]; print('phi',time.time()-t0)
def H(t,x):
    return mp.fsum(w*p*mp.e**(t*u*u)*mp.cos(x*u) for u,w,p in zip(us,ws,Ph))
for x in (0,10,28.27,60,100,140):
    a=H(0,x); b=mp.re(mp.mpf(1)/8*(lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))(mp.mpf(1)/2+1j*x/2))
    print(x, mp.nstr(a,12), mp.nstr(b,12), mp.nstr(abs(a-b)/abs(b),3))
print(time.time()-t0)
