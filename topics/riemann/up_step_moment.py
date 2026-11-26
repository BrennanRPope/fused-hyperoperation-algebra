import mpmath as mp, time
mp.mp.dps = 30
def xi(s): return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Xi(t): return xi(mp.mpf(1)/2+1j*t)
def Phi(u):
    u=mp.mpf(u); return mp.nsum(lambda n:(2*mp.pi**2*n**4*mp.e**(9*u/2)-3*mp.pi*n**2*mp.e**(5*u/2))*mp.e**(-mp.pi*n**2*mp.e**(2*u)),[1,mp.inf])
t0=time.time()
print('Phi(0)',Phi(0),'pi*Phi0',mp.pi*Phi(0))
M1 = mp.quad(lambda t: mp.re(Xi(t)), mp.linspace(0,120,61))
print('M1 = int_0^inf Xi', M1, time.time()-t0)
# F1 on real axis: min over grid
I=mp.mpf(0); prev=0; mn=10
xs=mp.linspace(0,100,201)
for a,b in zip(xs[:-1],xs[1:]):
    I+=mp.quad(lambda t: mp.re(Xi(t)),[a,b])
    mn=min(mn,I)
print('min of int_0^x Xi over x in (0,100]:', mn, 'final', I, time.time()-t0)
