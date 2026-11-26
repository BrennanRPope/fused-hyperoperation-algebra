import mpmath as mp, numpy as np, time, pickle
from multiprocessing import Pool
mp.mp.dps = 20
def xi(s): return mp.mpf(1)/2*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Xi(t): return xi(mp.mpf(1)/2+1j*t)
N=400
xg,wg=np.polynomial.legendre.leggauss(N)
tau=(xg+1)/2; w=wg/2
def F(z,kappa=1):
    # F_kappa(z) = int_0^1 tau^(kappa-1) Xi(z tau) dtau
    return sum(mp.mpf(wi)*mp.mpf(ti)**(kappa-1)*Xi(z*mp.mpf(ti)) for ti,wi in zip(tau,w))
def cell(args):
    x,y,k=args
    v=F(mp.mpc(x,y),k); return (x,y,k,complex(v))
if __name__=='__main__':
    xs=np.arange(0,60.01,1.0); ys=np.arange(0,20.01,1.0)
    jobs=[(x,y,k) for k in (1,) for x in xs for y in ys]
    t=time.time()
    with Pool(4) as p: res=p.map(cell,jobs,chunksize=8)
    pickle.dump(res,open('up_step_grid.pkl','wb')); print('done',time.time()-t)
