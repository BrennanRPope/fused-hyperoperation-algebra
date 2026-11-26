import numpy as np, math, time
from scipy.optimize import minimize
def solve(s, L, N=128, box=8.0, w=0.12, eps=1e-3, psi0=None):
    h = box/N
    x = (np.arange(N)-N/2)*h
    X, Y = np.meshgrid(x, x, indexing='ij')
    rho = (np.exp(-((X+L/2)**2+Y**2)/(2*w*w)) - np.exp(-((X-L/2)**2+Y**2)/(2*w*w)))/(2*math.pi*w*w)
    # particular solution D0 = grad phi0, lap phi0 = rho (periodic, spectral) -> use FD-consistent forward differences
    kx = 2*math.pi*np.fft.fftfreq(N, d=h)
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    # forward-difference symbols
    sx = (np.exp(1j*KX*h)-1)/h; sy = (np.exp(1j*KY*h)-1)/h
    den = -(np.abs(sx)**2+np.abs(sy)**2); den[0,0]=1
    # div D with backward diff of forward-diff grad = discrete laplacian
    ph = np.fft.fft2(rho)/den; ph[0,0]=0
    phi0 = np.real(np.fft.ifft2(ph))
    Dx0 = (np.roll(phi0,-1,0)-phi0)/h; Dy0 = (np.roll(phi0,-1,1)-phi0)/h
    # divergence-free addition: Dx = d_y^- psi... use D += (dpsi/dy (backward-staggered), -dpsi/dx)
    def D(psi):
        dx = Dx0 + (psi - np.roll(psi,1,1))/h
        dy = Dy0 - (psi - np.roll(psi,1,0))/h
        return dx, dy
    def E(v):
        psi = v.reshape(N,N)
        dx, dy = D(psi)
        m2 = dx*dx+dy*dy+eps*eps
        e = np.sum(m2**(s/2))/s*h*h
        g = m2**(s/2-1)  # dE/d(dx) = g*dx
        gx, gy = g*dx*h*h, g*dy*h*h
        # adjoint
        G = (gx - np.roll(gx,-1,1))/h - (gy - np.roll(gy,-1,0))/h
        return e, G.ravel()
    v0 = np.zeros(N*N) if psi0 is None else psi0
    res = minimize(E, v0, jac=True, method='L-BFGS-B', options={'maxiter':20000,'maxfun':40000,'gtol':1e-10,'ftol':1e-14})
    dx, dy = D(res.x.reshape(N,N))
    return res.fun, np.hypot(dx,dy), x, res
if __name__ == '__main__':
    for s in (2.0, 1.5, 1.2, 1.0):
        t=time.time()
        Es=[]
        for L in (1.0, 1.5, 2.0, 3.0):
            e, Dm, x, res = solve(s, L)
            Es.append(e)
        Ls=np.array([1,1.5,2,3.]); Es=np.array(Es)
        if s>1:
            sl=np.polyfit(np.log(Ls), np.log(Es-Es[0]+1e-30) if False else np.log(Es), 1)[0]
        print(s, Es.round(5), 'dE/dL', np.diff(Es)/np.diff(Ls), round(time.time()-t,1))
