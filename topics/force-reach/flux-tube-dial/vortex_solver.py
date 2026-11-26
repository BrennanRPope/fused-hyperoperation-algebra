import numpy as np, math
from scipy.integrate import solve_bvp, simpson
def vortex(n, lam, R=30.0, N=3000):
    r = np.linspace(1e-5, R, N)
    # initial guess
    f0 = np.tanh(r/ (1.0+1/ math.sqrt(lam)))**n
    a0 = np.tanh(r/1.5)**2
    y0 = np.vstack([f0, np.gradient(f0, r), a0, np.gradient(a0, r)])
    def ode(r, y):
        f, fp, a, ap = y
        fpp = -fp/r + n*n*(1-a)**2*f/r**2 + lam/2*f*(f*f-1)
        app = ap/r - (1-a)*f*f
        return np.vstack([fp, fpp, ap, app])
    def bc(ya, yb):
        return np.array([ya[0], ya[2], yb[0]-1, yb[2]-1])
    sol = solve_bvp(ode, bc, r, y0, tol=1e-8, max_nodes=200000)
    rr = np.linspace(1e-5, R, 40001)
    f, fp, a, ap = sol.sol(rr)
    dens = rr*(0.5*fp**2 + 0.5*n*n*(1-a)**2*f**2/rr**2 + 0.5*(n*ap/rr)**2 + lam/8*(f*f-1)**2)
    return 2*math.pi*simpson(dens, x=rr), sol.status
if __name__ == '__main__':
    for lam in (0.25, 0.5, 1.0, 2.0, 4.0):
        E = [vortex(n, lam) for n in (1,2,3,4)]
        Es = [e[0]/math.pi for e in E]
        k = (Es[1]-2*Es[0])/Es[0]**2
        p3 = Es[0]+Es[1]+k*Es[0]*Es[1]; p4 = 2*Es[1]+k*Es[1]**2
        print(lam, [round(x,6) for x in Es], [e[1] for e in E], 'k=%.5f'%k, 'E3 pred %.5f act %.5f'%(p3,Es[2]), 'E4 pred %.5f act %.5f'%(p4,Es[3]))
