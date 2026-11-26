"""Round 6 lattice engine: the lattice zeta of a face-power gauge,
    Z(s) = sum_{m in L, m != 0} g(m)^(-2s),   g(x) = M_p(|<n_i, x>|)  (power mean over the facet normals n_i),
continued into the critical strip by Riemann's Mellin split:
    pi^-s Gamma(s) Z(s) = int_1^inf (Theta-1) x^(s-1) dx + int_0^1 (Theta - A/x) x^(s-1) dx + A/(s-1) - 1/s,
    Theta(x) = sum_m exp(-pi x g(m)^2),  A = area{g <= 1} / (pi * covolume).
Theta is summed directly on [x0, inf); below x0 the remainder Theta - A/x is fitted by a few powers of x and integrated exactly.
Square lattice with the square's normals is the roundness dial (l_p ball); the hexagonal lattice with the hexagon's normals is
the face-power dial on the next polytope."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import loggamma

SQRT3 = np.sqrt(3.0)

def square_setup():
    B = np.array([[1.0, 0.0], [0.0, 1.0]])
    N = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    return B, N

def hex_setup():
    B = np.array([[1.0, 0.0], [0.5, SQRT3/2]])
    ang = np.deg2rad(30 + 60*np.arange(6))
    N = (2/SQRT3)*np.stack([np.cos(ang), np.sin(ang)], 1)
    return B, N

def gauge(X, N, p):
    """X: (n,2) points; returns g(x) = (mean_i |<n_i,x>|^p)^(1/p), p = inf gives the polytope's own norm."""
    proj = np.abs(X @ N.T)
    if np.isinf(p):
        return proj.max(1)
    mx = proj.max(1, keepdims=True); mx[mx == 0] = 1
    return mx[:, 0]*np.mean((proj/mx)**p, 1)**(1/p)

class LatticeZeta:
    def __init__(self, B, N, p, x0=2e-3, nq_low=1000, nq_high=500, cut=46.0, theta=1.25, fit_pows=None, fit_span=12.0):
        self.B, self.N, self.p = B, N, p
        covol = abs(np.linalg.det(B))
        th = np.linspace(0, 2*np.pi, 20001)[:-1]
        G = gauge(np.stack([np.cos(th), np.sin(th)], 1), N, p)
        area = 0.5*np.mean(G**-2.0)*2*np.pi
        self.A = area/(np.pi*covol)
        # lattice points
        self.theta = theta; ct = np.cos(theta); eith = np.exp(1j*theta)
        gmax2 = cut/(np.pi*x0*ct)
        # bound coefficients: |a b1 + c b2| >= lam_min * |(a,c)|; g >= |x| * min G
        lam = np.min(np.linalg.svd(B, compute_uv=False)); R = np.sqrt(gmax2)/(G.min()*lam) + 2
        M = int(np.ceil(R))
        a, c = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1))
        pts = np.stack([a.ravel(), c.ravel()], 1) @ B
        g2 = gauge(pts, N, p)**2
        g2 = g2[(g2 > 1e-12) & (g2 < gmax2)]
        self.g2 = np.sort(g2)
        self.gmin2 = self.g2[0]
        # quadrature nodes, low part: x in [x0, 1] in u = ln x
        z, w = leggauss(nq_low)
        u0, u1 = np.log(x0), 0.0
        self.uL = 0.5*(u1 - u0)*z + 0.5*(u1 + u0); self.wL = 0.5*(u1 - u0)*w
        xL = np.exp(self.uL)*eith
        thL = np.array([np.exp(-np.pi*x*self.g2).sum() for x in xL]) + 1.0
        self.RL = thL - self.A/xL                      # Theta - A/x
        # high part: x in [1, X]
        X = cut/(np.pi*self.gmin2*ct) + 1
        z, w = leggauss(nq_high)
        u0, u1 = 0.0, np.log(X)
        self.uH = 0.5*(u1 - u0)*z + 0.5*(u1 + u0); self.wH = 0.5*(u1 - u0)*w
        xH = np.exp(self.uH)*eith
        self.TH = np.array([np.exp(-np.pi*x*self.g2).sum() for x in xH])  # Theta - 1
        # small-x fit of Theta - A/x on [x0, 12 x0]
        xs = np.exp(np.linspace(np.log(x0), np.log(fit_span*x0), 80))*eith
        Rs = np.array([np.exp(-np.pi*x*self.g2).sum() for x in xs]) + 1.0 - self.A/xs
        self.fit_pows = np.array([0.0, 0.5, 1.0, 1.5, 2.0] if fit_pows is None else sorted(set(fit_pows)), float)
        Mx = xs[:, None]**self.fit_pows[None, :]
        self.fit_coef, *_ = np.linalg.lstsq(Mx.astype(complex), Rs, rcond=None)
        self.fit_resid = float(np.max(np.abs(Mx @ self.fit_coef - Rs)))
        self.x0 = x0*eith

    def completed(self, s):
        s = np.atleast_1d(np.asarray(s, complex))
        out = np.empty(s.shape, complex)
        for i0 in range(0, len(s), 2000):
            ss = s[i0:i0 + 2000][:, None]
            th = 1j*self.theta; s1 = ss[:, 0]
            IH = (np.exp(ss*(self.uH[None, :] + th))*(self.TH*self.wH)[None, :]).sum(1)
            IL = (np.exp(ss*(self.uL[None, :] + th))*(self.RL*self.wL)[None, :]).sum(1)
            lx0 = np.log(self.x0)
            tail = sum(c*np.exp((s1 + a)*lx0)/(s1 + a) for a, c in zip(self.fit_pows, self.fit_coef))
            out[i0:i0 + 2000] = IH + IL + tail + self.A*np.exp(th*(s1 - 1))/(s1 - 1) - np.exp(th*s1)/s1
        return out

    def __call__(self, s):
        s = np.atleast_1d(np.asarray(s, complex))
        return self.completed(s)*np.exp(s*np.log(np.pi) - loggamma(s))

def find_zeros(Zf, sig=(0.15, 1.3), tt=(2.0, 30.0), ds=0.01, dt=0.04, tol=1e-10):
    S = np.arange(sig[0], sig[1] + 1e-9, ds); T = np.arange(tt[0], tt[1] + 1e-9, dt)
    SS, TT = np.meshgrid(S, T)
    V = np.abs(Zf((SS + 1j*TT).ravel())).reshape(SS.shape)
    cands = []
    for i in range(1, V.shape[0] - 1):
        for j in range(1, V.shape[1] - 1):
            w = V[i-1:i+2, j-1:j+2]
            if V[i, j] == w.min() and V[i, j] < 0.3*np.median(V[max(i-25, 0):i+25, :]):
                cands.append(SS[i, j] + 1j*TT[i, j])
    zeros = []
    for z in cands:
        for _ in range(40):
            f = Zf(z)[0]; h = 1e-6
            d = (Zf(z + h)[0] - Zf(z - h)[0])/(2*h)
            step = f/d; z = z - step
            if abs(step) < tol: break
        if abs(Zf(z)[0]) < 1e-7 and sig[0] - 0.05 < z.real < sig[1] + 0.05 and tt[0] - 0.5 < z.imag < tt[1] + 0.5:
            if all(abs(z - y) > 1e-6 for y in zeros): zeros.append(z)
    zeros.sort(key=lambda z: z.imag)
    return zeros, (S, T, V)
