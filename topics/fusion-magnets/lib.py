"""Shared models for the fusion-magnets pass. Known engineering maths, plus the theory's operation."""
import numpy as np

MU0 = 4e-7*np.pi

def circ(a, b, k):            # the theory's operation a o_k b = a + b + k ab
    return a + b + k*a*b

def fold(xs, k):
    out = 0.0
    for x in xs: out = circ(out, x, k)
    return out

# REBCO Ic(B,T), field perpendicular to tape. Fit form of Succi, Ballarino et al. (IEEE TAS 34(3) 2024),
# constants as fitted to Faraday Factory tape in arXiv 2411.18124 Table III (a fit to measured data, not data itself).
P = dict(Ic00=6098.7, B000=0.5769, Birr00=200.0, Tc=85.0, alpha=0.48175, beta=1.4, gamma=1.45, q=2.0)
def Ic(B, T, p=P):
    t = 1 - T/p['Tc']
    Ic0 = p['Ic00']*t**p['gamma']; B0 = p['B000']*t; Birr = p['Birr00']*t**p['beta']
    return Ic0*(1 + B/B0)**(-p['alpha'])*np.clip(1 - B/Birr, 0, None)**p['q']

def corner_sum_log(f, B1, B2, T1, T2):
    """2x2 factorial interaction contrast of log f (Yates). Zero exactly when log f separates in B and T."""
    return np.log(f(B2,T2)) - np.log(f(B2,T1)) - np.log(f(B1,T2)) + np.log(f(B1,T1))

def KS(g, rho):                # Kreisselmeier-Steinhauser = log-sum-exp at temperature 1/rho
    m = np.max(g); return m + np.log(np.sum(np.exp(rho*(g-m))))/rho

def Lp(g, p):
    g = np.asarray(g, float); m = g.max(); return m*np.sum((g/m)**p)**(1/p)

# ---- coil geometry -----------------------------------------------------------
def modular_coil(n=4000, R0=5.0, a=2.0, d1=0.15, d2=0.06):
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    z = d1*np.sin(th) + d2*np.sin(2*th)
    R = R0 + a*np.cos(th)
    return np.stack([R*np.cos(z), R*np.sin(z), a*np.sin(th)], 1)

def planar_coil(n=4000, a=2.0, b=1.3):
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    return np.stack([a*np.cos(th), np.zeros_like(th), b*np.sin(th)], 1)

def frames(X):
    """Arc length, tangent, Frenet curvature/torsion (spectral derivatives), closed curve."""
    n = len(X); kf = np.fft.fftfreq(n, 1.0/n)*1j
    d = lambda f, m: np.real(np.fft.ifft((kf**m)[:,None]*np.fft.fft(f, axis=0), axis=0))
    X1, X2, X3 = d(X,1), d(X,2), d(X,3)
    sp = np.linalg.norm(X1, axis=1); ds = sp*(2*np.pi/n)
    c = np.cross(X1, X2); cn = np.linalg.norm(c, axis=1)
    kappa = cn/sp**3
    tau = np.einsum('ij,ij->i', c, X3)/cn**2
    T = X1/sp[:,None]
    Nf = np.cross(c, X1); Nf /= np.linalg.norm(Nf, axis=1)[:,None]
    return dict(ds=ds, L=ds.sum(), T=T, N=Nf, kappa=kappa, tau=tau)

def rodrigues(v, k, ang):
    return v*np.cos(ang) + np.cross(k, v)*np.sin(ang) + k*np.dot(k, v)*(1-np.cos(ang))

def bishop(T, N0):
    """Discrete parallel transport (rotation-minimising frame) of N0 along tangents T; returns normals and closing angle."""
    n = len(T); U = np.zeros_like(T); U[0] = N0
    for i in range(n):
        t0, t1 = T[i], T[(i+1) % n]
        ax = np.cross(t0, t1); s = np.linalg.norm(ax); c = np.clip(np.dot(t0, t1), -1, 1)
        u = U[i] if s < 1e-15 else rodrigues(U[i], ax/s, np.arctan2(s, c))
        if i+1 < n: U[i+1] = u
        else: uend = u
    # angle from N0 to uend about T[0]
    ang = np.arctan2(np.dot(np.cross(N0, uend), T[0]), np.dot(N0, uend))
    return U, ang

def spherical_area(T):
    """Signed solid angle enclosed by closed spherical polygon T (unit vectors), via triangle fan from a pole."""
    c = T.mean(0); c /= np.linalg.norm(c)
    if np.linalg.norm(T.mean(0)) < 1e-6: c = np.array([0.0, 1.0, 0.0])  # great-circle-like: any pole off the curve
    A = 0.0
    for i in range(len(T)):
        a, b = T[i], T[(i+1) % len(T)]
        num = np.dot(c, np.cross(a, b)); den = 1 + np.dot(c, a) + np.dot(a, b) + np.dot(b, c)
        A += 2*np.arctan2(num, den)
    return A

def wrap(x): return (x + np.pi) % (2*np.pi) - np.pi

def dented_coil(n=6000, s=2.0, e=0.4, h=0.02):
    """Coil with a dent (in-plane curvature crosses zero for e > 0.2) and a small out-of-plane wiggle."""
    th = np.linspace(0, 2*np.pi, n, endpoint=False); r = s*(1 + e*np.cos(2*th))
    return np.stack([r*np.cos(th), s*h*np.sin(2*th), r*np.sin(th)], 1)
