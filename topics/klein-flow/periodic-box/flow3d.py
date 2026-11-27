"""3D incompressible Navier-Stokes, pseudo-spectral, on a box cover [0,2pi)^2 x [0,2pi h).
Surfaces (quotients of the cover, all flat 3-manifolds):
  'T3'  : the cover itself (periodic box of sides 2pi, 2pi, 2pi h)
  'TR'  : rectangle torus x circle: cover / tau, tau(x,y,z) = (x+pi, y, z)
  'K'   : Klein bottle x circle: cover / g, g(x,y,z) = (x+pi, -y, z)  (non-orientable)
TR and K have the same volume and the same fundamental box [0,pi]x[0,2pi]x[0,2pi h].
Velocity form, rotational nonlinearity u x omega, Leray projection, 2/3 dealiasing,
integrating-factor RK4 for viscosity, CFL time step.
"""
import numpy as np
from scipy import fft as sfft

class Box:
    def __init__(self, N, h, Nz=None):
        self.N = N; self.h = h
        self.Nz = Nz if Nz is not None else max(16, int(round(N*h/4))*4)
        kx = sfft.fftfreq(N, 1.0/N)
        kz = sfft.rfftfreq(self.Nz, 1.0/self.Nz) / h
        self.KX, self.KY, self.KZ = np.meshgrid(kx, kx, kz, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        self.K2i = np.where(self.K2 > 0, 1.0/np.where(self.K2 > 0, self.K2, 1), 0.0)
        kxm = np.abs(self.KX) <= N/3
        kym = np.abs(self.KY) <= N/3
        kzm = np.abs(self.KZ*h) <= self.Nz/3
        self.mask = (kxm & kym & kzm).astype(float)
        self.vol = (2*np.pi)**3 * h
        self.shape = (N, N, self.Nz)

    def fwd(self, f):  return sfft.rfftn(f, axes=(0, 1, 2), workers=1)
    def inv(self, F):  return sfft.irfftn(F, s=self.shape, axes=(0, 1, 2), workers=1)

    def leray(self, U):
        dot = self.KX*U[0] + self.KY*U[1] + self.KZ*U[2]
        return np.array([U[0] - self.KX*dot*self.K2i,
                         U[1] - self.KY*dot*self.K2i,
                         U[2] - self.KZ*dot*self.K2i])

    def curl(self, U):
        return 1j*np.array([self.KY*U[2] - self.KZ*U[1],
                            self.KZ*U[0] - self.KX*U[2],
                            self.KX*U[1] - self.KY*U[0]])

def glide(f, N):
    """f(x+pi, -y, z) on the grid."""
    return np.roll(f, -N//2, axis=0)[:, (-np.arange(N)) % N, :]

def shift(f, N):
    return np.roll(f, -N//2, axis=0)

def symmetrize(u, surf, N):
    """Project a physical velocity field (3,N,N,Nz) onto fields that live on the surface."""
    if surf == 'T3':
        return u
    if surf == 'TR':
        return 0.5*(u + np.array([shift(c, N) for c in u]))
    if surf == 'K':
        s = np.array([1.0, -1.0, 1.0])   # Dg = diag(1,-1,1)
        return 0.5*(u + np.array([s[i]*glide(u[i], N) for i in range(3)]))
    raise ValueError(surf)

def glide_defect(u, N):
    s = np.array([1.0, -1.0, 1.0])
    g = np.array([s[i]*glide(u[i], N) for i in range(3)])
    return np.abs(u - g).max() / max(np.abs(u).max(), 1e-300)

class NS3:
    def __init__(self, box, nu, surf):
        self.b = box; self.nu = nu; self.surf = surf

    def rhs(self, U):
        b = self.b
        u = np.array([b.inv(c) for c in U])
        w = np.array([b.inv(c) for c in b.curl(U)])
        c = np.array([u[1]*w[2] - u[2]*w[1], u[2]*w[0] - u[0]*w[2], u[0]*w[1] - u[1]*w[0]])
        C = np.array([b.fwd(x) for x in c]) * b.mask
        return b.leray(C)

    def step(self, U, dt):
        L = np.exp(-self.nu*self.b.K2*dt/2)
        k1 = self.rhs(U)
        k2 = self.rhs(L*(U + dt/2*k1))
        k3 = self.rhs(L*U + dt/2*k2)
        k4 = self.rhs(L*L*U + dt*L*k3)
        return L*L*U + dt/6*(L*L*k1 + 2*L*(k2 + k3) + k4)

    def umax(self, U):
        return max(np.abs(self.b.inv(c)).max() for c in U)

    def dt_cfl(self, U, cfl=1.0):
        b = self.b
        u = [np.abs(b.inv(c)) for c in U]
        dx = 2*np.pi/b.N; dz = 2*np.pi*b.h/b.Nz
        rate = (u[0]/dx + u[1]/dx + u[2]/dz).max()
        return cfl/max(rate, 1e-8)

    # diagnostics (cover averages are surface averages: the quotient is a fixed fraction)
    def _wsum(self, A):
        """sum over full spectrum of |A|^2 using rfft halving, normalized to mean of square."""
        b = self.b
        w = np.full(A.shape[-1], 2.0); w[0] = 1.0
        if b.Nz % 2 == 0: w[-1] = 1.0
        n = b.N*b.N*b.Nz
        return (np.abs(A)**2 * w).sum() / n**2

    def energy(self, U):
        return 0.5*sum(self._wsum(c) for c in U)

    def enstrophy(self, U):     # 1/2 <|omega|^2>
        return 0.5*sum(self._wsum(c) for c in self.b.curl(U))

    def energy3d(self, U):      # energy in modes with kz != 0
        V = U.copy(); V[..., 0] = 0
        return 0.5*sum(self._wsum(c) for c in V)

    def wmax(self, U):
        w = np.array([self.b.inv(c) for c in self.b.curl(U)])
        return np.sqrt((w**2).sum(0)).max()

def random_field(box, kmin, kmax, seed, kz_index=None, rng=None):
    """Random divergence-free physical field with horizontal |k| in [kmin,kmax].
    kz_index=0 -> z-independent with u_z=0 (2D base); kz_index=n -> only kz index +-n (3D part)."""
    rng = rng or np.random.default_rng(seed)
    b = box
    U = (rng.standard_normal((3,) + b.KX.shape) + 1j*rng.standard_normal((3,) + b.KX.shape))
    kh = np.sqrt(b.KX**2 + b.KY**2)
    kzi = np.rint(b.KZ*b.h)
    sel = (kh >= kmin) & (kh <= kmax) & (kzi == kz_index)
    U *= sel
    if kz_index == 0:
        U[2] = 0
    U = b.leray(U) * b.mask
    u = np.array([b.inv(c) for c in U])
    return u
