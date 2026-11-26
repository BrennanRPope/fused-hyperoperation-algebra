"""2D incompressible Navier-Stokes on three flat surfaces, all built on the same square cover [0,2pi]^2.

  T2  : the square torus (the cover itself).
  Trect: the rectangle torus  T2 / tau,  tau(x,y) = (x+pi,  y)   straight gluing.
  K   : the Klein bottle       T2 / g,    g(x,y)   = (x+pi, -y)   flipped gluing (glide reflection).

Trect and K have the same area (2 pi^2) and the same fundamental rectangle [0,pi] x [0,2pi]; only the
gluing of the two short edges differs (straight vs flipped). Fields on a quotient are fields on the cover
that respect the symmetry. Velocity is a vector, so on K:  u(x+pi,-y) = u(x,y),  v(x+pi,-y) = -v(x,y).
Then the streamfunction (u = d_y psi, v = -d_x psi) and the vorticity w = -Lap psi are ODD under g:
they are pseudo-scalars (they flip sign when you go round the twist).

Method: pseudo-spectral vorticity-streamfunction, 2/3 dealiasing, RK4 with an exact viscous integrating
factor, adaptive CFL time step, plus a uniform (harmonic) mean flow (U, V) that is carried exactly.
"""
import numpy as np

TWO_PI = 2 * np.pi


class Grid:
    def __init__(self, N):
        self.N = N
        x = np.arange(N) * TWO_PI / N
        self.x = x
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")  # axis 0 = x, axis 1 = y
        k = np.fft.fftfreq(N, 1.0 / N)
        kr = np.fft.rfftfreq(N, 1.0 / N)
        self.KX, self.KY = np.meshgrid(k, kr, indexing="ij")
        self.K2 = self.KX**2 + self.KY**2
        self.K2inv = np.where(self.K2 > 0, 1.0 / np.where(self.K2 > 0, self.K2, 1), 0.0)
        kmax = N // 3
        self.dealias = (np.abs(self.KX) <= kmax) & (np.abs(self.KY) <= kmax)
        self.negy = (-np.arange(N)) % N

    def fft(self, f):
        return np.fft.rfft2(f)

    def ifft(self, fh):
        return np.fft.irfft2(fh, s=(self.N, self.N))

    # symmetry maps on physical fields
    def glide(self, f):
        """(g* f)(x,y) = f(x+pi, -y)"""
        return np.roll(f, -self.N // 2, axis=0)[:, self.negy]

    def shift(self, f):
        """(tau* f)(x,y) = f(x+pi, y)"""
        return np.roll(f, -self.N // 2, axis=0)


def project(grid, w, surface):
    """Project a vorticity field on the cover onto the fields allowed on `surface`."""
    if surface == "T2":
        return w
    if surface == "Trect":
        return 0.5 * (w + grid.shift(w))
    if surface == "K":
        return 0.5 * (w - grid.glide(w))  # vorticity is odd under the glide
    raise ValueError(surface)


def project_mean(U, V, surface):
    """Uniform (harmonic) flows allowed: T2, Trect keep (U,V); K keeps only U (dg = diag(1,-1))."""
    if surface == "K":
        return U, 0.0
    return U, V


class Flow:
    def __init__(self, N, nu, surface, w0, U=0.0, V=0.0, cfl=0.4, reproject=True):
        self.g = Grid(N)
        self.nu = nu
        self.surface = surface
        self.reproject = reproject
        w0 = project(self.g, w0, surface)
        self.wh = self.g.fft(w0) * self.g.dealias
        self.wh[0, 0] = 0.0
        self.U, self.V = project_mean(U, V, surface)
        self.t = 0.0
        self.cfl = cfl

    # fields
    def w(self):
        return self.g.ifft(self.wh)

    def psi(self):
        return self.g.ifft(self.wh * self.g.K2inv)

    def velocity(self, wh=None):
        wh = self.wh if wh is None else wh
        ph = wh * self.g.K2inv
        u = self.g.ifft(1j * self.g.KY * ph)
        v = self.g.ifft(-1j * self.g.KX * ph)
        return u, v

    def rhs(self, wh):
        u, v = self.velocity(wh)
        wx = self.g.ifft(1j * self.g.KX * wh)
        wy = self.g.ifft(1j * self.g.KY * wh)
        nl = -((u + self.U) * wx + (v + self.V) * wy)
        out = self.g.fft(nl) * self.g.dealias
        out[0, 0] = 0.0
        return out

    def step(self, dt):
        L = np.exp(-self.nu * self.g.K2 * dt / 2)  # half-step integrating factor
        w0 = self.wh
        k1 = self.rhs(w0)
        k2 = self.rhs(L * (w0 + dt / 2 * k1))
        k3 = self.rhs(L * w0 + dt / 2 * k2)
        k4 = self.rhs(L * L * w0 + dt * L * k3)
        self.wh = L * L * w0 + dt / 6 * (L * L * k1 + 2 * L * (k2 + k3) + k4)
        if self.reproject and self.surface != "T2":
            self.wh = self.g.fft(project(self.g, self.g.ifft(self.wh), self.surface)) * self.g.dealias
        self.t += dt

    def dt_cfl(self):
        u, v = self.velocity()
        umax = np.max(np.abs(u + self.U)) + np.max(np.abs(v + self.V)) + 1e-12
        dx = TWO_PI / self.g.N
        return min(self.cfl * dx / umax, 0.1)

    def advance_to(self, T):
        while self.t < T - 1e-12:
            dt = min(self.dt_cfl(), T - self.t)
            self.step(dt)

    # diagnostics (per unit area, averaged over the cover = averaged over the surface)
    def energy(self, include_mean=True):
        u, v = self.velocity()
        e = 0.5 * np.mean(u**2 + v**2)
        if include_mean:
            e += 0.5 * (self.U**2 + self.V**2)
        return e

    def enstrophy(self):
        return 0.5 * np.mean(self.w() ** 2)

    def spectrum(self):
        """Energy spectrum E(k) in integer shells, from the rfft (counting both halves)."""
        N = self.g.N
        ph = self.wh * self.g.K2inv
        e = 0.5 * self.g.K2 * np.abs(ph) ** 2 / N**4
        e[:, 1:] *= 2  # rfft half-plane weighting (last column at Nyquist negligible after dealias)
        kk = np.sqrt(self.g.K2)
        kbin = np.rint(kk).astype(int)
        return np.bincount(kbin.ravel(), weights=e.ravel())


def random_vorticity(N, k0=6.0, seed=0):
    """Random vorticity on the cover with energy spectrum ~ k^7 exp(-2 (k/k0)^2) style bump at k0, normalized later."""
    rng = np.random.default_rng(seed)
    g = Grid(N)
    kk = np.sqrt(g.K2)
    amp = np.where(kk > 0, kk ** 3 * np.exp(-((kk / k0) ** 2)), 0.0)  # |w_hat| ~ k E(k)^(1/2)/... bump
    phase = np.exp(2j * np.pi * rng.random(g.K2.shape))
    w = g.ifft(amp * phase * g.dealias)
    return w


def normalize_energy(N, w, E=0.5):
    g = Grid(N)
    wh = g.fft(w)
    ph = wh * g.K2inv
    u = g.ifft(1j * g.KY * ph)
    v = g.ifft(-1j * g.KX * ph)
    e = 0.5 * np.mean(u**2 + v**2)
    return w * np.sqrt(E / e)
