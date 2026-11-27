"""Fast float64 kernels for the pictures (5D parent and bath).

H_t(z) = int_0^U e^{t u^2} K(u) cos(z u) du, z complex.
Riemann: K = Phi (Rodgers-Tao normalisation, H_0(z) = xi(1/2+iz/2)/8, zeros at z = 2*gamma).
p-ball: K = (1-u^p)^a on [0,1].
Float64 is good for the Riemann kernel up to about x = 75 (H_0 ~ e^{-pi x/8}).
"""
import numpy as np
from scipy.special import roots_jacobi, roots_legendre

GAMMAS = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
          37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478]


def Phi(u):
    u = np.asarray(u, float)
    s = np.zeros_like(u)
    for n in range(1, 12):
        s += (2 * np.pi**2 * n**4 * np.exp(9 * u) - 3 * np.pi * n**2 * np.exp(5 * u)) * np.exp(-np.pi * n**2 * np.exp(4 * u))
    return s


class Kernel:
    def __init__(self, u, w, name):
        self.u, self.w, self.name = np.asarray(u), np.asarray(w), name

    def H(self, z, t=0.0, weight=None):
        z = np.asarray(z, complex)
        w = self.w * np.exp(t * self.u**2)
        if weight is not None:
            w = w * weight(self.u)
        sh = z.shape
        zf = z.reshape(-1)
        out = np.empty(zf.shape, complex)
        for i in range(0, zf.size, 4000):
            out[i:i + 4000] = np.cos(np.outer(zf[i:i + 4000], self.u)) @ w
        return out.reshape(sh)

    def dH(self, z, t=0.0, weight=None):
        z = np.asarray(z, complex)
        w = self.w * np.exp(t * self.u**2) * self.u
        if weight is not None:
            w = w * weight(self.u)
        return -(np.sin(np.multiply.outer(z, self.u)) @ w)

    def moment_gen(self, x):
        """M(x) = int_{-U}^{U} K(s) e^{x s} ds = 2 int_0^U K cosh(x s)."""
        x = np.asarray(x, float)
        return 2 * (np.cosh(np.multiply.outer(x, self.u)) @ self.w)


def riemann(n=400, U=0.95):
    x, w = roots_legendre(n)
    u = (x + 1) * U / 2
    return Kernel(u, w * U / 2 * Phi(u), 'Riemann Phi')


def pball(p, a, n=400):
    x, w = roots_jacobi(n, a, 0.0)          # weight (1-x)^a on [-1,1]
    u = (x + 1) / 2                          # (1-x)^a = 2^a (1-u)^a
    g = np.where(u < 1, ((1 - u**p) / np.maximum(1 - u, 1e-300))**a, p**a)
    return Kernel(u, w / 2**(a + 1) * g, f'p-ball p={p:g}, a={a:g}')


def real_zeros(k, xmax, dx=0.01, t=0.0, weight=None, x0=0.05):
    x = np.arange(x0, xmax, dx)
    h = k.H(x + 0j, t, weight).real
    i = np.where(np.sign(h[:-1]) * np.sign(h[1:]) < 0)[0]
    return x[i] - h[i] * dx / (h[i + 1] - h[i])


def newton(k, z0, t=0.0, weight=None, it=40):
    z = np.asarray(z0, complex).copy()
    for _ in range(it):
        f = k.H(z, t, weight)
        d = k.dH(z, t, weight)
        step = f / d
        z = z - step
        if np.all(np.abs(step) < 1e-11):
            break
    return z


def complex_zeros(k, xmax, ymax, t=0.0, weight=None, nx=None, ny=None):
    """Local minima of |H| on a grid, refined by Newton; upper half and real line."""
    nx = nx or int(xmax * 12)
    ny = ny or int(ymax * 12) + 1
    X, Y = np.meshgrid(np.linspace(0.05, xmax, nx), np.linspace(-0.02, ymax, ny))
    A = np.abs(k.H(X + 1j * Y, t, weight))
    m = np.ones_like(A, bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx or dy:
                m &= A <= np.roll(np.roll(A, dy, 0), dx, 1)
    m[0, :] = m[-1, :] = m[:, 0] = m[:, -1] = False
    z0 = (X + 1j * Y)[m]
    z = newton(k, z0, t, weight)
    ok = (np.abs(k.H(z, t, weight)) < 1e-6 * np.abs(k.H(z + 0.3, t, weight)) + 1e-300) & (z.real > 0) & (z.real < xmax) & (z.imag > -0.05) & (z.imag < ymax)
    z = z[ok]
    z = np.where(np.abs(z.imag) < 1e-7, z.real + 0j, z)
    z = np.unique(np.round(z, 6))
    return z
