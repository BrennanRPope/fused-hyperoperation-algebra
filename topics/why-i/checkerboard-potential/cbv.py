"""Checkerboard / telegraph chain with a potential. Units c = hbar = 1, lattice step eps.
One tick: turn (1 + k t), potential (1 + k v) at the current site, then move.
t = -eps*mu*X (swap directions), v = eps*V(x)*I (same on both directions).
linear=True : potential factor 1 + k eps V   (the plain chain)
linear=False: potential factor exp(k eps V)  (k = -i: a phase, unitary; k = -1: a killing weight e^(-eps V))
normalized  : k = -1 -> turn step (1 - a eps, a eps); k = -i -> (cos m eps, i sin m eps)
theta (optional): link phases for a vector potential; an R move from site s gets e^(i theta[s]), an L move from s gets e^(-i theta[s-1]).
"""
import numpy as np

def tick_factors(eps, mu, k, normalized):
    if normalized:
        if np.isclose(k, -1): return 1 - eps*mu, eps*mu
        return np.cos(eps*mu), 1j*np.sin(eps*mu)
    return 1.0, -k*eps*mu

def evolve_v(psi0, N, eps, mu, k, V, linear=True, normalized=False, theta=None, keep=False):
    """psi0: (n, 2) complex, sites on a ring; V: (n,) potential at each site. Returns final state (or all)."""
    psi = psi0.astype(complex).copy()
    c_, t_ = tick_factors(eps, mu, k, normalized)
    pf = (1 + k*eps*V) if linear else np.exp(k*eps*V)
    out = [psi.copy()] if keep else None
    for _ in range(N):
        psi = np.stack([c_*psi[:, 0] + t_*psi[:, 1], c_*psi[:, 1] + t_*psi[:, 0]], axis=1)
        psi = psi*pf[:, None]
        R, L = psi[:, 0], psi[:, 1]
        if theta is not None:
            R = R*np.exp(1j*theta); L = L*np.exp(-1j*np.roll(theta, 1))
        psi = np.stack([np.roll(R, 1), np.roll(L, -1)], axis=1)
        if keep: out.append(psi.copy())
    return np.array(out) if keep else psi

def reference(psi0, T, dx, mu, k, V, normalized_telegraph=False, steps=4000):
    """Independent continuum solver: Strang splitting with FFT on a periodic grid.
    d_t psi = (-sigma_z d_x - k mu sigma_x [- mu I if telegraph-normalized]) psi + k V psi."""
    n = psi0.shape[0]; p = 2*np.pi*np.fft.fftfreq(n, d=dx); dt = T/steps
    sz = np.array([1., -1.])
    # free propagator over dt for each p: exp(dt*G), G = [[-ip, -k mu], [-k mu, ip]] (- mu I)
    w = -k*mu
    Om = np.sqrt(w**2 - p**2 + 0j)
    ch = np.cosh(dt*Om); sh = np.where(np.abs(Om) > 1e-14, np.sinh(dt*Om)/np.where(np.abs(Om) > 1e-14, Om, 1), dt)
    U11 = ch - 1j*p*sh; U22 = ch + 1j*p*sh; U12 = w*sh
    damp = np.exp(-mu*dt) if normalized_telegraph else 1.0
    half = np.exp(k*V*dt/2)
    psi = psi0.astype(complex).copy()
    for _ in range(steps):
        psi = psi*half[:, None]
        F = np.fft.fft(psi, axis=0)
        F = np.stack([U11*F[:, 0] + U12*F[:, 1], U12*F[:, 0] + U22*F[:, 1]], axis=1)*damp
        psi = np.fft.ifft(F, axis=0)*half[:, None]
    return psi

def schrodinger(phi0, T, dx, m, V, steps=4000):
    """i phi_t = -(1/2m) phi'' + V phi, split-step FFT."""
    n = phi0.size; p = 2*np.pi*np.fft.fftfreq(n, d=dx); dt = T/steps
    kin = np.exp(-1j*p**2/(2*m)*dt); half = np.exp(-1j*V*dt/2); phi = phi0.astype(complex)
    for _ in range(steps): phi = np.fft.ifft(kin*np.fft.fft(phi*half))*half
    return phi

def heat_killing(u0, T, dx, D, V, steps=4000):
    """u_t = D u'' - V u, split-step FFT."""
    n = u0.size; p = 2*np.pi*np.fft.fftfreq(n, d=dx); dt = T/steps
    kin = np.exp(-D*p**2*dt); half = np.exp(-V*dt/2); u = u0.astype(complex)
    for _ in range(steps): u = np.fft.ifft(kin*np.fft.fft(u*half))*half
    return u.real
