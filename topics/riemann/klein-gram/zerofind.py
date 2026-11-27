"""Find every zero of Z_kind in a rectangle: winding of the phase round each grid cell, then Newton."""
import numpy as np
from zetas import Z

def cell_zeros(F, sig, ts):
    ph = np.angle(F)
    d = lambda a, b: ((b - a + np.pi) % (2*np.pi)) - np.pi
    w = (d(ph[:-1, :-1], ph[:-1, 1:]) + d(ph[:-1, 1:], ph[1:, 1:]) +
         d(ph[1:, 1:], ph[1:, :-1]) + d(ph[1:, :-1], ph[:-1, :-1]))
    i, j = np.nonzero(np.abs(w) > np.pi)
    return (sig[j] + sig[j+1])/2 + 1j*(ts[i] + ts[i+1])/2, np.rint(w[i, j]/(2*np.pi)).astype(int)

def newton(z, kind, it=30, h=1e-6):
    z = np.array(z, dtype=complex)
    for _ in range(it):
        f = Z(z, kind); fp = (Z(z + h, kind) - Z(z - h, kind))/(2*h)
        step = f/fp
        step = np.where(np.abs(step) > 0.05, 0.05*step/np.abs(step), step)
        z = z - step
    return z, np.abs(Z(z, kind))

def find_zeros(kind, s0, s1, t0, t1, ds=0.02, dt=0.02, band=20.0):
    sig = np.arange(s0, s1 + ds/2, ds)
    out = []
    t = t0
    while t < t1:
        tb = min(t + band, t1)
        ts = np.arange(t, tb + dt/2, dt)
        F = Z(sig[None, :] + 1j*ts[:, None], kind)
        zs, wn = cell_zeros(F, sig, ts)
        out.extend(zs[wn > 0])
        t = tb
    z, res = newton(np.array(out), kind)
    z = z[(z.real > s0) & (z.real < s1) & (z.imag > t0) & (z.imag < t1)]
    # merge duplicates
    z = np.array(sorted(z, key=lambda c: c.imag))
    keep = []
    for c in z:
        if not keep or min(abs(c - k) for k in keep[-6:]) > 1e-6: keep.append(c)
    return np.array(keep)
