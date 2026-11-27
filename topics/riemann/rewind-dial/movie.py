"""Zeros of F_theta in the box |Re z| <= 1, 0.3 <= Im z <= 36, found independently at each angle (phase winding on a
grid, then Newton). DRAFT."""
import numpy as np, json, time, sys
from rotate import setup, newton
def Fgrid(S, Z, th):
    X = np.cos(th)*S["a"]*S["M"] + np.sin(th)*S["H"]; out = np.empty(Z.size, complex); zf = Z.ravel()
    for i in range(0, zf.size, 4000): out[i:i+4000] = (S["w"][None, :]*np.exp(np.outer(zf[i:i+4000], X))).sum(1)
    return out.reshape(Z.shape)
def zeros_box(S, th, h=0.01):
    xr = np.arange(-1.0, 1.0+h/2, h); yr = np.arange(0.3, 36.0+h/2, h); Z = xr[None, :] + 1j*yr[:, None]
    ph = np.angle(Fgrid(S, Z, th)); d = lambda a, b: np.angle(np.exp(1j*(b-a)))
    w = (d(ph[:-1,:-1], ph[:-1,1:]) + d(ph[:-1,1:], ph[1:,1:]) + d(ph[1:,1:], ph[1:,:-1]) + d(ph[1:,:-1], ph[:-1,:-1]))/(2*np.pi)
    I, J = np.nonzero(np.rint(w)); z = np.array([xr[j]+h/2 + 1j*(yr[i]+h/2) for i, j in zip(I, J)])
    return newton(S, z, th) if z.size else z
if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    S = setup(n); out = {"n": n, "a": S["a"], "frames": {}}; t0 = time.time()
    for deg in list(range(0, 91, 5)) + [1, 2, 3]:
        z = zeros_box(S, np.radians(deg)); z = z[(np.abs(z.real) <= 1.0) & (z.imag > 0.3) & (z.imag <= 36)]
        z = np.unique(np.round(z, 8))
        out["frames"][deg] = [[float(v.real), float(v.imag)] for v in z]
        near = np.abs(z.real) < 0.05
        print("n=%d %2d deg: %3d zeros in box, %3d within 0.05 of the line, max|Re| of lowest five: %s  (t=%.0fs)" %
              (n, deg, z.size, near.sum(), np.round(np.abs(z[np.argsort(z.imag)][:5].real).max(), 3) if z.size else None, time.time()-t0), flush=True)
        json.dump(out, open("movie_n%d.json" % n, "w"))
