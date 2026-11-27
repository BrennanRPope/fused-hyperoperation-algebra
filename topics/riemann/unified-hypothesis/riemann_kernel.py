"""Riemann's kernel Phi and H_t(z) = int_0^inf e^{t u^2} Phi(u) cos(z u) du (Rodgers–Tao normalisation;
H_0(z) = (1/8) xi(1/2 + i z/2), so zeros of H_0 sit at z = 2*gamma)."""
import mpmath as mp
mp.mp.dps = 45
UMAX = mp.mpf('0.95')
NQ = 320

def Phi(u):
    s = mp.mpf(0)
    e4 = mp.e ** (4 * u)
    for n in range(1, 40):
        term = (2 * mp.pi ** 2 * n ** 4 * mp.e ** (9 * u) - 3 * mp.pi * n ** 2 * mp.e ** (5 * u)) * mp.e ** (-mp.pi * n * n * e4)
        s += term
        if abs(term) < mp.mpf(10) ** (-60):
            break
    return s

_x, _w = None, None
def nodes():
    global _x, _w
    if _x is None:
        xs, ws = zip(*mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec))  # degree 3*2^8 nodes
        _x = [(x + 1) * UMAX / 2 for x in xs]
        _w = [w * UMAX / 2 * Phi(u) for u, w in zip(_x, ws)]
    return _x, _w

def H(z, t=0, weight=None):
    xs, ws = nodes()
    z = mp.mpmathify(z)
    if weight is None:
        return mp.fsum(w * mp.e ** (t * u * u) * mp.cos(z * u) for u, w in zip(xs, ws))
    return mp.fsum(w * weight(u) * mp.e ** (t * u * u) * mp.cos(z * u) for u, w in zip(xs, ws))

if __name__ == '__main__':
    import time
    t0 = time.time()
    xs, ws = nodes(); print(len(xs), 'nodes', time.time() - t0)
    t0 = time.time()
    for g in [14.134725141734693, 21.022039638771555, 49.773832477672302, 74.920370]:
        print(g, mp.nstr(H(2 * g), 5), mp.nstr(H(2 * g + 0.5), 5))
    print('time per eval', (time.time() - t0) / 8)
    # compare with xi
    z = mp.mpf(30)
    print(mp.nstr(H(z), 12), mp.nstr(mp.re(mp.xi(mp.mpf(0.5) + 1j * z / 2)) / 8, 12) if hasattr(mp, 'xi') else '')
