"""Shared numerics for the diamond-lattice Potts map in the transmissivity t = v/(Q+v).
T_Q(t) = F_{Q-2,1-Q}(t^2, t^2): series (square in t) then parallel (the landscape law)."""
import numpy as np

def lens(p, q, a, b):
    return (a + b + p*a*b) / (1 - q*a*b)

def T(t, Q):
    s = t*t
    return lens(Q - 2, 1 - Q, s, s)

def fixed_ys(Q):
    """All real y with (y-1)^2 (y+1) = Q; the fixed points are t = 1/y."""
    r = np.roots([1, -1, -1, 1 - Q])
    return sorted(z.real for z in r if abs(z.imag) < 1e-7)

def ystar(Q):
    return max(fixed_ys(Q))

def lam(Q):
    y = ystar(Q)
    return (2*y/(y + 1))**2

def nu(Q):
    return np.log(2)/np.log(lam(Q))

def beraha(n):
    return 4*np.cos(np.pi/n)**2

def basin_grid(Q, xlim=(-2.2, 2.2), ylim=(-2.2, 2.2), n=400, iters=120, tol=1e-3):
    """Label: 0 hot (t -> 0), 1 cold (t -> 1), 2 third attractor (antiferro fixed point, or infinity at Q = 1),
    3 undecided. Also the step at which each point was decided (for shading)."""
    xs = np.linspace(*xlim, n); ys = np.linspace(*ylim, n)
    Z = xs[None, :] + 1j*ys[:, None]
    att = []
    for y in fixed_ys(Q):
        if -1/3 - 1e-4 < y < 1 - 1e-6 and abs(y) > 1e-9:
            att.append(1/y)
    lab = np.full(Z.shape, 3, dtype=np.int8); step = np.zeros(Z.shape)
    with np.errstate(all='ignore'):
        for k in range(iters):
            Z = T(Z, Q)
            und = lab == 3
            if abs(Q - 1) < 1e-12:
                big = und & (~np.isfinite(Z) | (np.abs(Z) > 1e4))
                lab[big] = 2
            else:
                Z = np.where(np.isfinite(Z), Z, (Q - 2)/(Q - 1))
            for a in att:
                lab[und & (np.abs(Z - a) < tol)] = 2
            lab[und & (np.abs(Z) < tol)] = 0
            lab[und & (np.abs(Z - 1) < tol)] = 1
            newly = und & (lab != 3)
            step[newly] = k + 1
    return xs, ys, lab, step
