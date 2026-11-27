"""DRAFT. A real four-dimensional return map of the planar three-body problem.

Section: the line-ups (oriented area zero) with body 3 between bodies 1 and 2, crossed with dA/dt < 0,
on the level set {centre of mass 0, momentum 0, angular momentum 0, energy E}.
Coordinates s = (p, q, u1, u2): p = |r2 - r1|, q = where body 3 sits from the midpoint of 1 and 2 (along the
line, |q| < p/2), u1, u2 = velocities of bodies 1 and 2 along the line (frame turned so r2 - r1 points along +x).
The velocity across the line is fixed by L = 0, P = 0 and the energy, and its sign by dA/dt < 0.
Every Li-Liao orbit starts on this section (start: p = 2, q = 0, u1 = u2 = v1)."""
import numpy as np
import cb3

def to_state(s, E):
    p, q, u1, u2 = s
    c = -q / 3
    x = np.array([-p / 2 + c, p / 2 + c, q + c])
    ux = np.array([u1, u2, -u1 - u2])
    n = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
    U = 1 / abs(x[1] - x[0]) + 1 / abs(x[2] - x[0]) + 1 / abs(x[2] - x[1])
    b2 = 2 * (E + U - 0.5 * (ux ** 2).sum()) / (n @ n)
    if b2 < 0 or abs(q) >= p / 2: return None
    beta = np.sqrt(b2)
    uy = beta * n
    # dA/dt = -(v2y - v1y)(x3 - x1) + p (v3y - v1y)  must be < 0
    dA = -(uy[1] - uy[0]) * (x[2] - x[0]) + p * (uy[2] - uy[0])
    if dA > 0: uy = -uy
    X = np.zeros(12)
    X[0:6:2] = x; X[6::2] = ux; X[7::2] = uy
    return X

def to_coords(X):
    r = cb3.pos(X); v = cb3.vel(X)
    a = r[1] - r[0]; p = np.linalg.norm(a); e = a / p; ep = np.array([-e[1], e[0]])
    xs = r @ e; q = xs[2] - 0.5 * (xs[0] + xs[1])
    u = v @ e
    return np.array([p, q, u[0], u[1]])

def is_hit(X):
    """body 3 strictly between 1 and 2 along the line, and dA/dt < 0."""
    r = cb3.pos(X); v = cb3.vel(X)
    a = r[1] - r[0]; b = r[2] - r[0]
    t = (b @ a) / (a @ a)
    dA = (v[1, 0] - v[0, 0]) * b[1] - (v[1, 1] - v[0, 1]) * b[0] + a[0] * (v[2, 1] - v[0, 1]) - a[1] * (v[2, 0] - v[0, 0])
    return 0 < t < 1 and dA < 0

def next_hit(X, tmax=400.0):
    return cb3.next_hit_c(X, tmax=tmax)

def next_hit_py(X, dt=0.002, tmax=400.0):
    """flow to the next section hit; returns (state at hit, time taken)."""
    t = 0.0
    X = X.copy()
    A0 = cb3.area(X)
    # leave the section first
    X1 = cb3.flow(X, 1e-6); t += 1e-6; A0 = cb3.area(X1); X = X1
    while t < tmax:
        Xn = cb3.flow(X, dt); A1 = cb3.area(Xn)
        if A0 > 0 and A1 <= 0:          # crossing with A decreasing
            lo, hi, Xlo = 0.0, dt, X
            for _ in range(60):
                mid = 0.5 * (lo + hi); Xm = cb3.flow(X, mid)
                if cb3.area(Xm) > 0: lo = mid
                else: hi = mid
                if hi - lo < 1e-15: break
            Xh = cb3.flow(X, hi); th = t + hi
            if is_hit(Xh): return Xh, th
        X, A0, t = Xn, A1, t + dt
    raise RuntimeError("no hit")

def F(s, E, n=1):
    X = to_state(s, E)
    if X is None: return None
    T = 0.0
    for _ in range(n):
        X, dt = next_hit(X); T += dt
    return to_coords(X), T

def jac(s, E, n=1, h=1e-7):
    J = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4); e[j] = h
        J[:, j] = (F(s + e, E, n)[0] - F(s - e, E, n)[0]) / (2 * h)
    return J

def returns_per_period(s, E, T, tol=1e-6):
    """number of section hits in one period T (the section period of the orbit)."""
    X = to_state(s, E); t = 0.0; k = 0
    while t < T - tol:
        X, dt = next_hit(X); t += dt; k += 1
        if abs(t - T) < 1e-4 * max(1, T): return k
    raise RuntimeError("period mismatch: hits overshoot T (t=%g, T=%g)" % (t, T))

def refine(s, E, k, iters=8, h=1e-7):
    """Newton on F^k(s) - s."""
    s = np.array(s, float)
    for _ in range(iters):
        g = F(s, E, k)[0] - s
        if np.abs(g).max() < 1e-12: break
        J = jac(s, E, k, h)
        s = s - np.linalg.solve(J - np.eye(4), g)
    return s, float(np.abs(F(s, E, k)[0] - s).max())
