"""Orbit tools for topics/entropy-topology (DRAFT): fixed-point indices, Conley-index bookkeeping,
the logistic windows from entropy-knots, and the area-preserving Henon map."""
import json, os
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- logistic map ----------

def logistic_iter(x, r, n):
    for _ in range(n):
        x = r * x * (1 - x)
    return x


def logistic_dmult(x, r, n):
    """Multiplier (f^n)'(x) and the orbit."""
    m = 1.0
    orb = []
    for _ in range(n):
        orb.append(x)
        m *= r * (1 - 2 * x)
        x = r * x * (1 - x)
    return m, orb


def interval_index(f, a, b):
    """Sum of fixed-point indices of f in [a, b], read from the endpoints only.

    In one dimension the index of an isolated fixed point is sign(1 - f'), the local degree of
    g(x) = x - f(x). Their sum over [a, b] is (sign g(b) - sign g(a)) / 2: no derivative and no
    knowledge of where the fixed points are is needed, only that none sits on the boundary.
    """
    ga, gb = a - f(a), b - f(b)
    return int((np.sign(gb) - np.sign(ga)) // 2)


def periodic_points_logistic(r, n, grid=2_000_000):
    """All points of period dividing n of the logistic map, by sign changes of f^n(x) - x."""
    xs = np.linspace(0.0, 1.0, grid)
    g = logistic_iter(xs, r, n) - xs
    idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    pts = []
    for i in idx:
        pts.append(brentq(lambda x: logistic_iter(x, r, n) - x, xs[i], xs[i + 1], xtol=1e-15))
    for i in np.nonzero(g == 0)[0]:
        pts.append(xs[i])
    return np.array(sorted(pts))


def orbits_of_exact_period(r, n):
    pts = periodic_points_logistic(r, n)
    used = np.zeros(len(pts), bool)
    orbits = []
    for i, p in enumerate(pts):
        if used[i]:
            continue
        m, orb = logistic_dmult(p, r, n)
        # exact period
        per = n
        for q in range(1, n):
            if n % q == 0 and abs(logistic_iter(p, r, q) - p) < 1e-9:
                per = q
                break
        for y in orb:
            j = np.argmin(np.abs(pts - y))
            if abs(pts[j] - y) < 1e-7:
                used[j] = True
        if per == n:
            orbits.append((np.array(orb), m))
    return orbits


def pattern(orb, r):
    """The orbit's pattern: perm[i] = rank of f(p_i), p sorted (the 1D knot type)."""
    s = np.sort(orb)
    img = r * s * (1 - s)
    return tuple(int(np.argmin(np.abs(s - v))) for v in img)


def load_windows():
    # the 65 windows from the entropy-knots topic (repo layout first, then the project-files layout)
    for p in [os.path.join(HERE, '..', 'entropy-knots', 'windows.json'),
              os.path.join(HERE, '..', '..', 'entropy-knots', 'folder', 'windows.json')]:
        if os.path.exists(p):
            return json.load(open(p))
    raise FileNotFoundError('windows.json from topics/entropy-knots/ not found')

# ---------- area-preserving Henon map (Henon 1969) ----------

def henon(z, ca, sa):
    x, y = z[..., 0], z[..., 1]
    w = y - x * x
    return np.stack([x * ca - w * sa, x * sa + w * ca], axis=-1)


def henon_n(z, ca, sa, n):
    for _ in range(n):
        z = henon(z, ca, sa)
    return z


def henon_jac_n(z, ca, sa, n):
    J = np.eye(2)
    for _ in range(n):
        x = z[0]
        Dz = np.array([[ca + 2 * x * sa, -sa], [sa - 2 * x * ca, ca]])
        J = Dz @ J
        z = henon(z, ca, sa)
    return J


def winding_index(F, z0, rho, m=720):
    """Fixed-point index of F at z0: winding number of F(z) - z round a small circle."""
    th = np.linspace(0, 2 * np.pi, m + 1)
    circ = z0[None, :] + rho * np.stack([np.cos(th), np.sin(th)], axis=1)
    v = F(circ) - circ
    ang = np.unwrap(np.arctan2(v[:, 1], v[:, 0]))
    return int(round((ang[-1] - ang[0]) / (2 * np.pi)))


def find_henon_orbits(ca, sa, nmax, seeds=60, box=1.0, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    found = []
    for n in range(1, nmax + 1):
        for _ in range(seeds * n):
            z = rng.uniform(-box, box, 2)
            ok = False
            for _ in range(60):
                Fz = henon_n(z, ca, sa, n) - z
                if not np.all(np.isfinite(Fz)) or np.linalg.norm(z) > 5:
                    break
                J = henon_jac_n(z, ca, sa, n) - np.eye(2)
                try:
                    step = np.linalg.solve(J, -Fz)
                except np.linalg.LinAlgError:
                    break
                z = z + step
                if np.linalg.norm(step) < 1e-13:
                    ok = True
                    break
            if not ok or np.linalg.norm(z) > 3:
                continue
            # exact period n
            if any(np.linalg.norm(henon_n(z, ca, sa, q) - z) < 1e-8 for q in range(1, n) if n % q == 0):
                continue
            orb = [z]
            for _ in range(n - 1):
                orb.append(henon(orb[-1], ca, sa))
            orb = np.array(orb)
            if any(o['n'] == n and np.min(np.linalg.norm(o['orbit'] - z, axis=1)) < 1e-7 for o in found):
                continue
            J = henon_jac_n(z, ca, sa, n)
            found.append({'n': n, 'orbit': orb, 'trace': float(np.trace(J)), 'det': float(np.linalg.det(J))})
    return found


def linear_type(tr):
    if abs(tr) < 2:
        return 'elliptic'
    return 'hyperbolic' if tr > 2 else 'inverse hyperbolic'

# ---------- Chirikov standard map (area-preserving twist map), on the lift ----------

def stdmap(z, K):
    x, p = z[..., 0], z[..., 1]
    p2 = p + K * np.sin(x)
    return np.stack([x + p2, p2], axis=-1)


def stdmap_n(z, K, n):
    for _ in range(n):
        z = stdmap(z, K)
    return z


def stdmap_jac_n(z, K, n):
    J = np.eye(2)
    for _ in range(n):
        c = K * np.cos(z[0])
        J = np.array([[1 + c, 1.0], [c, 1.0]]) @ J
        z = stdmap(z, K)
    return J


def _action_grad_hess(x, K, m, n):
    """Gradient and Hessian of W = sum (x_{i+1}-x_i)^2/2 - K cos x_i, with x_n = x_0 + 2 pi m."""
    xp = np.roll(x, -1); xp[-1] += 2 * np.pi * m
    xm = np.roll(x, 1); xm[0] -= 2 * np.pi * m
    g = 2 * x - xp - xm + K * np.sin(x)
    H = np.diag(2 + K * np.cos(x))
    for i in range(n):
        H[i, (i + 1) % n] -= 1
        H[i, (i - 1) % n] -= 1
    if n == 1:
        H = np.array([[K * np.cos(x[0])]])
        g = np.array([K * np.sin(x[0])])
    return g, H


def _newton_action(x, K, m, n, it=100):
    for _ in range(it):
        g, H = _action_grad_hess(x, K, m, n)
        step = np.linalg.lstsq(H, -g, rcond=1e-12)[0]
        x = x + step
        if np.max(np.abs(step)) < 1e-14:
            return x, True
    return x, False


def birkhoff_pair(K, m, n):
    """The minimizing and the minimax Birkhoff orbits of rotation number m/n (Aubry-Mather).

    Both visit the circle in the same cyclic order as the rotation by m/n, so they have the same
    braid type (Boyland). For each: trace of D(F^n), the fixed-point indices of F^n and F^{2n}
    (winding numbers, no derivatives), and the Morse index of the action at the orbit.
    """
    base = 2 * np.pi * m * np.arange(n) / n
    # the two symmetric Birkhoff orbits: seeds at phase 0 and phase pi/n (pi for n = 1)
    sols = []
    for phase in (0.0, np.pi / n):
        xx, ok = _newton_action(base + phase, K, m, n)
        _, H = _action_grad_hess(xx, K, m, n)
        sols.append((xx, ok, int(np.sum(np.linalg.eigvalsh(H) < 0))))
    sols.sort(key=lambda t: t[2])
    (xmin, ok1, _), (xmm, ok2, _) = sols
    shift = np.array([2 * np.pi * m, 0.0])
    out = []
    for name, xx, ok in (('minimizing', xmin, ok1), ('minimax', xmm, ok2)):
        p0 = xx[0] - (xx[-1] - 2 * np.pi * m) if n > 1 else 2 * np.pi * m
        z0 = np.array([xx[0], p0])
        orb = [z0]
        for _ in range(n - 1):
            orb.append(stdmap(orb[-1], K))
        orb = np.array(orb)
        Fn = lambda w: stdmap_n(w, K, n) - shift
        F2n = lambda w: stdmap_n(w, K, 2 * n) - 2 * shift
        _, H = _action_grad_hess(xx, K, m, n)
        morse = int(np.sum(np.linalg.eigvalsh(H) < 0))
        xs = np.mod(orb[:, 0], 2 * np.pi)
        ranks = np.argsort(np.argsort(xs))
        perm = tuple(int(ranks[(i + 1) % n]) for i in np.argsort(xs))
        tr = float(np.trace(stdmap_jac_n(orb[0], K, n)))
        D = np.linalg.norm(orb[:, None] - orb[None], axis=2) + 9 * np.eye(n)
        rho = min(1e-6, D.min() / 20)
        out.append({'kind': name, 'ok': ok and np.linalg.norm(Fn(z0) - z0) < 1e-8, 'orbit': orb, 'trace': tr,
                    'i1': winding_index(Fn, z0, rho), 'i2': winding_index(F2n, z0, rho / 10),
                    'morse': morse, 'perm': perm, 'x': xx})
    return out
