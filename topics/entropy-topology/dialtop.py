"""Shared tools for topics/entropy-topology (DRAFT).

The two-dial law  x o a = (x + a + k*a*x) / (1 - d*a*x)  acts on the real projective line RP^1,
which is a circle. Everything here works in projective coordinates [u : v] (x = u/v) so the point
at infinity needs no special case.
"""
import numpy as np

# ---------- the law on RP^1 ----------

def compose_h(u1, v1, u2, v2, k, d):
    """[u1:v1] o [u2:v2] in homogeneous form."""
    return (u1 * v2 + u2 * v1 + k * u1 * u2, v1 * v2 - d * u1 * u2)


def dial_power_h(u, v, k, d, m):
    """m-fold self-composition x o x o ... o x (m copies), homogeneous."""
    U, V = u, v
    for _ in range(m - 1):
        U, V = compose_h(U, V, u, v, k, d)
    return U, V


def dial_square(x, k, d):
    return (2 * x + k * x * x) / (1 - d * x * x)


def newton_quadratic(y, k, d):
    """Newton's method on q(y) = y^2 + k y + d."""
    return (y * y - d) / (2 * y + k)

# ---------- circle degree of a map of RP^1 ----------

def rp1_degree(fh, n=200001):
    """Degree of a map of RP^1 given in homogeneous form fh(u, v) -> (U, V).

    RP^1 is parametrised by the angle th in [0, pi) with [u:v] = [sin th : cos th]. The image angle
    is taken mod pi; we follow it continuously (doubling the angle turns 'mod pi' into a winding
    number around 0) and count how many half-turns it makes.
    """
    th = np.linspace(0.0, np.pi, n)
    U, V = fh(np.sin(th), np.cos(th))
    w = (V + 1j * U) ** 2          # exp(2 i phi) up to a positive factor
    w = w / np.abs(w)
    dphi = np.angle(w[1:] / w[:-1])
    if np.max(np.abs(dphi)) > 2.5:
        raise ValueError("grid too coarse for degree")
    return int(round(np.sum(dphi) / (2 * np.pi)))

# ---------- real roots ----------

def distinct_real_roots(coeffs, tol=1e-7):
    """Distinct real roots of a polynomial (highest power first)."""
    r = np.roots(coeffs)
    real = np.sort(r[np.abs(r.imag) < tol].real)
    out = []
    for x in real:
        if not out or abs(x - out[-1]) > 1e-5 * (1 + abs(x)):
            out.append(x)
    return np.array(out)


def newton_h(coeffs):
    """Homogeneous form of Newton's map N(y) = y - q/q' for a real polynomial q."""
    q = np.poly1d(coeffs)
    dq = q.deriv()
    n = q.order
    # N(y) = (y q'(y) - q(y)) / q'(y); homogenise with degree n.
    num = np.poly1d([1, 0]) * dq - q          # degree n
    den = dq                                  # degree n-1
    num_c = np.concatenate([np.zeros(n + 1 - len(num.c)), num.c])
    den_c = np.concatenate([np.zeros(n + 1 - len(den.c)), den.c])

    def fh(u, v):
        U = sum(num_c[i] * u ** (n - i) * v ** i for i in range(n + 1))
        V = sum(den_c[i] * u ** (n - i) * v ** i for i in range(n + 1))
        return U, V
    return fh, q

# ---------- topological entropy of a real rational map on RP^1 by lap counting ----------

def lap_entropy(fh, crit_points, solve_preimages, depth, cap=400000):
    """Estimate entropy as growth of the number of turning points of f^j.

    Turning points of f^j on the circle are the points that land on a turning point of f within
    j-1 steps. crit_points: turning points of f (finite reals or np.inf). solve_preimages(c):
    real preimages of c under f. Returns the list of counts laps_j for j = 1..depth.
    """
    level = list(crit_points)
    seen_total = len(level)
    counts = [seen_total]
    for _ in range(depth - 1):
        nxt = []
        for c in level:
            nxt.extend(solve_preimages(c))
            if len(nxt) > cap:
                return counts
        level = nxt
        seen_total += len(level)
        counts.append(seen_total)
    return counts


def rp1_degree_alg(num, den, value=0.7316, tol=1e-9):
    """Degree of the real rational map y -> num(y)/den(y) on RP^1, counted algebraically.

    Pick a regular value c, solve num - c*den = 0 for its real preimages, and add up the signs of
    the derivative there (the map's orientation at each preimage). If infinity is also a preimage
    (when deg num = deg den and their leading ratio is c), we perturb c. num, den: np.poly1d.
    """
    num, den = np.poly1d(num), np.poly1d(den)
    for c in [value, 1.913, -0.4142, 3.1416, -2.718]:
        p = num - c * den
        n_top = max(num.order, den.order)
        if p.order < n_top:          # infinity maps to c: choose another value
            continue
        roots = np.roots(p.c)
        real = roots[np.abs(roots.imag) < 1e-7 * (1 + np.abs(roots))].real
        dnum, dden = num.deriv(), den.deriv()
        s = 0
        ok = True
        for y in real:
            a1, a2 = dnum(y) * den(y), num(y) * dden(y)
            dv = a1 - a2                               # sign of (num/den)' times den^2
            if abs(dv) < tol * (abs(a1) + abs(a2) + 1e-300):
                ok = False
                break
            s += int(np.sign(dv))
        if ok:
            return s
    raise ValueError("no regular value found")


def newton_parts(coeffs):
    q = np.poly1d(coeffs)
    dq = q.deriv()
    return np.poly1d([1, 0]) * dq - q, dq, q


def real_preimages(num, den, c):
    """Real y with num(y)/den(y) = c (c may be np.inf)."""
    p = den if np.isinf(c) else num - c * den
    if p.order < 1:
        return []
    r = np.roots(p.c)
    return list(r[np.abs(r.imag) < 1e-9 * (1 + np.abs(r))].real)


def newton_lap_counts(coeffs, depth, cap=300000):
    """Number of turning points of N^j on RP^1, j = 1..depth, for Newton's map of q."""
    num, den, q = newton_parts(coeffs)
    crit = list(distinct_real_roots(q.c)) + list(distinct_real_roots((q.deriv(2)).c)) \
        if q.order >= 3 else list(distinct_real_roots(q.c))
    seen = set()

    def key(y):
        return round(float(y), 9)
    level = []
    for c in crit:
        if key(c) not in seen:
            seen.add(key(c)); level.append(c)
    counts = [len(seen)]
    for _ in range(depth - 1):
        nxt = []
        for c in level:
            for y in real_preimages(num, den, c):
                if key(y) not in seen:
                    seen.add(key(y)); nxt.append(y)
        level = nxt
        counts.append(len(seen))
        if len(seen) > cap:
            break
    return counts


def newton_preimage_counts(coeffs, depth, cap=400000):
    """|N^{-j}(inf) on RP^1| for j = 0..depth, counting the fixed point inf once.

    inf is a repelling fixed point of Newton's map, so it lies in the real Julia set; the growth rate
    of its backward orbit is the entropy of N on the real line.
    """
    num, den, q = newton_parts(coeffs)
    level = [np.inf]
    counts = [1]
    total = 1
    for _ in range(depth):
        nxt = []
        for c in level:
            pre = real_preimages(num, den, c)
            nxt.extend(pre)             # inf's own preimage inf is not re-added
        level = nxt
        total += len(level)
        counts.append(total)
        if total > cap:
            break
    return counts


def newton_turning_counts(coeffs, depth, cap=400000):
    """Turning points of N^j on RP^1 (j = 1..depth), by the backward tree of N's turning points.

    N' = q q'' / q'^2, so the turning points of N are the real roots of q (fixed) and of q''.
    A fixed root's own preimage (itself) is not re-added.
    """
    num, den, q = newton_parts(coeffs)
    roots = list(distinct_real_roots(q.c))
    infl = list(distinct_real_roots(q.deriv(2).c)) if q.order >= 3 else []
    level = roots + infl
    counts = [len(level)]
    total = len(level)
    for _ in range(depth - 1):
        nxt = []
        for c in level:
            for y in real_preimages(num, den, c):
                if any(abs(y - r0) < 1e-9 * (1 + abs(r0)) and abs(c - r0) < 1e-12 * (1 + abs(r0)) for r0 in roots):
                    continue
                nxt.append(y)
        level = nxt
        total += len(level)
        counts.append(total)
        if total > cap:
            break
    return counts
