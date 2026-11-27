"""DRAFT. Derivative-free fixed-point index in R^4: the degree of g = F^k - id on the boundary of a small 4-cube.

The boundary (8 three-dimensional faces) is cut into n^3 small cubes per face and each into 6 Kuhn tetrahedra.
For a generic direction e, deg(g/|g|) = sum over tetrahedra whose image cone contains e of
sign det[g(v0..v3)] * sign det[v0-c .. v3-c] (the second factor is the outward orientation).
Only values of F are used, never its derivative."""
import itertools, numpy as np

def cube_boundary_tets(n):
    """vertices (integer grid coords in {0..n}^4) of all boundary tetrahedra."""
    tets = []
    perms = list(itertools.permutations(range(3)))
    for ax in range(4):
        for side in (0, n):
            others = [a for a in range(4) if a != ax]
            for base in itertools.product(range(n), repeat=3):
                for pm in perms:
                    v = [0] * 4; v[ax] = side
                    for i, a in enumerate(others): v[a] = base[i]
                    vs = [tuple(v)]
                    cur = list(v)
                    for j in pm:
                        cur[others[j]] += 1; vs.append(tuple(cur))
                    tets.append(vs)
    return tets

def degree(g, c, r, n=4, dirs=None, B=None):
    """g: function R^4 -> R^4; c: centre; r: half-width (scalar or 4-vector); B: optional 4x4 shape matrix
    (the box is c + B @ (r * [-1,1]^4)). Returns (list of degrees, one per direction; smallest |g| on the boundary)."""
    c = np.asarray(c, float); r = np.broadcast_to(np.asarray(r, float), (4,))
    B = np.eye(4) if B is None else np.asarray(B, float)
    tets = cube_boundary_tets(n)
    pts = sorted({v for t in tets for v in t})
    vals = {}
    for v in pts:
        x = c + B @ (r * (2 * np.array(v) / n - 1))
        vals[v] = g(x)
    if dirs is None:
        rng = np.random.default_rng(1)
        dirs = [d / np.linalg.norm(d) for d in rng.normal(size=(5, 4))]
    out = []
    for e in dirs:
        deg = 0
        for t in tets:
            W = np.column_stack([vals[v] for v in t])
            V = np.column_stack([B @ (r * (2 * np.array(v) / n - 1)) for v in t])
            try:
                lam = np.linalg.solve(W, e)
            except np.linalg.LinAlgError:
                continue
            if np.all(lam >= 0):
                deg += int(np.sign(np.linalg.det(W)) * np.sign(np.linalg.det(V)))
        out.append(deg)
    minnorm = min(np.linalg.norm(w) for w in vals.values())
    return out, minnorm

if __name__ == '__main__':
    # self-tests: identity-like maps
    rng = np.random.default_rng(0)
    for M, want in [(np.eye(4), 1), (np.diag([-1, 1, 1, 1]), -1), (np.diag([-1, -1, 1, 1]), 1), (np.diag([2, 3, -1, 1]), -1)]:
        d, _ = degree(lambda x: M @ x, np.zeros(4), 1.0, n=2)
        print(M.diagonal(), d, 'expected', want)

def robust_degree(g, c, r0, shrink=10.0, tries=5, B=None):
    """shrink the cube until two successive sizes, two subdivisions (n = 2, 3) and five directions all agree."""
    prev = None; hist = []
    r = r0
    for _ in range(tries):
        a, m = degree(g, c, r, n=2, B=B); b, _ = degree(g, c, r, n=3, B=B)
        vals = set(a) | set(b); hist.append((r, sorted(vals), m))
        cur = vals.pop() if len(vals) == 1 else None
        if cur is not None and cur == prev: return cur, hist
        prev = cur; r /= shrink
    return None, hist
