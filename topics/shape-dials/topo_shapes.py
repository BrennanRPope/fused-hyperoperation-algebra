"""Topological shapes as first-class citizens (DRAFT, Nov 26 2026).

Brennan: "I wonder if we can start adding common topological shapes as well?" ... "I think all of these shapes are
just one type of shape" ... "And a hole is another one".

Every convex shape on the shape-dials page is, topologically, one shape: a ball (its surface a sphere). The next types
are counted by holes (handles) and twists (crosscaps). The classical construction makes them from shapes we already
have: a polygon with its edges glued in pairs by a word (Möbius 1860s, Dyck 1888, Dehn–Heegaard 1907; the fundamental
polygon). Every closed surface arises this way (classification of surfaces).

    Surface("a b a' b'")   -> torus        (a' means a run backwards)
    Surface("a a")         -> projective plane
"""
import itertools
import numpy as np


class Surface:
    """A closed surface given by a gluing word around a polygon: each letter appears twice."""

    def __init__(self, word, name=""):
        self.name = name
        toks = word.split() if isinstance(word, str) else list(word)
        self.edges = [(t.rstrip("'"), -1 if t.endswith("'") else 1) for t in toks]
        self.n = len(self.edges)
        letters = [e[0] for e in self.edges]
        assert all(letters.count(l) == 2 for l in set(letters)), "each letter must appear exactly twice"
        self.word = " ".join(toks)

    # vertices of the polygon: edge i runs from corner i to corner i+1
    def vertex_classes(self):
        n = self.n
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            parent[find(a)] = find(b)

        first = {}
        for i, (l, s) in enumerate(self.edges):
            tail, head = (i, (i + 1) % n) if s == 1 else ((i + 1) % n, i)
            if l in first:
                t0, h0 = first[l]
                union(t0, tail); union(h0, head)
            else:
                first[l] = (tail, head)
        return len({find(i) for i in range(n)})

    @property
    def V(self):
        return self.vertex_classes()

    @property
    def E(self):
        return self.n // 2

    @property
    def chi(self):
        """Euler characteristic V - E + F (one face)."""
        return self.V - self.E + 1

    @property
    def orientable(self):
        signs = {}
        for l, s in self.edges:
            signs.setdefault(l, []).append(s)
        return all(sorted(v) == [-1, 1] for v in signs.values())

    @property
    def euler_genus(self):
        """2 - chi: twice the holes if orientable, the number of twists (crosscaps) if not. Adds under #."""
        return 2 - self.chi

    @property
    def holes(self):
        return self.euler_genus // 2 if self.orientable else None

    @property
    def twists(self):
        return None if self.orientable else self.euler_genus

    def corner_angle(self):
        """Each corner's angle when all corners glued to one vertex class share 2π: α = 2π V / n."""
        return 2 * np.pi * self.V / self.n

    def curvature_sign(self):
        """Gauss–Bonnet on the polygon: κ · area = 2π χ, so the curvature that fits has the sign of χ."""
        return int(np.sign(self.chi))

    def describe(self):
        kind = (f"{self.holes} hole{'s' if self.holes != 1 else ''}" if self.orientable
                else f"{self.twists} twist{'s' if self.twists != 1 else ''}")
        return f"{self.name}: word {self.word}, {self.n}-gon, χ = {self.chi}, {kind}"


def relabel(word, suffix):
    out = []
    for t in word.split():
        b = t.rstrip("'")
        out.append(b + suffix + ("'" if t.endswith("'") else ""))
    return " ".join(out)


def connected_sum(S, T, name=""):
    """M # N: cut a disk from each and join along the cut. On words: write one after the other."""
    return Surface(relabel(S.word, "1") + " " + relabel(T.word, "2"), name or f"{S.name} # {T.name}")


def genus_word(g):
    return " ".join(f"a{i} b{i} a{i}' b{i}'" for i in range(1, g + 1))


def crosscap_word(c):
    return " ".join(f"c{i} c{i}" for i in range(1, c + 1))


CATALOG = {
    "sphere": Surface("a a'", "sphere"),
    "projective plane": Surface("a a", "projective plane"),
    "torus": Surface("a b a' b'", "torus"),
    "Klein bottle": Surface("a b a b'", "Klein bottle"),
    "genus-2 surface": Surface(genus_word(2), "genus-2 surface"),
    "genus-3 surface": Surface(genus_word(3), "genus-3 surface"),
    "Dyck's surface": Surface(crosscap_word(3), "Dyck's surface"),
}


# ------------------------------------------------------------------ regular polygons in constant curvature
def regular_polygon_radius(n, alpha, kappa_sign):
    """Circumradius of the regular n-gon with corner angle alpha at curvature +1, 0 or -1 (unit |κ|).
    Flat: needs alpha = (n-2)π/n (any size). Hyperbolic: cosh ρ = cot(π/n) cot(α/2). Spherical: cos ρ = same."""
    c = 1 / np.tan(np.pi / n) / np.tan(alpha / 2)
    if kappa_sign < 0:
        return np.arccosh(c)
    if kappa_sign > 0:
        return np.arccos(c)
    return None


def polygon_area(n, alpha, kappa_sign):
    """Area of the regular n-gon with corner angle alpha at |κ| = 1 (Gauss–Bonnet)."""
    if kappa_sign < 0:
        return (n - 2) * np.pi - n * alpha
    if kappa_sign > 0:
        return n * alpha - (n - 2) * np.pi
    return None


# ------------------------------------------------------------------ embedded pictures
def torus_of_revolution(R, r, nu=160, nv=80, p=2.0):
    """Torus with core circle of radius R and tube radius r; p sets both the core's and the tube's ℓp roundness."""
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, 2 * np.pi, nv)
    U, Vv = np.meshgrid(u, v)

    def lp_circle(t):
        c, s = np.cos(t), np.sin(t)
        if np.isinf(p):
            m = np.maximum(np.abs(c), np.abs(s))
        else:
            m = (np.abs(c) ** p + np.abs(s) ** p) ** (1 / p)
        return c / m, s / m

    cu, su = lp_circle(U)
    cv, sv = lp_circle(Vv)
    X = (R + r * cv) * cu
    Y = (R + r * cv) * su
    Z = r * sv
    return X, Y, Z


def mobius(nu=160, nv=24):
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(-0.45, 0.45, nv)
    U, Vv = np.meshgrid(u, v)
    X = (1 + Vv * np.cos(U / 2)) * np.cos(U)
    Y = (1 + Vv * np.cos(U / 2)) * np.sin(U)
    Z = Vv * np.sin(U / 2)
    return X, Y, Z


def klein_figure8(nu=160, nv=80, a=2.2):
    u = np.linspace(0, 2 * np.pi, nu)
    v = np.linspace(0, 2 * np.pi, nv)
    U, Vv = np.meshgrid(u, v)
    w = np.cos(U / 2) * np.sin(Vv) - np.sin(U / 2) * np.sin(2 * Vv)
    X = (a + w) * np.cos(U)
    Y = (a + w) * np.sin(U)
    Z = np.sin(U / 2) * np.sin(Vv) + np.cos(U / 2) * np.sin(2 * Vv)
    return X, Y, Z


def roman_surface(nu=120, nv=120):
    """Steiner's Roman surface: a picture of the projective plane in 3D (with pinch points)."""
    u = np.linspace(0, np.pi, nu)
    v = np.linspace(0, np.pi, nv)
    U, Vv = np.meshgrid(u, v)
    x, y, z = np.sin(U) * np.cos(Vv), np.sin(U) * np.sin(Vv), np.cos(U)
    return y * z, z * x, x * y


def genus_sdf(g, X, Y, Z, R=1.0, r=0.38, spacing=2.05, k=0.18):
    """Signed distance (smoothly joined) for g tori in a row: the boundary is a genus-g surface."""
    d = None
    for i in range(g):
        cx = (i - (g - 1) / 2) * spacing
        di = np.sqrt((np.sqrt((X - cx) ** 2 + Y ** 2) - R) ** 2 + Z ** 2) - r
        if d is None:
            d = di
        else:
            h = np.clip(0.5 + 0.5 * (di - d) / k, 0, 1)
            d = di * (1 - h) + d * h - k * h * (1 - h)
    return d


def euler_char_voxels(B):
    """Euler characteristic of a union of closed unit cubes (3D boolean array): V - E + F - C on the cubical complex."""
    B = np.pad(B.astype(bool), 1)
    C = B.sum()
    # faces, edges, vertices present if any adjacent cube is present
    def any_shift(shifts):
        out = np.zeros(tuple(s + 1 for s in B.shape), bool)
        for sh in shifts:
            sl = tuple(slice(1 - o, 1 - o + B.shape[i]) for i, o in enumerate(sh))
            tmp = np.zeros_like(out)
            tmp[sl] = B
            out |= tmp
        return out
    Vt = any_shift(list(itertools.product([0, 1], repeat=3))).sum()
    E = 0
    F = 0
    for axis in range(3):
        # edges parallel to axis: shared by up to 4 cubes, offsets in the other two axes
        others = [a for a in range(3) if a != axis]
        sh = []
        for o in itertools.product([0, 1], repeat=2):
            s = [0, 0, 0]; s[others[0]], s[others[1]] = o
            sh.append(tuple(s))
        e = any_shift(sh)
        # an edge parallel to axis exists at grid index along axis < B.shape[axis] (shape+1 in others)
        sl = [slice(None)] * 3; sl[axis] = slice(0, B.shape[axis])
        E += e[tuple(sl)].sum()
        # faces normal to axis: shared by 2 cubes along axis
        s = [0, 0, 0]; s[axis] = 1
        f = any_shift([(0, 0, 0), tuple(s)])
        sl = [slice(None)] * 3
        for a in others:
            sl[a] = slice(0, B.shape[a])
        F += f[tuple(sl)].sum()
    return int(Vt - E + F - C)
