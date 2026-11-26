"""Checks for the topological shapes round (DRAFT, Nov 26 2026).

T1  catalog: chi, orientability, holes and twists match the classification of surfaces
T2  connected sum: Euler genus adds; orientable iff both are; Dyck's T # RP2 = 3 RP2 (same chi and orientability)
T3  Gauss-Bonnet on the polygon: hyperbolic area = -2 pi chi, regular polygon exists exactly at the sign of chi
T4  meshes we draw have the holes we say (marching cubes chi = 2 - 2g for g = 1, 2, 3)
T5  every convex shape on the page is one type: solid chi = 1, surface chi = 2, at every p and depth tried
T6  p does not change topology (lp torus chi = 0 for p = 1.3 .. inf), and total Gauss curvature of the round torus is 0
T7  aspect ratio: ring torus chi = 0, spindle torus chi = 2 (solid 0 -> 1): genus changes only at the horn R = r
T8  chi as the t = -1 reading: products multiply chi, so 1 - chi composes as x + y - xy (the k = -1 fused step)
T9  Clifford torus (p = 2 product of circles in R4) lies on the unit 3-sphere and is flat, matching chi = 0
"""
import itertools
import json
import os
import time

import numpy as np
from skimage import measure

import topo_shapes as T
import facedial as FD

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
RES = {}
t0 = time.time()


def check(cid, desc, ok, detail=""):
    OUT.append(f"[{'PASS' if ok else 'FAIL'}] {cid}: {desc}" + (f"  ({detail})" if detail else ""))
    RES[cid] = {"desc": desc, "ok": bool(ok), "detail": detail}


def mesh_chi(F, level=0.0, spacing=(1, 1, 1)):
    verts, faces, _, _ = measure.marching_cubes(F, level, spacing=spacing)
    edges = set()
    for a, b, c in faces:
        for e in ((a, b), (b, c), (c, a)):
            edges.add((min(e), max(e)))
    return len(verts) - len(edges) + len(faces)


# ---------------------------------------------------------------- T1
expect = {"sphere": (2, True, 0), "projective plane": (1, False, 1), "torus": (0, True, 1),
          "Klein bottle": (0, False, 2), "genus-2 surface": (-2, True, 2), "genus-3 surface": (-4, True, 3),
          "Dyck's surface": (-1, False, 3)}
got = {k: (s.chi, s.orientable, s.holes if s.orientable else s.twists) for k, s in T.CATALOG.items()}
check("T1", "catalog chi, orientability, holes/twists match the classification",
      got == expect, "; ".join(f"{k} chi={v[0]}" for k, v in got.items()))

# ---------------------------------------------------------------- T2
ok = True
for (a, S), (b, U) in itertools.product(T.CATALOG.items(), repeat=2):
    M = T.connected_sum(S, U)
    ok &= (M.euler_genus == S.euler_genus + U.euler_genus) and (M.orientable == (S.orientable and U.orientable))
tor, rp = T.CATALOG["torus"], T.CATALOG["projective plane"]
dy = T.connected_sum(tor, rp)
d3 = T.CATALOG["Dyck's surface"]
ok_dyck = (dy.chi, dy.orientable) == (d3.chi, d3.orientable)
check("T2", "connected sum: Euler genus adds and orientable iff both; T # RP2 matches RP2#RP2#RP2 (Dyck)",
      ok and ok_dyck, f"49 pairs; T#RP2 chi={dy.chi}, orientable={dy.orientable}")

# ---------------------------------------------------------------- T3
rows = []
ok = True
words = {"sphere": "a b b' a'", "projective plane": "a b a b"}  # 4-gon versions avoid the degenerate bigon
for k, S in T.CATALOG.items():
    S4 = T.Surface(words[k], k) if k in words else S
    n, al = S4.n, S4.corner_angle()
    flat_gap = (n - 2) * np.pi - n * al           # hyperbolic area at kappa = -1
    ok &= np.isclose(flat_gap, -2 * np.pi * S4.chi)
    sgn = S4.curvature_sign()
    c = 1 / np.tan(np.pi / n) / np.tan(al / 2) if not np.isclose(al, np.pi * 2) else -np.inf
    exists = {-1: c > 1 + 1e-9, 0: np.isclose(c, 1), 1: -1 - 1e-9 <= c < 1 - 1e-9}
    ok &= exists[sgn] and sum(bool(v) for v in exists.values()) == 1
    rows.append(f"{k}: {n}-gon angle {np.degrees(al):.0f}deg, curvature {['negative','flat','positive'][sgn+1]}")
check("T3", "Gauss-Bonnet: each gluing polygon fits curvature of the sign of chi, and only that sign", ok, "; ".join(rows))
RES["T3"]["rows"] = rows

# ---------------------------------------------------------------- T4
g_chis = []
L = np.linspace(-3.6, 3.6, 150)
X, Y, Z = np.meshgrid(L, L * 0.55, L * 0.3, indexing="ij")
for g in (1, 2, 3):
    D = T.genus_sdf(g, X, Y, Z)
    g_chis.append(mesh_chi(D))
check("T4", "drawn genus-g surfaces have mesh chi = 2 - 2g", g_chis == [0, -2, -4], f"chi = {g_chis}")

# ---------------------------------------------------------------- T5
L = np.linspace(-1.9, 1.9, 44)
G = np.stack(np.meshgrid(L, L, L, indexing="ij"), -1).reshape(-1, 3)
shapes = {"cube": FD.cube(3), "octahedron": FD.cross_polytope(3), "tetrahedron": FD.simplex(3),
          "icosahedron": FD.icosahedron(), "dodecahedron": FD.dodecahedron()}
rows = []
ok = True
for nm, V in shapes.items():
    fd = FD.FaceDial(V, nm)
    ext = 1.25 * fd.radial(FD.sphere_dirs(3, 400, 1), 1.0).max() / 1.9
    for p in (1.0, 2.0, 5.0, np.inf):
        val = fd.gauge_value(G * ext, p).reshape(len(L), len(L), len(L)) - (fd.a + 1)
        s_chi = T.euler_char_voxels(val <= 0)
        b_chi = mesh_chi(val)
        ok &= (s_chi, b_chi) == (1, 2)
    rows.append(nm)
# depth: Minkowski mixes of two corners are convex too
from scipy.spatial import ConvexHull, Delaunay
c1, c2 = FD.cube(3), FD.cuboctahedron() if hasattr(FD, "cuboctahedron") else FD.cross_polytope(3)
for f in (0.25, 0.5, 0.75):
    P = np.array([(1 - f) * x + f * y for x in c1 for y in c2])
    hull = Delaunay(P[ConvexHull(P).vertices])
    inside = (hull.find_simplex(G * 1.0) >= 0).reshape(len(L), len(L), len(L))
    ok &= T.euler_char_voxels(inside) == 1
check("T5", "every convex shape is one type: solid chi 1, surface chi 2 (5 solids x p in 1,2,5,inf; 3 depth mixes)",
      ok, ", ".join(rows))

# ---------------------------------------------------------------- T6
L = np.linspace(-1.75, 1.75, 120)
X, Y, Z = np.meshgrid(L, L, L * 0.5, indexing="ij")
chis = {}
for p in (1.3, 2.0, 4.0, np.inf):
    def mp(a, b):
        return np.maximum(np.abs(a), np.abs(b)) if np.isinf(p) else (np.abs(a) ** p + np.abs(b) ** p) ** (1 / p)
    D = mp(mp(X, Y) - 1.1, Z) - 0.4
    chis[str(p)] = mesh_chi(D)
# total Gauss curvature of the round torus: K dA = cos v / (R + r cos v) * r (R + r cos v) du dv
R, r = 1.0, 0.4
v = np.linspace(0, 2 * np.pi, 4001)
totK = 2 * np.pi * np.trapezoid(r * np.cos(v), v) if hasattr(np, "trapezoid") else 2 * np.pi * np.trapz(r * np.cos(v), v)
check("T6", "p does not change topology: lp torus chi = 0 at every p; round torus total curvature = 0",
      all(c == 0 for c in chis.values()) and abs(totK) < 1e-9, f"chi {chis}; total K = {totK:.1e}")

# ---------------------------------------------------------------- T7
ratios = [2.0, 1.4, 1.1, 0.9, 0.6, 0.3]
L = np.linspace(-1.6, 1.6, 130)
X, Y, Z = np.meshgrid(L, L, L, indexing="ij")
rows = []
for q in ratios:
    r = 0.5; R = q * r
    D = np.sqrt((np.sqrt(X ** 2 + Y ** 2) - R) ** 2 + Z ** 2) - r
    rows.append((q, mesh_chi(D), T.euler_char_voxels(D <= 0)))
ok = all((b, s) == ((0, 0) if q > 1 else (2, 1)) for q, b, s in rows)
check("T7", "aspect ratio R/r: ring (R>r) has a hole, spindle (R<r) has none; the change is at the horn R = r",
      ok, "; ".join(f"R/r={q}: surface chi {b}, solid chi {s}" for q, b, s in rows))
RES["T7"]["rows"] = rows

# ---------------------------------------------------------------- T8
chi_list = [2, 1, 0, -1, -2, -4]
ok = True
for a, b in itertools.product(chi_list, repeat=2):
    x, y = 1 - a, 1 - b
    ok &= (1 - a * b) == x + y - x * y
# product cell complexes: f-vectors convolve, so the alternating sum (t = -1 reading) multiplies
def fpoly_chi(fv):
    return sum((-1) ** i * c for i, c in enumerate(fv))
circle = [4, 4]           # square loop: 4 vertices, 4 edges
sphere = [6, 12, 8]       # octahedron surface
torus_f = np.convolve(circle, circle)
s2s2 = np.convolve(sphere, sphere)
ok &= fpoly_chi(torus_f) == 0 and fpoly_chi(s2s2) == 4
check("T8", "chi multiplies on products (f-vectors convolve, t=-1 reading), so 1-chi composes as x+y-xy",
      ok, f"circle x circle chi = {fpoly_chi(torus_f)}, S2 x S2 chi = {fpoly_chi(s2s2)}")

# ---------------------------------------------------------------- T9
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, 2 * np.pi, 60))
P = np.stack([np.cos(u), np.sin(u), np.cos(v), np.sin(v)], -1) / np.sqrt(2)
on_sphere = np.abs(np.linalg.norm(P, axis=-1) - 1).max()
# first fundamental form: E = |Pu|^2, F = Pu.Pv, G = |Pv|^2 constant -> flat
Pu = np.stack([-np.sin(u), np.cos(u), 0 * u, 0 * u], -1) / np.sqrt(2)
Pv = np.stack([0 * v, 0 * v, -np.sin(v), np.cos(v)], -1) / np.sqrt(2)
E, Fm, Gm = (Pu * Pu).sum(-1), (Pu * Pv).sum(-1), (Pv * Pv).sum(-1)
flat = np.ptp(E) < 1e-12 and np.abs(Fm).max() < 1e-12 and np.ptp(Gm) < 1e-12
check("T9", "Clifford torus (p=2 circles x circles in R4) sits on the unit 3-sphere with a constant (flat) metric",
      on_sphere < 1e-12 and flat, f"max | |x|-1 | = {on_sphere:.1e}; E = G = {E.mean():.3f}, F = 0")

# ---------------------------------------------------------------- write
n_ok = sum(v["ok"] for v in RES.values())
OUT.append(f"\n{n_ok}/{len(RES)} pass in {time.time() - t0:.0f} s")
txt = "\n".join(OUT)
print(txt)
open(os.path.join(HERE, "topo_checks_output.txt"), "w").write(txt + "\n")
json.dump(RES, open(os.path.join(HERE, "topo_results.json"), "w"), indent=1, default=str)
