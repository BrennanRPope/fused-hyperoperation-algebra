"""Checks for the shape-dial and temperature tests (see shape_dial.py for the method).

Run: python3 shape_dial_checks.py > shape_dial_checks_output.txt
Khat = 3 x Gaussian curvature of the canonical (Vinberg) metric: disk = -1 (hyperbolic, relativity), triangle = 0
(flat, Bayes).  The leftover turn around a small loop is -Khat x area, so Khat is the turn per unit area.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.spatial import ConvexHull
from shape_dial import *

results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")


print("== 1. The method reproduces the two ends ==")
D = disk()
Yd = np.array([[0, 0], [0.3, 0.1], [-0.5, 0.4], [0.1, -0.75]])
k = khat(D, Yd)
check("V1 disk: Khat = -1 everywhere (Klein model of hyperbolic space)", np.allclose(k, -1, atol=2e-3), np.round(k, 4))
Yt = np.array([[0, 0], [0.1, 0.05], [-0.2, -0.2], [0, 0.6]])
k = poly_khat(TRI_V, Yt)
check("V2 triangle: Khat = 0 everywhere (flat)", np.allclose(k, 0, atol=2e-3), np.round(k, 5))
# chi formula: polar-body area of (disk - y) is pi/(1-|y|^2)^(3/2)
y = np.array([[0.4, -0.2]])
area = np.exp(D.log_chi(y))[0] / 2
r2 = (y ** 2).sum()
check("V3 chi = area of polar body: disk gives pi/(1-|y|^2)^(3/2)", abs(area - np.pi / (1 - r2) ** 1.5) < 1e-6 * area, f"{area:.6f}")
# metric at the disk centre: g = 3 x (Klein metric) = 3 I at the origin? Klein metric at 0 is I (curvature -1).
g0 = D.metric(np.array([[0.0, 0.0]]))[0]
check("V4 disk centre metric = 3 x identity (normalisation behind Khat = 3K)", np.allclose(g0, 3 * np.eye(2), atol=1e-6), np.round(g0, 6))
# invariance: a projective map of the triangle (diagonal scaling in homogeneous coords) keeps Khat = 0 and moves the
# centre anywhere; for the disk, curvature does not depend on the radius (the plane dial's scale).
for R in (0.2, 5.0):
    k = khat(disk(R), np.array([[0, 0], [0.4 * R, 0.3 * R]]))
    check(f"V5 disk of radius {R}: still Khat = -1 (size = plane dial, does not touch shape)", np.allclose(k, -1, atol=2e-3), np.round(k, 4))

print("\n== 2. Shape dial with the plane dial fixed: Minkowski path disk -> triangle ==")
cent = []
for t in (0, 0.25, 0.5, 0.75, 1.0):
    B = minkowski(t)
    kc = khat(B, np.array([[0.0, 0.0]]))[0]
    cent.append(kc)
    print(f"   t={t:4}: centre Khat = {kc:+.4f}")
check("M1 centre Khat runs from -1 to 0", abs(cent[0] + 1) < 5e-3 and abs(cent[-1]) < 1e-2)
check("M2 and moves monotonically in between", np.all(np.diff(cent) > 0), np.round(cent, 3))
B = minkowski(0.5)
ys = np.array([[0, 0.0], [0, 0.3], [0, -0.3], [0, 0.42], [0, -0.42]])  # +y toward a rounded corner, -y toward a flat edge
kk = khat(B, ys)
print("   t=0.5 along the axis (corner side +y, edge side -y):", np.round(kk, 3))
check("M3 at t=0.5 Khat is a field, not a number: spread across the body > 0.3", kk.max() - kk.min() > 0.3,
      f"range {kk.min():.3f} .. {kk.max():.3f}")
check("M4 no single plane-dial setting can copy this (any d gives a constant Khat on a round body)", kk.std() > 0.1)

print("\n== 3. Second path: triangle cut by a disk of radius s ==")
cc = []
for s in (0.5, 0.6, 0.7, 0.8, 0.9, 0.97, 1.0):
    B = cut_triangle(s) if s < 1 else None
    kc = khat(B, np.array([[0.0, 0.0]]))[0] if B is not None else poly_khat(TRI_V, np.array([[0.0, 0.0]]))[0]
    cc.append(kc)
    print(f"   s={s:4}: centre Khat = {kc:+.4f}")
check("C1 cut path also runs -1 -> 0 monotonically", abs(cc[0] + 1) < 5e-3 and abs(cc[-1]) < 1e-2 and np.all(np.diff(cc) > 0))

print("\n== 4. Temperature: free-energy body (1/beta) log mean exp(-beta x) <= 1, hot = disk, cold = triangle ==")
betas = [1e-3, 1e-2, 0.1, 1, 3, 10, 30, 100, 300, 1000]
kt = []
for b in betas:
    kc = khat(free_energy_body(b), np.array([[0.0, 0.0]]))[0]
    kt.append(kc)
    print(f"   beta={b:7}: centre Khat = {kc:+.6f}")
kt = np.array(kt)
check("T1 hot end is round (Khat -> -1)", abs(kt[0] + 1) < 2e-3)
check("T2 cold end is the triangle (Khat -> 0)", abs(kt[-1]) < 5e-3)
check("T3 the fade is smooth and monotone (no jump, no overshoot)", np.all(np.diff(kt) > 0) and np.all(kt < 0))
slope_cold = kt[-3:] * np.array(betas[-3:])
check("T4 cold tail ~ -c/beta (turn per area ~ 1/beta)", np.ptp(slope_cold) < 0.25 * abs(slope_cold.mean()), np.round(slope_cold, 3))
slope_hot = (kt[:3] + 1) / np.array(betas[:3])
check("T5 hot start ~ -1 + c*beta", np.ptp(slope_hot) < 0.15 * slope_hot.mean(), np.round(slope_hot, 3))

print("\n== 5. Regular polygons and l^p balls: only the triangle is flat at its centre ==")
for n in (3, 4, 5, 6, 8, 12):
    a = np.radians(90) + TWO_PI * np.arange(n) / n
    V = np.column_stack([np.cos(a), np.sin(a)])
    kc = poly_khat(V, np.array([[0.0, 0.0]]))[0]
    print(f"   regular {n:2}-gon: centre Khat = {kc:+.4f}")
    if n == 3:
        check("P1 triangle centre = 0", abs(kc) < 1e-3)
    else:
        check(f"P{n} {n}-gon centre = -1 (symmetry kills the cubic form there)", abs(kc + 1) < 2e-3)
for p in (1.5, 4, 10):
    t = np.linspace(0, TWO_PI, 8000, endpoint=False)
    c, s = np.cos(t), np.sin(t)
    bd = np.column_stack([np.sign(c) * np.abs(c) ** (2 / p), np.sign(s) * np.abs(s) ** (2 / p)])
    B = Body(bd)
    k = khat(B, np.array([[0, 0], [0.4, 0.2], [0.6, 0.6]]))
    check(f"L{p} l^{p} ball: centre -1, elsewhere inside [-1, 0]", abs(k[0] + 1) < 5e-3 and np.all((k > -1.01) & (k < 0.01)), np.round(k, 3))

print("\n== 6. Random bodies: Khat always lies between -1 (round) and 0 (triangle) ==")
rng = np.random.default_rng(3)
allk = []
for trial in range(30):  # random polygons with >= 4 sides, exact arcs
    while True:
        P = rng.normal(size=(rng.integers(4, 11), 2))
        V = P[ConvexHull(P).vertices]
        if len(V) >= 4:
            break
    c = V.mean(0)
    Y = np.vstack([c + (V - c) * rng.uniform(0.05, 0.8, (len(V), 1)), c])
    allk += list(poly_khat(V, Y))
for trial in range(15):  # random smooth bodies: radius 1 + small Fourier modes, made convex by hull
    t = np.linspace(0, TWO_PI, 6000, endpoint=False)
    r = 1 + sum(rng.normal(0, 0.06) * np.cos(m * t + rng.uniform(0, TWO_PI)) for m in range(2, 6))
    P = np.column_stack([r * np.cos(t), r * np.sin(t)])
    P = P[ConvexHull(P).vertices]
    B = Body(P)
    Y = rng.uniform(-0.6, 0.6, (6, 2))
    Y = Y[B.inside(Y, 0.2)]
    allk += list(khat(B, Y))
allk = np.array(allk)
check("B1 every sampled point has -1 <= Khat <= 0", allk.min() > -1.005 and allk.max() < 0.005,
      f"{len(allk)} points, min {allk.min():.4f}, max {allk.max():.4f}")

print(f"\n{sum(results)}/{len(results)} checks pass")
