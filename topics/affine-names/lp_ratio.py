"""Blaschke's affine isoperimetric ratio Omega^3/(8 pi^2 A) of the l^p circle (1 at the ellipse, 0 at polygons)."""
import numpy as np
from scipy.integrate import quad


def lp_ratio(p):
    e = 2 / p

    def derivs(tt):
        c, s = np.cos(tt), np.sin(tt)
        x1 = -e * c ** (e - 1) * s
        y1 = e * s ** (e - 1) * c
        x2 = -e * ((e - 1) * c ** (e - 2) * (-s) * s + c ** (e - 1) * c)
        y2 = e * ((e - 1) * s ** (e - 2) * c * c + s ** (e - 1) * (-s))
        return c ** e, s ** e, x1, y1, x2, y2

    def aff(tt):
        _, _, x1, y1, x2, y2 = derivs(tt)
        return abs(x1 * y2 - y1 * x2) ** (1 / 3)

    def area(tt):
        x0, y0, x1, y1, _, _ = derivs(tt)
        return 0.5 * (x0 * y1 - y0 * x1)

    om = 4 * quad(aff, 0, np.pi / 2, limit=400)[0]
    ar = 4 * quad(area, 0, np.pi / 2, limit=400)[0]
    return om ** 3 / (8 * np.pi ** 2 * ar)
