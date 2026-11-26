"""Every shape this page aims the dial at, by family (DRAFT)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "exceptional-shapes"))
import shapes as EX  # 24-cell, 600-cell, 120-cell constructions (exceptional-shapes page)
from facedial import (regular_polygon, simplex, cube, cross_polytope, icosahedron, dodecahedron, permutohedron,
                      loday_associahedron, zonotope, cyclic_polytope, birkhoff, hypersimplex)

REGULAR = {}
for n in range(3, 9):
    REGULAR[f"{n}-gon"] = regular_polygon(n)
REGULAR["12-gon"] = regular_polygon(12)
REGULAR.update({
    "tetrahedron": simplex(3), "cube": cube(3), "octahedron": cross_polytope(3),
    "dodecahedron": dodecahedron(), "icosahedron": icosahedron(),
    "5-cell": simplex(4), "tesseract": cube(4), "16-cell": cross_polytope(4),
    "24-cell": EX.binary_tetrahedral(), "120-cell": EX.onetwenty_cell(), "600-cell": EX.binary_icosahedral(),
    "5-simplex": simplex(5), "5-cube": cube(5), "5-orthoplex": cross_polytope(5),
    "6-simplex": simplex(6), "6-cube": cube(6), "6-orthoplex": cross_polytope(6),
})

rng = np.random.default_rng(7)
NONREGULAR = {
    "permutohedron (3D)": permutohedron(4),
    "permutohedron (4D)": permutohedron(5),
    "associahedron (3D, Loday)": loday_associahedron(4),
    "associahedron (4D, Loday)": loday_associahedron(5),
    "zonotope (3D, 5 random zones)": zonotope(rng.standard_normal((5, 3))),
    "cyclic polytope C(7,3)": cyclic_polytope(7, 3),
    "cyclic polytope C(8,4)": cyclic_polytope(8, 4),
    "Birkhoff B3 (4D)": birkhoff(3),
    "hypersimplex D(2,5) (4D)": hypersimplex(2, 5),
    "hypersimplex D(2,6) (5D)": hypersimplex(2, 6),
}
