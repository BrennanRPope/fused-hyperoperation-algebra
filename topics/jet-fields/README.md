# Fields of jets in mult-add (Nov 25, 2026)

Page: https://claude.ai/artifact/Lg84gxR34iuBqZwQNZWjF3 (slug 9f53fea2-9e95-4578-854c-316b2d06687e)
Rebuild: python3 figures.py && python3 page/build.py (writes page/index.html from page/template.html + SVGs).

Lean: lean/JetFields.lean, 46 theorems, core Lean 4.37.0-rc1, compiled against the shared library;
in a temp copy of the library all 1251 theorems pass Check (standard axioms, no sorry).
Library-ready, not yet folded in: copy to lean-library/FusedAlgebra/Topics/JetFields.lean and add
FusedAlgebra.Topics.JetFields to build.sh and FusedAlgebra.lean.
Note for the library: build.sh lists Topics.NSContinuum before Topics.Topology.CubeProofs, which it
imports; a clean build (empty build/) fails there until the order is swapped.

Findings (all known math; the cube reading is the framework's):
- k-jet = element of R[e]/e^(k+1); relative jets combine by x+y+xy (multiplicative formal group);
  jet order = cube depth; weight k <-> e^2 = k e, two sample points k apart, k->0 tangent.
- NS = surface in J^2; u.grad u = e-part of u(x + e u) (self-composition).
- det(I + h grad u) = 1 + hP + h^2 Q + h^3 R = mult-add cube of eigen-rates; Q = |w|^2/4 - |S|^2/2,
  R = det S + w.Sw/4; null Lagrangians -> <|w|^2> = 2<|S|^2>, Betchov <w.Sw> = -4<det S>.
- Restricted Euler: Q' = 3R, R' = -(2/3)Q^2 (R = det), 27R^2+4Q^3 conserved, tail |R|' ~ |R|^(4/3)
  fails Osgood (blow-up ~1.26 from |R|=1); 2D frozen.
- Open lead: is the far-corner exponent 4/3 the same quantity as the one-rung deficiency?
