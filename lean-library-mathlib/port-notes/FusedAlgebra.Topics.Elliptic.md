# FusedAlgebra.Topics.Elliptic

- New definitions `EllipticLead.nodalCurve` and `EllipticLead.freyCurve` are mathlib
  `WeierstrassCurve`s, so that the discriminant and `c₄` statements use mathlib's `Δ`, `c₄`, `b₂`.
- New theorems with no original: `nodalCurve_Δ` (the nodal cubic has `Δ = 0`), `nodalCurve_c₄`,
  `freyCurve_equation` (the curve is `y² = x(x - A)(x + B)`), `gap6_one`.
- All polynomial identities are stated with `^` and hold in any commutative ring (original: `ℤ`).
