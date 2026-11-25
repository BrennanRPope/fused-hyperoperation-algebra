# Group theory and Lie groups through the circle operation (Nov 25, 2026)

Files
- index.html: the page (self-contained; slide-rule ruler figure, landscape map, dial-as-dilation, finite clock).
- checks.py → checks_output.txt: 59 of 59 checks pass (sympy, mpmath, numpy).
- lean/GroupsLie.lean: 52 theorems, compiled clean (see lean/README.md).

Findings
- Every landscape law F_{p,q} is a one-parameter group of Möbius maps: mob a = I + aX, X = [[p,1],[-q,0]],
  mob a · mob b = (1 - qab) mob F(a,b). Generator field 1 + px + qx² (the landscape's Riccati rate). Known (Möbius/Lie).
- The landscape is the plane α = 1 in Lie's sl2 of fields α + βx + γx²; Killing form of the generator = 2(p² - 4q).
  Ours as a placement.
- Chart change x ↦ x/(1 - cx) sends (p,q) to (p+2c, q+pc+c²), keeping p² - 4q; scaling multiplies it by λ².
  So three classes: line (two walls), border (one double wall), circle. Classical (hyperbolic/parabolic/elliptic).
  Parity (-2,0) is velocity addition re-charted (c = 1).
- θ-lift: (1 + aθ)(1 + bθ) = (1 - qab) + (a + b + pab)θ with θ² = pθ - q; norm = rate, so ∂F/∂b · V(b) = V(F)
  (invariant differential). Same object as the linear-algebra page's dial matrix U.
- Closed-form logarithm (ruler) for every point; checked against numerical integration.
- Dial: T_a(x) = (1 + ka)x + a is a dilation about the wall -1/k in the ax+b group; k → 0 is the Inönü–Wigner
  style contraction to translations. In 1-D the Lie algebra is always the line, so k is units + wall position.
- Finite clocks: (F_p minus wall, ∘₁) cyclic (Gauss); 1 ∘₁ … ∘₁ 1 = 2ⁿ - 1, so 1 generates iff 2 is a primitive root
  (Artin's conjecture). Circle (tan law) clock has p + 1 elements for p ≡ 3 mod 4.
- The matrix circle operation ∘₁ (called mult-add before 1.0): A ∘₁ B - B ∘₁ A = AB - BA exactly; group commutator = t²[A,B] + O(t³). Magnus 1937.

Verdict: not a lead on its own. It names the landscape correctly (one-parameter Möbius groups, Lie's sl2) and
predicts nothing new, since 1-D groups are fully classified. The live direction is the non-commuting corner
(matrices / two boost directions), already being followed on the twisted boosts page.
