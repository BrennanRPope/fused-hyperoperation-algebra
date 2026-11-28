# Logic and set theory through the circle operation (Nov 25, 2026)

> **Names:** some words or symbols on this page (mult-add, fused, rung, wall) also mean something else in standard math, physics or computing, or inside this project. [NAMING-COLLISIONS.md](../../NAMING-COLLISIONS.md) says which meaning is in play.

Page title: "Only two notches on the dial keep true and false: OR and XOR" (source index.html, not yet published).

## Files
- index.html: the page (dial and truth table, wall fans, Hamacher line on the landscape, other t-norm families,
  Euler-characteristic wall, rung ladder, known/ours, verdict).
- checks.py -> checks_output.txt: 67/67 pass (sympy, mpmath, exact fractions).
- lean/LogicSets.lean: 35 theorems, compiles clean (see lean/README.md).

## Findings
1. For idempotents, (a ∘ₖ b)² − (a ∘ₖ b) = (k+1)(k+2)ab: only k = −1 (OR) and k = −2 (XOR) keep truth values,
   i.e. bases 0 and −1 in q-number terms ([n]_0 = 1, [n]_{−1} = parity). Classical algebra of idempotents; framing ours.
2. Distributivity and absorption leak by k·abc(1−a) and ab(1−a): Boolean algebra = the dial pinned to 0, 1.
3. Walls are truth values: OR's wall is 1 (true), XOR's is 1/2 (fair coin = Shannon's one-time pad).
4. Hamacher t-conorms (1978) are the landscape line p + d = −1 (statistics page's line), and
   odds(S_γ(a,b)) = odds a ∘_γ odds b: Hamacher's γ is the dial on odds (implicit in his generator).
   Valid as fuzzy "or" for γ ≥ 0; a copula only for 0 ≤ γ ≤ 2. Frank, Yager, Dombi, Schweizer–Sklar,
   Aczél–Alsina leave the plane except where they meet Hamacher (Dombi 1 = SS −1 = Hamacher 0; AA 1 = Frank→1 = prob. sum).
5. Sets: |A ⊔ B ⊔ A×B| = a ∘₁ b; no wall for finite sets, but with Euler characteristic (Schanuel 1991) the open
   interval is the wall (its lift [0,1) counts 0) and two open intervals are the other unit −2.
6. Rung 2 (fused exponent) alone defines the circle operation ∘₁ (called mult-add before 1.0), shifted add and 0, so its theory is undecidable; rung 0 and
   rung 1 alone are decidable (Presburger, Skolem/Mostowski). Folklore for x^y, proved here in fused form.

## Verdict
A small lead in fuzzy logic: the dial applied in other coordinates gives associative rules, and asking which
coordinates land back on rational rules could re-derive Hamacher's theorem in one line or expose a gap.
The Boolean and rung parts are known facts re-read through the lift; tidy, closed.
