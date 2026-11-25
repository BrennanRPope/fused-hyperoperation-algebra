# P vs NP: what can't work, and where the circle operation touches it

Nov 25, 2026. Lean module: `lean-library/FusedAlgebra/Topics/PvsNP.lean` (20 theorems,
core Lean 4.37, no mathlib, no `sorry`, standard axioms only; library now 108 modules,
2,735 theorems).

## The result in two sentences

The circle operation ∘₁ (called mult-add before 1.0) can't change P vs NP: on the shifted ruler it is ordinary multiplication, and
Lean now proves that any circuit built from ∘₁ and shifted add turns into an ordinary circuit
of exactly the same size, and back. Where it does say something real is the dial: folded
over a problem's witnesses, k = −1 asks "is there one?" (NP), k = 0 asks "how many?" (#P),
and k = −2 asks "odd or even?" (⊕P), which is a clean picture of known classes and a good
landscape page, but not a route to a proof.

Your hunch is the mainstream one: almost everyone expects P ≠ NP. And you're right that
there's no easy counterexample to look for. A counterexample would be an algorithm, and the
map below is mostly a list of whole families of algorithms and whole families of proofs that
are already known to fail.

---

## 1. Methods known not to prove P ≠ NP

Each barrier is a theorem. It doesn't say P = NP; it says a certain kind of argument can
never settle the question.

| Barrier | What it kills | Why, in one sentence |
|---|---|---|
| **Relativization** (Baker, Gill, Solovay 1975) | Diagonalization and simulation: the arguments behind the time and space hierarchy theorems | There is an oracle that makes P = NP and another that makes P ≠ NP, so any proof that would work unchanged with an oracle attached proves both, which is impossible. |
| **Natural proofs** (Razborov, Rudich 1994) | Almost every combinatorial circuit lower bound we have: random restrictions (parity not in AC⁰), the Razborov–Smolensky polynomial method (AC⁰[p]), and their relatives | A property of Boolean functions that is easy to check from the truth table and true of many functions, if it ruled out small circuits, would also tell pseudorandom functions apart from random ones, which would break cryptography we believe is secure. |
| **Algebrization** (Aaronson, Wigderson 2008) | Arithmetization plus diagonalization: the techniques that got past relativization (IP = PSPACE, MIP = NEXP) | Even when the oracle is replaced by its low-degree polynomial extension, there are worlds on both sides, so arithmetizing doesn't escape the relativization problem. |
| **Monotone bounds don't transfer** (Razborov 1985; Tardos 1988) | Proving a lower bound for circuits without NOT gates and hoping it carries over | Razborov proved CLIQUE needs huge monotone circuits, but he also proved perfect matching does, and matching is in P; Tardos found a problem in P whose monotone circuits are exponential, so monotone hardness says nothing about real hardness. |
| **Occurrence obstructions in GCT** (Bürgisser, Ikenmeyer, Panova, JAMS 2019) | The simplest form of Mulmuley–Sohoni geometric complexity theory for permanent vs determinant | The representations that were supposed to appear for the determinant and not the padded permanent always appear in both; the multiplicity version is still open. |
| **Proof-complexity limits** (Cook–Reckhow 1979; Razborov 1995) | Hoping a weak logical system can prove the lower bound, and hoping a weak proof system can certify every tautology | NP ≠ coNP is exactly "no proof system has short proofs of every tautology"; we have exponential lower bounds for resolution and bounded-depth Frege but none for full Frege, and Razborov showed weak bounded arithmetic can't prove strong circuit lower bounds if strong pseudorandom generators exist. |
| **Algebraic natural proofs** (Forbes–Shpilka–Volk; Grochow–Kumar–Saks–Saraf 2017) | The algebraic version of the natural-proofs strategy for VP vs VNP | If certain succinct hitting sets exist, the polynomial analogue of a natural proof can't separate VP from VNP; this one is conditional and less settled. |

What has worked despite the barriers: Ryan Williams's 2011 proof that NEXP is not in ACC⁰
turned fast algorithms into lower bounds and avoids all three big barriers. It is far weaker
than P ≠ NP but it shows the barriers can be stepped around.

## 2. Things that can't be the counterexample (routes to P = NP that are closed)

This is the counterexample-map spirit from the Fermat work: what can't work.

- **Any SAT solver whose runs are resolution proofs.** Modern CDCL solvers are equivalent to
  resolution (Pipatsrisawat–Darwiche 2011), and Haken (1985) proved resolution needs
  exponential size on pigeonhole formulas. So those solvers must take exponential time there.
- **Any polynomial-size linear program for TSP.** Fiorini, Massar, Pokutta, Tiwary and de Wolf
  (2012) proved the TSP polytope has exponential extension complexity. Rothvoß (2014) proved
  the same for perfect matching, even though matching is in P, so LPs aren't even the right
  measure of "easy." Many claimed proofs of P = NP were exactly these LPs.
- **Any pure dynamic program over min-plus.** Jerrum and Snir (1982) proved monotone arithmetic
  circuits for the permanent need exponential size, and for homogeneous problems like
  assignment that transfers to min-plus (tropical) circuits. Yet assignment is solved in
  polynomial time by the Hungarian method, which uses subtraction and branching. So tropical
  lower bounds, like monotone ones, are real but don't reach general algorithms.
- **A universal search.** Levin (1973) wrote down an algorithm that is polynomial for SAT
  search if P = NP at all. So if a counterexample exists, a version of it is already known; what
  nobody can do is prove its running time. The difficulty is proving, not finding.

## 3. What the quantum Fourier transform can and can't tell us

- **What it does:** the QFT finds hidden periods. That gives exponential speedups for factoring
  and discrete log (Shor), and in general for the hidden subgroup problem over abelian groups.
- **Why that doesn't reach NP-complete problems:** factoring isn't believed NP-complete. For
  black-box search, Bennett, Bernstein, Brassard and Vazirani (1997) proved quantum computers
  need about √N queries, so Grover's square-root speedup is the best possible; relative to a
  random oracle, NP is not inside BQP. Exponential speedups need hidden algebraic structure,
  and NP-complete problems in general don't offer it.
- **The nonabelian case fails too:** for the symmetric group (the graph-isomorphism case), the
  standard Fourier sampling approach provably can't succeed (Moore, Russell, Schulman 2005;
  Hallgren et al.).
- **And it wouldn't settle P vs NP anyway:** P vs NP is about classical computers. Even NP ⊆ BQP
  would not imply P = NP. Raz and Tal (2018) found an oracle where BQP is not inside the
  polynomial hierarchy, so quantum and classical classes sit sideways to each other.

## 4. Where the circle operation touches the question, judged on its own merits

### 4.1 The shifted ruler: no new power (proved)

With shifted add sadd(a, b) = a + b + 1 and the circle operation a ∘₁ b = a + b + ab, the map x ↦ 1 + x turns
sadd into + and ∘₁ into ×, gate for gate. Lean (`evalS_lift`, `eval_from_shift`, `size_mapC`):
every sadd/∘₁ circuit is a +/× circuit with the same wiring and size, and every +/× circuit is a
sadd/∘₁ circuit with the same wiring and size. So circuit size, VP vs VNP, and every algebraic
lower-bound question are identical in both worlds. **Not a lead.** This is the honest answer
to "does anything change?": nothing does, and that's now a theorem rather than a hunch.

The circle operation a + b + k·ab is the same story for every k ≠ 0 (it's ∘₁ seen through
x ↦ kx, `CircleOp.circK_scale`). At k = 0 it becomes plain addition. That's a contraction, not
a tropical limit: the dial doesn't turn + into min.

### 4.2 The dial over witnesses: NP, #P, ⊕P and ModₚP at different settings (proved)

Fold the dial over a list of yes/no witnesses and you get the q-integer of the number of yes
answers, in base k + 1 (`witness_fold`):

| Dial | k + 1 | The fold gives | The class |
|---|---|---|---|
| k = −1 | 0 | 1 if any witness, else 0 | NP ("is there one?"), `np_point` |
| k = 0 | 1 | the number of witnesses | #P ("how many?"), `sharpP_point` |
| k = −2 | −1 | 1 if odd, 0 if even | ⊕P ("odd or even?"), `parityP_point` |
| (k + 1)ᵖ = 1 | p-th root of unity | the count mod p | ModₚP, `modp_period` |
| k = 1 | 2 | 2ᵐ − 1, the whole count | #P again, `qnum_two_strictMono` |

The picture is correct and tidy: NP is the one point on the dial that throws information away
(`np_point_forgets`: one witness and five look the same at k = −1), and almost every other
point keeps the count. It fits your theme of operations that keep information, with a twist:
keeping information makes problems harder, not easier. Toda's theorem says counting (#P) is
at least as strong as the whole polynomial hierarchy.

**Judgment:** the classes are all standard (this is the GapP and counting-class picture in
new coordinates), and relativization already separates them in different ways in different
oracle worlds, so the dial gives no new leverage on P vs NP. It **is** worth one landscape
page for the project, because it puts NP, ⊕P, ModₚP and #P on one line as settings of the
same operation, and the known hardness jumps (permanent easy mod 2, hard mod 3; Valiant 1979)
become places where the landscape changes character.

### 4.3 Over the two-element field, the circle operation is the OR gate (proved)

Over GF(2), a + b + ab is a OR b, and x ↦ 1 + x is NOT, so the lift
1 + (a ∘₁ b) = (1 + a)(1 + b) is De Morgan's law (`gf2_circ_is_or`,
`gf2_lift_de_morgan`). Over the integers, the k = −1 version a + b − ab is the standard
arithmetized OR.

**Judgment:** this is where the circle operation lives most naturally in complexity theory, and it's the
most instructive negative finding here. The arithmetized OR is exactly what drives
IP = PSPACE and the Razborov–Smolensky lower bounds, and those two families are precisely
what the algebrization and natural-proofs barriers cover. So any P ≠ NP argument built mainly
from the circle operation's algebra on yes/no values sits inside known barrier territory. **Not a lead for
P vs NP**, but a useful "what can't work" entry.

### 4.4 Easy because it factors (proved)

The fold 1 ∘₁ 1 ∘₁ … ∘₁ 1 over n inputs has 2ⁿ − 1 monomials but takes n − 1 gates
(`foldCirc_ones`), because 1 + (fold) = (1 + x₁)⋯(1 + xₙ) factors. The hard counting problems
(the permanent, #SAT) are hard because their sums don't factor: the constraints tie the
variables together. The circle operation is the independent case. **Not a lead**; it explains why
circle-operation circuits are always easy.

### 4.5 The min-plus end

The dynamic-programming algebra (min, +) is a semiring (`minplus_distrib`), and the same
polynomial changes difficulty with the semiring: the permanent is perfect-matching existence
over Booleans (P), the assignment problem over min-plus (P, but tropical circuits need
exponential size), the determinant over GF(2) (P), and #P-complete over the integers. Tropical
circuit complexity is an active area where real lower bounds get proved (Jukna's work).
**Judgment:** a place where project-style algebra could prove real theorems, since lower bounds
there are provable. It won't reach P vs NP, because, like monotone bounds, they don't transfer
to general algorithms.

## 5. Leads worth following

- **For P vs NP itself:** none from the circle operation. The shifted ruler is provably the same
  computation, and the dial's natural home (arithmetized OR, counting) is covered by known
  barriers. Saying so plainly is the useful result.
- **For the project:** one, the dial landscape of counting classes (4.2), drawn as a page:
  NP at k = −1, ⊕P at k = −2, ModₚP at roots of unity, #P everywhere else, with the
  permanent's known jumps marked. It's a real organising picture of known math, not a new
  prediction.

## 6. The Lean theorems (FusedAlgebra/Topics/PvsNP.lean)

- `size_mapC`, `evalS_lift`, `eval_from_shift`: circuits transport across the shifted ruler at equal size
- `sAdd_lift`, `circ_sAdd_distrib`, `shifted_identities`: the shifted semiring (zero −1, one 0, −1 absorbs)
- `witness_fold`: the dial over witnesses is the q-integer of the witness count
- `np_point`, `sharpP_point`, `parityP_point`: exists, count and parity at k = −1, 0, −2
- `modp_period`: at a p-th root of unity only the count mod p matters
- `np_point_forgets`, `qnum_two_nonneg`, `qnum_two_strictMono`: k = −1 forgets, k = 1 keeps the count
- `gf2_circ_is_or`, `gf2_shift_is_not`, `gf2_lift_de_morgan`: over GF(2) circle operation is OR, the lift is De Morgan
- `foldCirc_ones`: n gates, 2ⁿ − 1 monomials
- `minplus_distrib`, `minplus_assoc`: min-plus semiring laws

## Sources

- Baker, Gill, Solovay, "Relativizations of the P =? NP question," SIAM J. Comput. 1975.
- Razborov, Rudich, "Natural proofs," JCSS 1997 (STOC 1994).
- Aaronson, Wigderson, "Algebrization: a new barrier in complexity theory," ACM TOCT 2009 (STOC 2008).
- Razborov, "Lower bounds on the monotone complexity of some Boolean functions," 1985; Tardos, Combinatorica 1988.
- Bürgisser, Ikenmeyer, Panova, "No occurrence obstructions in geometric complexity theory," J. Amer. Math. Soc. 32 (2019) 163–193.
- Cook, Reckhow, J. Symbolic Logic 1979; Haken, "The intractability of resolution," 1985; Razborov, "Unprovability of lower bounds on circuit size in certain fragments of bounded arithmetic," 1995.
- Grochow, Kumar, Saks, Saraf, "Towards an algebraic natural proofs barrier via polynomial identity testing," 2017.
- Williams, "Non-uniform ACC circuit lower bounds," 2011.
- Pipatsrisawat, Darwiche, "On the power of clause-learning SAT solvers as resolution engines," AIJ 2011.
- Fiorini, Massar, Pokutta, Tiwary, de Wolf, STOC 2012; Rothvoß, "The matching polytope has exponential extension complexity," STOC 2014 (2023 Gödel Prize).
- Jerrum, Snir, "Some exact complexity results for straight-line computations over semirings," J. ACM 29 (1982); Jukna, "Lower bounds for tropical circuits and dynamic programs," ToCS 2014.
- Levin, "Universal sequential search problems," 1973.
- Bennett, Bernstein, Brassard, Vazirani, "Strengths and weaknesses of quantum computing," SIAM J. Comput. 1997; Moore, Russell, Schulman, FOCS 2005; Raz, Tal, "Oracle separation of BQP and PH," STOC 2019.
- Valiant, "The complexity of computing the permanent," 1979; Toda, "PP is as hard as the polynomial-time hierarchy," 1991.
