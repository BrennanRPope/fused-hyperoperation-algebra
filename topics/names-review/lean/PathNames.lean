/-
  Is the circle operation a "path" construction? (1.0 checklist, naming check)
  ===========================================================================

  Part of Fused Hyperoperation Algebra. Companion: math-theory/names-review/
  NAMING-COLLISIONS.md ("path theory" section). Core Lean 4.37 only, no `sorry`.

  Over a NON-commutative ring (for example, the free algebra on the steps of a
  chain, or matrices), the circle operation keeps the order of a chain, and the
  composite of a chain is the sum of the products of all its ordered sub-chains:
      a ∘ b ∘ c ∘ d = Σ over nonempty ordered sub-chains S of the product of S.
  In the free (tensor) algebra this is the discrete, iterated-sums form of Chen's
  path signature, and the dial k weights a sub-chain of length m by k^(m−1),
  which is the graded dilation used in rough-path theory. Over a commutative ring
  (numbers) the order is forgotten, so the scalar composite is a lossy shadow of
  the free one.
-/

namespace PathNames

section
variable {R : Type} [Lean.Grind.Ring R]

/-- The circle operation with dial k, over a possibly non-commutative ring
    (k on the left; k = 1 is the classical circle operation). -/
def circ (k a b : R) : R := a + b + k * a * b

/-- k = 1: the circle operation is associative in any ring. -/
theorem circ_one_assoc (a b c : R) : circ 1 (circ 1 a b) c = circ 1 a (circ 1 b c) := by
  unfold circ; grind

/-- **The whole landscape of a chain.** Four steps composed: every ordered
    sub-chain (15 of them) appears exactly once, in its original order. -/
theorem chain4_landscape (a b c d : R) :
    circ 1 (circ 1 (circ 1 a b) c) d =
      (a + b + c + d) + (a*b + a*c + a*d + b*c + b*d + c*d)
      + (a*b*c + a*b*d + a*c*d + b*c*d) + a*b*c*d := by
  unfold circ; grind

/-- **The dial as dilation.** Three steps with dial k: a sub-chain of length m is
    weighted by k^(m−1). k = 0 keeps only the single steps (the plain sum). -/
theorem chain3_dial (k a b c : R) :
    circ k (circ k a b) c =
      (a + b + c) + k * (a*b + a*c + b*c) + k * k * (a*b*c) := by
  unfold circ; grind

/-- **Order is kept, scaled by k.** Swapping two steps changes the composite by
    k times their commutator, so nothing about order is lost unless k = 0 or the
    steps commute. -/
theorem circ_swap (k a b : R) : circ k a b - circ k b a = k * (a * b - b * a) := by
  unfold circ; grind

end

/-- Over numbers the order is lost: different chains, same composite.
    (1 ∘ 2 = 5 = 2 ∘ 1, and the chains 1,2 and 5 alone also give 5.) -/
theorem scalar_shadow :
    circ (1 : Int) 1 2 = 5 ∧ circ (1 : Int) 2 1 = 5 := by decide

end PathNames
