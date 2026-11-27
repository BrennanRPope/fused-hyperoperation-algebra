/-
  Fused Hyperoperation Algebra: entropy, knots and zeros (DRAFT, Nov 27 2026).
  Core Lean 4.37 only (no mathlib), no `sorry`.
  Library path: FusedAlgebra/Topics/EntropyKnots.lean. Write-up: topics/entropy-knots/README.md.

  1. The Cartier–Foata polynomial of the trace monoid whose letters commute unless neighbours
     (the Temperley–Lieb generators e₁, …, e_m, dependence graph a path) is the independence
     polynomial of the path at fugacity −t. Homogenised, it is the Jones–Wenzl / Chebyshev
     recursion: dᵐ⁺¹ · I_m(−u/d²) = H_m(d, u), with H_m = d·H_{m−1} − u·H_{m−2}.
     (Heaps of dimers and Chebyshev: Viennot 1986; Jones–Wenzl: Wenzl 1987, Kauffman–Lins 1994.)
  2. Its zeros in d at u = 1 are the Jones-index walls d = 2cos(jπ/(m+2)): small cases factor.
  3. The trace monoid of three affine maps, two sharing a centre, is the period-3 window of the
     logistic map squared (a Graeffe step): μ(t²) = (1 − t − t²)(1 + t − t²).
  4. Positive 3-braids: Deligne's Möbius polynomial 1 − 2t + t³ = (1 − t)(1 − t − t²).
  5. Maps acting on disjoint coordinates commute (why far-apart strand maps form a trace monoid).
-/

namespace EntropyKnots

set_option linter.unusedSectionVars false

/-- Independence polynomial of the path on m vertices, as a function of the fugacity x. -/
def indepPath : Nat → Int → Int
  | 0, _ => 1
  | 1, x => 1 + x
  | m + 2, x => indepPath (m + 1) x + x * indepPath m x

/-- The Jones–Wenzl recursion, homogenised: H₀ = d, H₁ = d² − u, H_{m+2} = d·H_{m+1} − u·H_m. -/
def jw : Nat → Int → Int → Int
  | 0, d, _ => d
  | 1, d, u => d * d - u
  | m + 2, d, u => d * jw (m + 1) d u - u * jw m d u

/-- Chebyshev polynomials of the second kind in the variable d = 2·cos θ: U₀ = 1, U₁ = d. -/
def cheb : Nat → Int → Int
  | 0, _ => 1
  | 1, d => d
  | n + 2, d => d * cheb (n + 1) d - cheb n d

theorem indepPath_two (x : Int) : indepPath 2 x = 1 + 2 * x := by
  simp only [indepPath]; grind

/-- Three letters in a row (two of them commute): independence polynomial 1 + 3x + x². -/
theorem indepPath_three (x : Int) : indepPath 3 x = 1 + 3 * x + x * x := by
  simp only [indepPath]; grind

theorem indepPath_four (x : Int) : indepPath 4 x = 1 + 4 * x + 3 * (x * x) := by
  simp only [indepPath]; grind

private theorem pair (d x : Int) (m : Nat) :
    jw m d (-(x * d * d)) = d ^ (m + 1) * indepPath m x ∧
    jw (m + 1) d (-(x * d * d)) = d ^ (m + 2) * indepPath (m + 1) x := by
  induction m with
  | zero =>
    constructor
    · simp only [jw, indepPath]; grind
    · simp only [jw, indepPath]; grind
  | succ m ih =>
    obtain ⟨h0, h1⟩ := ih
    refine ⟨h1, ?_⟩
    show d * jw (m + 1) d (-(x * d * d)) - (-(x * d * d)) * jw m d (-(x * d * d)) =
      d ^ (m + 1 + 2) * (indepPath (m + 1) x + x * indepPath m x)
    rw [h0, h1]
    simp only [Int.pow_succ]
    grind

/-- **Trace monoid = Jones–Wenzl.** With u = −x·d², the Jones–Wenzl recursion is dᵐ⁺¹ times the
    path's independence polynomial: so μ(t) = I_m(−t) vanishes at t = 1/d² exactly where
    H_m(d, 1) does. -/
theorem path_jones_wenzl (m : Nat) (d x : Int) :
    jw m d (-(x * d * d)) = d ^ (m + 1) * indepPath m x := (pair d x m).1

private theorem cheb_pair (d : Int) (m : Nat) :
    jw m d 1 = cheb (m + 1) d ∧ jw (m + 1) d 1 = cheb (m + 2) d := by
  induction m with
  | zero => constructor <;> simp only [jw, cheb] <;> grind
  | succ m ih =>
    obtain ⟨h0, h1⟩ := ih
    refine ⟨h1, ?_⟩
    show d * jw (m + 1) d 1 - 1 * jw m d 1 = d * cheb (m + 2) d - cheb (m + 1) d
    rw [h0, h1]; grind

/-- At u = 1 the recursion is Chebyshev U_{m+1}(d/2), the Jones–Wenzl trace Δ_{m+1}(d). -/
theorem jw_chebyshev (m : Nat) (d : Int) : jw m d 1 = cheb (m + 1) d := (cheb_pair d m).1

/-- Path P₂: zero at d² = 2 (Beraha B₄ = 2). -/
theorem jw_two (d : Int) : jw 2 d 1 = d * (d * d - 2) := by
  simp only [jw]; grind

/-- Path P₃: zeros at d² − d − 1 = 0, d = φ (Beraha B₅ = φ²), the growth of "three maps, two
    sharing a centre". -/
theorem jw_three (d : Int) : jw 3 d 1 = (d * d - d - 1) * (d * d + d - 1) := by
  simp only [jw]; grind

/-- Path P₄: zeros at d = 1 and d² = 3 (Beraha B₆ = 3). -/
theorem jw_four (d : Int) : jw 4 d 1 = d * (d * d - 1) * (d * d - 3) := by
  simp only [jw]; grind

/-- **Graeffe step.** The Cartier–Foata polynomial of three affine maps with two sharing a centre,
    taken at t², is the period-3 kneading polynomial times its mirror: that trace monoid's growth
    is the square of the period-3 window's (φ² versus φ). -/
theorem graeffe_period3 (t : Int) :
    (1 - t - t * t) * (1 + t - t * t) = indepPath 3 (-(t * t)) := by
  rw [indepPath_three]; grind

/-- The period-3 pattern's Markov matrix [[0,1],[1,1]]: det(I − tA) = 1 − t − t². -/
theorem markov3_det (t : Int) : (1 - 0 * t) * (1 - 1 * t) - (0 - 1 * t) * (0 - 1 * t) = 1 - t - t * t := by
  grind

/-- Positive 3-braids: Deligne's Möbius polynomial factors through the golden polynomial, so they
    grow like φ, the same rate as the period-3 window. -/
theorem deligne_b3 (t : Int) : 1 - 2 * t + t * t * t = (1 - t) * (1 - t - t * t) := by
  grind

/-- Maps on disjoint coordinates commute: f on the first coordinate, g on the last. -/
theorem disjoint_commute {α β γ : Type} (f : α → α) (g : γ → γ) (p : α × β × γ) :
    (fun q : α × β × γ => (f q.1, q.2.1, q.2.2)) ((fun q : α × β × γ => (q.1, q.2.1, g q.2.2)) p) =
    (fun q : α × β × γ => (q.1, q.2.1, g q.2.2)) ((fun q : α × β × γ => (f q.1, q.2.1, q.2.2)) p) := rfl

end EntropyKnots
