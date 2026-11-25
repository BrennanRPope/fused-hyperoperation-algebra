/-
  Entropy and quantum noise on the circle-operation cube (called mult-add cube before 1.0)
  ==============================================

  Part of Fused Hyperoperation Algebra. Companion page: index.html.
  Numbers: checks.py (45 checks). Pictures: figures.py.

  STATUS: core Lean 4.37.0-rc1, no mathlib, no `sorry`. Imports the shared
  library's Base, the quantum-noise module (Pauli channels, the 3-bit vote)
  and the information-theory module (Tsallis, Gini–Simpson, outer sums).

  LIBRARY NOTE: copied from entropy-noise/lean/EntropyNoise.lean (Nov 25).
  Six theorems that state exactly what a library theorem already states were
  removed rather than kept as twins; a comment stands in each place naming the
  library theorem to use. They are tsallis_two_is_union, weights_compose and
  linear_entropy_parallel (all InfoFused.gini_simpson), tsallis_weight
  (InfoFused.tsallis at k = 1 − q), carnot_weight (ThermoFused.carnot_series)
  and decimation (QuantumNoise.bit_flip). Nothing else changed.

  What this file proves, by part:
   1. The entropy weight.  Certainty absorbs under a ◯[k] b exactly when k = −1,
      so Tsallis' line (k = 1 − q) meets the "or" line only at q = 2.  The
      weight k = 1 − q composes by the union rule along a chain of escorts,
      the same law as Carnot efficiencies, k = 1 − T₀/T.  Escorts of escorts
      multiply their orders.
   2. The Pauli tetrahedron.  The collision sum Σp² is (1 + |λ|²)/4, so the
      collision (Rényi-2) entropy has spheres for level sets; Σp² ≤ ½ exactly
      inside the midsphere |λ| = 1.  The octahedron |λ|₁ ≤ 1 sits inside that
      ball, which sits inside the cube; all three kiss at the edge midpoints.
      Composing with any channel keeps a channel inside each of the three
      (the wall's absorption, fattened).  Series composition and parallel
      (tensor) composition of collision sums.
   3. Stacking.  The n-fold fault pattern's linear entropy is the k = −1 climb.
      The net channel is always purer than its fault pattern:
      ((1 + x)/2)ⁿ ≤ (1 + xⁿ)/2 for 0 ≤ x ≤ 1; at n = 2 the gap is a square.
   4. The 3-bit code.  Syndromes are the 3-cube's 4 body diagonals; their
      chances in the polarization; an alarm leaves exactly one fresh coin;
      the entropy chain rule on one diagonal, for any log-like L.
   5. Nishimori.  On the Nishimori line tanh K = 1 − 2p = λ, and the
      Kramers–Wannier lens sends λ to the error odds p/(1 − p).

  NOT checked here (core Lean has no real numbers or logarithms): Shannon and
  von Neumann entropy as numbers, the hashing bound 1 − H(p), its depolarizing
  point 0.18929, the plane bound p_max > 0.8107, King's capacity formula,
  that the octahedron is exactly the entanglement-breaking set, Baez's
  Rényi = free energy identity, and h(p) = ½ at p = 0.110.  All are in checks.py.
-/
import FusedAlgebra.Base
import FusedAlgebra.Topics.QuantumNoise
import FusedAlgebra.Topics.InformationTheory

namespace EntropyNoise

open CircleOp

/-! ## Small order facts on Rat (core has no `nlinarith`) -/

theorem sq_nn (s : Rat) : 0 ≤ s * s := by
  rcases (@Rat.le_total 0 s) with h | h
  · exact Rat.mul_nonneg h h
  · have h1 : 0 ≤ -s := by grind
    have := Rat.mul_nonneg h1 h1
    grind

theorem pow_nn {x : Rat} (hx : 0 ≤ x) : ∀ n : Nat, 0 ≤ x ^ n
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by rw [Rat.pow_succ]; exact Rat.mul_nonneg (pow_nn hx n) hx

theorem pow_le_one' {x : Rat} (h0 : 0 ≤ x) (h1 : x ≤ 1) : ∀ n : Nat, x ^ n ≤ 1
  | 0 => by rw [Rat.pow_zero]; grind
  | n + 1 => by
    rw [Rat.pow_succ]
    have ih := pow_le_one' h0 h1 n
    have := Rat.mul_le_mul_of_nonneg_left h1 (pow_nn h0 n)
    grind

theorem mul_ne_zero' {a b : Rat} (ha : a ≠ 0) (hb : b ≠ 0) : a * b ≠ 0 := by
  intro h; rcases Rat.mul_eq_zero.mp h with h' | h' <;> contradiction

/-- |a·b| ≤ u whenever |a| ≤ u and |b| ≤ 1, written without `abs`. -/
theorem mul_bound (a b u : Rat) (h1 : a ≤ u) (h2 : -a ≤ u) (hb1 : b ≤ 1) (hb2 : -1 ≤ b) :
    a * b ≤ u ∧ -(a * b) ≤ u := by
  have hu : 0 ≤ u := by grind
  rcases (@Rat.le_total 0 b) with hb | hb
  · have e1 := Rat.mul_le_mul_of_nonneg_left h1 hb
    have e2 := Rat.mul_le_mul_of_nonneg_left h2 hb
    have e3 := Rat.mul_le_mul_of_nonneg_left hb1 hu
    refine ⟨?_, ?_⟩ <;> grind
  · have hb' : 0 ≤ -b := by grind
    have e1 := Rat.mul_le_mul_of_nonneg_left h1 hb'
    have e2 := Rat.mul_le_mul_of_nonneg_left h2 hb'
    have e3 : -b ≤ 1 := by grind
    have e4 := Rat.mul_le_mul_of_nonneg_left e3 hu
    refine ⟨?_, ?_⟩ <;> grind

/-! ## 1. The entropy weight -/

/-- Certainty swallows everything, a ◯[k] 1 = 1 for every a, exactly when k = −1.
    The statistics page's "or" laws are the line k + q_L = −1 of lens laws; the
    Tsallis laws have no denominator (q_L = 0), so they meet it only here. -/
theorem certainty_absorbs_iff (k : Rat) : (∀ a : Rat, circK k a 1 = 1) ↔ k = -1 := by
  constructor
  · intro h; have := h 1; unfold circK at this; grind
  · intro hk a; subst hk; unfold circK; grind

/- Tsallis' index q = 2 is that point: Gini–Simpson (linear entropy) follows the
   union rule, because its lift is the collision sum. Library: InfoFused.gini_simpson
   (was tsallis_two_is_union).

   Tsallis pseudo-additivity with the weight written as 1 − q: InfoFused.tsallis
   with k = 1 − q (was tsallis_weight).

   Orders multiply along a chain, so weights k = 1 − q compose by the union rule:
   InfoFused.gini_simpson again (was weights_compose).

   k(T, T₀) = 1 − T₀/T is Carnot's efficiency, and it cascades by the union rule:
   ThermoFused.carnot_series (was carnot_weight). -/

/-- Raising every entry to a then to b is raising to a·b: escorts of escorts. -/
theorem escort_powers (l : List Rat) (a b : Nat) :
    (l.map (fun x => x ^ a)).map (fun y => y ^ b) = l.map (fun x => x ^ (a * b)) := by
  rw [List.map_map]
  apply List.map_congr_left
  intro x _
  simp only [Function.comp]
  rw [← pw_eq_pow, ← pw_eq_pow, ← pw_eq_pow, pw_mul]

theorem lsum_div (c : Rat) (l : List Rat) :
    InfoFused.lsum (l.map (fun y => y / c)) = InfoFused.lsum l / c := by
  induction l with
  | nil => simp only [List.map_nil, InfoFused.lsum]; grind
  | cons y ys ih =>
    simp only [List.map_cons, InfoFused.lsum]; rw [ih]
    rw [Rat.div_def, Rat.div_def, Rat.div_def]; grind

/-- The escort of order n of a list of weights. -/
def escort (n : Nat) (l : List Rat) : List Rat :=
  l.map (fun x => x ^ n / InfoFused.lsum (l.map (fun y => y ^ n)))

theorem div_pow' (x c : Rat) : ∀ n : Nat, (x / c) ^ n = x ^ n / c ^ n
  | 0 => by simp only [Rat.pow_zero]; rw [Rat.div_def]; grind
  | n + 1 => by
    rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_succ, div_pow' x c n]
    rw [Rat.div_def, Rat.div_def, Rat.div_def, Rat.inv_mul_rev]; grind

/-- Escort of an escort: P_b(P_a(p)) = P_{ab}(p), when the normalizers are nonzero. -/
theorem pow_ne_zero' {Z : Rat} (h : Z ≠ 0) : ∀ n : Nat, Z ^ n ≠ 0
  | 0 => by rw [Rat.pow_zero]; decide
  | n + 1 => by rw [Rat.pow_succ]; exact mul_ne_zero' (pow_ne_zero' h n) h

theorem escort_escort (l : List Rat) (a b : Nat)
    (ha : InfoFused.lsum (l.map (fun y => y ^ a)) ≠ 0) :
    escort b (escort a l) = escort (a * b) l := by
  unfold escort
  generalize hZ : InfoFused.lsum (l.map (fun y => y ^ a)) = Z at ha ⊢
  have hpow : ∀ x : Rat, (x ^ a / Z) ^ b = x ^ (a * b) / Z ^ b := by
    intro x; rw [div_pow', ← pw_eq_pow, ← pw_eq_pow, ← pw_mul, pw_eq_pow]
  have hsum : InfoFused.lsum ((l.map (fun x => x ^ a / Z)).map (fun y => y ^ b))
      = InfoFused.lsum (l.map (fun x => x ^ (a * b))) / Z ^ b := by
    rw [List.map_map]
    have : (fun y => y ^ b) ∘ (fun x => x ^ a / Z) = (fun y => y / Z ^ b) ∘ (fun x => x ^ (a * b)) := by
      funext x; simp only [Function.comp]; exact hpow x
    rw [this, ← List.map_map, lsum_div]
  rw [hsum, List.map_map]
  apply List.map_congr_left
  intro x _
  simp only [Function.comp]
  rw [hpow]
  have hZb : Z ^ b ≠ 0 := pow_ne_zero' ha b
  generalize InfoFused.lsum (List.map (fun x => x ^ (a * b)) l) = S
  generalize x ^ (a * b) = X
  generalize Z ^ b = W at hZb
  rw [Rat.div_def, Rat.div_def, Rat.div_def, Rat.div_def, Rat.inv_mul_rev, Rat.inv_inv]
  calc X * W⁻¹ * (W * S⁻¹) = X * S⁻¹ * (W⁻¹ * W) := by grind
    _ = X * S⁻¹ := by rw [Rat.inv_mul_cancel _ hZb]; grind

/-! ## 2. The Pauli tetrahedron: three norms of λ -/

open QuantumNoise

/-- The squared length of a channel's point λ in the cube. -/
def len2 (p : Pauli) : Rat := p.lx * p.lx + p.ly * p.ly + p.lz * p.lz
/-- The collision sum Σp², whose −log is the collision (Rényi-2) entropy. -/
def coll (p : Pauli) : Rat := p.i * p.i + p.x * p.x + p.y * p.y + p.z * p.z

/-- Parseval on the tetrahedron: Σp² = (1 + |λ|²)/4. -/
theorem collision_parseval (p : Pauli) (h : p.i + p.x + p.y + p.z = 1) :
    coll p = (1 + len2 p) / 4 := by
  simp only [coll, len2, Pauli.lx, Pauli.ly, Pauli.lz]; grind

/-- The Choi state's linear entropy 1 − Σp² is (3 − |λ|²)/4: zero at the corners,
    largest (3/4) at the centre, constant on spheres about the centre. -/
theorem linear_entropy_choi (p : Pauli) (h : p.i + p.x + p.y + p.z = 1) :
    1 - coll p = (3 - len2 p) / 4 := by
  rw [collision_parseval p h]; grind

/-- Σp² ≤ ½ (collision entropy at least one bit) exactly inside the midsphere. -/
theorem midsphere_iff (p : Pauli) (h : p.i + p.x + p.y + p.z = 1) :
    coll p ≤ 1 / 2 ↔ len2 p ≤ 1 := by
  rw [collision_parseval p h]; constructor <;> intro h' <;> grind

def corner : Pauli := ⟨1, 0, 0, 0⟩
def centre : Pauli := ⟨1/4, 1/4, 1/4, 1/4⟩
def edgeMid : Pauli := ⟨1/2, 1/2, 0, 0⟩

theorem corner_point : corner.lx = 1 ∧ corner.ly = 1 ∧ corner.lz = 1 ∧ coll corner = 1 := by
  simp only [corner, Pauli.lx, Pauli.ly, Pauli.lz, coll]; decide +kernel
theorem centre_point : centre.lx = 0 ∧ centre.ly = 0 ∧ centre.lz = 0 ∧ coll centre = 1/4 := by
  simp only [centre, Pauli.lx, Pauli.ly, Pauli.lz, coll]; decide +kernel
/-- The edge midpoint (X with chance ½) is (1, 0, 0): on the cube's face, the
    sphere and the octahedron at once, with collision sum exactly ½. -/
theorem edge_midpoint : edgeMid.lx = 1 ∧ edgeMid.ly = 0 ∧ edgeMid.lz = 0 ∧
    coll edgeMid = 1/2 ∧ len2 edgeMid = 1 := by
  simp only [edgeMid, Pauli.lx, Pauli.ly, Pauli.lz, coll, len2]; decide +kernel

/-- Octahedron inside the ball: |x| + |y| + |z| ≤ 1 gives x² + y² + z² ≤ 1. -/
theorem octahedron_in_ball (x y z u v w : Rat)
    (hx : x ≤ u) (hx' : -x ≤ u) (hy : y ≤ v) (hy' : -y ≤ v) (hz : z ≤ w) (hz' : -z ≤ w)
    (h : u + v + w ≤ 1) : x * x + y * y + z * z ≤ 1 := by
  have hu : u ≤ 1 := by grind
  have hv : v ≤ 1 := by grind
  have hw : w ≤ 1 := by grind
  have a := (mul_bound x x u hx hx' (by grind) (by grind)).1
  have b := (mul_bound y y v hy hy' (by grind) (by grind)).1
  have c := (mul_bound z z w hz hz' (by grind) (by grind)).1
  grind

/-- Ball inside the cube. -/
theorem ball_in_cube (x y z : Rat) (h : x * x + y * y + z * z ≤ 1) : x * x ≤ 1 := by
  have := sq_nn y; have := sq_nn z; grind

/-- Composition multiplies coordinates (QuantumNoise.comp_lx …), and a coordinate
    of a real channel lies in [−1, 1].  So composing never leaves the ball … -/
theorem compose_stays_in_ball (ax ay az bx b_y bz : Rat)
    (hx : bx * bx ≤ 1) (hy : b_y * b_y ≤ 1) (hz : bz * bz ≤ 1) :
    (ax * bx) * (ax * bx) + (ay * b_y) * (ay * b_y) + (az * bz) * (az * bz)
      ≤ ax * ax + ay * ay + az * az := by
  have e1 := Rat.mul_le_mul_of_nonneg_left hx (sq_nn ax)
  have e2 := Rat.mul_le_mul_of_nonneg_left hy (sq_nn ay)
  have e3 := Rat.mul_le_mul_of_nonneg_left hz (sq_nn az)
  grind

/-- … nor the octahedron … -/
theorem compose_stays_in_octahedron (ax ay az bx b_y bz u v w : Rat)
    (hx : ax ≤ u) (hx' : -ax ≤ u) (hy : ay ≤ v) (hy' : -ay ≤ v) (hz : az ≤ w) (hz' : -az ≤ w)
    (bx1 : bx ≤ 1) (bx2 : -1 ≤ bx) (by1 : b_y ≤ 1) (by2 : -1 ≤ b_y) (bz1 : bz ≤ 1) (bz2 : -1 ≤ bz) :
    (ax * bx ≤ u ∧ -(ax * bx) ≤ u) ∧ (ay * b_y ≤ v ∧ -(ay * b_y) ≤ v) ∧
      (az * bz ≤ w ∧ -(az * bz) ≤ w) :=
  ⟨mul_bound ax bx u hx hx' bx1 bx2, mul_bound ay b_y v hy hy' by1 by2, mul_bound az bz w hz hz' bz1 bz2⟩

/-- … nor any cube max|λᵢ| ≤ c about the centre (classical capacity cannot rise). -/
theorem compose_stays_in_cube (a b c : Rat) (h1 : a ≤ c) (h2 : -a ≤ c) (hb1 : b ≤ 1) (hb2 : -1 ≤ b) :
    a * b ≤ c ∧ -(a * b) ≤ c := mul_bound a b c h1 h2 hb1 hb2

/-- Series: the collision sum of a composed channel, through Parseval. -/
theorem collision_series (p q : Pauli) (hp : p.i + p.x + p.y + p.z = 1) (hq : q.i + q.x + q.y + q.z = 1) :
    coll (p.comp q) = (1 + (p.lx * q.lx) * (p.lx * q.lx) + (p.ly * q.ly) * (p.ly * q.ly)
      + (p.lz * q.lz) * (p.lz * q.lz)) / 4 := by
  have ht : (p.comp q).i + (p.comp q).x + (p.comp q).y + (p.comp q).z = 1 := by
    rw [comp_total, hp, hq]; grind
  rw [collision_parseval _ ht]; unfold len2; rw [comp_lx, comp_ly, comp_lz]; grind

theorem map_sq_outer (xs ys : List Rat) :
    (InfoFused.outer xs ys).map (fun t => t * t)
      = InfoFused.outer (xs.map (fun t => t * t)) (ys.map (fun t => t * t)) := by
  induction xs with
  | nil => simp only [InfoFused.outer, List.map_nil]
  | cons x xs ih =>
    simp only [InfoFused.outer, List.map_cons, List.map_append, ih, List.map_map]
    congr 1
    apply List.map_congr_left; intro y _; simp only [Function.comp]; grind

/-- Parallel (two qubits, independent faults): the 16 joint fault chances'
    collision sum is the product, so linear entropy composes by the union rule. -/
theorem collision_parallel (p q : Pauli) :
    InfoFused.lsum ((InfoFused.outer [p.i, p.x, p.y, p.z] [q.i, q.x, q.y, q.z]).map (fun t => t * t))
      = coll p * coll q := by
  rw [map_sq_outer, InfoFused.outer_sum]
  simp only [List.map, InfoFused.lsum, coll]; grind

/- Linear entropy of a parallel pair follows the union rule:
   InfoFused.gini_simpson (coll p) (coll q) (was linear_entropy_parallel). -/

/-- The depolarizing diagonal: total fault chance s spread evenly over X, Y, Z. -/
def depol (s : Rat) : Pauli := ⟨1 - s, s / 3, s / 3, s / 3⟩

theorem depol_point (s : Rat) : (depol s).lx = 1 - 4/3 * s ∧ (depol s).ly = 1 - 4/3 * s ∧
    (depol s).lz = 1 - 4/3 * s := by
  simp only [depol, Pauli.lx, Pauli.ly, Pauli.lz]; refine ⟨?_, ?_, ?_⟩ <;> grind

/-- Along the diagonal (λ, λ, λ) the octahedron's face is λ = 1/3 and the
    midsphere is λ² = 1/3, so the octahedron is strictly inside there. -/
theorem diagonal_crossings : (1/3 : Rat) + 1/3 + 1/3 = 1 ∧ 3 * ((1/3 : Rat) * (1/3)) < 1 := by
  decide +kernel

/-! ## 3. Stacking: what cancellation hides -/

/-- One bit-flip channel's collision sum in the polarization λ = 1 − 2p. -/
theorem bit_collision (p : Rat) : (1 - p) * (1 - p) + p * p = (1 + (1 - 2*p) * (1 - 2*p)) / 2 := by
  grind

/-- The n-fold fault pattern (a corner of the n-cube) has collision sum cⁿ, so its
    linear entropy climbs at weight −1: 1 − cⁿ = (1 − c)^{◯[−1] n}. -/
theorem pattern_linear_entropy (c : Rat) (n : Nat) : 1 - c ^ n = wpow (-1) (1 - c) n := by
  rw [wpow_neg_one]; grind

/-- The net channel after n stacked flips is purer than its fault pattern:
    ((1 + x)/2)ⁿ ≤ (1 + xⁿ)/2 for x = λ² in [0, 1]. -/
theorem cancellation_hides (x : Rat) (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    ∀ n : Nat, ((1 + x) / 2) ^ n ≤ (1 + x ^ n) / 2
  | 0 => by simp only [Rat.pow_zero]; grind
  | n + 1 => by
    have ih := cancellation_hides x h0 h1 n
    rw [Rat.pow_succ, Rat.pow_succ]
    have hc : 0 ≤ (1 + x) / 2 := by grind
    have e := Rat.mul_le_mul_of_nonneg_right ih hc
    have hxn := pow_nn h0 n
    have hxn1 := pow_le_one' h0 h1 n
    have key : 0 ≤ (1 - x) * (1 - x ^ n) := Rat.mul_nonneg (by grind) (by grind)
    grind

/-- At n = 2 the purity gap is a perfect square. -/
theorem gap_two (l : Rat) : (1 + l^4) / 2 - ((1 + l^2) / 2)^2 = (1 - l^2)^2 / 4 := by grind

/-! ## 4. The 3-bit code: syndromes are body diagonals -/

/-- The two parity checks of the 3-bit code on a corner (a, b, c) of the cube. -/
def synd (a b c : Bool) : Bool × Bool := (a != b, b != c)

/-- Opposite corners give the same syndrome … -/
theorem synd_antipodal : ∀ a b c : Bool, synd (!a) (!b) (!c) = synd a b c := by decide
/-- … and corners with the same syndrome are equal or opposite: each syndrome is one
    body diagonal of the 3-cube. -/
theorem synd_diagonal : ∀ a b c d e f : Bool, synd a b c = synd d e f →
    (a = d ∧ b = e ∧ c = f) ∨ (a = !d ∧ b = !e ∧ c = !f) := by decide

/-- No alarm: both ends of the main diagonal. Alarm: each of the three other diagonals. -/
theorem syndromes_total (p : Rat) :
    ((1 - p)^3 + p^3) + 3 * (p * (1 - p)^2 + p^2 * (1 - p)) = 1 := by grind

theorem no_alarm_polar (p : Rat) : (1 - p)^3 + p^3 = (1 + 3 * (1 - 2*p)^2) / 4 := by grind
theorem alarm_polar (p : Rat) : p * (1 - p)^2 + p^2 * (1 - p) = (1 - (1 - 2*p)^2) / 4 := by grind

/-- The vote picks the lighter end of each diagonal; the failures add up to the
    top two layers of the cube (QuantumNoise.maj3). -/
theorem decoder_failure (p : Rat) : p^3 + 3 * (p^2 * (1 - p)) = maj3 p := by unfold maj3; grind

/-- After an alarm, the leftover doubt is exactly one fresh coin of chance p. -/
theorem alarm_one_coin (p q : Rat) (hpq : p + q = 1) (hp : p ≠ 0) (hq : q ≠ 0) :
    p * q * q / (p * q * q + p * p * q) = q := by
  have e : p * q * q + p * p * q = p * q := by
    have : p * q * q + p * p * q = p * q * (q + p) := by grind
    rw [this]; grind
  rw [e]
  have hpq0 : p * q ≠ 0 := mul_ne_zero' hp hq
  rw [Rat.div_def]
  calc p * q * q * (p * q)⁻¹ = q * ((p * q) * (p * q)⁻¹) := by grind
    _ = q := by rw [Rat.mul_inv_cancel _ hpq0]; grind

/-- The chain rule on one diagonal, for any L with L(xy) = L(x) + L(y) (log is one):
    a·L(a) + b·L(b) = s·L(s) + [a·L(a/s) + b·L(b/s)] with s = a + b.
    Entropy of the faults = entropy of the syndrome + entropy left on the diagonal. -/
theorem chain_rule_pair (L : Rat → Rat) (hL : ∀ x y, L (x * y) = L x + L y)
    (a b : Rat) (hs : a + b ≠ 0) :
    a * L a + b * L b = (a + b) * L (a + b) + (a * L (a / (a + b)) + b * L (b / (a + b))) := by
  have ea : a = (a + b) * (a / (a + b)) := by
    rw [Rat.div_def]
    calc a = a * ((a + b) * (a + b)⁻¹) := by rw [Rat.mul_inv_cancel _ hs]; grind
      _ = (a + b) * (a * (a + b)⁻¹) := by grind
  have eb : b = (a + b) * (b / (a + b)) := by
    rw [Rat.div_def]
    calc b = b * ((a + b) * (a + b)⁻¹) := by rw [Rat.mul_inv_cancel _ hs]; grind
      _ = (a + b) * (b * (a + b)⁻¹) := by grind
  have la : L a = L (a + b) + L (a / (a + b)) := by rw [← hL]; exact congrArg L ea
  have lb : L b = L (a + b) + L (b / (a + b)) := by rw [← hL]; exact congrArg L eb
  rw [la, lb]; grind

/-! ## 5. Nishimori: the polarization is tanh of the coupling -/

/-- On the Nishimori line e^{−2K} = w = p/(1 − p).  Then tanh K = (1 − w)/(1 + w) = 1 − 2p. -/
theorem nishimori_tanh (p : Rat) (hp : 1 - p ≠ 0) :
    (1 - p / (1 - p)) / (1 + p / (1 - p)) = 1 - 2 * p := by
  have e1 : 1 - p / (1 - p) = (1 - 2 * p) * (1 - p)⁻¹ := by
    rw [Rat.div_def]
    have : (1 - p) * (1 - p)⁻¹ = 1 := Rat.mul_inv_cancel _ hp
    grind
  have e2 : 1 + p / (1 - p) = (1 - p)⁻¹ := by
    rw [Rat.div_def]
    have : (1 - p) * (1 - p)⁻¹ = 1 := Rat.mul_inv_cancel _ hp
    grind
  rw [e1, e2, Rat.div_def, Rat.inv_inv]
  calc (1 - 2 * p) * (1 - p)⁻¹ * (1 - p) = (1 - 2 * p) * ((1 - p) * (1 - p)⁻¹) := by grind
    _ = 1 - 2 * p := by rw [Rat.mul_inv_cancel _ hp]; grind

/-- The Kramers–Wannier lens M(t) = (1 − t)/(1 + t) sends λ = 1 − 2p to p/(1 − p). -/
theorem kw_lens_polar (p : Rat) :
    (1 - (1 - 2 * p)) / (1 + (1 - 2 * p)) = p / (1 - p) := by
  have a : 1 - (1 - 2 * p) = 2 * p := by grind
  have b : 1 + (1 - 2 * p) = 2 * (1 - p) := by grind
  rw [a, b, Rat.div_def, Rat.div_def, Rat.inv_mul_rev]
  have : (2 : Rat) * 2⁻¹ = 1 := by decide +kernel
  calc 2 * p * ((1 - p)⁻¹ * 2⁻¹) = (2 * 2⁻¹) * (p * (1 - p)⁻¹) := by grind
    _ = p * (1 - p)⁻¹ := by rw [this]; grind

/- Bonds in series (summing out a spin) multiply tanh K: the bit-flip law,
   QuantumNoise.bit_flip (was decimation). -/

end EntropyNoise
