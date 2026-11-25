module

public import FusedAlgebra.Basic

/-!
# The shape's directions: which dials are truly separate

Three levels of "the same law": units (`a ↦ l * a`, the size of the measuring stick), rulers (any
change of coordinate near `0`) and global (how many periods the log of the group has).

## Main results

* **A. The weight is units.** `weight_is_units`: `l * wAdd (l * k) a b = wAdd k (l * a) (l * b)`.
  On a graph (three sites in a row, `path3`) the weight can no longer be scaled away:
  `graph_weight_survives`.
* **B. The landscape bend.** A projective ruler `x ↦ x / (1 + c * x)` slides `(p, q)` along the
  parabola of constant discriminant `p ^ 2 - 4 * q` (`ruler_num`, `ruler_den`,
  `ruler_keeps_disc`); units multiply it by `l ^ 2` (`zoom_disc`).
* **C. Rung height.** The ladder `E^s ∘ E = E ∘ E^s` forces the jets of the fractional rung
  (`ladder_forces_c3`, `ladder_forces_c4`, `ladder_jets_commute`); the rung-`s` law has
  discriminant `s` through third order (`rung_disc_is_height`) and a fourth-order coefficient that
  no point of the landscape plane has (`rung_leaves_plane`, `rung_crossings`, `jet_rank`).
* **D. The nome.** Gluing `a ↦ a ✳ (q - 1)` commutes with the law and the law descends to every
  quotient (`glue_commutes`, `glue_twice`, `glue_descends`). Torsion through the 3-division
  polynomial: `cusp_psi3`, `node_psi3`, `node_three_points`.
* **E. Several waves.** `pair_weight_survives`, `periods_outnumber_curves`.
* **F. Whole numbers.** The `p`-fold weighted power of `x` is `p` times something when `p ∣ k`
  and is `k ^ (p - 1) * x ^ p` plus `p` times something otherwise, for `p = 2, 3, 5`.

These are elementary polynomial identities. The 3-division polynomial is mathlib's
`WeierstrassCurve.Ψ₃` (`psi3_eq_eval_Ψ₃`); the determinant in `jet_rank` is `Matrix.det`. The
core-Lean original stated the identities over `ℤ`; here they hold in any commutative ring. The
jet computations stay over `ℚ`.
-/

@[expose] public section

namespace ShapeDirections

open MultAdd

section CommRing

variable {R : Type*} [CommRing R]

/-! ### A. The weight is the size of the measuring stick -/

/-- Zooming the inputs by `l` and the weight by `1 / l` gives the same law: the weight is
units. -/
theorem weight_is_units (k l a b : R) : l * wAdd (l * k) a b = wAdd k (l * a) (l * b) := by
  unfold wAdd; ring

/-- The `q` of the `q`-integers (and of the Potts colours) is the weight plus one. -/
theorem qcount_is_weight (q a b : R) : wAdd (q - 1) a b = a + b + q * a * b - a * b := by
  unfold wAdd; ring

/-- Three sites in a row, each site with fugacity `z` and each occupied bond with weight `k`:
`Z = 1 + 3 * z + (2 * k + 1) * z ^ 2 + k ^ 2 * z ^ 3`. -/
def path3 (k z : R) : R := 1 + 3 * z + (2 * k + 1) * z ^ 2 + k ^ 2 * z ^ 3

/-- The list of all eight occupation patterns adds up to the closed form. -/
theorem path3_sum (k z : R) :
    path3 k z = 1 + z + z + z + (z * z * k) + (z * z * k) + (z * z) + (z * z * z * k * k) := by
  unfold path3; ring

end CommRing

/-- On the graph the weight survives: if a rescaling of the sites (`z ↦ l * z`) and of the whole
sum (by `m`) turns the coefficients for weight `k` into those for weight `k'`, then `k = k'`. -/
theorem graph_weight_survives (k k' l m : ℚ)
    (h0 : m * 1 = 1) (h1 : m * (3 * l) = 3) (h2 : m * ((2 * k + 1) * l ^ 2) = 2 * k' + 1) :
    k = k' := by
  obtain rfl : m = 1 := by linarith
  obtain rfl : l = 1 := by linarith
  linarith

section CommRing

variable {R : Type*} [CommRing R]

/-! ### B. The landscape bend: rulers slide along parabolas, units zoom them -/

/-- Cross-multiplied numerator of `F_{p,q} (N a) (N b)` with `N x = x / (1 + c * x)`: it is the
numerator of `F_{p + 2 * c, ·} a b`. -/
theorem ruler_num (p c a b : R) :
    a * (1 + c * b) + b * (1 + c * a) + p * a * b = wAdd (p + 2 * c) a b := by
  unfold wAdd; ring

/-- Cross-multiplied denominator: `(1 + c * a) * (1 + c * b) - q * a * b` is the denominator of
`N` applied to `F_{p + 2 * c, q + p * c + c ^ 2} a b`. Together with `ruler_num`:
`N⁻¹ (F_{p,q} (N a) (N b)) = F_{p + 2 * c, q + p * c + c ^ 2} a b`. -/
theorem ruler_den (p q c a b : R) :
    (1 + c * a) * (1 + c * b) - q * a * b
      = (1 - (q + p * c + c ^ 2) * a * b) + c * wAdd (p + 2 * c) a b := by
  unfold wAdd; ring

/-- The ruler keeps the discriminant: `(p + 2 * c) ^ 2 - 4 * (q + p * c + c ^ 2) = p ^ 2 - 4 * q`.
The ruler's orbits are the parabolas `q = (p ^ 2 - D) / 4`, one for each `D`. -/
theorem ruler_keeps_disc (p q c : R) :
    (p + 2 * c) ^ 2 - 4 * (q + p * c + c ^ 2) = p ^ 2 - 4 * q := by
  ring

/-- Units zoom the discriminant: `(l * p) ^ 2 - 4 * (l ^ 2 * q) = l ^ 2 * (p ^ 2 - 4 * q)`. Only
its sign is free of units. -/
theorem zoom_disc (p q l : R) : (l * p) ^ 2 - 4 * (l ^ 2 * q) = l ^ 2 * (p ^ 2 - 4 * q) := by
  ring

end CommRing

/-! ### C. Rung height: the ladder forces the jets -/

/-- A jet `x + c₂ * x ^ 2 + c₃ * x ^ 3 + c₄ * x ^ 4 + c₅ * x ^ 5`, stored as its four
coefficients. -/
structure Jet where
  /-- The coefficient of `x ^ 2`. -/
  c2 : ℚ
  /-- The coefficient of `x ^ 3`. -/
  c3 : ℚ
  /-- The coefficient of `x ^ 4`. -/
  c4 : ℚ
  /-- The coefficient of `x ^ 5`. -/
  c5 : ℚ

/-- The coefficient of `x ^ 4` in the composition `f ∘ g` of two jets. -/
def comp4 (f g : Jet) : ℚ := f.c2 * (g.c2 * g.c2 + 2 * g.c3) + 3 * f.c3 * g.c2 + f.c4 + g.c4

/-- The coefficient of `x ^ 5` in the composition `f ∘ g` of two jets. -/
def comp5 (f g : Jet) : ℚ :=
  f.c2 * (2 * g.c2 * g.c3 + 2 * g.c4) + f.c3 * (3 * g.c2 * g.c2 + 3 * g.c3) + 4 * f.c4 * g.c2
    + f.c5 + g.c5

/-- `E x = exp x - 1` to fifth order. -/
def E : Jet := ⟨1 / 2, 1 / 6, 1 / 24, 1 / 120⟩

/-- The ladder at fourth order forces `c₃` of any rung with `c₂ = s / 2`. -/
theorem ladder_forces_c3 (s c3 c4 c5 : ℚ)
    (h4 : comp4 E ⟨s / 2, c3, c4, c5⟩ = comp4 ⟨s / 2, c3, c4, c5⟩ E) :
    c3 = s * (3 * s - 1) / 12 := by
  simp only [comp4, E] at h4
  linarith

/-- The ladder at fifth order then forces `c₄`. -/
theorem ladder_forces_c4 (s c3 c4 c5 : ℚ)
    (h4 : comp4 E ⟨s / 2, c3, c4, c5⟩ = comp4 ⟨s / 2, c3, c4, c5⟩ E)
    (h5 : comp5 E ⟨s / 2, c3, c4, c5⟩ = comp5 ⟨s / 2, c3, c4, c5⟩ E) :
    c4 = s * (2 * s - 1) * (3 * s - 1) / 48 := by
  obtain rfl := ladder_forces_c3 s c3 c4 c5 h4
  simp only [comp5, E] at h5
  linarith

/-- The forced jets do satisfy the ladder (so the rung exists to this order). -/
theorem ladder_jets_commute (s c5 : ℚ) :
    comp4 E ⟨s / 2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩
      = comp4 ⟨s / 2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩ E ∧
    comp5 E ⟨s / 2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩
      = comp5 ⟨s / 2, s * (3 * s - 1) / 12, s * (2 * s - 1) * (3 * s - 1) / 48, c5⟩ E := by
  simp only [comp4, comp5, E]
  constructor <;> ring

/-- The rung-`s` law is plain addition in the ruler `ℓ = E^(-s)`, whose jets are the forced ones
at `-s`: `ℓ₂ = -s / 2`, `ℓ₃ = s * (3 * s + 1) / 12`, `ℓ₄ = -s * (2 * s + 1) * (3 * s + 1) / 48`.
Its rate `R = 1 / ℓ'` has coefficients `r₁ = s`, `r₂ = s * (s - 1) / 4`,
`r₃ = -s * (s - 1) / 12`: these satisfy `ℓ' * R = 1` through `a ^ 3`. -/
theorem rung_rate_inverts (s : ℚ) :
    s + 2 * (-s / 2) = 0 ∧
    s * (s - 1) / 4 + 2 * (-s / 2) * s + 3 * (s * (3 * s + 1) / 12) = 0 ∧
    -(s * (s - 1)) / 12 + 2 * (-s / 2) * (s * (s - 1) / 4) + 3 * (s * (3 * s + 1) / 12) * s
      + 4 * (-(s * (2 * s + 1) * (3 * s + 1)) / 48) = 0 := by
  refine ⟨?_, ?_, ?_⟩ <;> ring

/-- The jets of `E^(-s)` are the forced jets with `s` replaced by `-s`. -/
theorem lens_jets_at_minus (s : ℚ) :
    (-s) * (3 * (-s) - 1) / 12 = s * (3 * s + 1) / 12 ∧
    (-s) * (2 * (-s) - 1) * (3 * (-s) - 1) / 48 = -(s * (2 * s + 1) * (3 * s + 1)) / 48 := by
  constructor <;> ring

/-- Through third order the rung's discriminant is its height: `r₁ ^ 2 - 4 * r₂ = s`. So rung `s`
looks, to third order, like the law with one period `2 * π * i / √s` (rung `1`: mult-add; rung
`0`: no period; rung `-1`: real period `2 * π`, a circle). -/
theorem rung_disc_is_height (s : ℚ) : s ^ 2 - 4 * (s * (s - 1) / 4) = s := by ring

/-- Every point of the landscape plane has rate `1 + p * a + q * a ^ 2` exactly, so `r₃ = 0`
there; the rung's `r₃ = -s * (s - 1) / 12` differs from `0` at `s = 2`, `s = -1` and
`s = 1 / 2`. -/
theorem rung_leaves_plane :
    -((2 : ℚ) * (2 - 1)) / 12 ≠ 0 ∧ -((-1 : ℚ) * (-1 - 1)) / 12 ≠ 0 ∧
    -((1 / 2 : ℚ) * (1 / 2 - 1)) / 12 ≠ 0 := by
  norm_num

/-- The rung crosses the plane at rungs `0` and `1` (`r₃ = 0` there), crossing, not touching:
the derivative `-(2 * s - 1) / 12` is `1 / 12` at `s = 0` and `-1 / 12` at `s = 1`. -/
theorem rung_crossings :
    -((0 : ℚ) * (0 - 1)) / 12 = 0 ∧ -((1 : ℚ) * (1 - 1)) / 12 = 0 ∧
    -(2 * (0 : ℚ) - 1) / 12 = 1 / 12 ∧ -(2 * (1 : ℚ) - 1) / 12 = -1 / 12 := by
  norm_num

/-- Jet-space rank at mult-add, coordinates `(r₁, r₂, r₃)`, scaled by `12`: zoom `(12, 0, 0)`,
weight `(12, 0, 0)`, bend `(0, 12, 0)`, rung `(12, 3, -1)`. The weight is the zoom (no new
direction); zoom, bend, rung have determinant `-144 ≠ 0`. -/
theorem jet_rank :
    (!![12, 0, 0; 0, 12, 0; 12, 3, -1] : Matrix (Fin 3) (Fin 3) ℤ).det = -144 ∧
    (!![12, 0, 0; 12, 0, 0; 0, 12, 0] : Matrix (Fin 3) (Fin 3) ℤ).det = 0 := by
  simp [Matrix.det_fin_three]

section CommRing

variable {R : Type*} [CommRing R]

/-! ### D. The nome: a gluing the law never sees -/

/-- Gluing `a ↦ a ✳ (q - 1)` commutes with the law, for every `q`. -/
theorem glue_commutes (q a b : R) : a ✳ (q - 1) ✳ b = a ✳ b ✳ (q - 1) :=
  multAdd_right_comm a (q - 1) b

/-- Two gluings are one gluing by the product: the gluings form the group `q ^ ℤ` (in mult-add
form). -/
theorem glue_twice (q q' a : R) : a ✳ (q - 1) ✳ (q' - 1) = a ✳ (q * q' - 1) := by
  unfold multAdd; ring

/-- The law descends to the quotient: glued inputs give a glued output, whatever the gluings. So
the torus `ℂˣ / q ^ ℤ` carries mult-add for every nome, and the nome is invisible near `0`. -/
theorem glue_descends (x y a b : R) :
    (a ✳ (x - 1)) ✳ (b ✳ (y - 1)) = a ✳ b ✳ (x * y - 1) := by
  unfold multAdd; ring

/-- The 3-division polynomial of `y ^ 2 + a₁ * x * y = x ^ 3 + a₂ * x ^ 2`
(`a₃ = a₄ = a₆ = 0`): `ψ₃ = 3 * x ^ 4 + b₂ * x ^ 3` with `b₂ = a₁ ^ 2 + 4 * a₂`. -/
def psi3 (a1 a2 x : R) : R := 3 * x ^ 4 + (a1 ^ 2 + 4 * a2) * x ^ 3

/-- `psi3` is mathlib's 3-division polynomial `WeierstrassCurve.Ψ₃` of the curve
`y ^ 2 + a₁ * x * y = x ^ 3 + a₂ * x ^ 2`, evaluated at `x`. -/
theorem psi3_eq_eval_Ψ₃ (a1 a2 x : R) :
    psi3 a1 a2 x = (WeierstrassCurve.Ψ₃ ⟨a1, a2, 0, 0, 0⟩ : Polynomial R).eval x := by
  simp [psi3, WeierstrassCurve.Ψ₃, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈]

/-- Cusp `y ^ 2 = x ^ 3` (addition): `ψ₃ = 3 * x ^ 4`, whose only root is the cusp itself: no
point other than the start comes back in three steps. -/
theorem cusp_psi3 (x : R) : psi3 0 0 x = 3 * x ^ 4 := by unfold psi3; ring

/-- Node `y ^ 2 + x * y = x ^ 3` (mult-add): `ψ₃ = x ^ 3 * (3 * x + 1)`. Besides the node
(`x = 0`) the root `x = -1 / 3`. -/
theorem node_psi3 (x : R) : psi3 1 0 x = x ^ 3 * (3 * x + 1) := by unfold psi3; ring

end CommRing

/-- At `x = -1 / 3` the node's equation `y ^ 2 + x * y - x ^ 3 = 0` reads
`y ^ 2 - y / 3 + 1 / 27 = 0`, whose discriminant `1 / 9 - 4 / 27 = -1 / 27` is not `0`: two
points. With the start, three: mult-add's three cube roots of `1`. -/
theorem node_three_points :
    discrim (1 : ℚ) (-1 / 3) (-(-1 / 3 : ℚ) ^ 3) = -1 / 27 ∧ (-1 / 27 : ℚ) ≠ 0 := by
  norm_num [discrim]

section CommRing

variable {R : Type*} [CommRing R]

/-! ### E. Several waves -/

/-- Two waves with pair weight `A`: `τ = 1 + w₁ + w₂ + A * w₁ * w₂`. -/
def tau2 (A w1 w2 : R) : R := 1 + w1 + w2 + A * w1 * w2

/-- `tau2` is the lift of the weighted mult-add. -/
theorem tau2_eq_one_add_wAdd (A w1 w2 : R) : tau2 A w1 w2 = 1 + wAdd A w1 w2 := by
  unfold tau2 wAdd; ring

/-- No separate rescaling of the waves changes the pair weight: if
`τ_A (l₁ * w₁) (l₂ * w₂) = τ_A' w₁ w₂` for all `w₁`, `w₂`, then `A = A'`. (A single law's weight
is units; a pair weight is not.) -/
theorem pair_weight_survives (A A' l1 l2 : R)
    (h : ∀ w1 w2 : R, tau2 A (l1 * w1) (l2 * w2) = tau2 A' w1 w2) : A = A' := by
  have h10 := h 1 0
  have h01 := h 0 1
  have h11 := h 1 1
  unfold tau2 at h10 h01 h11
  obtain rfl : l1 = 1 := by linear_combination h10
  obtain rfl : l2 = 1 := by linear_combination h01
  linear_combination h11

end CommRing

/-- Period matrices for `g` waves have `g * (g + 1) / 2` entries; curves of genus `g` have
`3 * g - 3` moduli. The first exceeds the second exactly from `g = 4` on: Schottky's problem. -/
theorem periods_outnumber_curves (g : ℕ) (hg : 2 ≤ g) : 6 * g - 6 < g * (g + 1) ↔ 4 ≤ g := by
  constructor
  · intro h
    by_contra hc
    obtain rfl | rfl : g = 2 ∨ g = 3 := by omega
    all_goals omega
  · intro h
    obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le h
    have : 6 * (4 + k) - 6 = 18 + 6 * k := by omega
    rw [this]
    nlinarith

section CommRing

variable {R : Type*} [CommRing R]

/-! ### F. Whole numbers: the weight keeps its primes -/

theorem pseries_two (k x : R) : wpow k x 2 = 2 * x + k * x ^ 2 := by
  simp only [wpow, wAdd]; ring

theorem pseries_three (k x : R) : wpow k x 3 = 3 * x + 3 * k * x ^ 2 + k ^ 2 * x ^ 3 := by
  simp only [wpow, wAdd]; ring

theorem pseries_five (k x : R) :
    wpow k x 5 = 5 * x + 10 * k * x ^ 2 + 10 * k ^ 2 * x ^ 3 + 5 * k ^ 3 * x ^ 4
      + k ^ 4 * x ^ 5 := by
  simp only [wpow, wAdd]; ring

/-- `2 ∣ k`: the 2-fold power is `2` times a whole number (the law is additive at `2`). -/
theorem additive_at_two (m x : R) : wpow (2 * m) x 2 = 2 * (x + m * x ^ 2) := by
  rw [pseries_two]; ring

/-- `3 ∣ k`: the 3-fold power is `3` times a whole number. -/
theorem additive_at_three (m x : R) :
    wpow (3 * m) x 3 = 3 * (x + 3 * m * x ^ 2 + 3 * m ^ 2 * x ^ 3) := by
  rw [pseries_three]; ring

/-- `5 ∣ k`: the 5-fold power is `5` times a whole number. -/
theorem additive_at_five (m x : R) :
    wpow (5 * m) x 5 = 5 * (x + 10 * m * x ^ 2 + 50 * m ^ 2 * x ^ 3
      + 125 * m ^ 3 * x ^ 4 + 125 * m ^ 4 * x ^ 5) := by
  rw [pseries_five]; ring

/-- The 3-fold power is `k ^ 2 * x ^ 3` plus `3` times a whole number, and `k ^ 2` is `1` more
than a multiple of `3` when `3 ∤ k` (`fermat_three_a`, `fermat_three_b`; in general Fermat's
little theorem, `ZMod.pow_card_sub_one_eq_one`): the law is multiplicative at `3`. -/
theorem mult_at_three (k x : R) : wpow k x 3 - k ^ 2 * x ^ 3 = 3 * (x + k * x ^ 2) := by
  rw [pseries_three]; ring

theorem fermat_three_a (m : R) : (3 * m + 1) ^ 2 - 1 = 3 * (3 * m ^ 2 + 2 * m) := by ring

theorem fermat_three_b (m : R) : (3 * m + 2) ^ 2 - 1 = 3 * (3 * m ^ 2 + 4 * m + 1) := by ring

/-- The 2-fold power is `k * x ^ 2` plus `2` times a whole number. -/
theorem mult_at_two (k x : R) : wpow k x 2 - k * x ^ 2 = 2 * x := by
  rw [pseries_two]; ring

end CommRing

end ShapeDirections
