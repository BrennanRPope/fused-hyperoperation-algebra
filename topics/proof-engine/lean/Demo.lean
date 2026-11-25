import FusedPictures
open CircleOp FusedPictures Lean.Grind
set_option trace.fused.pictures true

namespace FusedPictures.Demo
variable {R : Type} [CommRing R]

/-- circ_comm, again: lift, open the box, flip it on its diagonal. -/
theorem comm_by_moves (a b : R) : circ a b = circ b a := by
  lift
  multadd
  swap

/-- circ_assoc, again: slabs and columns are one box. -/
theorem assoc_by_moves (a b c : R) : circ (circ a b) c = circ a (circ b c) := by
  lift
  multadd
  regroup

/-- circ_zero, again: a side of length 0 flattens away. -/
theorem zero_by_moves (a : R) : circ a 0 = a := by
  lift
  multadd
  collapse

/-- neg_one_absorbs, again: the wall. A side of length −1 makes the lifted rod 0. -/
theorem wall_by_moves (b : R) : circ (-1) b = -1 := by
  lift
  multadd
  collapse

/-- mpow_add, again: a cube of m + n sides is a cube of m by a cube of n. -/
theorem mpow_add_by_moves (a : R) (m n : Nat) : mpow a (m + n) = circ (mpow a m) (mpow a n) := by
  lift
  multadd
  flatten

/-- circK_comm (scaled), again: dial the corner, then flip. -/
theorem wcomm_by_moves (k a b : R) : k * circK k a b = k * circK k b a := by
  dial
  dial
  lift
  multadd
  swap

/-- fexp_iter, again (Nat, rung 2): one more slab extrudes the fused exponent. -/
theorem fexp_iter_by_moves (a b : Nat) : fexp a (b + 1) = circ a (fexp a b) := by
  lift
  extrude

/-- circ_circInv, again (fields): squash the second rod to 1/(1+u); the rectangle
    has area 1, so the mult-add is 0. -/
theorem inv_by_moves {K : Type} [Field K] (u : K) (h : 1 + u ≠ 0) :
    circ u (circInv u) = 0 := by
  lift
  multadd
  squash u h
  collapse

end FusedPictures.Demo

#shape CircleOp.circ_assoc
#shape CircleOp.circK_lift
#shape CircleOp.circ_circInv
#shape CircleOp.mpow_lift
#print axioms FusedPictures.Demo.inv_by_moves
#print axioms FusedPictures.Demo.mpow_add_by_moves
#print axioms FusedPictures.Demo.assoc_by_moves

-- goal pictures for the proof films (read by film.py)
#shape FusedPictures.Demo.comm_by_moves
#shape FusedPictures.Demo.assoc_by_moves
#shape FusedPictures.Demo.zero_by_moves
#shape FusedPictures.Demo.wall_by_moves
#shape FusedPictures.Demo.mpow_add_by_moves
#shape FusedPictures.Demo.wcomm_by_moves
#shape FusedPictures.Demo.fexp_iter_by_moves
#shape FusedPictures.Demo.inv_by_moves
