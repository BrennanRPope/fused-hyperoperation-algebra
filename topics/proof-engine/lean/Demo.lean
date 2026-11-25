import FusedPictures
open MultAdd FusedPictures Lean.Grind
set_option trace.fused.pictures true

namespace FusedPictures.Demo
variable {R : Type} [CommRing R]

/-- multAdd_comm, again: lift, open the box, flip it on its diagonal. -/
theorem comm_by_moves (a b : R) : multAdd a b = multAdd b a := by
  lift
  multadd
  swap

/-- multAdd_assoc, again: slabs and columns are one box. -/
theorem assoc_by_moves (a b c : R) : multAdd (multAdd a b) c = multAdd a (multAdd b c) := by
  lift
  multadd
  regroup

/-- multAdd_zero, again: a side of length 0 flattens away. -/
theorem zero_by_moves (a : R) : multAdd a 0 = a := by
  lift
  multadd
  collapse

/-- neg_one_absorbs, again: the wall. A side of length −1 makes the lifted rod 0. -/
theorem wall_by_moves (b : R) : multAdd (-1) b = -1 := by
  lift
  multadd
  collapse

/-- mpow_add, again: a cube of m + n sides is a cube of m by a cube of n. -/
theorem mpow_add_by_moves (a : R) (m n : Nat) : mpow a (m + n) = multAdd (mpow a m) (mpow a n) := by
  lift
  multadd
  flatten

/-- wAdd_comm (scaled), again: dial the corner, then flip. -/
theorem wcomm_by_moves (k a b : R) : k * wAdd k a b = k * wAdd k b a := by
  dial
  dial
  lift
  multadd
  swap

/-- fexp_iter, again (Nat, rung 2): one more slab extrudes the fused exponent. -/
theorem fexp_iter_by_moves (a b : Nat) : fexp a (b + 1) = multAdd a (fexp a b) := by
  lift
  extrude

/-- multAdd_multAddInv, again (fields): squash the second rod to 1/(1+u); the rectangle
    has area 1, so the mult-add is 0. -/
theorem inv_by_moves {K : Type} [Field K] (u : K) (h : 1 + u ≠ 0) :
    multAdd u (multAddInv u) = 0 := by
  lift
  multadd
  squash u h
  collapse

end FusedPictures.Demo

#shape MultAdd.multAdd_assoc
#shape MultAdd.wAdd_lift
#shape MultAdd.multAdd_multAddInv
#shape MultAdd.mpow_lift
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
