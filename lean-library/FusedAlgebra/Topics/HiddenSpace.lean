/-
  Hidden space built from all the primes (DRAFT, Nov 27 2026). Core Lean 4 only, no imports, no `sorry`.
  NOT COMPILED where it was written (no Lean toolchain there). Built with Lean 4.37 in the repo on Nov 27,
  no proof changed; in the library as FusedAlgebra.Topics.HiddenSpace (namespace CircleOp.HiddenSpace).

  The finite piece of the hidden space seen through 2, 3, 5, 7 is the residues mod 210.
  1. The grid (n mod 14, n mod 15) gives each residue its own cell (Chinese remainder).
  2. The prime window there, "coprime to 210", has 48 cells = (2-1)(3-1)(5-1)(7-1),
     the window share 48/210 = product of (1 - 1/p).
-/
namespace CircleOp.HiddenSpace

theorem crt_grid_210 : ((List.range 210).map fun r => (r % 14, r % 15)).Nodup := by decide +kernel

theorem window_210 :
    ((List.range 210).filter fun r => Nat.gcd r 210 == 1).length = (2-1)*(3-1)*(5-1)*(7-1) := by decide +kernel

end CircleOp.HiddenSpace
