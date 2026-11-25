/-
  Axioms.Rederive: the library's theorems, re-derived from the axioms
  ====================================================================

  1. Every library theorem stated for any commutative semiring, ring or field
     is re-derived for an arbitrary MA-ring (or MA-field, for the field ones):
     its hypotheses are now consequences of M1–M4 (+ F1–F2). The command
     `#rederive` below makes one new theorem `MAFound.Rederived.<name>` for
     each, and Lean's kernel checks it. The ones stated for commutative
     semirings are also re-derived for any model of the fused Peano axioms.

  2. In both settings the library's `multAdd` (a + b + ab with the derived
     + and ×) is the primitive mult-add (`multAdd_eq_madd`, `multAdd_eq_madd'`),
     so these are theorems about the primitive.

  3. The rung laws of `FusedAlgebra.Rungs` are restated for an arbitrary model
     of the fused Peano axioms, with each rung given by the one ladder schema
     (`MAFound.Rederived.Rungs`).
-/
import FusedAlgebra
import Axioms.MultAddRing
import Axioms.FusedPeano
import Lean

open Lean Elab Command Meta

namespace MAFound

theorem multAdd_eq_madd {R : Type} [MARing R] (a b : R) :
    MultAdd.multAdd a b = MARing.madd a b := MARing.multAdd_eq a b

theorem multAdd_eq_madd' {N : Type} [FusedPeano N] (a b : N) :
    MultAdd.multAdd a b = FusedPeano.madd a b := FusedPeano.multAdd_eq a b

end MAFound

/-- Instantiate theorem `thm` at a fresh type `R` carrying class `cls`. -/
def rederiveOne (thm : Name) (cls : Name) (newName : Name) : MetaM Bool := do
  let some ci := (← getEnv).find? thm | return false
  withLocalDecl `R .implicit (mkSort Level.one) fun R => do
  withLocalDecl `inst .instImplicit (mkApp (mkConst cls) R) fun inst => do
    let c ← mkConstWithFreshMVarLevels thm
    let mut e := c
    let mut ty ← inferType c
    let mut usedR := false
    repeat
      ty ← whnfR ty
      match ty with
      | .forallE _ d b bi =>
        if d.isSort && !usedR then
          unless ← isDefEq d (mkSort Level.one) do return false
          e := mkApp e R; ty := b.instantiate1 R; usedR := true
        else if bi.isInstImplicit && usedR then
          let i? ← (try some <$> synthInstance d catch _ => pure none)
          let i? ← match i? with
            | some i => pure (some i)
            | none =>
              -- the toric files ask for decidable equality; classical logic supplies it
              if d.getAppFn.constName? == some ``DecidableEq then
                pure (some (← mkAppM ``Classical.typeDecidableEq #[R]))
              else pure none
          let some i := i? | return false
          e := mkApp e i; ty := b.instantiate1 i
        else break
      | _ => break
    unless usedR do return false
    let e' ← instantiateMVars e
    let ety ← instantiateMVars (← inferType e')
    if e'.hasMVar || ety.hasMVar then return false
    let stmt ← mkForallFVars #[R, inst] ety
    let val ← mkLambdaFVars #[R, inst] e'
    let _ := ci
    addDecl (Declaration.thmDecl
      { name := newName
        levelParams := []
        type := stmt
        value := val })
    return true

elab "#rederive" : command => do
  let env ← getEnv
  let mut okMA : Array Name := #[]
  let mut okPeano : Array Name := #[]
  let mut failed : Array Name := #[]
  let mut thms : Array (Name × String) := #[]
  for h : i in [0:env.header.moduleNames.size] do
    let mod := env.header.moduleNames[i]
    unless (`FusedAlgebra).isPrefixOf mod do continue
    for n in env.header.moduleData[i]!.constNames do
      let some ci := env.find? n | continue
      unless ci matches ConstantInfo.thmInfo _ do continue
      if n.isInternalDetail then continue
      let some _ ← liftCoreM (findDeclarationRanges? n) | continue
      let cls ← liftTermElabM do
        forallTelescope ci.type fun xs _ => do
          let mut c := ""
          for x in xs do
            match (← inferType x).getAppFn.constName? with
            | some `Lean.Grind.Field => c := "Field"
            | some `Lean.Grind.CommRing | some `Lean.Grind.Ring =>
                if c != "Field" then c := "Ring"
            | some `Lean.Grind.CommSemiring | some `Lean.Grind.Semiring =>
                if c == "" then c := "Semiring"
            | _ => pure ()
          return c
      unless cls == "" do thms := thms.push (n, cls)
  for (n, cls) in thms do
    let target := if cls == "Field" then ``MAFound.MAField else ``MAFound.MARing
    let ok ← liftTermElabM <| rederiveOne n target (`MAFound.Rederived ++ n)
    if ok then okMA := okMA.push n else failed := failed.push n
    if cls == "Semiring" then
      let ok2 ← liftTermElabM <| rederiveOne n ``MAFound.FusedPeano (`MAFound.RederivedPeano ++ n)
      if ok2 then okPeano := okPeano.push n
  let rows := okMA.map (fun n => s!"{n}\tMA-ring/field") ++ okPeano.map (fun n => s!"{n}\tfused Peano")
  IO.FS.writeFile "build/rederived.tsv" ("\n".intercalate rows.toList ++ "\n")
  logInfo m!"generic library theorems: {thms.size}; re-derived from the MA axioms: {okMA.size}; \
    also from the fused Peano axioms: {okPeano.size}; failed: {failed}"

#rederive

/-! ## The rung laws, on any model of the fused Peano axioms -/

namespace MAFound.Rederived.Rungs

open MAFound FusedPeano

theorem lad_eq_sadd (a b : Nat) : lad 1 a b = _root_.MultAdd.sadd a b := by rw [lad_one]; rfl
theorem lad_eq_multAdd (a b : Nat) : lad 2 a b = _root_.MultAdd.multAdd a b := by rw [lad_two]; rfl
theorem lad_eq_fexp (a b : Nat) : lad 3 a b = _root_.MultAdd.fexp a b := by
  induction b with
  | zero => rw [_root_.MultAdd.fexp_zero_right]; rfl
  | succ b ih => rw [_root_.MultAdd.fexp_iter, ← ih, ← lad_eq_multAdd]; rfl
theorem lad_eq_ftet (a b : Nat) : lad 4 a b = _root_.MultAdd.ftet a b := by
  induction b with
  | zero => rw [_root_.MultAdd.ftet_zero_right]; rfl
  | succ b ih => rw [_root_.MultAdd.ftet_iter, ← ih, ← lad_eq_fexp]; rfl
theorem lad_eq_fpent (a b : Nat) : lad 5 a b = _root_.MultAdd.fpent a b := by
  induction b with
  | zero => rw [_root_.MultAdd.fpent_zero_right]; rfl
  | succ b ih => rw [_root_.MultAdd.fpent_iter, ← ih, ← lad_eq_ftet]; rfl

variable {N : Type} [FusedPeano N]

theorem ladN_toN (k a b : Nat) : ladN k (toN (N := N) a) (toN b) = toN (lad k a b) :=
  (toN_lad k a b).symm

theorem toN_zero : (zero : N) = toN 0 := rfl
theorem toN_succ (n : N) : succ n = toN (fromN n + 1) := by
  show succ n = succ (toN (fromN n)); rw [toN_fromN]

macro "lift_rungs" : tactic => `(tactic| (
  simp only [← ladN_one, ← ladN_two, toN_zero, ladN_toN,
    lad_eq_sadd, lad_eq_multAdd, lad_eq_fexp, lad_eq_ftet, lad_eq_fpent]))

/-- Each rung is the rung below, iterated (library: `multAdd_iter`). -/
theorem multAdd_iter (x y : N) : madd x (succ y) = sadd x (madd x y) := by
  rw [← ladN_one, ← ladN_two, ← ladN_two, ladN_succ]
/-- library: `fexp_iter` -/
theorem fexp_iter (x y : N) : ladN 3 x (succ y) = madd x (ladN 3 x y) := by
  rw [← ladN_two, ladN_succ]
/-- library: `ftet_iter` -/
theorem ftet_iter (x y : N) : ladN 4 x (succ y) = ladN 3 x (ladN 4 x y) := ladN_succ 3 x y
/-- library: `fpent_iter` -/
theorem fpent_iter (x y : N) : ladN 5 x (succ y) = ladN 4 x (ladN 5 x y) := ladN_succ 4 x y

/-- 0 is the right identity at rungs 3, 4, 5 (library: `fexp_zero_right`, …). -/
theorem fexp_zero_right (x : N) : ladN 3 x zero = x := by rw [ladN_zero_right]; rfl
theorem ftet_zero_right (x : N) : ladN 4 x zero = x := by rw [ladN_zero_right]; rfl
theorem fpent_zero_right (x : N) : ladN 5 x zero = x := by rw [ladN_zero_right]; rfl

/-- 0 absorbs from the left (library: `fexp_zero_left`, `ftet_zero_left`, `fpent_zero_left`). -/
theorem fexp_zero_left (y : N) : ladN 3 zero y = zero := by
  rw [← toN_fromN y]; lift_rungs; rw [_root_.MultAdd.fexp_zero_left]
theorem ftet_zero_left (y : N) : ladN 4 zero y = zero := by
  rw [← toN_fromN y]; lift_rungs; rw [_root_.MultAdd.ftet_zero_left]
theorem fpent_zero_left (y : N) : ladN 5 zero y = zero := by
  rw [← toN_fromN y]; lift_rungs; rw [_root_.MultAdd.fpent_zero_left]

/-- Rung 3 distributes over mult-add (library: `fexp_multAdd_right`). -/
theorem fexp_multAdd_right (x y z : N) :
    ladN 3 (madd x y) z = madd (ladN 3 x z) (ladN 3 y z) := by
  rw [← toN_fromN x, ← toN_fromN y, ← toN_fromN z]; lift_rungs; rw [_root_.MultAdd.fexp_multAdd_right]

/-- The tower law (library: `fexp_tower`). -/
theorem fexp_tower (x y z : N) : ladN 3 (ladN 3 x y) z = ladN 3 x (madd y z) := by
  rw [← toN_fromN x, ← toN_fromN y, ← toN_fromN z]; lift_rungs; rw [_root_.MultAdd.fexp_tower]

/-- The sum law (library: `fexp_sadd`). -/
theorem fexp_sadd (x y z : N) : ladN 3 x (sadd y z) = madd (ladN 3 x y) (ladN 3 x z) := by
  rw [← toN_fromN x, ← toN_fromN y, ← toN_fromN z]; lift_rungs; rw [_root_.MultAdd.fexp_sadd]

theorem toN_ne {a b : Nat} (h : a ≠ b) : toN (N := N) a ≠ toN b := fun e => h (toN_injective e)

/-- Rung 3 is neither associative nor commutative (library: `fexp_not_assoc`, `fexp_not_comm`). -/
theorem fexp_not_assoc :
    ladN 3 (ladN 3 (toN 1 : N) (toN 2)) (toN 1) ≠ ladN 3 (toN 1) (ladN 3 (toN 2) (toN 1)) := by
  lift_rungs; exact toN_ne _root_.MultAdd.fexp_not_assoc
theorem fexp_not_comm : ladN 3 (toN 1 : N) (toN 2) ≠ ladN 3 (toN 2) (toN 1) := by
  lift_rungs; exact toN_ne _root_.MultAdd.fexp_not_comm
theorem fexp_example :
    ladN 3 (ladN 3 (toN 1 : N) (toN 2)) (toN 1) = toN 63 ∧
      ladN 3 (toN 1 : N) (ladN 3 (toN 2) (toN 1)) = toN 511 := by
  lift_rungs; rw [_root_.MultAdd.fexp_example.1, _root_.MultAdd.fexp_example.2]; exact ⟨rfl, rfl⟩

/-- Rung 4 values and failures (library: `ftet_values`, `fpent_values`, `ftet_not_assoc`, …). -/
theorem ftet_values :
    ladN 4 (toN 1 : N) (toN 1) = toN 3 ∧ ladN 4 (toN 1 : N) (toN 2) = toN 15 ∧
      ladN 4 (toN 2 : N) (toN 1) = toN 26 ∧ ladN 4 (toN 3 : N) (toN 1) = toN 255 ∧
      ladN 4 (toN 1 : N) (toN 3) = toN 65535 := by
  lift_rungs
  obtain ⟨h1, h2, h3, h4, h5⟩ := _root_.MultAdd.ftet_values
  rw [h1, h2, h3, h4, h5]; exact ⟨rfl, rfl, rfl, rfl, rfl⟩
theorem fpent_values : ladN 5 (toN 1 : N) (toN 1) = toN 3 ∧ ladN 5 (toN 1 : N) (toN 2) = toN 65535 := by
  lift_rungs; rw [_root_.MultAdd.fpent_values.1, _root_.MultAdd.fpent_values.2]; exact ⟨rfl, rfl⟩
theorem ftet_not_assoc :
    ladN 4 (ladN 4 (toN 1 : N) (toN 1)) (toN 1) ≠ ladN 4 (toN 1) (ladN 4 (toN 1) (toN 1)) := by
  lift_rungs; exact toN_ne _root_.MultAdd.ftet_not_assoc
theorem ftet_not_comm : ladN 4 (toN 1 : N) (toN 2) ≠ ladN 4 (toN 2) (toN 1) := by
  lift_rungs; exact toN_ne _root_.MultAdd.ftet_not_comm
/-- The rung-3 laws stop at rung 4 (library: `ftet_no_multAdd_dist`, `ftet_no_tower`, `ftet_no_sum`). -/
theorem ftet_no_multAdd_dist :
    ladN 4 (madd (toN 1 : N) (toN 1)) (toN 1) ≠ madd (ladN 4 (toN 1) (toN 1)) (ladN 4 (toN 1) (toN 1)) := by
  lift_rungs; exact toN_ne _root_.MultAdd.ftet_no_multAdd_dist
theorem ftet_no_tower :
    ladN 4 (ladN 4 (toN 1 : N) (toN 1)) (toN 1) ≠ ladN 4 (toN 1) (ladN 3 (toN 1) (toN 1)) := by
  lift_rungs; exact toN_ne _root_.MultAdd.ftet_no_tower
theorem ftet_no_sum :
    ladN 4 (toN 1 : N) (sadd (toN 1) (toN 1)) ≠ ladN 3 (ladN 4 (toN 1) (toN 1)) (ladN 4 (toN 1) (toN 1)) := by
  lift_rungs; exact toN_ne _root_.MultAdd.ftet_no_sum

end MAFound.Rederived.Rungs
