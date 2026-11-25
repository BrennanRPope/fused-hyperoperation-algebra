import ProofComplexity.Measure

open Lean Elab Command Meta

def glue (n : Name) : Bool :=
  let s := n.toString
  s.startsWith "Mathlib.Tactic" || s.startsWith "Mathlib.Meta" || s.startsWith "Lean." ||
  s.startsWith "Eq." || s.startsWith "Iff." || s.startsWith "And." || s.startsWith "Or." ||
  s.startsWith "Decidable." || s.startsWith "Classical." || s.startsWith "Bool." ||
  ["congrArg", "congr", "congrFun", "congrFun'", "of_eq_true", "eq_self", "id", "funext", "propext",
   "Not.intro", "True.intro", "absurd", "Eq.mpr", "rfl", "trivial", "eq_true", "implies_congr",
   "forall_congr", "forall_congr'", "ne_eq", "not_false_eq_true", "imp_self", "iff_self",
   "eq_comm", "HEq.refl", "cast_eq", "Subtype.property", "Exists.intro", "of_eq_false",
   "eq_false", "Nat.le_refl", "le_refl", "dite_congr", "ite_congr"].contains s

elab "#tree" : command => liftTermElabM do
  let env ← getEnv
  let roots : List Name := [``Finset.sum_powerset_neg_one_pow_card_of_nonempty,
    ``HCT.Results.op_sum_powerset_neg_one_pow_card_of_nonempty,
    ``HCT.Control.ctl_sum_powerset_neg_one_pow_card_of_nonempty]
  let mut out := ""
  for r in roots do
    let mut level : List (Name × Nat) := [(r, 0)]
    let mut seen : NameSet := {}
    while !level.isEmpty do
      match level with
      | [] => break
      | (n, d) :: rest =>
        level := rest
        if seen.contains n then continue
        seen := seen.insert n
        let kids := (Yard.directThms env n).toList.filter (fun k => !glue k)
        out := out ++ s!"{r}\t{d}\t{n}\t{← Yard.declLines n}\t{kids}\n"
        if d < 3 then
          level := level ++ kids.map (fun k => (k, d + 1))
  IO.FS.writeFile "tree.tsv" out
  logInfo "ok"
#tree
