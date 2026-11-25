/-
  FusedPictures: a geometric proof engine for Fused Hyperoperation Algebra.

  STATUS: compiled with core Lean 4.37.0-rc1 against the shared library (no mathlib,
  no ProofWidgets, no `sorry`).

  What it is
  ----------
  1. Shape     : a small datatype for pictures (rods, boxes with corner weights, squashes).
  2. Reader    : reads any expression or goal and finds the boxes in it: lifted products
                 (1+x)(1+y)..., mult-adds a ◯ b ◯ c (a box with its unit removed), the
                 weighted family (a dialed corner), mult-add powers (a cube of n equal
                 sides), inverses (a squash onto the hyperbola). The picture is computed
                 from the goal, never drawn by hand per theorem.
  3. Widget    : one inline JavaScript module that draws a Shape as SVG in the infoview
                 (VS Code). Outside an editor, the same shape is written as JSON to the
                 `fused.pictures` trace, which `lean --json` prints.
  4. Moves     : named tactics (lift, multadd, swap, regroup, collapse, squash, dial,
                 flatten, extrude) that each run an ordinary rewrite with a library lemma,
                 then attach the picture of the new goal.
  The kernel checks every rewrite. The pictures are a view of the goal and carry no trust.

  Dependency direction (one-way on purpose): the engine imports the shared library;
  the library never imports the engine. The move table lives here, so no library
  lemma needs a picture attribute.
-/
import FusedAlgebra
import Lean

open Lean Elab Tactic Meta
open Lean.Grind

namespace FusedPictures

/-! ## 1. Shapes -/

/-- One side of a box: the variable as printed and its colour key (a, b, c, ...). -/
structure Side where
  name : String
  key  : String
  deriving ToJson, FromJson, Inhabited, BEq, Repr

/-- A box (1+x₁)(1+x₂)…(1+xₙ): its 2ⁿ pieces are the corners of the n-cube.
    * `weight`: the far corner's weight k (the a + b + k·ab family), if dialed.
    * `power`: the box has this many equal sides (a symbolic count, e.g. "n").
    * `dropUnit`: the unit cube is removed (a mult-add, as opposed to a product).
    * `cut`: how the box was built, e.g. "((0 1) 2)" = slabs, "(0 (1 2))" = columns.
    * `any`: a list of rods, so any number of sides. -/
structure Box where
  sides    : Array Side
  weight   : Option String := none
  power    : Option String := none
  dropUnit : Bool := true
  cut      : String := ""
  any      : Bool := false
  deriving ToJson, FromJson, Inhabited, BEq, Repr

inductive Shape where
  /-- a rod: one length; `lifted` means 1 + x (the unit attached) -/
  | rod    (side : Side) (lifted : Bool)
  /-- a box, with or without its unit -/
  | box    (b : Box)
  /-- one corner piece only: a plain product x·y·z, no unit, no faces -/
  | piece  (sides : Array Side)
  /-- a rod rescaled to 1/(1+a): area exactly 1 on the hyperbola xy = 1 -/
  | squash (inner : Shape)
  /-- a mirror image: minus sign -/
  | flip   (inner : Shape)
  /-- a plain power x^n that is not a lifted box: a tower one rung up -/
  | tower  (base : Shape) (exp : String)
  /-- a sum of pieces laid end to end -/
  | sum    (parts : Array Shape)
  /-- a relation between two pictures: "=", "≤", "<", "≠", "∣", "↔" -/
  | rel    (r : String) (lhs rhs : Shape)
  /-- several facts at once (∧) -/
  | both   (parts : Array Shape)
  | note   (text : String)
  deriving Inhabited, Repr

partial def Shape.toJson : Shape → Json
  | .rod s l     => Json.mkObj [("kind", "rod"), ("side", Lean.toJson s), ("lifted", l)]
  | .box b       => Json.mkObj [("kind", "box"), ("box", Lean.toJson b)]
  | .piece ss    => Json.mkObj [("kind", "piece"), ("sides", Lean.toJson ss)]
  | .squash s    => Json.mkObj [("kind", "squash"), ("inner", s.toJson)]
  | .flip s      => Json.mkObj [("kind", "flip"), ("inner", s.toJson)]
  | .tower b e   => Json.mkObj [("kind", "tower"), ("base", b.toJson), ("exp", e)]
  | .sum ps      => Json.mkObj [("kind", "sum"), ("parts", Json.arr (ps.map Shape.toJson))]
  | .rel r l rr  => Json.mkObj [("kind", "rel"), ("rel", r), ("lhs", l.toJson), ("rhs", rr.toJson)]
  | .both ps     => Json.mkObj [("kind", "both"), ("parts", Json.arr (ps.map Shape.toJson))]
  | .note t      => Json.mkObj [("kind", "note"), ("text", t)]

/-- The size of the largest box in a picture (a symbolic power counts as 99). -/
partial def Shape.maxSides : Shape → Nat
  | .rod _ _    => 1
  | .box b      => if b.power.isSome || b.any then 99 else b.sides.size
  | .piece ss   => ss.size
  | .squash s   => s.maxSides
  | .flip s     => s.maxSides
  | .tower b _  => max b.maxSides 2
  | .sum ps | .both ps => ps.foldl (fun m s => max m s.maxSides) 0
  | .rel _ l r  => max l.maxSides r.maxSides
  | .note _     => 0

/-- Marks a picture carries: dialed corners, squashes, flips, towers, unit removed. -/
partial def Shape.marks : Shape → List String
  | .rod _ l    => if l then ["lifted"] else []
  | .box b      => (if b.weight.isSome then ["dial"] else []) ++
                   (if b.power.isSome then ["cube"] else []) ++
                   (if b.dropUnit then ["unit-removed"] else ["lifted"]) ++
                   (if b.any then ["lists"] else [])
  | .piece _    => ["piece"]
  | .squash s   => "squash" :: s.marks
  | .flip s     => "flip" :: s.marks
  | .tower b _  => "tower" :: b.marks
  | .sum ps | .both ps => ps.toList.flatMap Shape.marks
  | .rel _ l r  => l.marks ++ r.marks
  | .note _     => []

/-! ## 2. Reader: from an expression to a shape -/

/-- Colour keys are given out in order of first appearance, so one variable keeps one
    colour on both sides of an equation. -/
abbrev ReadM := StateT (Array String) MetaM

def keyOf (i : Nat) : String :=
  if i < 26 then String.singleton (Char.ofNat ('a'.toNat + i)) else s!"v{i}"

def sideOf (e : Expr) : ReadM Side := do
  let nm := toString (← ppExpr e)
  let ks ← get
  match ks.findIdx? (· == nm) with
  | some i => return { name := nm, key := keyOf i }
  | none   => set (ks.push nm); return { name := nm, key := keyOf ks.size }

def lastArgs (e : Expr) (k : Nat) : Array Expr :=
  let as := e.getAppArgs
  as.extract (as.size - k) as.size

def constIs (e : Expr) (n : Name) : Bool := e.getAppFn.isConstOf n

/-- Is `e` the numeral `n` (in any semiring)? -/
def isLit (e : Expr) (n : Nat) : Bool :=
  match e.consumeMData with
  | .lit (.natVal m) => m == n
  | e => if constIs e ``OfNat.ofNat && e.getAppNumArgs ≥ 2 then
            match (e.getArg! 1).consumeMData with
            | .lit (.natVal m) => m == n
            | _ => false
         else if n == 1 then constIs e ``One.one else if n == 0 then constIs e ``Zero.zero else false

/-- If `e` is `1 + x` or `x + 1`, return `x`. -/
def asLifted? (e : Expr) : Option Expr :=
  if constIs e ``HAdd.hAdd && e.getAppNumArgs == 6 then
    let x := e.getArg! 4; let y := e.getArg! 5
    if isLit x 1 then some y else if isLit y 1 then some x else none
  else none

/-- Flatten a tree of one binary operation, remembering how it was bracketed. -/
partial def flattenOp (isOp : Expr → Option (Expr × Expr)) (e : Expr) (off : Nat) :
    Array Expr × String :=
  match isOp e with
  | some (l, r) =>
      let (fl, sl) := flattenOp isOp l off
      let (fr, sr) := flattenOp isOp r (off + fl.size)
      (fl ++ fr, "(" ++ sl ++ " " ++ sr ++ ")")
  | none => (#[e], toString off)

def maOp? (e : Expr) : Option (Expr × Expr) :=
  if constIs e ``CircleOp.circ then let a := lastArgs e 2; some (a[0]!, a[1]!) else none

def mulOp? (e : Expr) : Option (Expr × Expr) :=
  if constIs e ``HMul.hMul && e.getAppNumArgs == 6 then some (e.getArg! 4, e.getArg! 5) else none

def addOp? (e : Expr) : Option (Expr × Expr) :=
  if constIs e ``HAdd.hAdd && e.getAppNumArgs == 6 then some (e.getArg! 4, e.getArg! 5) else none

mutual
/-- Read one expression as a picture. -/
partial def read (e0 : Expr) : ReadM Shape := do
  let e := (← instantiateMVars e0).consumeMData
  -- mult-add a ◯ b ◯ ...: a box with its unit removed
  if (maOp? e).isSome then
    let (ls, cut) := flattenOp maOp? e 0
    return .box { sides := ← ls.mapM sideOf, cut }
  -- weighted family a ◯[k] b: a box whose far corner is dialed to weight k
  if constIs e ``CircleOp.circK then
    let as := lastArgs e 3
    let k := as[0]!
    let wOp : Expr → Option (Expr × Expr) := fun x =>
      if constIs x ``CircleOp.circK && (lastArgs x 3)[0]! == k then
        let bs := lastArgs x 3; some (bs[1]!, bs[2]!) else none
    let (ls, cut) := flattenOp wOp e 0
    return .box { sides := ← ls.mapM sideOf, cut, weight := some (toString (← ppExpr k)) }
  -- mult-add power a^(◯n): a cube of n equal sides, unit removed
  if constIs e ``CircleOp.mpow then
    let as := lastArgs e 2
    return .box { sides := #[← sideOf as[0]!], power := some (toString (← ppExpr as[1]!)) }
  if constIs e ``CircleOp.wpow then
    let as := lastArgs e 3
    return .box { sides := #[← sideOf as[1]!], power := some (toString (← ppExpr as[2]!)),
                  weight := some (toString (← ppExpr as[0]!)) }
  -- a fold over a list: a box with any number of sides
  if constIs e ``CircleOp.foldCirc then
    return .box { sides := #[← sideOf (lastArgs e 1)[0]!], any := true }
  -- the inverse and mult-add subtraction: a squash onto the hyperbola
  if constIs e ``CircleOp.circInv then
    return .squash (← read (lastArgs e 1)[0]!)
  if constIs e ``CircleOp.msub then
    let as := lastArgs e 2
    return .squash (.rel "⊟" (← read as[0]!) (← read as[1]!))
  if constIs e ``Inv.inv then
    return .squash (← read (lastArgs e 1)[0]!)
  if constIs e ``HDiv.hDiv && e.getAppNumArgs == 6 then
    return .squash (.rel "/" (← read (e.getArg! 4)) (← read (e.getArg! 5)))
  if constIs e ``Neg.neg then
    return .flip (← read (lastArgs e 1)[0]!)
  -- 1 + x: attach the unit
  if let some x := asLifted? e then
    let saved ← get
    match ← read x with
    | .box b => return .box { b with dropUnit := false }
    | _      => set saved; return .rod (← sideOf x) true
  -- a product: a box when every factor is lifted, else a single corner piece
  if (mulOp? e).isSome then
    let (fs, cut) := flattenOp mulOp? e 0
    if fs.all (fun f => (asLifted? f).isSome || liftedPow? f) then
      let mut sides : Array Side := #[]
      let mut powers : Array String := #[]
      let mut weight : Option String := none
      for f in fs do
        match ← readFactor f with
        | .box b =>
          for s in b.sides do
            -- equal sides that come as blocks of a power merge into one cube
            unless b.power.isSome && sides.contains s do sides := sides.push s
          if let some p := b.power then powers := powers.push p
          if b.weight.isSome then weight := b.weight
        | .rod s _ => sides := sides.push s
        | _ => pure ()
      let power := if powers.isEmpty then none else some (" + ".intercalate powers.toList)
      return .box { sides, cut, dropUnit := false, power, weight }
    return .piece (← fs.mapM sideOf)
  -- a power: a lifted base is a cube of n sides; otherwise a tower one rung up
  if constIs e ``HPow.hPow && e.getAppNumArgs == 6 then
    let b := e.getArg! 4; let n := e.getArg! 5
    if let some x := asLifted? b then
      return .box { sides := #[← sideOf x], power := some (toString (← ppExpr n)), dropUnit := false }
    return .tower (← read b) (toString (← ppExpr n))
  if (addOp? e).isSome then
    let (ps, _) := flattenOp addOp? e 0
    return .sum (← ps.mapM read)
  return .rod (← sideOf e) false

partial def liftedPow? (e : Expr) : Bool :=
  constIs e ``HPow.hPow && e.getAppNumArgs == 6 && (asLifted? (e.getArg! 4)).isSome

/-- A factor of a lifted product: (1 + x) gives one side; (1 + a ◯ b) gives all of the
    box's sides, so nested mult-adds merge into one box. -/
partial def readFactor (f : Expr) : ReadM Shape := do
  match asLifted? f with
  | some x =>
    let saved ← get
    match ← read x with
    | .box b => return .box b
    | _ => set saved; return .rod (← sideOf x) true
  | none => read f
end

/-- Read a proposition (the statement of a theorem, or a goal). -/
partial def readProp (p0 : Expr) : ReadM Shape := do
  let p := (← instantiateMVars p0).consumeMData
  if p.isConstOf ``False then return .note "no picture fits"
  if let some (_, l, r) := p.eq? then return .rel "=" (← read l) (← read r)
  if let some (l, r) := p.iff? then return .rel "↔" (← readProp l) (← readProp r)
  if let some (l, r) := p.and? then
    return .both #[← readProp l, ← readProp r]
  if let some q := p.not? then
    match ← readProp q with
    | .rel "=" l r => return .rel "≠" l r
    | s => return .rel "¬" s (.note "")
  if constIs p ``Ne then let as := lastArgs p 2; return .rel "≠" (← read as[0]!) (← read as[1]!)
  if constIs p ``LE.le then let as := lastArgs p 2; return .rel "≤" (← read as[0]!) (← read as[1]!)
  if constIs p ``LT.lt then let as := lastArgs p 2; return .rel "<" (← read as[0]!) (← read as[1]!)
  if constIs p ``GE.ge then let as := lastArgs p 2; return .rel "≥" (← read as[0]!) (← read as[1]!)
  if constIs p ``GT.gt then let as := lastArgs p 2; return .rel ">" (← read as[0]!) (← read as[1]!)
  if constIs p ``Dvd.dvd then let as := lastArgs p 2; return .rel "∣" (← read as[0]!) (← read as[1]!)
  if constIs p ``Exists then
    let body := (lastArgs p 1)[0]!
    match body with
    | .lam n ty b bi =>
      return ← withLocalDecl n bi ty fun x => do
        return .rel "∃" (.note n.toString) (← readProp (b.instantiate1 x))
    | _ => return .note "a piece is found that fits"
  if p.isForall then
    return ← forallTelescope p fun _ b => readProp b
  return .note (toString (← ppExpr p))

/-- Read a theorem's statement, under its binders. -/
def readStatement (ty : Expr) : MetaM Shape :=
  forallTelescope ty fun _ b => (readProp b).run' #[]

/-- Read the main goal. -/
def readGoal : TacticM Shape := do
  let g ← getMainGoal
  g.withContext do
    let t ← instantiateMVars (← g.getType)
    (readProp t).run' #[]

/-! ## 3. Widget: one inline module that draws a Shape in the infoview -/

@[widget_module]
def FusedBoxWidget : Widget.Module where
  javascript := "
import * as React from 'react';
const e = React.createElement;
const BASE = ['#1d5aa6','#d0312d','#e8b417','#3f8f4f','#6a3d9a','#e0701f'];
const MIX = {'':'#d9d6cf', a:'#1d5aa6', b:'#d0312d', c:'#e8b417', ab:'#6a3d9a', ac:'#3f8f4f',
             bc:'#e0701f', abc:'#2b2b33'};
function col(keys){ const k = keys.slice().sort().join(''); return MIX[k] || BASE[keys.length % BASE.length]; }
function Box({b}) {
  const n = b.sides.length, S = 34, pieces = [];
  if (n <= 3 && !b.power && !b.any) {
    // corners of the n-cube drawn as an oblique box: unit 18, side S
    const len = [18, S];
    for (let m = 0; m < (1 << n); m++) {
      const bits = [0,1,2].map(i => (m >> i) & 1);
      const keys = b.sides.filter((_, i) => bits[i]).map(s => s.key);
      const x = 10 + (bits[0] ? 18 : 0) + (n > 2 && bits[2] ? 10 : 0);
      const y = 10 + (n > 1 && bits[1] ? 18 : 0) - (n > 2 && bits[2] ? 10 : 0) + 10;
      const w = len[bits[0]], h = n > 1 ? len[bits[1]] : 10;
      const full = m === (1 << n) - 1;
      pieces.push(e('rect', {key: m, x, y, width: w - 2, height: h - 2, fill: col(keys),
        opacity: (m === 0 && b.dropUnit) ? 0.2 : 1, stroke: '#10131a',
        strokeDasharray: (m === 0 && b.dropUnit) ? '3 2' : 'none'}));
      if (full && b.weight) pieces.push(e('text', {key: 't', x: x + 4, y: y + 14, fill: '#fff', fontSize: 11}, b.weight));
    }
  } else {
    pieces.push(e('text', {key: 'p', x: 10, y: 30, fontSize: 12},
      (b.any ? 'any number of sides' : (b.power ? b.power + ' equal sides' : n + ' sides')) +
      (b.weight ? ', corner weight ' + b.weight : '') + (b.dropUnit ? ', unit removed' : '')));
  }
  return e('svg', {width: 110, height: 80}, pieces);
}
function Draw({s}) {
  if (!s) return null;
  switch (s.kind) {
    case 'box': return e(Box, {b: s.box});
    case 'rel': return e('div', {style: {display: 'flex', gap: '0.6em', alignItems: 'center'}},
                         e(Draw, {s: s.lhs}), e('b', null, s.rel), e(Draw, {s: s.rhs}));
    case 'both': return e('div', null, s.parts.map((p, i) => e(Draw, {s: p, key: i})));
    case 'rod': return e('code', {style: {borderBottom: '4px solid ' + col([s.side.key])}},
                         (s.lifted ? '1+' : '') + s.side.name);
    case 'piece': return e('code', null, s.sides.map(x => x.name).join('·'));
    case 'squash': return e('span', null, '1/(', e(Draw, {s: s.inner}), ')');
    case 'flip': return e('span', null, '−', e(Draw, {s: s.inner}));
    case 'tower': return e('span', null, e(Draw, {s: s.base}), e('sup', null, s.exp));
    case 'sum': return e('span', null, s.parts.map((p, i) => e('span', {key: i}, i ? ' + ' : '', e(Draw, {s: p}))));
    default: return e('i', null, s.text || '');
  }
}
export default function(props) {
  return e('div', null, e('div', {style: {fontSize: '0.8em', opacity: 0.7}}, 'move: ' + props.move),
           e(Draw, {s: props.shape}));
}
"

/-! ## 4. Moves: named rewrites that draw themselves -/

initialize registerTraceClass `fused.pictures

/-- After a move: draw the new goal in the infoview and write its shape to the trace. -/
def showGoal (stx : Syntax) (move : String) : TacticM Unit := do
  let gs ← getUnsolvedGoals
  let sh ← if gs.isEmpty then pure (Shape.note "proved") else readGoal
  let props := Json.mkObj [("shape", sh.toJson), ("move", move)]
  Widget.savePanelWidgetInfo FusedBoxWidget.javascriptHash (pure props) stx
  trace[fused.pictures] "{props.compress}"

/-- Lift: attach the unit cube to both sides (x = y becomes 1 + x = 1 + y). -/
theorem lift_back_ring {R : Type} [CommRing R] {x y : R} (h : 1 + x = 1 + y) : x = y := by
  have h2 : (-1) + (1 + x) = (-1) + (1 + y) := by rw [h]
  grind

theorem lift_back_nat {x y : Nat} (h : 1 + x = 1 + y) : x = y := by omega

theorem one_add_neg_one {R : Type} [CommRing R] : (1 : R) + -1 = 0 := by grind

/-- Run a list of rewrites: the first one that applies wins. -/
def tryRewrites (rws : List (TSyntax `tactic)) : TacticM Unit := do
  for r in rws do
    let ok ← (do evalTactic r; pure true) <|> pure false
    if ok then return
  throwError "no rewrite of this move applies to the goal"

/-- `lift`: attach the unit cube to both sides of an equation. -/
syntax (name := liftMove) "lift" : tactic
@[tactic liftMove] def evalLift : Tactic := fun stx => do
  tryRewrites [← `(tactic| refine lift_back_ring ?_), ← `(tactic| refine lift_back_nat ?_),
               ← `(tactic| apply CircleOp.lift_inj)]
  showGoal stx "lift"

/-- `multadd`: open every lifted mult-add into its box, (1 + a ◯ b) ↦ (1 + a)(1 + b).
    Uses the library's lift lemmas only. -/
syntax (name := multaddMove) "multadd" : tactic
@[tactic multaddMove] def evalMultadd : Tactic := fun stx => do
  evalTactic (← `(tactic| simp only [CircleOp.circ_lift, CircleOp.circ_lift',
      CircleOp.mpow_lift, CircleOp.circK_lift, CircleOp.wpow_lift, CircleOp.foldCirc_lift,
      CircleOp.fexp_lift, CircleOp.ftet_lift, CircleOp.sadd_lift]))
  showGoal stx "multadd"

/-- `swap`: reflect a box across its diagonal (mul_comm). Optional arguments pick the factors. -/
syntax (name := swapMove) "swap" (ppSpace colGt term:max)* : tactic
@[tactic swapMove] def evalSwap : Tactic := fun stx => do
  let args := stx[1].getArgs.map (⟨·⟩ : Syntax → Term)
  let t ← if args.isEmpty then `(CommSemiring.mul_comm) else `(CommSemiring.mul_comm $args*)
  evalTactic (← `(tactic| rw [$t:term]))
  showGoal stx "swap"

/-- `regroup`: the same box cut as columns instead of slabs (mul_assoc), or back. -/
syntax (name := regroupMove) "regroup" (ppSpace colGt term:max)* : tactic
@[tactic regroupMove] def evalRegroup : Tactic := fun stx => do
  let args := stx[1].getArgs.map (⟨·⟩ : Syntax → Term)
  let t ← if args.isEmpty then `(Semiring.mul_assoc) else `(Semiring.mul_assoc $args*)
  tryRewrites [← `(tactic| rw [$t:term]), ← `(tactic| rw [← $t:term])]
  showGoal stx "regroup"

/-- `collapse`: a side of length 0, or a factor of 1, flattens a dimension away. -/
syntax (name := collapseMove) "collapse" : tactic
@[tactic collapseMove] def evalCollapse : Tactic := fun stx => do
  evalTactic (← `(tactic| simp only [Semiring.mul_one, Semiring.one_mul, Semiring.add_zero,
      Semiring.zero_mul, Semiring.mul_zero, AddCommMonoid.zero_add, one_add_neg_one,
      CircleOp.circ_zero, CircleOp.zero_circ, CircleOp.circK_zero, CircleOp.circK_weight_zero,
      Nat.mul_one, Nat.one_mul, Nat.add_zero, Nat.zero_add]))
  showGoal stx "collapse"

/-- `squash`: rescale a rod to 1/(1+a), the way back (the inverse). Takes the hypothesis
    that 1 + a is not the wall. -/
syntax (name := squashMove) "squash" (ppSpace colGt term:max)* : tactic
@[tactic squashMove] def evalSquash : Tactic := fun stx => do
  let args := stx[1].getArgs.map (⟨·⟩ : Syntax → Term)
  let t ← `(CircleOp.inv_one_add $args*)
  tryRewrites [← `(tactic| rw [← $t:term]), ← `(tactic| rw [$t:term])]
  -- the squashed rectangle has area exactly 1
  evalTactic (← `(tactic| try rw [Lean.Grind.Field.mul_inv_cancel (by assumption)]))
  showGoal stx "squash"

/-- `dial`: weight the far corner. k·(a ◯[k] b) is the plain mult-add of ka and kb;
    weight 1 is mult-add, weight 0 is plain addition. -/
syntax (name := dialMove) "dial" : tactic
@[tactic dialMove] def evalDial : Tactic := fun stx => do
  tryRewrites [← `(tactic| rw [CircleOp.circK_scale]), ← `(tactic| rw [CircleOp.circK_one]),
               ← `(tactic| rw [CircleOp.circK_weight_zero]), ← `(tactic| rw [CircleOp.wpow_one])]
  showGoal stx "dial"

/-- `flatten`: area becomes length, a product of powers becomes a sum of exponents. -/
syntax (name := flattenMove) "flatten" : tactic
@[tactic flattenMove] def evalFlatten : Tactic := fun stx => do
  tryRewrites [← `(tactic| rw [CircleOp.pow_add']), ← `(tactic| rw [CircleOp.pw_add]),
               ← `(tactic| rw [CircleOp.mpow_add]), ← `(tactic| rw [Nat.pow_add])]
  showGoal stx "flatten"

/-- `extrude`: one more dimension. The next rung is the rung below iterated. -/
syntax (name := extrudeMove) "extrude" : tactic
@[tactic extrudeMove] def evalExtrude : Tactic := fun stx => do
  tryRewrites [← `(tactic| rw [CircleOp.fexp_iter]), ← `(tactic| rw [CircleOp.ftet_iter]),
               ← `(tactic| rw [CircleOp.fpent_iter]), ← `(tactic| rw [Semiring.pow_succ]),
               ← `(tactic| rw [Nat.pow_succ])]
  showGoal stx "extrude"

/-- `#shape name`: print the picture of a theorem's statement (as JSON). -/
elab "#shape " id:ident : command => Command.liftTermElabM do
  let n ← realizeGlobalConstNoOverloadWithInfo id
  let ci ← getConstInfo n
  let sh ← readStatement ci.type
  logInfo m!"{sh.toJson.compress}"

end FusedPictures
