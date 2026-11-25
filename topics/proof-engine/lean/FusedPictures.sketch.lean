-- SUPERSEDED: the compiled engine is FusedPictures.lean in this folder. Kept for history.
/-
  FusedPictures (SKETCH): a geometric proof engine for Fused Hyperoperation Algebra.

  STATUS: design sketch. NOT COMPILED. No Lean toolchain in the session that wrote it,
  and installing one needs Brennan's yes. Core Lean 4 only (no mathlib, no ProofWidgets).
  Lean's user-widget API is marked unstable, so names below may need small changes
  against the toolchain the library settles on (other files here were checked with
  Lean 4.37.0-rc1).

  Shape of the engine
  -------------------
  1. Shape     : a small datatype for pictures (rod, rectangle, box with weighted corners).
  2. Reader    : looks at the current goal and finds lifted products (x+1)*(y+1)*...
                 The picture is computed from the goal, never drawn by hand per theorem.
  3. Widget    : one inline JavaScript module that draws a Shape as SVG in the infoview.
  4. Moves     : named tactics (lift, swap, regroup, collapse, ...) that each run an
                 ordinary rewrite with a library lemma and then attach the picture.
  The kernel checks every rewrite. The pictures are a view of the goal and carry no trust.

  Dependency direction (kept one-way on purpose): the engine imports the shared library;
  the library never imports the engine. Moves are listed in a table on the engine side,
  so no library lemma needs a picture attribute.
-/
import Lean
open Lean Elab Tactic Meta

namespace FusedPictures

/-! ## 1. Shapes -/

/-- One side of a box: a name for the variable (as printed) and a colour key. -/
structure Side where
  name : String
  key  : String          -- "a", "b", "c", ... colour of the rod
  deriving ToJson, Inhabited

/-- A box (1+x₁)(1+x₂)…(1+xₙ). Its 2ⁿ pieces are the corners of the n-cube.
    `weights` gives a corner weight for a subset (as a bitmask), default 1.
    `cut` records how the box was built (which sides were grouped first), so the
    widget can draw slabs vs columns; it does not change the box. -/
structure Box where
  sides   : Array Side
  weights : Array (Nat × String) := #[]   -- (mask, weight as text, e.g. "A12")
  cut     : String := ""                  -- bracketing, e.g. "((0 1) 2)" = slabs (a⊕b)⊕c
  dropUnit : Bool := true                 -- draw the unit cube dashed (removed)
  deriving ToJson, Inhabited

inductive Shape where
  | rod  (parts : Array Side)
  | box  (b : Box)
  | eq   (lhs rhs : Shape)
  | note (text : String)
  deriving Inhabited

partial def Shape.toJson : Shape → Json
  | .rod ps  => Json.mkObj [("kind", "rod"), ("parts", toJson ps)]
  | .box b   => Json.mkObj [("kind", "box"), ("box", toJson b)]
  | .eq l r  => Json.mkObj [("kind", "eq"), ("lhs", l.toJson), ("rhs", r.toJson)]
  | .note t  => Json.mkObj [("kind", "note"), ("text", t)]

/-! ## 2. Reader: from a goal to a shape -/

/-- If `e` is `x + 1` (as the library's lift lemmas state it), return `x`. -/
def asLifted? (e : Expr) : Option Expr :=
  match e.getAppFnArgs with
  | (``HAdd.hAdd, #[_, _, _, _, x, one]) =>
      if one.nat? == some 1 then some x else none
  | _ => none

/-- Split a product into its factors, remembering how it was bracketed, e.g. "((0 1) 2)"
    for slabs (a⊕b)⊕c and "(0 (1 2))" for columns a⊕(b⊕c). -/
partial def factorsAux (e : Expr) (off : Nat) : Array Expr × String :=
  match e.getAppFnArgs with
  | (``HMul.hMul, #[_, _, _, _, l, r]) =>
      let (fl, sl) := factorsAux l off
      let (fr, sr) := factorsAux r (off + fl.size)
      (fl ++ fr, "(" ++ sl ++ " " ++ sr ++ ")")
  | _ => (#[e], toString off)

def factors (e : Expr) : Array Expr × String := factorsAux e 0

/-- Colour key for the i-th distinct variable: a, b, c, then d, e, … -/
def keyOf (i : Nat) : String := String.singleton (Char.ofNat ('a'.toNat + i))

/-- Read one side of an equation. A product of lifted rods is a box; anything else a rod. -/
def readSide (e : Expr) : MetaM Shape := do
  let (fs, cut) := factors e
  let lifted := fs.filterMap asLifted?
  if fs.size ≥ 2 && lifted.size == fs.size then
    let sides ← lifted.mapIdxM fun i x => do
      return { name := toString (← ppExpr x), key := keyOf i : Side }
    return .box { sides, cut }
  else
    return .rod #[{ name := toString (← ppExpr e), key := "?" }]

/-- Read the main goal `lhs = rhs`. -/
def readGoal : TacticM Shape := do
  let g ← getMainTarget
  match (← instantiateMVars g).eq? with
  | some (_, l, r) => return .eq (← readSide l) (← readSide r)
  | none => return .note "goal is not an equation"

/-! ## 3. Widget: one inline module that draws a Shape -/

/-- The JavaScript draws the same vocabulary as shapes.py (rods, boxes, corner weights,
    Byrne colours where a corner's colour mixes its sides' colours). Kept short here;
    the full renderer would be a straight port of shapes.py's draw_box. -/
@[widget_module]
def FusedBoxWidget : Widget.Module where
  javascript := "
    import * as React from 'react';
    const C = {a:'#1d5aa6', b:'#d0312d', c:'#e8b417', ab:'#6a3d9a', ac:'#3f8f4f',
               bc:'#e0701f', abc:'#2b2b33', '':'#d9d6cf'};
    function Box({box}) {
      const n = box.sides.length, S = 40, pieces = [];
      for (let m = 0; m < (1 << n); m++) {
        const key = box.sides.filter((_, i) => m >> i & 1).map(s => s.key).join('');
        pieces.push(React.createElement('rect', {key: m,
          x: 10 + (m & 1) * S, y: 10 + (m >> 1 & 1) * S, width: S - 4, height: S - 4,
          fill: C[key] || '#777', opacity: m === 0 && box.dropUnit ? 0.25 : 1}));
      }
      return React.createElement('svg', {width: 200, height: 120}, pieces);
    }
    function Draw({shape}) {
      if (shape.kind === 'box') return React.createElement(Box, {box: shape.box});
      if (shape.kind === 'eq') return React.createElement('div', {style:{display:'flex',gap:'1em'}},
        React.createElement(Draw, {shape: shape.lhs}), '=', React.createElement(Draw, {shape: shape.rhs}));
      if (shape.kind === 'rod') return React.createElement('code', null,
        shape.parts.map(p => p.name).join(' '));
      return React.createElement('i', null, shape.text || '');
    }
    export default function(props) { return React.createElement(Draw, {shape: props.shape}); }
  "

/-- Attach the current goal's picture to a piece of syntax (shown when the cursor is on it). -/
def showGoal (stx : Syntax) (move : String) : TacticM Unit := do
  let sh ← readGoal
  Widget.savePanelWidgetInfo FusedBoxWidget.javascriptHash
    (pure <| Json.mkObj [("shape", sh.toJson), ("move", move)]) stx

/-! ## 4. Moves: named rewrites that draw themselves

  Each move is an ordinary rewrite with a library lemma. The table is on the engine side.

    move       lemma (library or core)     what the picture does
    lift       lift_inj, madd_lift         attach the unit; a ⊕ b becomes a rectangle / box
    swap       Nat.mul_comm                reflect the box across a diagonal
    regroup    Nat.mul_assoc               same box, cut as slabs instead of columns
    collapse   Nat.mul_one, Nat.add_zero   a zero-length rod flattens a dimension
    squash     (inverse, needs Int/Rat)    rescale a rod to 1/(1+a); area exactly 1
    dial       weighted lift               corner pieces shrink to their weight
    flatten    Nat.pow_add / log           area becomes length (exponents add)
    extrude    fexp_lift, Nat.pow_succ     one more dimension: rung n+1 iterates rung n
-/

syntax (name := liftTac) "lift" : tactic
syntax (name := regroupTac) "regroup" : tactic
syntax (name := swapTac) "swap" : tactic

/- `lift` and the others name lemmas that live in the shared library. In this sketch the
   library's `lift_inj` and `madd_lift` are the ones in ns-rungs/lean/FusedRungs.lean. -/

@[tactic regroupTac] def evalRegroup : Tactic := fun stx => do
  evalTactic (← `(tactic| rw [Nat.mul_assoc]))
  showGoal stx "regroup"

@[tactic swapTac] def evalSwap : Tactic := fun stx => do
  evalTactic (← `(tactic| rw [Nat.mul_comm]))
  showGoal stx "swap"

end FusedPictures

/-! ## Worked example as it would read with the engine

  theorem madd_assoc' (a b c : Nat) : madd (madd a b) c = madd a (madd b c) := by
    apply lift_inj            -- attach the unit cube to both sides
    rw [madd_lift]            -- left: rod 1+(a⊕b) by rod 1+c, a rectangle
    rw [madd_lift]            -- left: the rod opens into (1+a)(1+b): a box of slabs
    rw [madd_lift, madd_lift] -- right: (1+a) by the (1+b)(1+c) face: a box of columns
    regroup                   -- Nat.mul_assoc: the slabs and the columns are one box

  Plain-Lean version of the same proof (each lemma it uses is checked in FusedRungs.lean;
  this exact chain of rewrites has not been compiled):

    apply lift_inj; rw [madd_lift, madd_lift, madd_lift, madd_lift, Nat.mul_assoc]
-/
