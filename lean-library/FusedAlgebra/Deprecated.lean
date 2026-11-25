/-
  FusedAlgebra.Deprecated: the pre-1.0 names, kept for one minor release
  =====================================================================

  The 1.0 rename (names-review/RENAME-PLAN.md) changed mult-add's Lean names to the
  ring-theory circle operation: multAdd → circ, wAdd → circK, namespace MultAdd →
  CircleOp, module Foundations/MultAddRing → Foundations/CircleRing, and so on. Every
  theorem and definition whose full name changed keeps its old name here as a
  deprecated alias of the new one (deprecated since 2026-11-25), so code written
  against v0.1.0 still compiles, with a warning that names the new constant.

  Core Lean has no `alias` command, so `deprecated_alias Old := New` below adds the
  old name with the same type and universe parameters, whose value is the new
  constant: a theorem for a theorem, a reducible definition otherwise. The aliases
  are added with `addDecl`, so they have no source range and Check.lean does not
  count them as theorems.

  Not kept: the old notation a ⊞ b (now a ◯ b, `open scoped CircleOp`) and the
  local ✳ infixes, which were never exported.
-/

import FusedAlgebra.Base
import FusedAlgebra.Rungs
import FusedAlgebra.Cube
import FusedAlgebra.Foundations.CircleRing
import FusedAlgebra.Foundations.Independence
import FusedAlgebra.Foundations.FusedPeano
import FusedAlgebra.Topics.ClosingRate
import FusedAlgebra.Topics.Consensus
import FusedAlgebra.Topics.DarkSector
import FusedAlgebra.Topics.DataCenter
import FusedAlgebra.Topics.E
import FusedAlgebra.Topics.Finance
import FusedAlgebra.Topics.FusedConstants
import FusedAlgebra.Topics.FractionalRungs
import FusedAlgebra.Topics.Infinities
import FusedAlgebra.Topics.InformationTheory
import FusedAlgebra.Topics.Landscape
import FusedAlgebra.Topics.RungShape
import FusedAlgebra.Topics.Elliptic
import FusedAlgebra.Topics.EllipticS3
import FusedAlgebra.Topics.Thermodynamics
import FusedAlgebra.Topics.Relativity
import FusedAlgebra.Topics.LeeYang
import FusedAlgebra.Topics.LeeYangCluster
import FusedAlgebra.Topics.LeeYangSpill
import FusedAlgebra.Topics.TwistedBoosts
import FusedAlgebra.Topics.GravityProbeB
import FusedAlgebra.Topics.BinaryPulsar
import FusedAlgebra.Topics.Chemistry.ChemQuantum
import FusedAlgebra.Topics.Chemistry.ChemConservation
import FusedAlgebra.Topics.Chemistry.ChemPeriodic
import FusedAlgebra.Topics.Chemistry.ChemBonding
import FusedAlgebra.Topics.Chemistry.ChemStoichiometry
import FusedAlgebra.Topics.Chemistry.ChemGas
import FusedAlgebra.Topics.Chemistry.ChemGibbs
import FusedAlgebra.Topics.Chemistry.ChemEquilibrium
import FusedAlgebra.Topics.Chemistry.ChemRates
import FusedAlgebra.Topics.Chemistry.ChemElectro
import FusedAlgebra.Topics.Chemistry.ChemAcids
import FusedAlgebra.Topics.Chemistry.ChemOrganic
import FusedAlgebra.Topics.Chemistry.ChemBiochem
import FusedAlgebra.Topics.ShapeDirections
import FusedAlgebra.Topics.LogicSets
import FusedAlgebra.Topics.CategoryTheory
import FusedAlgebra.Topics.Combinatorics
import FusedAlgebra.Topics.GroupsLie
import FusedAlgebra.Topics.Clifford
import FusedAlgebra.Topics.GraphTheory
import FusedAlgebra.Topics.MeasureProbability
import FusedAlgebra.Topics.ComplexAnalysis
import FusedAlgebra.Topics.FunctionalAnalysis
import FusedAlgebra.Topics.DiffGeometry
import FusedAlgebra.Topics.KnotTheory
import FusedAlgebra.Topics.Dynamics
import FusedAlgebra.Topics.InformationFull
import FusedAlgebra.Topics.Numerical
import FusedAlgebra.Topics.GameTheory
import FusedAlgebra.Topics.DiamondLattice
import FusedAlgebra.Topics.DialPointProcesses
import FusedAlgebra.Topics.SeveralWaves
import FusedAlgebra.Topics.RuleBreakers
import FusedAlgebra.Topics.Biology.MolecularBiology
import FusedAlgebra.Topics.Biology.Metabolism
import FusedAlgebra.Topics.Biology.CellBiology
import FusedAlgebra.Topics.Biology.Genetics
import FusedAlgebra.Topics.Biology.Development
import FusedAlgebra.Topics.Biology.Physiology
import FusedAlgebra.Topics.Biology.NeuroScience
import FusedAlgebra.Topics.Biology.Evolution
import FusedAlgebra.Topics.Biology.Population
import FusedAlgebra.Topics.Biology.Ecology
import FusedAlgebra.Topics.Biology.Behaviour
import FusedAlgebra.Topics.FractionalPowerset
import FusedAlgebra.Topics.NumberPlane
import FusedAlgebra.Topics.GeometryTrig
import FusedAlgebra.Topics.EarlyScience
import FusedAlgebra.Topics.Patterns
import FusedAlgebra.Topics.EllipticCounting
import FusedAlgebra.Topics.Pi
import FusedAlgebra.Topics.Calculus
import FusedAlgebra.Topics.AbstractAlgebra
import FusedAlgebra.Topics.Statistics
import FusedAlgebra.Topics.NSContinuum
import FusedAlgebra.Topics.JetFields
import FusedAlgebra.Topics.OsgoodInverse
import FusedAlgebra.Topics.LinearAlgebra
import FusedAlgebra.Topics.DialSweep
import FusedAlgebra.Topics.EntropyNoise
import FusedAlgebra.Topics.QAcross.QDial
import FusedAlgebra.Topics.QAcross.Knots
import FusedAlgebra.Topics.QAcross.Partitions
import FusedAlgebra.Topics.QAcross.QCounting
import FusedAlgebra.Topics.SolitonsRungs
import FusedAlgebra.Topics.NavierStokes
import FusedAlgebra.Topics.QuantumNoise
import FusedAlgebra.Topics.Rulecraft
import FusedAlgebra.Topics.Toric.ToricCharts
import FusedAlgebra.Topics.Toric.ProjectivePlane
import FusedAlgebra.Topics.Topology.CubeProofs
import FusedAlgebra.Topics.Topology.KummerWilson
import FusedAlgebra.Topics.Topology.HopfInvariantOne
import FusedAlgebra.Topics.Fermat.Fermat4
import FusedAlgebra.Topics.Fermat.Fermat3
import FusedAlgebra.Topics.Fermat.FreyPieces
import FusedAlgebra.Topics.Fermat.FermatRoute
import FusedAlgebra.Topics.Fermat.CounterexampleMap
import FusedAlgebra.Topics.GeneralRelativity
import FusedAlgebra.Topics.GRBridges
import FusedAlgebra.Topics.Electromagnetism
import FusedAlgebra.Topics.PvsNP
import FusedAlgebra.Topics.PowersetDials
import FusedAlgebra.Topics.DialSurvey
import FusedAlgebra.Topics.SizeDialResidual
import FusedAlgebra.Topics.DialWall
import FusedAlgebra.Topics.FractalBoxes
import FusedAlgebra.Topics.FractalDials
import FusedAlgebra.Topics.DialPotential
import FusedAlgebra.Topics.RepetitionDial
import FusedAlgebra.Topics.ZoomDial
import FusedAlgebra.Topics.Conventions
import FusedAlgebra.Topics.OrderDial
import Lean

open Lean Elab Command Meta

/-- `deprecated_alias Old.name := New.name` adds `Old.name` with the type and universe
    parameters of `New.name`, whose value is `New.name` itself (a theorem for a theorem,
    a reducible definition otherwise), and marks it deprecated in favour of `New.name`. -/
elab "deprecated_alias " old:ident " := " new:ident : command => liftTermElabM do
  let oldN := old.getId
  let newN := new.getId
  let some ci := (← getEnv).find? newN | throwError "deprecated_alias: unknown constant {newN}"
  if (← getEnv).contains oldN then throwError "deprecated_alias: {oldN} already exists"
  let val := mkConst newN (ci.levelParams.map mkLevelParam)
  match ci with
  | .thmInfo _ =>
    addDecl <| .thmDecl { name := oldN, levelParams := ci.levelParams, type := ci.type, value := val }
  | _ =>
    let safety := if ci.isUnsafe then DefinitionSafety.unsafe else DefinitionSafety.safe
    let decl := Declaration.defnDecl
      { name := oldN, levelParams := ci.levelParams, type := ci.type, value := val,
        hints := .abbrev, safety := safety }
    addDecl decl
    setReducibilityStatus oldN .reducible
    compileDecl decl (logErrors := false)
  match Linter.deprecatedAttr.setParam (← getEnv) oldN
      { newName? := some newN, since? := some "2026-11-25" } with
  | .ok env => setEnv env
  | .error e => throwError e


-- FusedAlgebra.Base
deprecated_alias MultAdd.foldMA := CircleOp.foldCirc
deprecated_alias MultAdd.foldMA_append := CircleOp.foldCirc_append
deprecated_alias MultAdd.foldMA_lift := CircleOp.foldCirc_lift
deprecated_alias MultAdd.foldMA_perm := CircleOp.foldCirc_perm
deprecated_alias MultAdd.groups := CircleOp.groups
deprecated_alias MultAdd.inv_one_add := CircleOp.inv_one_add
deprecated_alias MultAdd.mpow := CircleOp.mpow
deprecated_alias MultAdd.mpow_add := CircleOp.mpow_add
deprecated_alias MultAdd.mpow_eq := CircleOp.mpow_eq
deprecated_alias MultAdd.mpow_lift := CircleOp.mpow_lift
deprecated_alias MultAdd.msub := CircleOp.msub
deprecated_alias MultAdd.msub_multAdd := CircleOp.msub_circ
deprecated_alias MultAdd.mul_pow' := CircleOp.mul_pow'
deprecated_alias MultAdd.multAdd := CircleOp.circ
deprecated_alias MultAdd.multAddInv := CircleOp.circInv
deprecated_alias MultAdd.multAddInv_multAdd := CircleOp.circInv_circ
deprecated_alias MultAdd.multAdd_assoc := CircleOp.circ_assoc
deprecated_alias MultAdd.multAdd_comm := CircleOp.circ_comm
deprecated_alias MultAdd.multAdd_def := CircleOp.circ_def
deprecated_alias MultAdd.multAdd_excess := CircleOp.circ_excess
deprecated_alias MultAdd.multAdd_lift := CircleOp.circ_lift
deprecated_alias MultAdd.multAdd_lift' := CircleOp.circ_lift'
deprecated_alias MultAdd.multAdd_multAddInv := CircleOp.circ_circInv
deprecated_alias MultAdd.multAdd_zero := CircleOp.circ_zero
deprecated_alias MultAdd.neg_one_absorbs := CircleOp.neg_one_absorbs
deprecated_alias MultAdd.no_multAddInv_at_wall := CircleOp.no_circInv_at_wall
deprecated_alias MultAdd.one_add_foldMA_eq_corners := CircleOp.one_add_foldCirc_eq_corners
deprecated_alias MultAdd.one_pow' := CircleOp.one_pow'
deprecated_alias MultAdd.pow_add' := CircleOp.pow_add'
deprecated_alias MultAdd.prodL := CircleOp.prodL
deprecated_alias MultAdd.pw := CircleOp.pw
deprecated_alias MultAdd.pw_add := CircleOp.pw_add
deprecated_alias MultAdd.pw_eq_pow := CircleOp.pw_eq_pow
deprecated_alias MultAdd.pw_mul := CircleOp.pw_mul
deprecated_alias MultAdd.pw_mul_base := CircleOp.pw_mul_base
deprecated_alias MultAdd.pw_one := CircleOp.pw_one
deprecated_alias MultAdd.sumL := CircleOp.sumL
deprecated_alias MultAdd.sumL_append := CircleOp.sumL_append
deprecated_alias MultAdd.sumL_map_cons := CircleOp.sumL_map_cons
deprecated_alias MultAdd.wAdd := CircleOp.circK
deprecated_alias MultAdd.wAdd_assoc := CircleOp.circK_assoc
deprecated_alias MultAdd.wAdd_comm := CircleOp.circK_comm
deprecated_alias MultAdd.wAdd_def := CircleOp.circK_def
deprecated_alias MultAdd.wAdd_either := CircleOp.circK_either
deprecated_alias MultAdd.wAdd_flip := CircleOp.circK_flip
deprecated_alias MultAdd.wAdd_lift := CircleOp.circK_lift
deprecated_alias MultAdd.wAdd_neg_one_mirror := CircleOp.circK_neg_one_mirror
deprecated_alias MultAdd.wAdd_one := CircleOp.circK_one
deprecated_alias MultAdd.wAdd_scale := CircleOp.circK_scale
deprecated_alias MultAdd.wAdd_wall := CircleOp.circK_wall
deprecated_alias MultAdd.wAdd_weight_zero := CircleOp.circK_weight_zero
deprecated_alias MultAdd.wAdd_zero := CircleOp.circK_zero
deprecated_alias MultAdd.wpow := CircleOp.wpow
deprecated_alias MultAdd.wpow_lift := CircleOp.wpow_lift
deprecated_alias MultAdd.wpow_neg_one := CircleOp.wpow_neg_one
deprecated_alias MultAdd.wpow_one := CircleOp.wpow_one
deprecated_alias MultAdd.zero_multAdd := CircleOp.zero_circ

-- FusedAlgebra.Cube
deprecated_alias MultAdd.Cube.absorb := CircleOp.Cube.absorb
deprecated_alias MultAdd.Cube.fact := CircleOp.Cube.fact
deprecated_alias MultAdd.Cube.fact_pos := CircleOp.Cube.fact_pos
deprecated_alias MultAdd.Cube.layer := CircleOp.Cube.layer
deprecated_alias MultAdd.Cube.layer_eq_zero := CircleOp.Cube.layer_eq_zero
deprecated_alias MultAdd.Cube.layer_self := CircleOp.Cube.layer_self
deprecated_alias MultAdd.Cube.layer_step := CircleOp.Cube.layer_step
deprecated_alias MultAdd.Cube.layer_succ_succ := CircleOp.Cube.layer_succ_succ
deprecated_alias MultAdd.Cube.layer_zero_right := CircleOp.Cube.layer_zero_right

-- FusedAlgebra.Foundations.FusedPeano
deprecated_alias MAFound.FusedPeano := CircFound.FusedPeano
deprecated_alias MAFound.FusedPeano.add := CircFound.FusedPeano.add
deprecated_alias MAFound.FusedPeano.add_assoc' := CircFound.FusedPeano.add_assoc'
deprecated_alias MAFound.FusedPeano.add_comm' := CircFound.FusedPeano.add_comm'
deprecated_alias MAFound.FusedPeano.add_definable := CircFound.FusedPeano.add_definable
deprecated_alias MAFound.FusedPeano.add_right_cancel' := CircFound.FusedPeano.add_right_cancel'
deprecated_alias MAFound.FusedPeano.add_succ' := CircFound.FusedPeano.add_succ'
deprecated_alias MAFound.FusedPeano.add_toN_fromN := CircFound.FusedPeano.add_toN_fromN
deprecated_alias MAFound.FusedPeano.add_zero' := CircFound.FusedPeano.add_zero'
deprecated_alias MAFound.FusedPeano.casesOn := CircFound.FusedPeano.casesOn
deprecated_alias MAFound.FusedPeano.fromN := CircFound.FusedPeano.fromN
deprecated_alias MAFound.FusedPeano.fromN_toN := CircFound.FusedPeano.fromN_toN
deprecated_alias MAFound.FusedPeano.hyp := CircFound.FusedPeano.hyp
deprecated_alias MAFound.FusedPeano.hypBase := CircFound.FusedPeano.hypBase
deprecated_alias MAFound.FusedPeano.hyp_add := CircFound.FusedPeano.hyp_add
deprecated_alias MAFound.FusedPeano.hyp_mul := CircFound.FusedPeano.hyp_mul
deprecated_alias MAFound.FusedPeano.hyp_one := CircFound.FusedPeano.hyp_one
deprecated_alias MAFound.FusedPeano.hyp_pow := CircFound.FusedPeano.hyp_pow
deprecated_alias MAFound.FusedPeano.ind := CircFound.FusedPeano.ind
deprecated_alias MAFound.FusedPeano.instAdd := CircFound.FusedPeano.instAdd
deprecated_alias MAFound.FusedPeano.instCommSemiring := CircFound.FusedPeano.instCommSemiring
deprecated_alias MAFound.FusedPeano.instHPowNat := CircFound.FusedPeano.instHPowNat
deprecated_alias MAFound.FusedPeano.instMul := CircFound.FusedPeano.instMul
deprecated_alias MAFound.FusedPeano.instNatCast := CircFound.FusedPeano.instNatCast
deprecated_alias MAFound.FusedPeano.instOfNat := CircFound.FusedPeano.instOfNat
deprecated_alias MAFound.FusedPeano.instSMulNat := CircFound.FusedPeano.instSMulNat
deprecated_alias MAFound.FusedPeano.lad := CircFound.FusedPeano.lad
deprecated_alias MAFound.FusedPeano.ladBase := CircFound.FusedPeano.ladBase
deprecated_alias MAFound.FusedPeano.ladN := CircFound.FusedPeano.ladN
deprecated_alias MAFound.FusedPeano.ladN_one := CircFound.FusedPeano.ladN_one
deprecated_alias MAFound.FusedPeano.ladN_succ := CircFound.FusedPeano.ladN_succ
deprecated_alias MAFound.FusedPeano.ladN_two := CircFound.FusedPeano.ladN_two
deprecated_alias MAFound.FusedPeano.ladN_zero_right := CircFound.FusedPeano.ladN_zero_right
deprecated_alias MAFound.FusedPeano.lad_lift := CircFound.FusedPeano.lad_lift
deprecated_alias MAFound.FusedPeano.lad_one := CircFound.FusedPeano.lad_one
deprecated_alias MAFound.FusedPeano.lad_three := CircFound.FusedPeano.lad_three
deprecated_alias MAFound.FusedPeano.lad_two := CircFound.FusedPeano.lad_two
deprecated_alias MAFound.FusedPeano.ma := CircFound.FusedPeano.ma
deprecated_alias MAFound.FusedPeano.madd := CircFound.FusedPeano.circ
deprecated_alias MAFound.FusedPeano.madd_assoc := CircFound.FusedPeano.circ_assoc
deprecated_alias MAFound.FusedPeano.madd_comm := CircFound.FusedPeano.circ_comm
deprecated_alias MAFound.FusedPeano.madd_eq' := CircFound.FusedPeano.circ_eq'
deprecated_alias MAFound.FusedPeano.madd_is := CircFound.FusedPeano.circ_is
deprecated_alias MAFound.FusedPeano.madd_sadd := CircFound.FusedPeano.circ_sadd
deprecated_alias MAFound.FusedPeano.madd_succ := CircFound.FusedPeano.circ_succ
deprecated_alias MAFound.FusedPeano.madd_zero := CircFound.FusedPeano.circ_zero
deprecated_alias MAFound.FusedPeano.mk := CircFound.FusedPeano.mk
deprecated_alias MAFound.FusedPeano.mk.noConfusion := CircFound.FusedPeano.mk.noConfusion
deprecated_alias MAFound.FusedPeano.mul := CircFound.FusedPeano.mul
deprecated_alias MAFound.FusedPeano.mul_add' := CircFound.FusedPeano.mul_add'
deprecated_alias MAFound.FusedPeano.mul_assoc' := CircFound.FusedPeano.mul_assoc'
deprecated_alias MAFound.FusedPeano.mul_comm' := CircFound.FusedPeano.mul_comm'
deprecated_alias MAFound.FusedPeano.mul_definable := CircFound.FusedPeano.mul_definable
deprecated_alias MAFound.FusedPeano.mul_one' := CircFound.FusedPeano.mul_one'
deprecated_alias MAFound.FusedPeano.mul_succ' := CircFound.FusedPeano.mul_succ'
deprecated_alias MAFound.FusedPeano.mul_zero' := CircFound.FusedPeano.mul_zero'
deprecated_alias MAFound.FusedPeano.multAdd_eq := CircFound.FusedPeano.circ_eq
deprecated_alias MAFound.FusedPeano.noConfusion := CircFound.FusedPeano.noConfusion
deprecated_alias MAFound.FusedPeano.npow' := CircFound.FusedPeano.npow'
deprecated_alias MAFound.FusedPeano.rec := CircFound.FusedPeano.rec
deprecated_alias MAFound.FusedPeano.recN := CircFound.FusedPeano.recN
deprecated_alias MAFound.FusedPeano.recN_succ := CircFound.FusedPeano.recN_succ
deprecated_alias MAFound.FusedPeano.recN_zero := CircFound.FusedPeano.recN_zero
deprecated_alias MAFound.FusedPeano.recOn := CircFound.FusedPeano.recOn
deprecated_alias MAFound.FusedPeano.sadd := CircFound.FusedPeano.sadd
deprecated_alias MAFound.FusedPeano.sadd_assoc := CircFound.FusedPeano.sadd_assoc
deprecated_alias MAFound.FusedPeano.sadd_comm := CircFound.FusedPeano.sadd_comm
deprecated_alias MAFound.FusedPeano.sadd_is_succ := CircFound.FusedPeano.sadd_is_succ
deprecated_alias MAFound.FusedPeano.sadd_ne_zero := CircFound.FusedPeano.sadd_ne_zero
deprecated_alias MAFound.FusedPeano.sadd_right_cancel := CircFound.FusedPeano.sadd_right_cancel
deprecated_alias MAFound.FusedPeano.sadd_succ := CircFound.FusedPeano.sadd_succ
deprecated_alias MAFound.FusedPeano.sadd_zero := CircFound.FusedPeano.sadd_zero
deprecated_alias MAFound.FusedPeano.succ := CircFound.FusedPeano.succ
deprecated_alias MAFound.FusedPeano.succ_add' := CircFound.FusedPeano.succ_add'
deprecated_alias MAFound.FusedPeano.succ_add_eq := CircFound.FusedPeano.succ_add_eq
deprecated_alias MAFound.FusedPeano.succ_inj := CircFound.FusedPeano.succ_inj
deprecated_alias MAFound.FusedPeano.succ_madd := CircFound.FusedPeano.succ_circ
deprecated_alias MAFound.FusedPeano.succ_ne_zero := CircFound.FusedPeano.succ_ne_zero
deprecated_alias MAFound.FusedPeano.succ_sadd := CircFound.FusedPeano.succ_sadd
deprecated_alias MAFound.FusedPeano.toN := CircFound.FusedPeano.toN
deprecated_alias MAFound.FusedPeano.toN_add := CircFound.FusedPeano.toN_add
deprecated_alias MAFound.FusedPeano.toN_bijective := CircFound.FusedPeano.toN_bijective
deprecated_alias MAFound.FusedPeano.toN_fromN := CircFound.FusedPeano.toN_fromN
deprecated_alias MAFound.FusedPeano.toN_injective := CircFound.FusedPeano.toN_injective
deprecated_alias MAFound.FusedPeano.toN_lad := CircFound.FusedPeano.toN_lad
deprecated_alias MAFound.FusedPeano.toN_madd := CircFound.FusedPeano.toN_circ
deprecated_alias MAFound.FusedPeano.toN_mul := CircFound.FusedPeano.toN_mul
deprecated_alias MAFound.FusedPeano.toN_sadd := CircFound.FusedPeano.toN_sadd
deprecated_alias MAFound.FusedPeano.toN_surjective := CircFound.FusedPeano.toN_surjective
deprecated_alias MAFound.FusedPeano.transfer2 := CircFound.FusedPeano.transfer2
deprecated_alias MAFound.FusedPeano.transfer3 := CircFound.FusedPeano.transfer3
deprecated_alias MAFound.FusedPeano.zero := CircFound.FusedPeano.zero
deprecated_alias MAFound.FusedPeano.zero_add' := CircFound.FusedPeano.zero_add'
deprecated_alias MAFound.FusedPeano.zero_madd := CircFound.FusedPeano.zero_circ
deprecated_alias MAFound.FusedPeano.zero_sadd := CircFound.FusedPeano.zero_sadd

-- FusedAlgebra.Foundations.Independence
deprecated_alias MAFound.Independence.A1 := CircFound.Independence.A1
deprecated_alias MAFound.Independence.A2 := CircFound.Independence.A2
deprecated_alias MAFound.Independence.A3 := CircFound.Independence.A3
deprecated_alias MAFound.Independence.A4 := CircFound.Independence.A4
deprecated_alias MAFound.Independence.M1 := CircFound.Independence.M1
deprecated_alias MAFound.Independence.M2 := CircFound.Independence.M2
deprecated_alias MAFound.Independence.M3 := CircFound.Independence.M3
deprecated_alias MAFound.Independence.M4 := CircFound.Independence.M4
deprecated_alias MAFound.Independence.Mat := CircFound.Independence.Mat
deprecated_alias MAFound.Independence.Mat.a := CircFound.Independence.Mat.a
deprecated_alias MAFound.Independence.Mat.add := CircFound.Independence.Mat.add
deprecated_alias MAFound.Independence.Mat.b := CircFound.Independence.Mat.b
deprecated_alias MAFound.Independence.Mat.c := CircFound.Independence.Mat.c
deprecated_alias MAFound.Independence.Mat.casesOn := CircFound.Independence.Mat.casesOn
deprecated_alias MAFound.Independence.Mat.d := CircFound.Independence.Mat.d
deprecated_alias MAFound.Independence.Mat.ext' := CircFound.Independence.Mat.ext'
deprecated_alias MAFound.Independence.Mat.mk := CircFound.Independence.Mat.mk
deprecated_alias MAFound.Independence.Mat.mk.noConfusion := CircFound.Independence.Mat.mk.noConfusion
deprecated_alias MAFound.Independence.Mat.mul := CircFound.Independence.Mat.mul
deprecated_alias MAFound.Independence.Mat.neg := CircFound.Independence.Mat.neg
deprecated_alias MAFound.Independence.Mat.noConfusion := CircFound.Independence.Mat.noConfusion
deprecated_alias MAFound.Independence.Mat.rec := CircFound.Independence.Mat.rec
deprecated_alias MAFound.Independence.Mat.recOn := CircFound.Independence.Mat.recOn
deprecated_alias MAFound.Independence.Ops := CircFound.Independence.Ops
deprecated_alias MAFound.Independence.Ops.add := CircFound.Independence.Ops.add
deprecated_alias MAFound.Independence.Ops.casesOn := CircFound.Independence.Ops.casesOn
deprecated_alias MAFound.Independence.Ops.madd := CircFound.Independence.Ops.circ
deprecated_alias MAFound.Independence.Ops.mk := CircFound.Independence.Ops.mk
deprecated_alias MAFound.Independence.Ops.mk.noConfusion := CircFound.Independence.Ops.mk.noConfusion
deprecated_alias MAFound.Independence.Ops.neg := CircFound.Independence.Ops.neg
deprecated_alias MAFound.Independence.Ops.noConfusion := CircFound.Independence.Ops.noConfusion
deprecated_alias MAFound.Independence.Ops.rec := CircFound.Independence.Ops.rec
deprecated_alias MAFound.Independence.Ops.recOn := CircFound.Independence.Ops.recOn
deprecated_alias MAFound.Independence.Ops.wall := CircFound.Independence.Ops.wall
deprecated_alias MAFound.Independence.Ops.zero := CircFound.Independence.Ops.zero
deprecated_alias MAFound.Independence.Tri := CircFound.Independence.Tri
deprecated_alias MAFound.Independence.Tri.add := CircFound.Independence.Tri.add
deprecated_alias MAFound.Independence.Tri.casesOn := CircFound.Independence.Tri.casesOn
deprecated_alias MAFound.Independence.Tri.ext' := CircFound.Independence.Tri.ext'
deprecated_alias MAFound.Independence.Tri.mk := CircFound.Independence.Tri.mk
deprecated_alias MAFound.Independence.Tri.mk.noConfusion := CircFound.Independence.Tri.mk.noConfusion
deprecated_alias MAFound.Independence.Tri.mul := CircFound.Independence.Tri.mul
deprecated_alias MAFound.Independence.Tri.neg := CircFound.Independence.Tri.neg
deprecated_alias MAFound.Independence.Tri.noConfusion := CircFound.Independence.Tri.noConfusion
deprecated_alias MAFound.Independence.Tri.p := CircFound.Independence.Tri.p
deprecated_alias MAFound.Independence.Tri.q := CircFound.Independence.Tri.q
deprecated_alias MAFound.Independence.Tri.r := CircFound.Independence.Tri.r
deprecated_alias MAFound.Independence.Tri.rec := CircFound.Independence.Tri.rec
deprecated_alias MAFound.Independence.Tri.recOn := CircFound.Independence.Tri.recOn
deprecated_alias MAFound.Independence.instDecidableEqMat := CircFound.Independence.instDecidableEqMat
deprecated_alias MAFound.Independence.instDecidableEqMat.decEq := CircFound.Independence.instDecidableEqMat.decEq
deprecated_alias MAFound.Independence.instDecidableEqTri := CircFound.Independence.instDecidableEqTri
deprecated_alias MAFound.Independence.instDecidableEqTri.decEq := CircFound.Independence.instDecidableEqTri.decEq
deprecated_alias MAFound.Independence.intOps := CircFound.Independence.intOps
deprecated_alias MAFound.Independence.int_group := CircFound.Independence.int_group
deprecated_alias MAFound.Independence.matOps := CircFound.Independence.matOps
deprecated_alias MAFound.Independence.no_M1 := CircFound.Independence.no_M1
deprecated_alias MAFound.Independence.no_M2 := CircFound.Independence.no_M2
deprecated_alias MAFound.Independence.no_M3 := CircFound.Independence.no_M3
deprecated_alias MAFound.Independence.no_M4 := CircFound.Independence.no_M4
deprecated_alias MAFound.Independence.triOps := CircFound.Independence.triOps

-- FusedAlgebra.Foundations.CircleRing (was FusedAlgebra.Foundations.MultAddRing)
deprecated_alias MAFound.MAField := CircFound.CircField
deprecated_alias MAFound.MAField.add_wall_ne := CircFound.CircField.add_wall_ne
deprecated_alias MAFound.MAField.casesOn := CircFound.CircField.casesOn
deprecated_alias MAFound.MAField.instField := CircFound.CircField.instField
deprecated_alias MAFound.MAField.inv := CircFound.CircField.inv
deprecated_alias MAFound.MAField.inv_inv' := CircFound.CircField.inv_inv'
deprecated_alias MAFound.MAField.inv_ne_zero := CircFound.CircField.inv_ne_zero
deprecated_alias MAFound.MAField.inv_one' := CircFound.CircField.inv_one'
deprecated_alias MAFound.MAField.inv_zero' := CircFound.CircField.inv_zero'
deprecated_alias MAFound.MAField.mk := CircFound.CircField.mk
deprecated_alias MAFound.MAField.mk.noConfusion := CircFound.CircField.mk.noConfusion
deprecated_alias MAFound.MAField.mul_inv := CircFound.CircField.mul_inv
deprecated_alias MAFound.MAField.noConfusion := CircFound.CircField.noConfusion
deprecated_alias MAFound.MAField.rec := CircFound.CircField.rec
deprecated_alias MAFound.MAField.recOn := CircFound.CircField.recOn
deprecated_alias MAFound.MAField.toMARing := CircFound.CircField.toCircRing
deprecated_alias MAFound.MAField.way_back := CircFound.CircField.way_back
deprecated_alias MAFound.MAField.way_back_formula := CircFound.CircField.way_back_formula
deprecated_alias MAFound.MAField.zero_ne_one' := CircFound.CircField.zero_ne_one'
deprecated_alias MAFound.MAField.zero_ne_wall := CircFound.CircField.zero_ne_wall
deprecated_alias MAFound.MAField.zpow := CircFound.CircField.zpow
deprecated_alias MAFound.MARing := CircFound.CircRing
deprecated_alias MAFound.MARing.add := CircFound.CircRing.add
deprecated_alias MAFound.MARing.add_assoc := CircFound.CircRing.add_assoc
deprecated_alias MAFound.MARing.add_comm := CircFound.CircRing.add_comm
deprecated_alias MAFound.MARing.add_left_cancel := CircFound.CircRing.add_left_cancel
deprecated_alias MAFound.MARing.add_left_comm := CircFound.CircRing.add_left_comm
deprecated_alias MAFound.MARing.add_mul := CircFound.CircRing.add_mul
deprecated_alias MAFound.MARing.add_neg := CircFound.CircRing.add_neg
deprecated_alias MAFound.MARing.add_right_cancel := CircFound.CircRing.add_right_cancel
deprecated_alias MAFound.MARing.add_zero := CircFound.CircRing.add_zero
deprecated_alias MAFound.MARing.casesOn := CircFound.CircRing.casesOn
deprecated_alias MAFound.MARing.instAdd := CircFound.CircRing.instAdd
deprecated_alias MAFound.MARing.instAssociativeAdd := CircFound.CircRing.instAssociativeAdd
deprecated_alias MAFound.MARing.instCommRing := CircFound.CircRing.instCommRing
deprecated_alias MAFound.MARing.instCommutativeAdd := CircFound.CircRing.instCommutativeAdd
deprecated_alias MAFound.MARing.instHPowNat := CircFound.CircRing.instHPowNat
deprecated_alias MAFound.MARing.instIntCast := CircFound.CircRing.instIntCast
deprecated_alias MAFound.MARing.instMul := CircFound.CircRing.instMul
deprecated_alias MAFound.MARing.instNatCast := CircFound.CircRing.instNatCast
deprecated_alias MAFound.MARing.instNeg := CircFound.CircRing.instNeg
deprecated_alias MAFound.MARing.instOfNat := CircFound.CircRing.instOfNat
deprecated_alias MAFound.MARing.instSMulInt := CircFound.CircRing.instSMulInt
deprecated_alias MAFound.MARing.instSMulNat := CircFound.CircRing.instSMulNat
deprecated_alias MAFound.MARing.instSub := CircFound.CircRing.instSub
deprecated_alias MAFound.MARing.intCast := CircFound.CircRing.intCast
deprecated_alias MAFound.MARing.intCast_neg' := CircFound.CircRing.intCast_neg'
deprecated_alias MAFound.MARing.madd := CircFound.CircRing.circ
deprecated_alias MAFound.MARing.madd_add := CircFound.CircRing.circ_add
deprecated_alias MAFound.MARing.madd_assoc := CircFound.CircRing.circ_assoc
deprecated_alias MAFound.MARing.madd_comm := CircFound.CircRing.circ_comm
deprecated_alias MAFound.MARing.madd_eq := CircFound.CircRing.circ_eq
deprecated_alias MAFound.MARing.madd_wall := CircFound.CircRing.circ_wall
deprecated_alias MAFound.MARing.madd_zero := CircFound.CircRing.circ_zero
deprecated_alias MAFound.MARing.mk := CircFound.CircRing.mk
deprecated_alias MAFound.MARing.mk.noConfusion := CircFound.CircRing.mk.noConfusion
deprecated_alias MAFound.MARing.mul := CircFound.CircRing.mul
deprecated_alias MAFound.MARing.mul_add := CircFound.CircRing.mul_add
deprecated_alias MAFound.MARing.mul_assoc := CircFound.CircRing.mul_assoc
deprecated_alias MAFound.MARing.mul_comm := CircFound.CircRing.mul_comm
deprecated_alias MAFound.MARing.mul_neg := CircFound.CircRing.mul_neg
deprecated_alias MAFound.MARing.mul_one := CircFound.CircRing.mul_one
deprecated_alias MAFound.MARing.mul_wall := CircFound.CircRing.mul_wall
deprecated_alias MAFound.MARing.mul_zero := CircFound.CircRing.mul_zero
deprecated_alias MAFound.MARing.multAdd_eq := CircFound.CircRing.add_add_mul_eq_circ
deprecated_alias MAFound.MARing.natCast := CircFound.CircRing.natCast
deprecated_alias MAFound.MARing.natCast_succ := CircFound.CircRing.natCast_succ
deprecated_alias MAFound.MARing.neg := CircFound.CircRing.neg
deprecated_alias MAFound.MARing.neg_add := CircFound.CircRing.neg_add
deprecated_alias MAFound.MARing.neg_add_distrib := CircFound.CircRing.neg_add_distrib
deprecated_alias MAFound.MARing.neg_neg := CircFound.CircRing.neg_neg
deprecated_alias MAFound.MARing.neg_unique := CircFound.CircRing.neg_unique
deprecated_alias MAFound.MARing.neg_zero := CircFound.CircRing.neg_zero
deprecated_alias MAFound.MARing.noConfusion := CircFound.CircRing.noConfusion
deprecated_alias MAFound.MARing.npow := CircFound.CircRing.npow
deprecated_alias MAFound.MARing.one := CircFound.CircRing.one
deprecated_alias MAFound.MARing.one_eq := CircFound.CircRing.one_eq
deprecated_alias MAFound.MARing.rec := CircFound.CircRing.rec
deprecated_alias MAFound.MARing.recOn := CircFound.CircRing.recOn
deprecated_alias MAFound.MARing.wall := CircFound.CircRing.wall
deprecated_alias MAFound.MARing.wall_eq := CircFound.CircRing.wall_eq
deprecated_alias MAFound.MARing.wall_unique := CircFound.CircRing.wall_unique
deprecated_alias MAFound.MARing.zero := CircFound.CircRing.zero
deprecated_alias MAFound.MARing.zero_add := CircFound.CircRing.zero_add
deprecated_alias MAFound.MARing.zero_eq := CircFound.CircRing.zero_eq
deprecated_alias MAFound.MARing.zero_madd := CircFound.CircRing.zero_circ
deprecated_alias MAFound.MARing.zero_mul := CircFound.CircRing.zero_mul
deprecated_alias MAFound.ofCommRing := CircFound.ofCommRing
deprecated_alias MAFound.roundtrip_madd := CircFound.roundtrip_circ
deprecated_alias MAFound.roundtrip_mul := CircFound.roundtrip_mul
deprecated_alias MAFound.roundtrip_wall := CircFound.roundtrip_wall

-- FusedAlgebra.Rungs
deprecated_alias MultAdd.fexp := CircleOp.fexp
deprecated_alias MultAdd.fexp_base_iter := CircleOp.fexp_base_iter
deprecated_alias MultAdd.fexp_carries_multAdd := CircleOp.fexp_carries_circ
deprecated_alias MultAdd.fexp_example := CircleOp.fexp_example
deprecated_alias MultAdd.fexp_iter := CircleOp.fexp_iter
deprecated_alias MultAdd.fexp_lift := CircleOp.fexp_lift
deprecated_alias MultAdd.fexp_multAdd_right := CircleOp.fexp_circ_right
deprecated_alias MultAdd.fexp_not_assoc := CircleOp.fexp_not_assoc
deprecated_alias MultAdd.fexp_not_comm := CircleOp.fexp_not_comm
deprecated_alias MultAdd.fexp_sadd := CircleOp.fexp_sadd
deprecated_alias MultAdd.fexp_tower := CircleOp.fexp_tower
deprecated_alias MultAdd.fexp_zero_left := CircleOp.fexp_zero_left
deprecated_alias MultAdd.fexp_zero_right := CircleOp.fexp_zero_right
deprecated_alias MultAdd.fpent := CircleOp.fpent
deprecated_alias MultAdd.fpent_iter := CircleOp.fpent_iter
deprecated_alias MultAdd.fpent_lift := CircleOp.fpent_lift
deprecated_alias MultAdd.fpent_values := CircleOp.fpent_values
deprecated_alias MultAdd.fpent_zero_left := CircleOp.fpent_zero_left
deprecated_alias MultAdd.fpent_zero_right := CircleOp.fpent_zero_right
deprecated_alias MultAdd.ftet := CircleOp.ftet
deprecated_alias MultAdd.ftet_iter := CircleOp.ftet_iter
deprecated_alias MultAdd.ftet_lift := CircleOp.ftet_lift
deprecated_alias MultAdd.ftet_no_multAdd_dist := CircleOp.ftet_no_circ_dist
deprecated_alias MultAdd.ftet_no_sum := CircleOp.ftet_no_sum
deprecated_alias MultAdd.ftet_no_tower := CircleOp.ftet_no_tower
deprecated_alias MultAdd.ftet_not_assoc := CircleOp.ftet_not_assoc
deprecated_alias MultAdd.ftet_not_comm := CircleOp.ftet_not_comm
deprecated_alias MultAdd.ftet_values := CircleOp.ftet_values
deprecated_alias MultAdd.ftet_zero_left := CircleOp.ftet_zero_left
deprecated_alias MultAdd.ftet_zero_right := CircleOp.ftet_zero_right
deprecated_alias MultAdd.lift_inj := CircleOp.lift_inj
deprecated_alias MultAdd.multAdd_iter := CircleOp.circ_iter
deprecated_alias MultAdd.pent := CircleOp.pent
deprecated_alias MultAdd.pent_one_base := CircleOp.pent_one_base
deprecated_alias MultAdd.pent_pos := CircleOp.pent_pos
deprecated_alias MultAdd.pow_pos' := CircleOp.pow_pos'
deprecated_alias MultAdd.sadd := CircleOp.sadd
deprecated_alias MultAdd.sadd_lift := CircleOp.sadd_lift
deprecated_alias MultAdd.tet := CircleOp.tet
deprecated_alias MultAdd.tet_one := CircleOp.tet_one
deprecated_alias MultAdd.tet_one_base := CircleOp.tet_one_base
deprecated_alias MultAdd.tet_pos := CircleOp.tet_pos

-- FusedAlgebra.Topics.AbstractAlgebra
deprecated_alias AbstractAlgebra.int_multAdd_units := AbstractAlgebra.int_circ_units
deprecated_alias AbstractAlgebra.mpow_multAdd := AbstractAlgebra.mpow_circ
deprecated_alias AbstractAlgebra.multAdd_commMonoid := AbstractAlgebra.circ_commMonoid
deprecated_alias AbstractAlgebra.multAdd_fadd_distrib := AbstractAlgebra.circ_fadd_distrib
deprecated_alias AbstractAlgebra.multAdd_not_distrib_add := AbstractAlgebra.circ_not_distrib_add
deprecated_alias AbstractAlgebra.multAdd_one_wall_mod7 := AbstractAlgebra.circ_one_wall_mod7
deprecated_alias AbstractAlgebra.multAdd_unit_iff := AbstractAlgebra.circ_unit_iff
deprecated_alias AbstractAlgebra.tateRel_multAdd := AbstractAlgebra.tateRel_circ
deprecated_alias AbstractAlgebra.wAdd_commMonoid := AbstractAlgebra.circK_commMonoid
deprecated_alias AbstractAlgebra.wAdd_unit_iff := AbstractAlgebra.circK_unit_iff

-- FusedAlgebra.Topics.Biology.CellBiology
deprecated_alias CellBiology.ci_is_wAdd := CellBiology.ci_is_circK

-- FusedAlgebra.Topics.Biology.Metabolism
deprecated_alias Metabolism.bisub_is_wAdd := Metabolism.bisub_is_circK
deprecated_alias Metabolism.wAdd_comm' := Metabolism.circK_comm'

-- FusedAlgebra.Topics.Biology.MolecularBiology
deprecated_alias MolecularBiology.bind2_is_wAdd := MolecularBiology.bind2_is_circK

-- FusedAlgebra.Topics.Calculus
deprecated_alias MultAdd.Calculus.cubes_example := CircleOp.Calculus.cubes_example
deprecated_alias MultAdd.Calculus.dif := CircleOp.Calculus.dif
deprecated_alias MultAdd.Calculus.dif_one := CircleOp.Calculus.dif_one
deprecated_alias MultAdd.Calculus.dif_shift := CircleOp.Calculus.dif_shift
deprecated_alias MultAdd.Calculus.diff := CircleOp.Calculus.diff
deprecated_alias MultAdd.Calculus.diff_multAdd := CircleOp.Calculus.diff_circ
deprecated_alias MultAdd.Calculus.diff_pow := CircleOp.Calculus.diff_pow
deprecated_alias MultAdd.Calculus.ftc := CircleOp.Calculus.ftc
deprecated_alias MultAdd.Calculus.ftc_add := CircleOp.Calculus.ftc_add
deprecated_alias MultAdd.Calculus.ftc_box := CircleOp.Calculus.ftc_box
deprecated_alias MultAdd.Calculus.ftc_mslope := CircleOp.Calculus.ftc_mslope
deprecated_alias MultAdd.Calculus.ftc_step := CircleOp.Calculus.ftc_step
deprecated_alias MultAdd.Calculus.isum := CircleOp.Calculus.isum
deprecated_alias MultAdd.Calculus.isum_add := CircleOp.Calculus.isum_add
deprecated_alias MultAdd.Calculus.isum_congr := CircleOp.Calculus.isum_congr
deprecated_alias MultAdd.Calculus.jackson_gap := CircleOp.Calculus.jackson_gap
deprecated_alias MultAdd.Calculus.leibniz_corner := CircleOp.Calculus.leibniz_corner
deprecated_alias MultAdd.Calculus.liftProd := CircleOp.Calculus.liftProd
deprecated_alias MultAdd.Calculus.mint := CircleOp.Calculus.mint
deprecated_alias MultAdd.Calculus.mint_add := CircleOp.Calculus.mint_add
deprecated_alias MultAdd.Calculus.mint_const := CircleOp.Calculus.mint_const
deprecated_alias MultAdd.Calculus.mint_lift := CircleOp.Calculus.mint_lift
deprecated_alias MultAdd.Calculus.mint_multAdd := CircleOp.Calculus.mint_circ
deprecated_alias MultAdd.Calculus.mpow_mul := CircleOp.Calculus.mpow_mul
deprecated_alias MultAdd.Calculus.mslope := CircleOp.Calculus.mslope
deprecated_alias MultAdd.Calculus.mslope_step := CircleOp.Calculus.mslope_step
deprecated_alias MultAdd.Calculus.newton := CircleOp.Calculus.newton
deprecated_alias MultAdd.Calculus.pascal_sum := CircleOp.Calculus.pascal_sum
deprecated_alias MultAdd.Calculus.product_rule_step := CircleOp.Calculus.product_rule_step
deprecated_alias MultAdd.Calculus.rate_product := CircleOp.Calculus.rate_product
deprecated_alias MultAdd.Calculus.rate_product_limit := CircleOp.Calculus.rate_product_limit
deprecated_alias MultAdd.Calculus.step_iff := CircleOp.Calculus.step_iff
deprecated_alias MultAdd.Calculus.sumTo := CircleOp.Calculus.sumTo
deprecated_alias MultAdd.Calculus.ts_exp := CircleOp.Calculus.ts_exp
deprecated_alias MultAdd.Calculus.ts_exp_mul := CircleOp.Calculus.ts_exp_mul
deprecated_alias MultAdd.Calculus.ts_inverse := CircleOp.Calculus.ts_inverse

-- FusedAlgebra.Topics.Chemistry.ChemConservation
deprecated_alias ChemConservation.fugacity_multAdd := ChemConservation.fugacity_circ

-- FusedAlgebra.Topics.Chemistry.ChemGibbs
deprecated_alias ChemGibbs.independent_multAdd := ChemGibbs.independent_circ

-- FusedAlgebra.Topics.Chemistry.ChemQuantum
deprecated_alias ChemQuantum.fragments_multAdd := ChemQuantum.fragments_circ

-- FusedAlgebra.Topics.Clifford
deprecated_alias CliffordW.B4.lift_multAdd := CliffordW.B4.lift_circ
deprecated_alias CliffordW.dual_multAdd := CliffordW.dual_circ
deprecated_alias CliffordW.multAdd_dial_null := CliffordW.circ_dial_null
deprecated_alias CliffordW.multAdd_form_split := CliffordW.circ_form_split
deprecated_alias CliffordW.multAdd_lift_tconj := CliffordW.circ_lift_tconj

-- FusedAlgebra.Topics.Combinatorics
deprecated_alias Combinatorics.foldMA_replicate := Combinatorics.foldCirc_replicate

-- FusedAlgebra.Topics.ComplexAnalysis
deprecated_alias ComplexAnalysis.two_walls_multAdd := ComplexAnalysis.two_walls_circ

-- FusedAlgebra.Topics.DarkSector
deprecated_alias DarkSector.foldMA_ge_sum := DarkSector.foldCirc_ge_sum
deprecated_alias DarkSector.multAdd_ge_add := DarkSector.circ_ge_add
deprecated_alias DarkSector.multAdd_growth_is_exponential := DarkSector.circ_growth_is_exponential

-- FusedAlgebra.Topics.DataCenter
deprecated_alias DataCenterFused.foldMA_lift_nat := DataCenterFused.foldCirc_lift_nat
deprecated_alias DataCenterFused.foldMA_split := DataCenterFused.foldCirc_split
deprecated_alias DataCenterFused.foldMA_swap := DataCenterFused.foldCirc_swap

-- FusedAlgebra.Topics.DialPotential
deprecated_alias DialPotential.multAdd_mirror := DialPotential.circ_mirror

-- FusedAlgebra.Topics.DialSurvey
deprecated_alias MultAdd.DialSurvey.M1 := CircleOp.DialSurvey.M1
deprecated_alias MultAdd.DialSurvey.M2 := CircleOp.DialSurvey.M2
deprecated_alias MultAdd.DialSurvey.M3 := CircleOp.DialSurvey.M3
deprecated_alias MultAdd.DialSurvey.bose_exponents := CircleOp.DialSurvey.bose_exponents
deprecated_alias MultAdd.DialSurvey.cyclotomic_z2 := CircleOp.DialSurvey.cyclotomic_z2
deprecated_alias MultAdd.DialSurvey.cyclotomic_z3 := CircleOp.DialSurvey.cyclotomic_z3
deprecated_alias MultAdd.DialSurvey.exponent_ruler := CircleOp.DialSurvey.exponent_ruler
deprecated_alias MultAdd.DialSurvey.fall := CircleOp.DialSurvey.fall
deprecated_alias MultAdd.DialSurvey.fall_full := CircleOp.DialSurvey.fall_full
deprecated_alias MultAdd.DialSurvey.fall_vanish := CircleOp.DialSurvey.fall_vanish
deprecated_alias MultAdd.DialSurvey.fct := CircleOp.DialSurvey.fct
deprecated_alias MultAdd.DialSurvey.fct_succ := CircleOp.DialSurvey.fct_succ
deprecated_alias MultAdd.DialSurvey.fermi_exponents := CircleOp.DialSurvey.fermi_exponents
deprecated_alias MultAdd.DialSurvey.kRise := CircleOp.DialSurvey.kRise
deprecated_alias MultAdd.DialSurvey.kRise_bose := CircleOp.DialSurvey.kRise_bose
deprecated_alias MultAdd.DialSurvey.kRise_fermi := CircleOp.DialSurvey.kRise_fermi
deprecated_alias MultAdd.DialSurvey.kRise_parafermi := CircleOp.DialSurvey.kRise_parafermi
deprecated_alias MultAdd.DialSurvey.kRise_poisson := CircleOp.DialSurvey.kRise_poisson
deprecated_alias MultAdd.DialSurvey.kRise_succ := CircleOp.DialSurvey.kRise_succ
deprecated_alias MultAdd.DialSurvey.necklace_whole := CircleOp.DialSurvey.necklace_whole
deprecated_alias MultAdd.DialSurvey.odds_complement := CircleOp.DialSurvey.odds_complement
deprecated_alias MultAdd.DialSurvey.odds_union := CircleOp.DialSurvey.odds_union
deprecated_alias MultAdd.DialSurvey.units_change_weight := CircleOp.DialSurvey.units_change_weight
deprecated_alias MultAdd.DialSurvey.weight_to_one := CircleOp.DialSurvey.weight_to_one

-- FusedAlgebra.Topics.DialSweep
deprecated_alias DialSweep.kmadd := DialSweep.kcirc
deprecated_alias DialSweep.kmadd_lift := DialSweep.kcirc_lift

-- FusedAlgebra.Topics.DialWall
deprecated_alias MultAdd.DialWall.M2_zero_iff := CircleOp.DialWall.M2_zero_iff
deprecated_alias MultAdd.DialWall.M3_zero_iff := CircleOp.DialWall.M3_zero_iff
deprecated_alias MultAdd.DialWall.M4 := CircleOp.DialWall.M4
deprecated_alias MultAdd.DialWall.M5 := CircleOp.DialWall.M5
deprecated_alias MultAdd.DialWall.M5_zero_iff := CircleOp.DialWall.M5_zero_iff
deprecated_alias MultAdd.DialWall.M6 := CircleOp.DialWall.M6
deprecated_alias MultAdd.DialWall.bose_all_zero := CircleOp.DialWall.bose_all_zero
deprecated_alias MultAdd.DialWall.cancel_half := CircleOp.DialWall.cancel_half
deprecated_alias MultAdd.DialWall.cancel_half' := CircleOp.DialWall.cancel_half'
deprecated_alias MultAdd.DialWall.copies_len2 := CircleOp.DialWall.copies_len2
deprecated_alias MultAdd.DialWall.copies_len3 := CircleOp.DialWall.copies_len3
deprecated_alias MultAdd.DialWall.fermi_from_three := CircleOp.DialWall.fermi_from_three
deprecated_alias MultAdd.DialWall.half_exponent := CircleOp.DialWall.half_exponent
deprecated_alias MultAdd.DialWall.one_contributor := CircleOp.DialWall.one_contributor
deprecated_alias MultAdd.DialWall.unit_plus_multiple_ne_zero := CircleOp.DialWall.unit_plus_multiple_ne_zero
deprecated_alias MultAdd.DialWall.words_4 := CircleOp.DialWall.words_4
deprecated_alias MultAdd.DialWall.words_5 := CircleOp.DialWall.words_5
deprecated_alias MultAdd.DialWall.words_6 := CircleOp.DialWall.words_6

-- FusedAlgebra.Topics.DiffGeometry
deprecated_alias DiffGeometry.multAdd_curv := DiffGeometry.circ_curv

-- FusedAlgebra.Topics.Dynamics
deprecated_alias Dynamics.joukowski_multAdd := Dynamics.joukowski_circ

-- FusedAlgebra.Topics.E
deprecated_alias MultAdd.E.bernoulli := CircleOp.E.bernoulli
deprecated_alias MultAdd.E.bernoulli_strict := CircleOp.E.bernoulli_strict
deprecated_alias MultAdd.E.binomial := CircleOp.E.binomial
deprecated_alias MultAdd.E.collisions := CircleOp.E.collisions
deprecated_alias MultAdd.E.fact_cast_pos := CircleOp.E.fact_cast_pos
deprecated_alias MultAdd.E.fall := CircleOp.E.fall
deprecated_alias MultAdd.E.fall_le_pow := CircleOp.E.fall_le_pow
deprecated_alias MultAdd.E.inv_pow_mul := CircleOp.E.inv_pow_mul
deprecated_alias MultAdd.E.layer_fact_le := CircleOp.E.layer_fact_le
deprecated_alias MultAdd.E.layer_mul_fact := CircleOp.E.layer_mul_fact
deprecated_alias MultAdd.E.le_div_of_mul_le := CircleOp.E.le_div_of_mul_le
deprecated_alias MultAdd.E.natcast_ne := CircleOp.E.natcast_ne
deprecated_alias MultAdd.E.one_div_le_one_div := CircleOp.E.one_div_le_one_div
deprecated_alias MultAdd.E.one_div_nonneg := CircleOp.E.one_div_nonneg
deprecated_alias MultAdd.E.pow2_le_fact := CircleOp.E.pow2_le_fact
deprecated_alias MultAdd.E.rsum := CircleOp.E.rsum
deprecated_alias MultAdd.E.rsum_add := CircleOp.E.rsum_add
deprecated_alias MultAdd.E.rsum_le := CircleOp.E.rsum_le
deprecated_alias MultAdd.E.rsum_mono := CircleOp.E.rsum_mono
deprecated_alias MultAdd.E.rsum_mul := CircleOp.E.rsum_mul
deprecated_alias MultAdd.E.rsum_nonneg := CircleOp.E.rsum_nonneg
deprecated_alias MultAdd.E.rsum_shift := CircleOp.E.rsum_shift
deprecated_alias MultAdd.E.s := CircleOp.E.s
deprecated_alias MultAdd.E.s_bound := CircleOp.E.s_bound
deprecated_alias MultAdd.E.s_le_stair := CircleOp.E.s_le_stair
deprecated_alias MultAdd.E.s_le_u := CircleOp.E.s_le_u
deprecated_alias MultAdd.E.s_lt_three := CircleOp.E.s_lt_three
deprecated_alias MultAdd.E.s_mono := CircleOp.E.s_mono
deprecated_alias MultAdd.E.stair := CircleOp.E.stair
deprecated_alias MultAdd.E.stair_eq_cube := CircleOp.E.stair_eq_cube
deprecated_alias MultAdd.E.stair_eq_mpow := CircleOp.E.stair_eq_mpow
deprecated_alias MultAdd.E.stair_le_s := CircleOp.E.stair_le_s
deprecated_alias MultAdd.E.stair_le_u := CircleOp.E.stair_le_u
deprecated_alias MultAdd.E.stair_lt_succ := CircleOp.E.stair_lt_succ
deprecated_alias MultAdd.E.stair_nat_lt := CircleOp.E.stair_nat_lt
deprecated_alias MultAdd.E.tri := CircleOp.E.tri
deprecated_alias MultAdd.E.tri_mono := CircleOp.E.tri_mono
deprecated_alias MultAdd.E.u := CircleOp.E.u
deprecated_alias MultAdd.E.u_anti := CircleOp.E.u_anti
deprecated_alias MultAdd.E.u_step := CircleOp.E.u_step
deprecated_alias MultAdd.E.wt := CircleOp.E.wt
deprecated_alias MultAdd.E.wt_ge := CircleOp.E.wt_ge
deprecated_alias MultAdd.E.wt_le := CircleOp.E.wt_le
deprecated_alias MultAdd.E.wt_nonneg := CircleOp.E.wt_nonneg

-- FusedAlgebra.Topics.Electromagnetism
deprecated_alias Electromagnetism.doppler_multAdd := Electromagnetism.doppler_circ
deprecated_alias Electromagnetism.stokes_multAdd := Electromagnetism.stokes_circ

-- FusedAlgebra.Topics.Elliptic
deprecated_alias EllipticLead.multAdd_curve_d := EllipticLead.circ_curve_d
deprecated_alias EllipticLead.tate_q0_is_multAdd := EllipticLead.tate_q0_is_circ

-- FusedAlgebra.Topics.EllipticS3
deprecated_alias EllipticS3.multAdd_is_1x1 := EllipticS3.circ_is_1x1

-- FusedAlgebra.Topics.Fermat.Fermat3
deprecated_alias Fermat3.flt3_multAdd := Fermat3.flt3_circ
deprecated_alias Fermat3.flt3_multAdd_iff := Fermat3.flt3_circ_iff

-- FusedAlgebra.Topics.Fermat.FreyPieces
deprecated_alias FreyPieces.flt_multAdd := FreyPieces.flt_circ
deprecated_alias FreyPieces.node_chord_multAdd := FreyPieces.node_chord_circ

-- FusedAlgebra.Topics.Finance
deprecated_alias FinanceMultAdd.am_hm := FinanceCircle.am_hm
deprecated_alias FinanceMultAdd.average_factors := FinanceCircle.average_factors
deprecated_alias FinanceMultAdd.cheb := FinanceCircle.cheb
deprecated_alias FinanceMultAdd.drag := FinanceCircle.drag
deprecated_alias FinanceMultAdd.drag_le := FinanceCircle.drag_le
deprecated_alias FinanceMultAdd.fisher := FinanceCircle.fisher
deprecated_alias FinanceMultAdd.gordon := FinanceCircle.gordon
deprecated_alias FinanceMultAdd.growth := FinanceCircle.growth
deprecated_alias FinanceMultAdd.gsum := FinanceCircle.gsum
deprecated_alias FinanceMultAdd.gsum_le := FinanceCircle.gsum_le
deprecated_alias FinanceMultAdd.gsum_pole := FinanceCircle.gsum_pole
deprecated_alias FinanceMultAdd.lowest_rate_wins := FinanceCircle.lowest_rate_wins
deprecated_alias FinanceMultAdd.mall := FinanceCircle.mall
deprecated_alias FinanceMultAdd.mall_growth := FinanceCircle.mall_growth
deprecated_alias FinanceMultAdd.one_le_pow := FinanceCircle.one_le_pow
deprecated_alias FinanceMultAdd.pow_le_pow' := FinanceCircle.pow_le_pow'
deprecated_alias FinanceMultAdd.pow_nonneg' := FinanceCircle.pow_nonneg'
deprecated_alias FinanceMultAdd.power_mean := FinanceCircle.power_mean
deprecated_alias FinanceMultAdd.sq_nonneg' := FinanceCircle.sq_nonneg'

-- FusedAlgebra.Topics.FractionalPowerset
deprecated_alias FractionalPowerset.flat_madd_iff_flat_mul := FractionalPowerset.flat_circ_iff_flat_mul

-- FusedAlgebra.Topics.FunctionalAnalysis
deprecated_alias FunctionalAnalysis.M2.madd := FunctionalAnalysis.M2.circ
deprecated_alias FunctionalAnalysis.det_madd_swap := FunctionalAnalysis.det_circ_swap
deprecated_alias FunctionalAnalysis.opMA := FunctionalAnalysis.opCirc
deprecated_alias FunctionalAnalysis.opMA_comm_eq := FunctionalAnalysis.opCirc_comm_eq
deprecated_alias FunctionalAnalysis.tr_madd_swap := FunctionalAnalysis.tr_circ_swap

-- FusedAlgebra.Topics.GRBridges
deprecated_alias GRBridges.backward_part_multAdd := GRBridges.backward_part_circ
deprecated_alias GRBridges.forward_part_multAdd := GRBridges.forward_part_circ

-- FusedAlgebra.Topics.GameTheory
deprecated_alias GameTheory.odds_multAdd := GameTheory.odds_circ

-- FusedAlgebra.Topics.GroupsLie
deprecated_alias GroupsLie.T_eq_wAdd := GroupsLie.T_eq_circK
deprecated_alias GroupsLie.mob_multAdd := GroupsLie.mob_circ
deprecated_alias GroupsLie.theta_multAdd := GroupsLie.theta_circ

-- FusedAlgebra.Topics.Infinities
deprecated_alias Infinities.maddOnes := Infinities.circOnes
deprecated_alias Infinities.maddOnes_succ := Infinities.circOnes_succ

-- FusedAlgebra.Topics.InformationFull
deprecated_alias InformationFull.graft_multAdd := InformationFull.graft_circ

-- FusedAlgebra.Topics.InformationTheory
deprecated_alias InfoFused.union_is_flipped_madd := InfoFused.union_is_flipped_circ

-- FusedAlgebra.Topics.JetFields
deprecated_alias MultAdd.JetFields.Jet2 := CircleOp.JetFields.Jet2
deprecated_alias MultAdd.JetFields.Jet2.c0 := CircleOp.JetFields.Jet2.c0
deprecated_alias MultAdd.JetFields.Jet2.c1 := CircleOp.JetFields.Jet2.c1
deprecated_alias MultAdd.JetFields.Jet2.c2 := CircleOp.JetFields.Jet2.c2
deprecated_alias MultAdd.JetFields.Jet2.casesOn := CircleOp.JetFields.Jet2.casesOn
deprecated_alias MultAdd.JetFields.Jet2.mk := CircleOp.JetFields.Jet2.mk
deprecated_alias MultAdd.JetFields.Jet2.mk.noConfusion := CircleOp.JetFields.Jet2.mk.noConfusion
deprecated_alias MultAdd.JetFields.Jet2.noConfusion := CircleOp.JetFields.Jet2.noConfusion
deprecated_alias MultAdd.JetFields.Jet2.rec := CircleOp.JetFields.Jet2.rec
deprecated_alias MultAdd.JetFields.Jet2.recOn := CircleOp.JetFields.Jet2.recOn
deprecated_alias MultAdd.JetFields.M2 := CircleOp.JetFields.M2
deprecated_alias MultAdd.JetFields.M2.casesOn := CircleOp.JetFields.M2.casesOn
deprecated_alias MultAdd.JetFields.M2.det := CircleOp.JetFields.M2.det
deprecated_alias MultAdd.JetFields.M2.m11 := CircleOp.JetFields.M2.m11
deprecated_alias MultAdd.JetFields.M2.m12 := CircleOp.JetFields.M2.m12
deprecated_alias MultAdd.JetFields.M2.m21 := CircleOp.JetFields.M2.m21
deprecated_alias MultAdd.JetFields.M2.m22 := CircleOp.JetFields.M2.m22
deprecated_alias MultAdd.JetFields.M2.mk := CircleOp.JetFields.M2.mk
deprecated_alias MultAdd.JetFields.M2.mk.noConfusion := CircleOp.JetFields.M2.mk.noConfusion
deprecated_alias MultAdd.JetFields.M2.mul := CircleOp.JetFields.M2.mul
deprecated_alias MultAdd.JetFields.M2.noConfusion := CircleOp.JetFields.M2.noConfusion
deprecated_alias MultAdd.JetFields.M2.rec := CircleOp.JetFields.M2.rec
deprecated_alias MultAdd.JetFields.M2.recOn := CircleOp.JetFields.M2.recOn
deprecated_alias MultAdd.JetFields.M2.tr := CircleOp.JetFields.M2.tr
deprecated_alias MultAdd.JetFields.M3 := CircleOp.JetFields.M3
deprecated_alias MultAdd.JetFields.M3.addS := CircleOp.JetFields.M3.addS
deprecated_alias MultAdd.JetFields.M3.casesOn := CircleOp.JetFields.M3.casesOn
deprecated_alias MultAdd.JetFields.M3.det := CircleOp.JetFields.M3.det
deprecated_alias MultAdd.JetFields.M3.e2 := CircleOp.JetFields.M3.e2
deprecated_alias MultAdd.JetFields.M3.fsq := CircleOp.JetFields.M3.fsq
deprecated_alias MultAdd.JetFields.M3.m11 := CircleOp.JetFields.M3.m11
deprecated_alias MultAdd.JetFields.M3.m12 := CircleOp.JetFields.M3.m12
deprecated_alias MultAdd.JetFields.M3.m13 := CircleOp.JetFields.M3.m13
deprecated_alias MultAdd.JetFields.M3.m21 := CircleOp.JetFields.M3.m21
deprecated_alias MultAdd.JetFields.M3.m22 := CircleOp.JetFields.M3.m22
deprecated_alias MultAdd.JetFields.M3.m23 := CircleOp.JetFields.M3.m23
deprecated_alias MultAdd.JetFields.M3.m31 := CircleOp.JetFields.M3.m31
deprecated_alias MultAdd.JetFields.M3.m32 := CircleOp.JetFields.M3.m32
deprecated_alias MultAdd.JetFields.M3.m33 := CircleOp.JetFields.M3.m33
deprecated_alias MultAdd.JetFields.M3.mk := CircleOp.JetFields.M3.mk
deprecated_alias MultAdd.JetFields.M3.mk.noConfusion := CircleOp.JetFields.M3.mk.noConfusion
deprecated_alias MultAdd.JetFields.M3.mul := CircleOp.JetFields.M3.mul
deprecated_alias MultAdd.JetFields.M3.noConfusion := CircleOp.JetFields.M3.noConfusion
deprecated_alias MultAdd.JetFields.M3.rec := CircleOp.JetFields.M3.rec
deprecated_alias MultAdd.JetFields.M3.recOn := CircleOp.JetFields.M3.recOn
deprecated_alias MultAdd.JetFields.M3.step := CircleOp.JetFields.M3.step
deprecated_alias MultAdd.JetFields.M3.symm2 := CircleOp.JetFields.M3.symm2
deprecated_alias MultAdd.JetFields.M3.tr := CircleOp.JetFields.M3.tr
deprecated_alias MultAdd.JetFields.M3.vort := CircleOp.JetFields.M3.vort
deprecated_alias MultAdd.JetFields.WDual := CircleOp.JetFields.WDual
deprecated_alias MultAdd.JetFields.WDual.casesOn := CircleOp.JetFields.WDual.casesOn
deprecated_alias MultAdd.JetFields.WDual.du := CircleOp.JetFields.WDual.du
deprecated_alias MultAdd.JetFields.WDual.mk := CircleOp.JetFields.WDual.mk
deprecated_alias MultAdd.JetFields.WDual.mk.noConfusion := CircleOp.JetFields.WDual.mk.noConfusion
deprecated_alias MultAdd.JetFields.WDual.noConfusion := CircleOp.JetFields.WDual.noConfusion
deprecated_alias MultAdd.JetFields.WDual.re := CircleOp.JetFields.WDual.re
deprecated_alias MultAdd.JetFields.WDual.rec := CircleOp.JetFields.WDual.rec
deprecated_alias MultAdd.JetFields.WDual.recOn := CircleOp.JetFields.WDual.recOn
deprecated_alias MultAdd.JetFields.box_minor_zero := CircleOp.JetFields.box_minor_zero
deprecated_alias MultAdd.JetFields.cayley_hamilton := CircleOp.JetFields.cayley_hamilton
deprecated_alias MultAdd.JetFields.cube_depth_one := CircleOp.JetFields.cube_depth_one
deprecated_alias MultAdd.JetFields.cube_depth_two := CircleOp.JetFields.cube_depth_two
deprecated_alias MultAdd.JetFields.d1 := CircleOp.JetFields.d1
deprecated_alias MultAdd.JetFields.d2 := CircleOp.JetFields.d2
deprecated_alias MultAdd.JetFields.dDet := CircleOp.JetFields.dDet
deprecated_alias MultAdd.JetFields.dE2 := CircleOp.JetFields.dE2
deprecated_alias MultAdd.JetFields.det2_step := CircleOp.JetFields.det2_step
deprecated_alias MultAdd.JetFields.det9 := CircleOp.JetFields.det9
deprecated_alias MultAdd.JetFields.det_diag_cube := CircleOp.JetFields.det_diag_cube
deprecated_alias MultAdd.JetFields.det_expand := CircleOp.JetFields.det_expand
deprecated_alias MultAdd.JetFields.det_one_add_smul := CircleOp.JetFields.det_one_add_smul
deprecated_alias MultAdd.JetFields.det_vorticity_strain := CircleOp.JetFields.det_vorticity_strain
deprecated_alias MultAdd.JetFields.diag3 := CircleOp.JetFields.diag3
deprecated_alias MultAdd.JetFields.diag_layers := CircleOp.JetFields.diag_layers
deprecated_alias MultAdd.JetFields.diff_leibniz := CircleOp.JetFields.diff_leibniz
deprecated_alias MultAdd.JetFields.dual_inverse := CircleOp.JetFields.dual_inverse
deprecated_alias MultAdd.JetFields.dual_units_add := CircleOp.JetFields.dual_units_add
deprecated_alias MultAdd.JetFields.e2_expand := CircleOp.JetFields.e2_expand
deprecated_alias MultAdd.JetFields.e2_vorticity_strain := CircleOp.JetFields.e2_vorticity_strain
deprecated_alias MultAdd.JetFields.eps_sq := CircleOp.JetFields.eps_sq
deprecated_alias MultAdd.JetFields.jma := CircleOp.JetFields.jma
deprecated_alias MultAdd.JetFields.jma_lift := CircleOp.JetFields.jma_lift
deprecated_alias MultAdd.JetFields.jmul := CircleOp.JetFields.jmul
deprecated_alias MultAdd.JetFields.leibniz2 := CircleOp.JetFields.leibniz2
deprecated_alias MultAdd.JetFields.minor := CircleOp.JetFields.minor
deprecated_alias MultAdd.JetFields.minor_is_divergence := CircleOp.JetFields.minor_is_divergence
deprecated_alias MultAdd.JetFields.nilpotent_inverse := CircleOp.JetFields.nilpotent_inverse
deprecated_alias MultAdd.JetFields.ofTwoPt := CircleOp.JetFields.ofTwoPt
deprecated_alias MultAdd.JetFields.ofTwoPt_twoPt := CircleOp.JetFields.ofTwoPt_twoPt
deprecated_alias MultAdd.JetFields.off_wall_inverse := CircleOp.JetFields.off_wall_inverse
deprecated_alias MultAdd.JetFields.quad := CircleOp.JetFields.quad
deprecated_alias MultAdd.JetFields.re2_frozen := CircleOp.JetFields.re2_frozen
deprecated_alias MultAdd.JetFields.reStep := CircleOp.JetFields.reStep
deprecated_alias MultAdd.JetFields.re_conserved := CircleOp.JetFields.re_conserved
deprecated_alias MultAdd.JetFields.re_dQ := CircleOp.JetFields.re_dQ
deprecated_alias MultAdd.JetFields.re_dR := CircleOp.JetFields.re_dR
deprecated_alias MultAdd.JetFields.re_traceless := CircleOp.JetFields.re_traceless
deprecated_alias MultAdd.JetFields.shift_is_corner := CircleOp.JetFields.shift_is_corner
deprecated_alias MultAdd.JetFields.sq_nonneg_int := CircleOp.JetFields.sq_nonneg_int
deprecated_alias MultAdd.JetFields.sumTo := CircleOp.JetFields.sumTo
deprecated_alias MultAdd.JetFields.sumTo_congr := CircleOp.JetFields.sumTo_congr
deprecated_alias MultAdd.JetFields.sumTo_sub := CircleOp.JetFields.sumTo_sub
deprecated_alias MultAdd.JetFields.sumTo_telescope := CircleOp.JetFields.sumTo_telescope
deprecated_alias MultAdd.JetFields.sym_e2_nonpos := CircleOp.JetFields.sym_e2_nonpos
deprecated_alias MultAdd.JetFields.sym_e2_sumsq := CircleOp.JetFields.sym_e2_sumsq
deprecated_alias MultAdd.JetFields.sym_tr_cube := CircleOp.JetFields.sym_tr_cube
deprecated_alias MultAdd.JetFields.tr_cube_2d := CircleOp.JetFields.tr_cube_2d
deprecated_alias MultAdd.JetFields.tr_cube_traceless := CircleOp.JetFields.tr_cube_traceless
deprecated_alias MultAdd.JetFields.tr_sq_traceless := CircleOp.JetFields.tr_sq_traceless
deprecated_alias MultAdd.JetFields.twoPt := CircleOp.JetFields.twoPt
deprecated_alias MultAdd.JetFields.twoPt_mul := CircleOp.JetFields.twoPt_mul
deprecated_alias MultAdd.JetFields.twoPt_ofTwoPt := CircleOp.JetFields.twoPt_ofTwoPt
deprecated_alias MultAdd.JetFields.twoPt_weight_zero := CircleOp.JetFields.twoPt_weight_zero
deprecated_alias MultAdd.JetFields.units_multAdd := CircleOp.JetFields.units_circ
deprecated_alias MultAdd.JetFields.units_wAdd := CircleOp.JetFields.units_circK
deprecated_alias MultAdd.JetFields.vsq := CircleOp.JetFields.vsq
deprecated_alias MultAdd.JetFields.wall_no_inverse := CircleOp.JetFields.wall_no_inverse
deprecated_alias MultAdd.JetFields.wmul := CircleOp.JetFields.wmul
deprecated_alias MultAdd.JetFields.wmul_assoc := CircleOp.JetFields.wmul_assoc
deprecated_alias MultAdd.JetFields.wmul_comm := CircleOp.JetFields.wmul_comm

-- FusedAlgebra.Topics.KnotTheory
deprecated_alias KnotTheory.closeD_multAdd := KnotTheory.closeD_circ
deprecated_alias KnotTheory.poly_multAdd := KnotTheory.poly_circ

-- FusedAlgebra.Topics.LeeYangCluster
deprecated_alias LeeYangCluster.wAdd_units := LeeYangCluster.circK_units

-- FusedAlgebra.Topics.LinearAlgebra
deprecated_alias LinearAlgebraMA.M2 := LinearAlgebraCirc.M2
deprecated_alias LinearAlgebraMA.M2.Dk := LinearAlgebraCirc.M2.Dk
deprecated_alias LinearAlgebraMA.M2.Dk_hom := LinearAlgebraCirc.M2.Dk_hom
deprecated_alias LinearAlgebraMA.M2.a := LinearAlgebraCirc.M2.a
deprecated_alias LinearAlgebraMA.M2.add := LinearAlgebraCirc.M2.add
deprecated_alias LinearAlgebraMA.M2.add_one_eq_one := LinearAlgebraCirc.M2.add_one_eq_one
deprecated_alias LinearAlgebraMA.M2.b := LinearAlgebraCirc.M2.b
deprecated_alias LinearAlgebraMA.M2.c := LinearAlgebraCirc.M2.c
deprecated_alias LinearAlgebraMA.M2.casesOn := LinearAlgebraCirc.M2.casesOn
deprecated_alias LinearAlgebraMA.M2.cayley := LinearAlgebraCirc.M2.cayley
deprecated_alias LinearAlgebraMA.M2.cayley_eq := LinearAlgebraCirc.M2.cayley_eq
deprecated_alias LinearAlgebraMA.M2.cayley_is_msub := LinearAlgebraCirc.M2.cayley_is_msub
deprecated_alias LinearAlgebraMA.M2.cayley_on_circle := LinearAlgebraCirc.M2.cayley_on_circle
deprecated_alias LinearAlgebraMA.M2.cayley_on_wall_circle := LinearAlgebraCirc.M2.cayley_on_wall_circle
deprecated_alias LinearAlgebraMA.M2.coupling_sets_weight := LinearAlgebraCirc.M2.coupling_sets_weight
deprecated_alias LinearAlgebraMA.M2.d := LinearAlgebraCirc.M2.d
deprecated_alias LinearAlgebraMA.M2.det := LinearAlgebraCirc.M2.det
deprecated_alias LinearAlgebraMA.M2.det_mul := LinearAlgebraCirc.M2.det_mul
deprecated_alias LinearAlgebraMA.M2.det_one_add_skew := LinearAlgebraCirc.M2.det_one_add_skew
deprecated_alias LinearAlgebraMA.M2.det_one_add_smul := LinearAlgebraCirc.M2.det_one_add_smul
deprecated_alias LinearAlgebraMA.M2.diagonal_is_multAdd := LinearAlgebraCirc.M2.diagonal_is_circ
deprecated_alias LinearAlgebraMA.M2.ext' := LinearAlgebraCirc.M2.ext'
deprecated_alias LinearAlgebraMA.M2.inv := LinearAlgebraCirc.M2.inv
deprecated_alias LinearAlgebraMA.M2.lift := LinearAlgebraCirc.M2.lift
deprecated_alias LinearAlgebraMA.M2.madd := LinearAlgebraCirc.M2.circ
deprecated_alias LinearAlgebraMA.M2.madd_assoc := LinearAlgebraCirc.M2.circ_assoc
deprecated_alias LinearAlgebraMA.M2.madd_inverse := LinearAlgebraCirc.M2.circ_inverse
deprecated_alias LinearAlgebraMA.M2.madd_not_comm := LinearAlgebraCirc.M2.circ_not_comm
deprecated_alias LinearAlgebraMA.M2.madd_zero := LinearAlgebraCirc.M2.circ_zero
deprecated_alias LinearAlgebraMA.M2.mk := LinearAlgebraCirc.M2.mk
deprecated_alias LinearAlgebraMA.M2.mk.noConfusion := LinearAlgebraCirc.M2.mk.noConfusion
deprecated_alias LinearAlgebraMA.M2.msubM := LinearAlgebraCirc.M2.msubM
deprecated_alias LinearAlgebraMA.M2.msub_undo := LinearAlgebraCirc.M2.msub_undo
deprecated_alias LinearAlgebraMA.M2.msub_undo_W := LinearAlgebraCirc.M2.msub_undo_W
deprecated_alias LinearAlgebraMA.M2.mul := LinearAlgebraCirc.M2.mul
deprecated_alias LinearAlgebraMA.M2.mul_inv := LinearAlgebraCirc.M2.mul_inv
deprecated_alias LinearAlgebraMA.M2.noConfusion := LinearAlgebraCirc.M2.noConfusion
deprecated_alias LinearAlgebraMA.M2.one := LinearAlgebraCirc.M2.one
deprecated_alias LinearAlgebraMA.M2.orth_lift := LinearAlgebraCirc.M2.orth_lift
deprecated_alias LinearAlgebraMA.M2.rank_one_is_add := LinearAlgebraCirc.M2.rank_one_is_add
deprecated_alias LinearAlgebraMA.M2.rec := LinearAlgebraCirc.M2.rec
deprecated_alias LinearAlgebraMA.M2.recOn := LinearAlgebraCirc.M2.recOn
deprecated_alias LinearAlgebraMA.M2.rotation_condition := LinearAlgebraCirc.M2.rotation_condition
deprecated_alias LinearAlgebraMA.M2.rotation_condition_flat := LinearAlgebraCirc.M2.rotation_condition_flat
deprecated_alias LinearAlgebraMA.M2.shifted_det_hom := LinearAlgebraCirc.M2.shifted_det_hom
deprecated_alias LinearAlgebraMA.M2.skew := LinearAlgebraCirc.M2.skew
deprecated_alias LinearAlgebraMA.M2.smul := LinearAlgebraCirc.M2.smul
deprecated_alias LinearAlgebraMA.M2.sub := LinearAlgebraCirc.M2.sub
deprecated_alias LinearAlgebraMA.M2.tr := LinearAlgebraCirc.M2.tr
deprecated_alias LinearAlgebraMA.M2.tr_add := LinearAlgebraCirc.M2.tr_add
deprecated_alias LinearAlgebraMA.M2.transpose := LinearAlgebraCirc.M2.transpose
deprecated_alias LinearAlgebraMA.M2.wlift := LinearAlgebraCirc.M2.wlift
deprecated_alias LinearAlgebraMA.M2.wmadd := LinearAlgebraCirc.M2.circK
deprecated_alias LinearAlgebraMA.M2.wmadd_zero_weight := LinearAlgebraCirc.M2.circK_zero_weight
deprecated_alias LinearAlgebraMA.M2.zero := LinearAlgebraCirc.M2.zero
deprecated_alias LinearAlgebraMA.M2.zero_madd := LinearAlgebraCirc.M2.zero_circ
deprecated_alias LinearAlgebraMA.M3 := LinearAlgebraCirc.M3
deprecated_alias LinearAlgebraMA.M3.Dk := LinearAlgebraCirc.M3.Dk
deprecated_alias LinearAlgebraMA.M3.Dk_hom := LinearAlgebraCirc.M3.Dk_hom
deprecated_alias LinearAlgebraMA.M3.a11 := LinearAlgebraCirc.M3.a11
deprecated_alias LinearAlgebraMA.M3.a12 := LinearAlgebraCirc.M3.a12
deprecated_alias LinearAlgebraMA.M3.a13 := LinearAlgebraCirc.M3.a13
deprecated_alias LinearAlgebraMA.M3.a21 := LinearAlgebraCirc.M3.a21
deprecated_alias LinearAlgebraMA.M3.a22 := LinearAlgebraCirc.M3.a22
deprecated_alias LinearAlgebraMA.M3.a23 := LinearAlgebraCirc.M3.a23
deprecated_alias LinearAlgebraMA.M3.a31 := LinearAlgebraCirc.M3.a31
deprecated_alias LinearAlgebraMA.M3.a32 := LinearAlgebraCirc.M3.a32
deprecated_alias LinearAlgebraMA.M3.a33 := LinearAlgebraCirc.M3.a33
deprecated_alias LinearAlgebraMA.M3.add := LinearAlgebraCirc.M3.add
deprecated_alias LinearAlgebraMA.M3.box_cube := LinearAlgebraCirc.M3.box_cube
deprecated_alias LinearAlgebraMA.M3.casesOn := LinearAlgebraCirc.M3.casesOn
deprecated_alias LinearAlgebraMA.M3.det := LinearAlgebraCirc.M3.det
deprecated_alias LinearAlgebraMA.M3.det_mul := LinearAlgebraCirc.M3.det_mul
deprecated_alias LinearAlgebraMA.M3.det_one_add := LinearAlgebraCirc.M3.det_one_add
deprecated_alias LinearAlgebraMA.M3.det_one_add_smul := LinearAlgebraCirc.M3.det_one_add_smul
deprecated_alias LinearAlgebraMA.M3.e2 := LinearAlgebraCirc.M3.e2
deprecated_alias LinearAlgebraMA.M3.ext' := LinearAlgebraCirc.M3.ext'
deprecated_alias LinearAlgebraMA.M3.lift := LinearAlgebraCirc.M3.lift
deprecated_alias LinearAlgebraMA.M3.madd := LinearAlgebraCirc.M3.circ
deprecated_alias LinearAlgebraMA.M3.madd_assoc := LinearAlgebraCirc.M3.circ_assoc
deprecated_alias LinearAlgebraMA.M3.madd_zero := LinearAlgebraCirc.M3.circ_zero
deprecated_alias LinearAlgebraMA.M3.mk := LinearAlgebraCirc.M3.mk
deprecated_alias LinearAlgebraMA.M3.mk.noConfusion := LinearAlgebraCirc.M3.mk.noConfusion
deprecated_alias LinearAlgebraMA.M3.mul := LinearAlgebraCirc.M3.mul
deprecated_alias LinearAlgebraMA.M3.noConfusion := LinearAlgebraCirc.M3.noConfusion
deprecated_alias LinearAlgebraMA.M3.one := LinearAlgebraCirc.M3.one
deprecated_alias LinearAlgebraMA.M3.orth_lift := LinearAlgebraCirc.M3.orth_lift
deprecated_alias LinearAlgebraMA.M3.rec := LinearAlgebraCirc.M3.rec
deprecated_alias LinearAlgebraMA.M3.recOn := LinearAlgebraCirc.M3.recOn
deprecated_alias LinearAlgebraMA.M3.smul := LinearAlgebraCirc.M3.smul
deprecated_alias LinearAlgebraMA.M3.tr := LinearAlgebraCirc.M3.tr
deprecated_alias LinearAlgebraMA.M3.tr_add := LinearAlgebraCirc.M3.tr_add
deprecated_alias LinearAlgebraMA.M3.transpose := LinearAlgebraCirc.M3.transpose
deprecated_alias LinearAlgebraMA.M3.triangular_eigen := LinearAlgebraCirc.M3.triangular_eigen
deprecated_alias LinearAlgebraMA.M3.wlift := LinearAlgebraCirc.M3.wlift
deprecated_alias LinearAlgebraMA.M3.wmadd := LinearAlgebraCirc.M3.circK
deprecated_alias LinearAlgebraMA.M3.zero := LinearAlgebraCirc.M3.zero
deprecated_alias LinearAlgebraMA.fl1 := LinearAlgebraCirc.fl1
deprecated_alias LinearAlgebraMA.fl2 := LinearAlgebraCirc.fl2
deprecated_alias LinearAlgebraMA.fl_fixes_wall := LinearAlgebraCirc.fl_fixes_wall
deprecated_alias LinearAlgebraMA.fl_sadd := LinearAlgebraCirc.fl_sadd
deprecated_alias LinearAlgebraMA.fl_scale := LinearAlgebraCirc.fl_scale
deprecated_alias LinearAlgebraMA.kdv_two_soliton := LinearAlgebraCirc.kdv_two_soliton
deprecated_alias LinearAlgebraMA.mono := LinearAlgebraCirc.mono
deprecated_alias LinearAlgebraMA.mono_compose := LinearAlgebraCirc.mono_compose
deprecated_alias LinearAlgebraMA.mono_hom := LinearAlgebraCirc.mono_hom
deprecated_alias LinearAlgebraMA.multAdd_is_zoom := LinearAlgebraCirc.circ_is_zoom
deprecated_alias LinearAlgebraMA.pow_mul'' := LinearAlgebraCirc.pow_mul''
deprecated_alias LinearAlgebraMA.rotation_circle := LinearAlgebraCirc.rotation_circle
deprecated_alias LinearAlgebraMA.sadd := LinearAlgebraCirc.sadd

-- FusedAlgebra.Topics.LogicSets
deprecated_alias LogicSets.multAdd_definable := LogicSets.circ_definable

-- FusedAlgebra.Topics.NSContinuum
deprecated_alias MultAdd.NSContinuum.annulus_cells := CircleOp.NSContinuum.annulus_cells
deprecated_alias MultAdd.NSContinuum.annulus_graded := CircleOp.NSContinuum.annulus_graded
deprecated_alias MultAdd.NSContinuum.box_doubles := CircleOp.NSContinuum.box_doubles
deprecated_alias MultAdd.NSContinuum.energy_telescopes := CircleOp.NSContinuum.energy_telescopes
deprecated_alias MultAdd.NSContinuum.fan_pivot := CircleOp.NSContinuum.fan_pivot
deprecated_alias MultAdd.NSContinuum.flux := CircleOp.NSContinuum.flux
deprecated_alias MultAdd.NSContinuum.gain := CircleOp.NSContinuum.gain
deprecated_alias MultAdd.NSContinuum.gain_is_model := CircleOp.NSContinuum.gain_is_model
deprecated_alias MultAdd.NSContinuum.k41Margin := CircleOp.NSContinuum.k41Margin
deprecated_alias MultAdd.NSContinuum.k41_3d_positive := CircleOp.NSContinuum.k41_3d_positive
deprecated_alias MultAdd.NSContinuum.k41_5d_negative := CircleOp.NSContinuum.k41_5d_negative
deprecated_alias MultAdd.NSContinuum.k41_edge_at_four := CircleOp.NSContinuum.k41_edge_at_four
deprecated_alias MultAdd.NSContinuum.lions_on_2d_ray := CircleOp.NSContinuum.lions_on_2d_ray
deprecated_alias MultAdd.NSContinuum.mpow_one_nat := CircleOp.NSContinuum.mpow_one_nat
deprecated_alias MultAdd.NSContinuum.ns3d_on_bmr_line := CircleOp.NSContinuum.ns3d_on_bmr_line
deprecated_alias MultAdd.NSContinuum.ns3d_on_ray_of_2d_four_fifths := CircleOp.NSContinuum.ns3d_on_ray_of_2d_four_fifths
deprecated_alias MultAdd.NSContinuum.ns3d_supercritical := CircleOp.NSContinuum.ns3d_supercritical
deprecated_alias MultAdd.NSContinuum.onRay := CircleOp.NSContinuum.onRay
deprecated_alias MultAdd.NSContinuum.onRay_refl := CircleOp.NSContinuum.onRay_refl
deprecated_alias MultAdd.NSContinuum.onRay_symm := CircleOp.NSContinuum.onRay_symm
deprecated_alias MultAdd.NSContinuum.onRay_trans := CircleOp.NSContinuum.onRay_trans
deprecated_alias MultAdd.NSContinuum.planar_stretching_zero := CircleOp.NSContinuum.planar_stretching_zero
deprecated_alias MultAdd.NSContinuum.ssum := CircleOp.NSContinuum.ssum
deprecated_alias MultAdd.NSContinuum.stretch := CircleOp.NSContinuum.stretch

-- FusedAlgebra.Topics.Numerical
deprecated_alias Numerical.hornerMA := Numerical.hornerCirc
deprecated_alias Numerical.hornerMA_eq := Numerical.hornerCirc_eq
deprecated_alias Numerical.multAdd_bound := Numerical.circ_bound
deprecated_alias Numerical.multAdd_sensitivity := Numerical.circ_sensitivity

-- FusedAlgebra.Topics.OrderDial
deprecated_alias OrderDial.one_multAdd_one := OrderDial.one_circ_one
deprecated_alias OrderDial.wAdd_idem_iff := OrderDial.circK_idem_iff

-- FusedAlgebra.Topics.OsgoodInverse
deprecated_alias MultAdd.OsgoodInverse.Antitone1 := CircleOp.OsgoodInverse.Antitone1
deprecated_alias MultAdd.OsgoodInverse.blk := CircleOp.OsgoodInverse.blk
deprecated_alias MultAdd.OsgoodInverse.blk_le := CircleOp.OsgoodInverse.blk_le
deprecated_alias MultAdd.OsgoodInverse.block_start := CircleOp.OsgoodInverse.block_start
deprecated_alias MultAdd.OsgoodInverse.condensation := CircleOp.OsgoodInverse.condensation
deprecated_alias MultAdd.OsgoodInverse.const_block := CircleOp.OsgoodInverse.const_block
deprecated_alias MultAdd.OsgoodInverse.csum := CircleOp.OsgoodInverse.csum
deprecated_alias MultAdd.OsgoodInverse.csum_harm := CircleOp.OsgoodInverse.csum_harm
deprecated_alias MultAdd.OsgoodInverse.csum_isq := CircleOp.OsgoodInverse.csum_isq
deprecated_alias MultAdd.OsgoodInverse.csum_le_dsum := CircleOp.OsgoodInverse.csum_le_dsum
deprecated_alias MultAdd.OsgoodInverse.dsum := CircleOp.OsgoodInverse.dsum
deprecated_alias MultAdd.OsgoodInverse.dsum_le_csum := CircleOp.OsgoodInverse.dsum_le_csum
deprecated_alias MultAdd.OsgoodInverse.dsum_succ := CircleOp.OsgoodInverse.dsum_succ
deprecated_alias MultAdd.OsgoodInverse.fan := CircleOp.OsgoodInverse.fan
deprecated_alias MultAdd.OsgoodInverse.fan_distinct := CircleOp.OsgoodInverse.fan_distinct
deprecated_alias MultAdd.OsgoodInverse.fan_step := CircleOp.OsgoodInverse.fan_step
deprecated_alias MultAdd.OsgoodInverse.fan_zero := CircleOp.OsgoodInverse.fan_zero
deprecated_alias MultAdd.OsgoodInverse.fsum := CircleOp.OsgoodInverse.fsum
deprecated_alias MultAdd.OsgoodInverse.fsum_append := CircleOp.OsgoodInverse.fsum_append
deprecated_alias MultAdd.OsgoodInverse.fsum_le_dsum := CircleOp.OsgoodInverse.fsum_le_dsum
deprecated_alias MultAdd.OsgoodInverse.fsum_mono := CircleOp.OsgoodInverse.fsum_mono
deprecated_alias MultAdd.OsgoodInverse.fsum_nonneg := CircleOp.OsgoodInverse.fsum_nonneg
deprecated_alias MultAdd.OsgoodInverse.harm := CircleOp.OsgoodInverse.harm
deprecated_alias MultAdd.OsgoodInverse.harm_anti := CircleOp.OsgoodInverse.harm_anti
deprecated_alias MultAdd.OsgoodInverse.harm_nonneg := CircleOp.OsgoodInverse.harm_nonneg
deprecated_alias MultAdd.OsgoodInverse.harmonic_grows := CircleOp.OsgoodInverse.harmonic_grows
deprecated_alias MultAdd.OsgoodInverse.inverse_squares_bounded := CircleOp.OsgoodInverse.inverse_squares_bounded
deprecated_alias MultAdd.OsgoodInverse.isq := CircleOp.OsgoodInverse.isq
deprecated_alias MultAdd.OsgoodInverse.isq_anti := CircleOp.OsgoodInverse.isq_anti
deprecated_alias MultAdd.OsgoodInverse.le_blk := CircleOp.OsgoodInverse.le_blk
deprecated_alias MultAdd.OsgoodInverse.lipschitz_unique := CircleOp.OsgoodInverse.lipschitz_unique
deprecated_alias MultAdd.OsgoodInverse.lt_pow2 := CircleOp.OsgoodInverse.lt_pow2
deprecated_alias MultAdd.OsgoodInverse.mirror_lift := CircleOp.OsgoodInverse.mirror_lift
deprecated_alias MultAdd.OsgoodInverse.mirror_step := CircleOp.OsgoodInverse.mirror_step
deprecated_alias MultAdd.OsgoodInverse.nr := CircleOp.OsgoodInverse.nr
deprecated_alias MultAdd.OsgoodInverse.nr_add := CircleOp.OsgoodInverse.nr_add
deprecated_alias MultAdd.OsgoodInverse.nr_mono := CircleOp.OsgoodInverse.nr_mono
deprecated_alias MultAdd.OsgoodInverse.nr_nonneg := CircleOp.OsgoodInverse.nr_nonneg
deprecated_alias MultAdd.OsgoodInverse.nr_pos := CircleOp.OsgoodInverse.nr_pos
deprecated_alias MultAdd.OsgoodInverse.nr_pow2 := CircleOp.OsgoodInverse.nr_pow2
deprecated_alias MultAdd.OsgoodInverse.one_div_pos := CircleOp.OsgoodInverse.one_div_pos
deprecated_alias MultAdd.OsgoodInverse.one_le_pow2 := CircleOp.OsgoodInverse.one_le_pow2
deprecated_alias MultAdd.OsgoodInverse.pow2_succ := CircleOp.OsgoodInverse.pow2_succ
deprecated_alias MultAdd.OsgoodInverse.rad1 := CircleOp.OsgoodInverse.rad1
deprecated_alias MultAdd.OsgoodInverse.rad2 := CircleOp.OsgoodInverse.rad2
deprecated_alias MultAdd.OsgoodInverse.rad_out := CircleOp.OsgoodInverse.rad_out
deprecated_alias MultAdd.OsgoodInverse.rad_speed := CircleOp.OsgoodInverse.rad_speed
deprecated_alias MultAdd.OsgoodInverse.recip_diff := CircleOp.OsgoodInverse.recip_diff
deprecated_alias MultAdd.OsgoodInverse.ric := CircleOp.OsgoodInverse.ric
deprecated_alias MultAdd.OsgoodInverse.ric_lift := CircleOp.OsgoodInverse.ric_lift
deprecated_alias MultAdd.OsgoodInverse.ric_near_one := CircleOp.OsgoodInverse.ric_near_one
deprecated_alias MultAdd.OsgoodInverse.ric_step := CircleOp.OsgoodInverse.ric_step
deprecated_alias MultAdd.OsgoodInverse.ric_step_multAdd := CircleOp.OsgoodInverse.ric_step_circ
deprecated_alias MultAdd.OsgoodInverse.riccati_clock := CircleOp.OsgoodInverse.riccati_clock
deprecated_alias MultAdd.OsgoodInverse.riccati_clock_step := CircleOp.OsgoodInverse.riccati_clock_step
deprecated_alias MultAdd.OsgoodInverse.riccati_hits_wall := CircleOp.OsgoodInverse.riccati_hits_wall
deprecated_alias MultAdd.OsgoodInverse.riccati_quarter_blows_up := CircleOp.OsgoodInverse.riccati_quarter_blows_up
deprecated_alias MultAdd.OsgoodInverse.root := CircleOp.OsgoodInverse.root
deprecated_alias MultAdd.OsgoodInverse.root_sq := CircleOp.OsgoodInverse.root_sq
deprecated_alias MultAdd.OsgoodInverse.rot1 := CircleOp.OsgoodInverse.rot1
deprecated_alias MultAdd.OsgoodInverse.rot2 := CircleOp.OsgoodInverse.rot2
deprecated_alias MultAdd.OsgoodInverse.rot_perp := CircleOp.OsgoodInverse.rot_perp
deprecated_alias MultAdd.OsgoodInverse.rot_speed := CircleOp.OsgoodInverse.rot_speed
deprecated_alias MultAdd.OsgoodInverse.rung1_lift := CircleOp.OsgoodInverse.rung1_lift
deprecated_alias MultAdd.OsgoodInverse.rung1_never_wall := CircleOp.OsgoodInverse.rung1_never_wall
deprecated_alias MultAdd.OsgoodInverse.steps := CircleOp.OsgoodInverse.steps
deprecated_alias MultAdd.OsgoodInverse.steps_add := CircleOp.OsgoodInverse.steps_add
deprecated_alias MultAdd.OsgoodInverse.steps_quarter := CircleOp.OsgoodInverse.steps_quarter
deprecated_alias MultAdd.OsgoodInverse.time_left_pos := CircleOp.OsgoodInverse.time_left_pos

-- FusedAlgebra.Topics.Pi
deprecated_alias PiMultAdd.G := PiCircle.G
deprecated_alias PiMultAdd.G.add := PiCircle.G.add
deprecated_alias PiMultAdd.G.casesOn := PiCircle.G.casesOn
deprecated_alias PiMultAdd.G.conj := PiCircle.G.conj
deprecated_alias PiMultAdd.G.ext' := PiCircle.G.ext'
deprecated_alias PiMultAdd.G.ext'_iff := PiCircle.G.ext'_iff
deprecated_alias PiMultAdd.G.im := PiCircle.G.im
deprecated_alias PiMultAdd.G.mk := PiCircle.G.mk
deprecated_alias PiMultAdd.G.mk.noConfusion := PiCircle.G.mk.noConfusion
deprecated_alias PiMultAdd.G.mul := PiCircle.G.mul
deprecated_alias PiMultAdd.G.noConfusion := PiCircle.G.noConfusion
deprecated_alias PiMultAdd.G.norm := PiCircle.G.norm
deprecated_alias PiMultAdd.G.one := PiCircle.G.one
deprecated_alias PiMultAdd.G.pow := PiCircle.G.pow
deprecated_alias PiMultAdd.G.re := PiCircle.G.re
deprecated_alias PiMultAdd.G.rec := PiCircle.G.rec
deprecated_alias PiMultAdd.G.recOn := PiCircle.G.recOn
deprecated_alias PiMultAdd.G.zero := PiCircle.G.zero
deprecated_alias PiMultAdd.Slope := PiCircle.Slope
deprecated_alias PiMultAdd.Slope.casesOn := PiCircle.Slope.casesOn
deprecated_alias PiMultAdd.Slope.den := PiCircle.Slope.den
deprecated_alias PiMultAdd.Slope.mk := PiCircle.Slope.mk
deprecated_alias PiMultAdd.Slope.mk.noConfusion := PiCircle.Slope.mk.noConfusion
deprecated_alias PiMultAdd.Slope.noConfusion := PiCircle.Slope.noConfusion
deprecated_alias PiMultAdd.Slope.num := PiCircle.Slope.num
deprecated_alias PiMultAdd.Slope.rec := PiCircle.Slope.rec
deprecated_alias PiMultAdd.Slope.recOn := PiCircle.Slope.recOn
deprecated_alias PiMultAdd.circle_closed := PiCircle.circle_closed
deprecated_alias PiMultAdd.circle_real_points := PiCircle.circle_real_points
deprecated_alias PiMultAdd.cross_10 := PiCircle.cross_10
deprecated_alias PiMultAdd.cross_100 := PiCircle.cross_100
deprecated_alias PiMultAdd.cross_1000 := PiCircle.cross_1000
deprecated_alias PiMultAdd.div_le_div_of := PiCircle.div_le_div_of
deprecated_alias PiMultAdd.euler := PiCircle.euler
deprecated_alias PiMultAdd.firstCross := PiCircle.firstCross
deprecated_alias PiMultAdd.gma := PiCircle.gma
deprecated_alias PiMultAdd.gma_assoc := PiCircle.gma_assoc
deprecated_alias PiMultAdd.gma_comm := PiCircle.gma_comm
deprecated_alias PiMultAdd.gma_lift := PiCircle.gma_lift
deprecated_alias PiMultAdd.gma_zero := PiCircle.gma_zero
deprecated_alias PiMultAdd.gpow := PiCircle.gpow
deprecated_alias PiMultAdd.gpow_lift := PiCircle.gpow_lift
deprecated_alias PiMultAdd.greg := PiCircle.greg
deprecated_alias PiMultAdd.greg_even_le_odd := PiCircle.greg_even_le_odd
deprecated_alias PiMultAdd.greg_even_up := PiCircle.greg_even_up
deprecated_alias PiMultAdd.greg_odd_above := PiCircle.greg_odd_above
deprecated_alias PiMultAdd.greg_odd_down := PiCircle.greg_odd_down
deprecated_alias PiMultAdd.greg_two := PiCircle.greg_two
deprecated_alias PiMultAdd.gt := PiCircle.gt
deprecated_alias PiMultAdd.gt_nonneg := PiCircle.gt_nonneg
deprecated_alias PiMultAdd.hermann := PiCircle.hermann
deprecated_alias PiMultAdd.hutton := PiCircle.hutton
deprecated_alias PiMultAdd.instDecidableEqG := PiCircle.instDecidableEqG
deprecated_alias PiMultAdd.instDecidableEqG.decEq := PiCircle.instDecidableEqG.decEq
deprecated_alias PiMultAdd.instDecidableEqSlope := PiCircle.instDecidableEqSlope
deprecated_alias PiMultAdd.instDecidableEqSlope.decEq := PiCircle.instDecidableEqSlope.decEq
deprecated_alias PiMultAdd.instReprG := PiCircle.instReprG
deprecated_alias PiMultAdd.instReprG.repr := PiCircle.instReprG.repr
deprecated_alias PiMultAdd.instReprSlope := PiCircle.instReprSlope
deprecated_alias PiMultAdd.instReprSlope.repr := PiCircle.instReprSlope.repr
deprecated_alias PiMultAdd.lens := PiCircle.lens
deprecated_alias PiMultAdd.lens_conj := PiCircle.lens_conj
deprecated_alias PiMultAdd.lens_infinity := PiCircle.lens_infinity
deprecated_alias PiMultAdd.lens_on_circle := PiCircle.lens_on_circle
deprecated_alias PiMultAdd.lens_zero := PiCircle.lens_zero
deprecated_alias PiMultAdd.machin := PiCircle.machin
deprecated_alias PiMultAdd.machinHi := PiCircle.machinHi
deprecated_alias PiMultAdd.machinLo := PiCircle.machinLo
deprecated_alias PiMultAdd.machin_50_digits := PiCircle.machin_50_digits
deprecated_alias PiMultAdd.machin_lo_le_hi := PiCircle.machin_lo_le_hi
deprecated_alias PiMultAdd.machin_nested := PiCircle.machin_nested
deprecated_alias PiMultAdd.machin_slopes := PiCircle.machin_slopes
deprecated_alias PiMultAdd.mul_assoc' := PiCircle.mul_assoc'
deprecated_alias PiMultAdd.multAdd_self_zero := PiCircle.circ_self_zero
deprecated_alias PiMultAdd.neg_one_pow_even := PiCircle.neg_one_pow_even
deprecated_alias PiMultAdd.neg_one_pow_odd := PiCircle.neg_one_pow_odd
deprecated_alias PiMultAdd.norm_conj := PiCircle.norm_conj
deprecated_alias PiMultAdd.norm_mul := PiCircle.norm_mul
deprecated_alias PiMultAdd.parity_twice := PiCircle.parity_twice
deprecated_alias PiMultAdd.q4 := PiCircle.q4
deprecated_alias PiMultAdd.quarter_four := PiCircle.quarter_four
deprecated_alias PiMultAdd.quarter_on_circle := PiCircle.quarter_on_circle
deprecated_alias PiMultAdd.quarter_twice := PiCircle.quarter_twice
deprecated_alias PiMultAdd.slope_one_four := PiCircle.slope_one_four
deprecated_alias PiMultAdd.slope_one_twice := PiCircle.slope_one_twice
deprecated_alias PiMultAdd.term_shrink := PiCircle.term_shrink
deprecated_alias PiMultAdd.tl := PiCircle.tl
deprecated_alias PiMultAdd.toG := PiCircle.toG

-- FusedAlgebra.Topics.PvsNP
deprecated_alias PvsNP.foldMA_ones := PvsNP.foldCirc_ones
deprecated_alias PvsNP.gf2MultAdd := PvsNP.gf2Circle
deprecated_alias PvsNP.gf2_multAdd_is_or := PvsNP.gf2_circ_is_or
deprecated_alias PvsNP.multAdd_sAdd_distrib := PvsNP.circ_sAdd_distrib

-- FusedAlgebra.Topics.QAcross.Partitions
deprecated_alias QPartitions.eulerMA := QPartitions.eulerCirc
deprecated_alias QPartitions.euler_is_multAdd_30 := QPartitions.euler_is_circ_30

-- FusedAlgebra.Topics.QAcross.QCounting
deprecated_alias QCounting.ising_multAdd := QCounting.ising_circ

-- FusedAlgebra.Topics.QuantumNoise
deprecated_alias QuantumNoise.madd_inverse := QuantumNoise.circ_inverse

-- FusedAlgebra.Topics.Rulecraft
deprecated_alias MultAdd.Rulecraft.balance := CircleOp.Rulecraft.balance
deprecated_alias MultAdd.Rulecraft.balance_le := CircleOp.Rulecraft.balance_le
deprecated_alias MultAdd.Rulecraft.boundary := CircleOp.Rulecraft.boundary
deprecated_alias MultAdd.Rulecraft.fexp_not_comm_values := CircleOp.Rulecraft.fexp_not_comm_values
deprecated_alias MultAdd.Rulecraft.fused_price_expand := CircleOp.Rulecraft.fused_price_expand
deprecated_alias MultAdd.Rulecraft.fused_price_step := CircleOp.Rulecraft.fused_price_step
deprecated_alias MultAdd.Rulecraft.panel_grid := CircleOp.Rulecraft.panel_grid
deprecated_alias MultAdd.Rulecraft.recovery := CircleOp.Rulecraft.recovery
deprecated_alias MultAdd.Rulecraft.sq_nonneg' := CircleOp.Rulecraft.sq_nonneg'
deprecated_alias MultAdd.Rulecraft.thermostat_residue := CircleOp.Rulecraft.thermostat_residue

-- FusedAlgebra.Topics.SolitonsRungs
deprecated_alias SolitonsRungs.glue_is_multAdd := SolitonsRungs.glue_is_circ
deprecated_alias SolitonsRungs.soliton_vs_wAdd_corner := SolitonsRungs.soliton_vs_circK_corner
deprecated_alias SolitonsRungs.tau2_is_wAdd := SolitonsRungs.tau2_is_circK
deprecated_alias SolitonsRungs.wAdd_three_corner := SolitonsRungs.circK_three_corner

-- FusedAlgebra.Topics.Statistics
deprecated_alias StatFused.evidence_multAdd := StatFused.evidence_circ
deprecated_alias StatFused.union_neg_foldMA := StatFused.union_neg_foldCirc

-- FusedAlgebra.Topics.Thermodynamics
deprecated_alias ThermoFused.chain3_multAdd := ThermoFused.chain3_circ

-- FusedAlgebra.Topics.Topology.CubeProofs
deprecated_alias MultAdd.Cube.IsPrime := CircleOp.Cube.IsPrime
deprecated_alias MultAdd.Cube.MEq := CircleOp.Cube.MEq
deprecated_alias MultAdd.Cube.MEq.add := CircleOp.Cube.MEq.add
deprecated_alias MultAdd.Cube.MEq.mul := CircleOp.Cube.MEq.mul
deprecated_alias MultAdd.Cube.MEq.refl := CircleOp.Cube.MEq.refl
deprecated_alias MultAdd.Cube.MEq.symm := CircleOp.Cube.MEq.symm
deprecated_alias MultAdd.Cube.MEq.trans := CircleOp.Cube.MEq.trans
deprecated_alias MultAdd.Cube.MEq.zero_of_dvd := CircleOp.Cube.MEq.zero_of_dvd
deprecated_alias MultAdd.Cube.aks_middle := CircleOp.Cube.aks_middle
deprecated_alias MultAdd.Cube.binomial := CircleOp.Cube.binomial
deprecated_alias MultAdd.Cube.blockCount := CircleOp.Cube.blockCount
deprecated_alias MultAdd.Cube.blockCount_eq := CircleOp.Cube.blockCount_eq
deprecated_alias MultAdd.Cube.bsum := CircleOp.Cube.bsum
deprecated_alias MultAdd.Cube.bsum_add := CircleOp.Cube.bsum_add
deprecated_alias MultAdd.Cube.bsum_congr := CircleOp.Cube.bsum_congr
deprecated_alias MultAdd.Cube.bsum_dvd := CircleOp.Cube.bsum_dvd
deprecated_alias MultAdd.Cube.bsum_mul := CircleOp.Cube.bsum_mul
deprecated_alias MultAdd.Cube.bsum_shift := CircleOp.Cube.bsum_shift
deprecated_alias MultAdd.Cube.coprime_of_lt := CircleOp.Cube.coprime_of_lt
deprecated_alias MultAdd.Cube.cubeSum := CircleOp.Cube.cubeSum
deprecated_alias MultAdd.Cube.digitProd := CircleOp.Cube.digitProd
deprecated_alias MultAdd.Cube.divmod_of_eq := CircleOp.Cube.divmod_of_eq
deprecated_alias MultAdd.Cube.exists_prime_factor := CircleOp.Cube.exists_prime_factor
deprecated_alias MultAdd.Cube.fermat := CircleOp.Cube.fermat
deprecated_alias MultAdd.Cube.frobenius := CircleOp.Cube.frobenius
deprecated_alias MultAdd.Cube.layer_zero_left := CircleOp.Cube.layer_zero_left
deprecated_alias MultAdd.Cube.lucas := CircleOp.Cube.lucas
deprecated_alias MultAdd.Cube.lucas_step := CircleOp.Cube.lucas_step
deprecated_alias MultAdd.Cube.prime_dvd_middle := CircleOp.Cube.prime_dvd_middle

-- FusedAlgebra.Topics.Topology.HopfInvariantOne
deprecated_alias MultAdd.HopfOne.Elt := CircleOp.HopfOne.Elt
deprecated_alias MultAdd.HopfOne.Elt.casesOn := CircleOp.HopfOne.Elt.casesOn
deprecated_alias MultAdd.HopfOne.Elt.mk := CircleOp.HopfOne.Elt.mk
deprecated_alias MultAdd.HopfOne.Elt.mk.noConfusion := CircleOp.HopfOne.Elt.mk.noConfusion
deprecated_alias MultAdd.HopfOne.Elt.noConfusion := CircleOp.HopfOne.Elt.noConfusion
deprecated_alias MultAdd.HopfOne.Elt.rec := CircleOp.HopfOne.Elt.rec
deprecated_alias MultAdd.HopfOne.Elt.recOn := CircleOp.HopfOne.Elt.recOn
deprecated_alias MultAdd.HopfOne.Elt.«α» := CircleOp.HopfOne.Elt.«α»
deprecated_alias MultAdd.HopfOne.Elt.«β» := CircleOp.HopfOne.Elt.«β»
deprecated_alias MultAdd.HopfOne.IsOdd := CircleOp.HopfOne.IsOdd
deprecated_alias MultAdd.HopfOne.OneMod := CircleOp.HopfOne.OneMod
deprecated_alias MultAdd.HopfOne.Trunc := CircleOp.HopfOne.Trunc
deprecated_alias MultAdd.HopfOne.Trunc.casesOn := CircleOp.HopfOne.Trunc.casesOn
deprecated_alias MultAdd.HopfOne.Trunc.mk := CircleOp.HopfOne.Trunc.mk
deprecated_alias MultAdd.HopfOne.Trunc.mk.noConfusion := CircleOp.HopfOne.Trunc.mk.noConfusion
deprecated_alias MultAdd.HopfOne.Trunc.noConfusion := CircleOp.HopfOne.Trunc.noConfusion
deprecated_alias MultAdd.HopfOne.Trunc.rec := CircleOp.HopfOne.Trunc.rec
deprecated_alias MultAdd.HopfOne.Trunc.recOn := CircleOp.HopfOne.Trunc.recOn
deprecated_alias MultAdd.HopfOne.Trunc.«α» := CircleOp.HopfOne.Trunc.«α»
deprecated_alias MultAdd.HopfOne.Trunc.«β» := CircleOp.HopfOne.Trunc.«β»
deprecated_alias MultAdd.HopfOne.commute_equation := CircleOp.HopfOne.commute_equation
deprecated_alias MultAdd.HopfOne.cp2_case := CircleOp.HopfOne.cp2_case
deprecated_alias MultAdd.HopfOne.hopf_invariant_one := CircleOp.HopfOne.hopf_invariant_one
deprecated_alias MultAdd.HopfOne.instDecidableEqTrunc := CircleOp.HopfOne.instDecidableEqTrunc
deprecated_alias MultAdd.HopfOne.instDecidableEqTrunc.decEq := CircleOp.HopfOne.instDecidableEqTrunc.decEq
deprecated_alias MultAdd.HopfOne.le_of_dvd_two_pow_odd := CircleOp.HopfOne.le_of_dvd_two_pow_odd
deprecated_alias MultAdd.HopfOne.level_three_pow := CircleOp.HopfOne.level_three_pow
deprecated_alias MultAdd.HopfOne.madd := CircleOp.HopfOne.circ
deprecated_alias MultAdd.HopfOne.mpow := CircleOp.HopfOne.mpow
deprecated_alias MultAdd.HopfOne.mpow_formula := CircleOp.HopfOne.mpow_formula
deprecated_alias MultAdd.HopfOne.odd_mul := CircleOp.HopfOne.odd_mul
deprecated_alias MultAdd.HopfOne.odd_three_pow := CircleOp.HopfOne.odd_three_pow
deprecated_alias MultAdd.HopfOne.oneMod_oddPow := CircleOp.HopfOne.oneMod_oddPow
deprecated_alias MultAdd.HopfOne.oneMod_sq := CircleOp.HopfOne.oneMod_sq
deprecated_alias MultAdd.HopfOne.one_two_four_work := CircleOp.HopfOne.one_two_four_work
deprecated_alias MultAdd.HopfOne.only_one_two_four := CircleOp.HopfOne.only_one_two_four
deprecated_alias MultAdd.HopfOne.pow_ge := CircleOp.HopfOne.pow_ge
deprecated_alias MultAdd.HopfOne.psi := CircleOp.HopfOne.psi
deprecated_alias MultAdd.HopfOne.psiT := CircleOp.HopfOne.psiT
deprecated_alias MultAdd.HopfOne.psi_compose := CircleOp.HopfOne.psi_compose
deprecated_alias MultAdd.HopfOne.three_pow_two_pow := CircleOp.HopfOne.three_pow_two_pow
deprecated_alias MultAdd.HopfOne.two_adic_split := CircleOp.HopfOne.two_adic_split
deprecated_alias MultAdd.HopfOne.two_pow_dvd_of_mul_odd := CircleOp.HopfOne.two_pow_dvd_of_mul_odd

-- FusedAlgebra.Topics.Topology.KummerWilson
deprecated_alias MultAdd.Cube.S := CircleOp.Cube.S
deprecated_alias MultAdd.Cube.S_add := CircleOp.Cube.S_add
deprecated_alias MultAdd.Cube.S_bsum := CircleOp.Cube.S_bsum
deprecated_alias MultAdd.Cube.S_mul := CircleOp.Cube.S_mul
deprecated_alias MultAdd.Cube.S_pow := CircleOp.Cube.S_pow
deprecated_alias MultAdd.Cube.S_sub := CircleOp.Cube.S_sub
deprecated_alias MultAdd.Cube.S_succ := CircleOp.Cube.S_succ
deprecated_alias MultAdd.Cube.S_zero := CircleOp.Cube.S_zero
deprecated_alias MultAdd.Cube.Val := CircleOp.Cube.Val
deprecated_alias MultAdd.Cube.Val.mul := CircleOp.Cube.Val.mul
deprecated_alias MultAdd.Cube.Val.unique := CircleOp.Cube.Val.unique
deprecated_alias MultAdd.Cube.alt_layer := CircleOp.Cube.alt_layer
deprecated_alias MultAdd.Cube.bsum_meq := CircleOp.Cube.bsum_meq
deprecated_alias MultAdd.Cube.bsum_one := CircleOp.Cube.bsum_one
deprecated_alias MultAdd.Cube.carries := CircleOp.Cube.carries
deprecated_alias MultAdd.Cube.carries_dsum := CircleOp.Cube.carries_dsum
deprecated_alias MultAdd.Cube.carries_eq := CircleOp.Cube.carries_eq
deprecated_alias MultAdd.Cube.carries_zero := CircleOp.Cube.carries_zero
deprecated_alias MultAdd.Cube.dsum := CircleOp.Cube.dsum
deprecated_alias MultAdd.Cube.dsum_eq := CircleOp.Cube.dsum_eq
deprecated_alias MultAdd.Cube.dsum_succ := CircleOp.Cube.dsum_succ
deprecated_alias MultAdd.Cube.dsum_zero := CircleOp.Cube.dsum_zero
deprecated_alias MultAdd.Cube.dvd_fact := CircleOp.Cube.dvd_fact
deprecated_alias MultAdd.Cube.euclid := CircleOp.Cube.euclid
deprecated_alias MultAdd.Cube.fermat_unit := CircleOp.Cube.fermat_unit
deprecated_alias MultAdd.Cube.kummer := CircleOp.Cube.kummer
deprecated_alias MultAdd.Cube.layer_pos := CircleOp.Cube.layer_pos
deprecated_alias MultAdd.Cube.layer_succ_self := CircleOp.Cube.layer_succ_self
deprecated_alias MultAdd.Cube.layer_val_dsum := CircleOp.Cube.layer_val_dsum
deprecated_alias MultAdd.Cube.val := CircleOp.Cube.val
deprecated_alias MultAdd.Cube.val_of_mod_eq := CircleOp.Cube.val_of_mod_eq
deprecated_alias MultAdd.Cube.val_of_mod_ne := CircleOp.Cube.val_of_mod_ne
deprecated_alias MultAdd.Cube.val_spec := CircleOp.Cube.val_spec
deprecated_alias MultAdd.Cube.wilson := CircleOp.Cube.wilson
deprecated_alias MultAdd.Cube.wilson_iff := CircleOp.Cube.wilson_iff

-- FusedAlgebra.Topics.Toric.ProjectivePlane
deprecated_alias MultAdd.Toric.P2.Cone := CircleOp.Toric.P2.Cone
deprecated_alias MultAdd.Toric.P2.Cone.c0 := CircleOp.Toric.P2.Cone.c0
deprecated_alias MultAdd.Toric.P2.Cone.c0.elim := CircleOp.Toric.P2.Cone.c0.elim
deprecated_alias MultAdd.Toric.P2.Cone.c1 := CircleOp.Toric.P2.Cone.c1
deprecated_alias MultAdd.Toric.P2.Cone.c1.elim := CircleOp.Toric.P2.Cone.c1.elim
deprecated_alias MultAdd.Toric.P2.Cone.c2 := CircleOp.Toric.P2.Cone.c2
deprecated_alias MultAdd.Toric.P2.Cone.c2.elim := CircleOp.Toric.P2.Cone.c2.elim
deprecated_alias MultAdd.Toric.P2.Cone.casesOn := CircleOp.Toric.P2.Cone.casesOn
deprecated_alias MultAdd.Toric.P2.Cone.ctorElim := CircleOp.Toric.P2.Cone.ctorElim
deprecated_alias MultAdd.Toric.P2.Cone.noConfusion := CircleOp.Toric.P2.Cone.noConfusion
deprecated_alias MultAdd.Toric.P2.Cone.rec := CircleOp.Toric.P2.Cone.rec
deprecated_alias MultAdd.Toric.P2.Cone.recOn := CircleOp.Toric.P2.Cone.recOn
deprecated_alias MultAdd.Toric.P2.Dom := CircleOp.Toric.P2.Dom
deprecated_alias MultAdd.Toric.P2.InChart := CircleOp.Toric.P2.InChart
deprecated_alias MultAdd.Toric.P2.P2 := CircleOp.Toric.P2.P2
deprecated_alias MultAdd.Toric.P2.chart := CircleOp.Toric.P2.chart
deprecated_alias MultAdd.Toric.P2.chart_cover := CircleOp.Toric.P2.chart_cover
deprecated_alias MultAdd.Toric.P2.chart_injective := CircleOp.Toric.P2.chart_injective
deprecated_alias MultAdd.Toric.P2.chart_trans := CircleOp.Toric.P2.chart_trans
deprecated_alias MultAdd.Toric.P2.chart_wind := CircleOp.Toric.P2.chart_wind
deprecated_alias MultAdd.Toric.P2.fan_count := CircleOp.Toric.P2.fan_count
deprecated_alias MultAdd.Toric.P2.inChart_iff := CircleOp.Toric.P2.inChart_iff
deprecated_alias MultAdd.Toric.P2.instDecidableEqCone := CircleOp.Toric.P2.instDecidableEqCone
deprecated_alias MultAdd.Toric.P2.instReprCone := CircleOp.Toric.P2.instReprCone
deprecated_alias MultAdd.Toric.P2.instReprCone.repr := CircleOp.Toric.P2.instReprCone.repr
deprecated_alias MultAdd.Toric.P2.inv_inv' := CircleOp.Toric.P2.inv_inv'
deprecated_alias MultAdd.Toric.P2.inv_ne := CircleOp.Toric.P2.inv_ne
deprecated_alias MultAdd.Toric.P2.inv_pow' := CircleOp.Toric.P2.inv_pow'
deprecated_alias MultAdd.Toric.P2.origin_in_one_chart := CircleOp.Toric.P2.origin_in_one_chart
deprecated_alias MultAdd.Toric.P2.pow_succ_eq_zero := CircleOp.Toric.P2.pow_succ_eq_zero
deprecated_alias MultAdd.Toric.P2.pw := CircleOp.Toric.P2.pw
deprecated_alias MultAdd.Toric.P2.torus_in_all := CircleOp.Toric.P2.torus_in_all
deprecated_alias MultAdd.Toric.P2.trans := CircleOp.Toric.P2.trans
deprecated_alias MultAdd.Toric.P2.trans01_multAdd := CircleOp.Toric.P2.trans01_circ
deprecated_alias MultAdd.Toric.P2.trans02_multAdd := CircleOp.Toric.P2.trans02_circ
deprecated_alias MultAdd.Toric.P2.trans_cocycle := CircleOp.Toric.P2.trans_cocycle
deprecated_alias MultAdd.Toric.P2.trans_wind := CircleOp.Toric.P2.trans_wind
deprecated_alias MultAdd.Toric.P2.wind := CircleOp.Toric.P2.wind

-- FusedAlgebra.Topics.Toric.ToricCharts
deprecated_alias MultAdd.Toric.InChart1 := CircleOp.Toric.InChart1
deprecated_alias MultAdd.Toric.InChart2 := CircleOp.Toric.InChart2
deprecated_alias MultAdd.Toric.P1 := CircleOp.Toric.P1
deprecated_alias MultAdd.Toric.chart1 := CircleOp.Toric.chart1
deprecated_alias MultAdd.Toric.chart1_cover := CircleOp.Toric.chart1_cover
deprecated_alias MultAdd.Toric.chart1_injective := CircleOp.Toric.chart1_injective
deprecated_alias MultAdd.Toric.chart1_trans := CircleOp.Toric.chart1_trans
deprecated_alias MultAdd.Toric.chart2 := CircleOp.Toric.chart2
deprecated_alias MultAdd.Toric.chart2_cover := CircleOp.Toric.chart2_cover
deprecated_alias MultAdd.Toric.chart2_injective := CircleOp.Toric.chart2_injective
deprecated_alias MultAdd.Toric.chart2_trans := CircleOp.Toric.chart2_trans
deprecated_alias MultAdd.Toric.corner := CircleOp.Toric.corner
deprecated_alias MultAdd.Toric.corner_all_on := CircleOp.Toric.corner_all_on
deprecated_alias MultAdd.Toric.corner_scale := CircleOp.Toric.corner_scale
deprecated_alias MultAdd.Toric.fan_count := CircleOp.Toric.fan_count
deprecated_alias MultAdd.Toric.fixed_point_in_one_chart := CircleOp.Toric.fixed_point_in_one_chart
deprecated_alias MultAdd.Toric.inChart2_iff := CircleOp.Toric.inChart2_iff
deprecated_alias MultAdd.Toric.origin_in_one_chart := CircleOp.Toric.origin_in_one_chart
deprecated_alias MultAdd.Toric.torus_in_all := CircleOp.Toric.torus_in_all
deprecated_alias MultAdd.Toric.torus_in_both := CircleOp.Toric.torus_in_both
deprecated_alias MultAdd.Toric.trans1 := CircleOp.Toric.trans1
deprecated_alias MultAdd.Toric.trans1_cocycle := CircleOp.Toric.trans1_cocycle
deprecated_alias MultAdd.Toric.trans2 := CircleOp.Toric.trans2
deprecated_alias MultAdd.Toric.trans2_cocycle := CircleOp.Toric.trans2_cocycle

/-! ## Test: an old name still elaborates, and says what replaced it -/

/--
warning: `MultAdd.multAdd_comm` has been deprecated: Use `CircleOp.circ_comm` instead
-/
#guard_msgs (warning, drop info, substring := true) in
example (a b : Nat) : CircleOp.circ a b = CircleOp.circ b a := MultAdd.multAdd_comm a b

