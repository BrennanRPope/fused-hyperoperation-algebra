#!/bin/bash
# Compile the whole library with core Lean (no mathlib, no lake needed).
# Usage: ./build.sh            (expects `lean` on PATH, version in lean-toolchain)
set -e
cd "$(dirname "$0")"
mods=(
  FusedAlgebra.Base FusedAlgebra.Rungs FusedAlgebra.Cube
  FusedAlgebra.Foundations.CircleRing FusedAlgebra.Foundations.Independence FusedAlgebra.Foundations.FusedPeano FusedAlgebra.Foundations.RungOneSemiring
  FusedAlgebra.Topics.ClosingRate FusedAlgebra.Topics.Consensus FusedAlgebra.Topics.DarkSector
  FusedAlgebra.Topics.DataCenter FusedAlgebra.Topics.E FusedAlgebra.Topics.Finance
  FusedAlgebra.Topics.FusedConstants FusedAlgebra.Topics.FractionalRungs FusedAlgebra.Topics.Infinities FusedAlgebra.Topics.InformationTheory
  FusedAlgebra.Topics.Landscape FusedAlgebra.Topics.RungShape FusedAlgebra.Topics.Elliptic FusedAlgebra.Topics.EllipticS3
  FusedAlgebra.Topics.Thermodynamics FusedAlgebra.Topics.Relativity FusedAlgebra.Topics.LeeYang FusedAlgebra.Topics.LeeYangCluster FusedAlgebra.Topics.LeeYangSpill FusedAlgebra.Topics.TwistedBoosts FusedAlgebra.Topics.GravityProbeB FusedAlgebra.Topics.BinaryPulsar
  FusedAlgebra.Topics.Chemistry.ChemQuantum FusedAlgebra.Topics.Chemistry.ChemConservation FusedAlgebra.Topics.Chemistry.ChemPeriodic FusedAlgebra.Topics.Chemistry.ChemBonding FusedAlgebra.Topics.Chemistry.ChemStoichiometry FusedAlgebra.Topics.Chemistry.ChemGas FusedAlgebra.Topics.Chemistry.ChemGibbs FusedAlgebra.Topics.Chemistry.ChemEquilibrium FusedAlgebra.Topics.Chemistry.ChemRates FusedAlgebra.Topics.Chemistry.ChemElectro FusedAlgebra.Topics.Chemistry.ChemAcids FusedAlgebra.Topics.Chemistry.ChemOrganic FusedAlgebra.Topics.Chemistry.ChemBiochem FusedAlgebra.Topics.ShapeDirections
  FusedAlgebra.Topics.NumberPlane FusedAlgebra.Topics.GeometryTrig FusedAlgebra.Topics.EarlyScience FusedAlgebra.Topics.Patterns FusedAlgebra.Topics.EllipticCounting FusedAlgebra.Topics.Pi FusedAlgebra.Topics.Calculus FusedAlgebra.Topics.AbstractAlgebra FusedAlgebra.Topics.Statistics
  FusedAlgebra.Topics.NavierStokes FusedAlgebra.Topics.QuantumNoise FusedAlgebra.Topics.Rulecraft
  FusedAlgebra.Topics.Toric.ToricCharts FusedAlgebra.Topics.Toric.ProjectivePlane
  FusedAlgebra.Topics.Topology.CubeProofs FusedAlgebra.Topics.Topology.KummerWilson
  FusedAlgebra.Topics.Topology.HopfInvariantOne
  FusedAlgebra.Topics.NSContinuum FusedAlgebra.Topics.JetFields FusedAlgebra.Topics.JetRiver FusedAlgebra.Topics.SpinningRiver FusedAlgebra.Topics.OsgoodInverse FusedAlgebra.Topics.LinearAlgebra FusedAlgebra.Topics.DialSweep FusedAlgebra.Topics.EntropyNoise
  FusedAlgebra.Topics.QAcross.QDial FusedAlgebra.Topics.QAcross.Knots FusedAlgebra.Topics.QAcross.Partitions FusedAlgebra.Topics.QAcross.QCounting
  FusedAlgebra.Topics.SolitonsRungs FusedAlgebra.Topics.PvsNP
  FusedAlgebra.Topics.Fermat.Fermat4 FusedAlgebra.Topics.Fermat.Fermat3 FusedAlgebra.Topics.Fermat.FreyPieces FusedAlgebra.Topics.Fermat.FermatRoute FusedAlgebra.Topics.Fermat.CounterexampleMap
  FusedAlgebra.Topics.LogicSets FusedAlgebra.Topics.CategoryTheory FusedAlgebra.Topics.Combinatorics FusedAlgebra.Topics.GroupsLie FusedAlgebra.Topics.Clifford FusedAlgebra.Topics.GraphTheory
  FusedAlgebra.Topics.MeasureProbability FusedAlgebra.Topics.ComplexAnalysis FusedAlgebra.Topics.FunctionalAnalysis FusedAlgebra.Topics.DiffGeometry FusedAlgebra.Topics.KnotTheory
  FusedAlgebra.Topics.Dynamics FusedAlgebra.Topics.InformationFull FusedAlgebra.Topics.Numerical FusedAlgebra.Topics.GameTheory FusedAlgebra.Topics.DiamondLattice FusedAlgebra.Topics.DialPointProcesses FusedAlgebra.Topics.SeveralWaves FusedAlgebra.Topics.RuleBreakers \
  FusedAlgebra.Topics.Biology.MolecularBiology FusedAlgebra.Topics.Biology.Metabolism FusedAlgebra.Topics.Biology.CellBiology FusedAlgebra.Topics.Biology.Genetics FusedAlgebra.Topics.Biology.Development FusedAlgebra.Topics.Biology.Physiology FusedAlgebra.Topics.Biology.NeuroScience FusedAlgebra.Topics.Biology.Evolution FusedAlgebra.Topics.Biology.Population FusedAlgebra.Topics.Biology.Ecology FusedAlgebra.Topics.Biology.Behaviour FusedAlgebra.Topics.FractionalPowerset
  FusedAlgebra.Topics.GeneralRelativity FusedAlgebra.Topics.GRBridges FusedAlgebra.Topics.WarpWormhole FusedAlgebra.Topics.Electromagnetism FusedAlgebra.Topics.PowersetDials FusedAlgebra.Topics.DialSurvey FusedAlgebra.Topics.SizeDialResidual FusedAlgebra.Topics.DialWall FusedAlgebra.Topics.FractalBoxes FusedAlgebra.Topics.FractalDials FusedAlgebra.Topics.DialPotential FusedAlgebra.Topics.RepetitionDial FusedAlgebra.Topics.ZoomDial FusedAlgebra.Topics.Conventions FusedAlgebra.Topics.PathNames FusedAlgebra.Topics.OrderDial FusedAlgebra.Topics.WhyI FusedAlgebra.Topics.Checkerboard FusedAlgebra.Topics.CheckerboardPotential FusedAlgebra.Topics.EIPi FusedAlgebra.Foundations.RungOneLineage FusedAlgebra.Topics.ShapeDependencies FusedAlgebra.Topics.BooleanLattice FusedAlgebra.Topics.BeyondOneLoop FusedAlgebra.Topics.CornerWeights FusedAlgebra.Topics.Precession FusedAlgebra.Topics.Hypersphere FusedAlgebra.Topics.QuantumAction FusedAlgebra.Topics.HoloGluing FusedAlgebra.Topics.HoloGluingField FusedAlgebra.Topics.FuzzyFrame FusedAlgebra.Topics.FreshLook FusedAlgebra.Topics.CutIndependence FusedAlgebra.Topics.YangMillsGap FusedAlgebra.Topics.HyperLadder FusedAlgebra.Topics.Riemann FusedAlgebra.Topics.RiemannElliptic FusedAlgebra.Topics.ForceReach FusedAlgebra.Topics.FluxTubeDial FusedAlgebra.Topics.Tunnelling FusedAlgebra.Topics.EllipticModulus FusedAlgebra.Topics.NewForce FusedAlgebra.Topics.RiemannRound3 FusedAlgebra.Topics.RiemannRound4 FusedAlgebra.Topics.RiemannRound5 FusedAlgebra.Topics.DialAudit FusedAlgebra.Topics.RiemannRound6 FusedAlgebra.Topics.RiemannNeighbours FusedAlgebra.Topics.RiemannRound7 FusedAlgebra.Topics.RiemannRound8 FusedAlgebra.Topics.RiemannRound9 FusedAlgebra.Topics.EllipticPendulum FusedAlgebra.Topics.ChaosComposition FusedAlgebra.Topics.RiemannUnified FusedAlgebra.Topics.RiemannRound10 FusedAlgebra.Topics.EntropyKnots FusedAlgebra.Deprecated
  FusedAlgebra Check
)
export LEAN_PATH="$PWD/build"
for m in "${mods[@]}"; do
  f="$(echo "$m" | tr . /).lean"; o="build/$(echo "$m" | tr . /).olean"
  mkdir -p "$(dirname "$o")"
  echo "== $m"
  lean -o "$o" "$f"
done
echo "All modules compiled."
