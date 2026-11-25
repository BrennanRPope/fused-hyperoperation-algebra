#!/bin/bash
# Compile the whole library with core Lean (no mathlib, no lake needed).
# Usage: ./build.sh            (expects `lean` on PATH, version in lean-toolchain)
set -e
cd "$(dirname "$0")"
mods=(
  FusedAlgebra.Base FusedAlgebra.Rungs FusedAlgebra.Cube
  FusedAlgebra.Foundations.MultAddRing FusedAlgebra.Foundations.Independence FusedAlgebra.Foundations.FusedPeano
  FusedAlgebra.Topics.ClosingRate FusedAlgebra.Topics.Consensus FusedAlgebra.Topics.DarkSector
  FusedAlgebra.Topics.DataCenter FusedAlgebra.Topics.E FusedAlgebra.Topics.Finance
  FusedAlgebra.Topics.FusedConstants FusedAlgebra.Topics.FractionalRungs FusedAlgebra.Topics.Infinities FusedAlgebra.Topics.InformationTheory
  FusedAlgebra.Topics.Landscape FusedAlgebra.Topics.RungShape FusedAlgebra.Topics.Elliptic FusedAlgebra.Topics.EllipticS3
  FusedAlgebra.Topics.Thermodynamics
  FusedAlgebra.Topics.NumberPlane FusedAlgebra.Topics.GeometryTrig FusedAlgebra.Topics.EarlyScience FusedAlgebra.Topics.Pi FusedAlgebra.Topics.Calculus FusedAlgebra.Topics.AbstractAlgebra FusedAlgebra.Topics.Statistics
  FusedAlgebra.Topics.NavierStokes FusedAlgebra.Topics.QuantumNoise FusedAlgebra.Topics.Rulecraft
  FusedAlgebra.Topics.Toric.ToricCharts FusedAlgebra.Topics.Toric.ProjectivePlane
  FusedAlgebra.Topics.Topology.CubeProofs FusedAlgebra.Topics.Topology.KummerWilson
  FusedAlgebra.Topics.Topology.HopfInvariantOne
  FusedAlgebra.Topics.NSContinuum FusedAlgebra.Topics.JetFields FusedAlgebra.Topics.OsgoodInverse FusedAlgebra.Topics.LinearAlgebra FusedAlgebra.Topics.DialSweep FusedAlgebra.Topics.EntropyNoise
  FusedAlgebra.Topics.QAcross.QDial FusedAlgebra.Topics.QAcross.Knots FusedAlgebra.Topics.QAcross.Partitions FusedAlgebra.Topics.QAcross.QCounting
  FusedAlgebra.Topics.Fermat.Fermat4 FusedAlgebra.Topics.Fermat.Fermat3 FusedAlgebra.Topics.Fermat.FreyPieces FusedAlgebra.Topics.Fermat.FermatRoute
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
