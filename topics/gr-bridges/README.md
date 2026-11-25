# Curved space from the circle operation: four bridges to Einstein
Page: https://claude.ai/artifact/LhAbnuPMJnNoq9cRYBz3xc (source index.html, pictures in figures/)

Question (Brennan, Nov 25, voice): can we crack Einstein's field equations? Lead question added on the call:
can the circle operation's own rules build a metric with freely varying curvature (pure math first)? Then the direction-dependent
dial ("the promising lead"), k as a Clifford-style fused number, and four known bridges: Jacobson, teleparallel,
Lovelock, Deser (with Brans-Dicke, Cassini and the sombrero).

## The answer in plain words

- **The circle operation ∘₁ (called mult-add before 1.0) on its own builds one fixed curvature, not a varying one.** Its two-direction rule (boosts close on a
  twist = q x area) makes relativity's velocity space, a curved space with one fixed curvature.
- **A dial that varies only by place works in a 2D world and is blocked in 3D and 4D.** In 2D any curvature
  q(x, y) is reachable. From 3D on, Schur's theorem says a dial that is the same in every direction at each point
  must be the same everywhere. That is the exact step where "its own rules" stop.
- **So the dial has to depend on direction, which confirms the "add to the circle operation" call.** The smallest addition that
  works is a dial on each direction pair: the square rule of a Clifford (geometric) algebra, e_i e_j + e_j e_i = 2 k_ij.
  That k_ij is exactly a metric. It is a renaming of g_ij, plus a reason the metric has to be there.
- **k as a Clifford number works only in the sandwich form a + b + a·k·b.** That law is associative for every
  multivector k and keeps the exact lift 1 + k(a∘₁b) = (1 + k a)(1 + k b). A weight field k(x) used alone makes no
  curvature (it is pure gauge). Add a twist field and you get gauge-theory gravity (Lasenby, Doran and Gull), which is
  general relativity, and our river is its natural ("Newtonian") gauge.
- **The bridges:** every bridge reaches Einstein's equations once the direction dial is in. The circle operation supplies the most
  in Jacobson's (the area bookkeeping and the Unruh temperature). None gives a new prediction.

## Part S. Can the circle operation build curved geometry? (pure math)

| Step | What it gives | Check |
|---|---|---|
| One direction (S1) | In Cl(1,1), boosts are 1 + βσ (σ² = +1). Products give Einstein's velocity addition; the two light-cone parts multiply, so D² − 1 combines by the circle operation; the norm 1 − β² is the wall (light speed). | sympy, Lean |
| Two directions (S2) | [e0e1, e0e2] = −2 k0 e1e2: the twisted-boost law is the Clifford commutator, and q is the square of the shared direction. k0 = 0 (Galilean) gives no twist. | sympy |
| The circle operation's own rules (S3) | Velocity space has one fixed curvature; the Thomas–Wigner twist of a velocity loop equals its angle defect (area x curvature) to 1e-15. | numeric |
| Dial by place, 2D (S4) | Metric e^{2u}: curvature q = −e^{−2u}Δu (Liouville), any q reachable locally; loop twist / area matched q at three places with three values. | sympy, numeric |
| Dial by place, 3D+ (S5) | Isotropic curvature gives G = −(n−1)(n−2)/2 · q g, so the Bianchi identity forces dq = 0 for n ≥ 3 (Schur). Explicit 3D example: making all planes agree forces a constant curvature. | sympy, Lean |
| Multivector weight (S6) | a + b + a·k·b: exact lift and associativity for any k in Cl(3,0). a + b + k·a·b fails unless k is central (scalars, plus the volume element I in Cl(3,0); scalars only in Cl(1,3)). A constant invertible k is a change of coordinates. | sympy (all 8 parts), Lean (2×2 matrices) |
| Weight field k(x) (S7) | Changing weight from place to place is the connection A = k⁻¹dk, whose curvature is zero identically. Curvature needs a second field. | sympy |
| Direction dial (S8) | k_ij(x) in the Clifford square rule is a metric. Sectional curvatures differ by plane at one point (0.00325, 0.01138, −0.00345 in the example), and a numeric loop matches. | sympy, numeric |
| Which k_ij does matter make? (S9) | Weak field: in vacuum light focuses in some direction unless γ = 1. So "no focusing without energy crossing" (Jacobson's Clausius rule) forces γ = 1; the dial-only theory (γ = 0) breaks it. | sympy |

**Is the direction dial "just renaming g_ij"?** Yes, as an object. What it adds: (1) the reason the metric must appear
(Schur blocks every scalar dial), (2) a Clifford frame (the square root of the metric), which is what spinors, tetrads,
teleparallel gravity and gauge-theory gravity need, (3) the twisted-boost law falls out as the bivector commutator.
**Can k itself be a fused number?** Yes: with k = 1 + K, weights compose by the same multivector circle operation, so the
structure is self-similar, and the geometric product keeps both the dot (k-weighted) and the wedge (twist). But a weight
field alone is flat (S7), so self-similarity does not create curvature by itself.

## Part J. Jacobson 1995: Einstein's equations from heat flow (δQ = T dS) across every small horizon

| Ingredient the bridge needs | Supplied by the circle operation? |
|---|---|
| Area change of a light bundle (Raychaudhuri) | **Yes.** Area = (1+a)(1+b) = 1 + (a∘₁b). Exactly: θ′ + θ²/2 + σ² = a″/(1+a) + b″/(1+b) = −R_kk. With straight rays the area's second difference is 2αβ, the circle-operation cross term; that cross term is the focusing (J1). |
| Unruh temperature T = a/2π | **Yes.** Constant acceleration runs the Doppler factor along the circle-operation flow e^{aτ}; its Fourier transform is thermal, |β/α|² = e^{−2πΩ/a}. The 2π is the rotation period after flipping the dial's sign (J2). |
| "Heat flux = focusing in every light direction ⇒ Einstein" | Algebra, checked: X(k,k) = 0 for all light-like k forces X = φg; Bianchi fixes φ (J3). |
| Entropy proportional to area, with η = 1/(4Għ) | **No.** This is the one outside input (quantum entanglement / Bekenstein). |
| Matter's energy flux T_ab | No (matter). |

Extra links found:
- Our earlier "clock deficit ∘₁ ruler deficit = 0" is exactly Jacobson's 2007 condition g_tt g_rr = −1, which holds iff
  no energy flows along radial light rays (J4).
- If horizon patches combine entropy by the circle operation with dial k (that is Tsallis's rule), the entropy that adds is
  log(1+kS)/k (Rényi). Small patches can't see k, so the local Einstein equations are unchanged. Globally black holes
  change: Rényi temperature (1 + 4πλM²)/(8πM), stable above M = 1/√(4πλ) (Biró & Czinner 2013) (J5).
- Sombrero role: a sombrero field coupled as ξφ²R could set η (Newton's G) through its vev (induced gravity, Zee 1979),
  but the Higgs vev would need ξ ≈ 1e32 (J6).

**Verdict:** Jacobson's bridge is where the circle operation carries the most: all three geometric pieces (area bookkeeping, focusing,
temperature) are the circle operation or dial statements. Gap: the area law for entropy. This is a derivation chain worth writing up;
it gives GR exactly, no new prediction. One lead worth following on its own merits: whether the circle operation's
information-keeping can produce the area law (entropy per area). We did not find it, and the Tsallis check shows the
combination rule alone can't change the local equations.

## Part T. Teleparallel (torsion instead of curvature)

- The river is natively teleparallel. Frame e⁰ = dt, eⁱ = dxⁱ + vⁱ dt: the teleparallel connection is flat (no twist
  at all), and the torsion is the river's velocity gradient, Tⁱ_jt = ∂_j vⁱ (T1).
- On **any** river (not only spherical), checked symbolically with general v(t,x,y,z) (T2, T3):
  - energy constraint G_nn = e2(K) = 8πρ, where e2 is the sum of pairwise products of the river's three strain rates:
    the cross-term layer of a three-way the circle operation, (1+k1)(1+k2)(1+k3) = 1 + e1 + e2 + e3;
  - momentum constraint ∂_i(K_ij − δ_ij K) = −½ curl curl v, so moving mass makes the river swirl (an exact
    gravitomagnetic Ampère law);
  - e2(K) = (a divergence) − |curl v|²/4: only the swirl carries bulk action;
  - the teleparallel Lagrangian (torsion scalar) on the river is exactly 2·e2(K), and R = −T − 2∂T: TEGR holds.
- Twist versus torsion (T4): torsion = v′; Newton's pull = v v′; the curvature twist rate (v v′)′ = 2M/r³ is the
  holonomy found before. Twist (rotation defect) and torsion (translation defect) are the two halves of one loop
  defect; TEGR moves all of it into torsion. So twist is not torsion, but the river's torsion carries the same information.
- **Where it stops:** a river on flat slices has no transverse-traceless part (D1), so it can't carry gravitational
  waves. The general teleparallel frame (any tetrad, Hamilton–Lisle's "twisting river" for Kerr) can, and that is TEGR = GR.

**Verdict:** a clean, exact rewrite: gravity's energy rule = the circle-operation cross terms of the river's strain. Known in
substance (ADM constraints in flat slicing); the circle-operation reading is ours. No new prediction.

## Part L. Lovelock 1971: in 4D the only local, second-order, conserved metric law is Einstein + Λ

| Hypothesis | Does our construction meet it? |
|---|---|
| 4 dimensions | Yes. Gauss–Bonnet's field tensor vanishes identically in 4D and not in 5D (H_tt = 6(f−1)f f′/r³), checked (L1). |
| Conserved source | Yes: matter carried by the river is conserved; Bianchi holds on every river. |
| Second order | An input (a simplicity assumption), not from the circle operation. |
| Metric only | **Yes if the river is a gauge choice** (then Lovelock forces Einstein + Λ). **No if the river is a real field**: then the general theory is Einstein–aether. Exact Schwarzschild with a river aether needs c₂ = −c₁₃; GW170817 pins c₁₃ ≲ 1e-15, so the field is pinned to the gauge; c₁₄ stays open, bounded by solar-system tests at about 1e-5 (Oost, Mukohyama, Wang 2018) (L2). |
| Λ | Left free. The river carries it as a second escape term, v² = 2M/r + Λr²/3, exact Schwarzschild–de Sitter (L3). |

**Verdict:** uniqueness closes the gap once you accept "metric only, second order". With the direction dial, the
metric is ours, so Lovelock's hypotheses are met and Einstein + Λ is forced. The river-as-field option (Einstein–aether)
is a real, testable lead, but it is existing research and data already squeeze it.

## Part D. Deser 1970 bootstrap, and dial + river as Brans–Dicke

- Deser (and Weinberg 1964, Feynman) need a massless spin-2 field. The river is a vector, the place dial a scalar: no
  spin-2 (D1). Self-coupling a scalar gives Nordström's theory (Freund & Nambu 1968): the conformal dial, light bending 0,
  already ruled out. The direction dial's transverse-traceless part has exactly two polarizations (h+, hx): the spin-2 field.
  With it, the bootstrap gives Einstein's equations. So the missing rule can be traded for one observed fact: gravity's
  waves are massless spin-2 (LIGO polarizations, GW170817 speed).
- Dial touching curvature = Brans–Dicke: γ = (1+ω)/(2+ω). Cassini needs ω > about 40,000 (coupling α² < 1.2e-5) (D2).
- **The sombrero's role:** it gives the dial a mass. Above about 1e-17 eV (range under ~1/10 AU) the dial hides from
  Cassini at full strength; the Higgs-mass sombrero (range 1.6e-18 m) hides it completely. The light quintessence dial
  that the dial-potential search keeps alive (mass ~ H0 ~ 1e-33 eV) gets no shelter and must couple below α² ~ 1e-5.
  The sombrero does not supply the matter-bends-space rule in any bridge.

## Leads, judged on their own merits

- **Direction dial = metric, forced by Schur.** Worth writing up as pure math: "The circle operation's twisted-boost rule, made
  local, needs a direction-dependent dial, and that dial is a metric." True and geometric; not new physics.
- **Jacobson chain in circle-operation language.** Worth following: every geometric piece is ours; the one open piece is the
  area law. If information-keeping could give entropy ∝ area, that would be a real step. Nothing here produces it yet.
- **River = teleparallel frame; gravity's energy = strain cross terms; bulk action = swirl.** Exact and pretty; a good
  picture for Kerr (twisting river). Not a new prediction.
- **River as a physical field (Einstein–aether, c₁₄).** Testable, but known research, and GW170817 already removes most room.
- **No new prediction anywhere.** Every route with the direction dial gives GR exactly; every route without it is ruled out.

## Files (for the repo thread to push)
- `checks.py` → `checks_output.txt`, `results.json`: 65 of 65 checks pass (sympy, mpmath, scipy, numpy; about 2 min).
- `ga.py`: small Clifford-algebra helper (diagonal metric, symbolic).
- `lean/GRBridges.lean`: 22 theorems, core Lean 4.37.0-rc1, no sorry; merged as `FusedAlgebra/Topics/GRBridges`
  (library builds clean: 112 modules, 2,823 theorems; build.sh, FusedAlgebra.lean, README.md, theorems.tsv, THEOREMS.md updated).
- `figures/*.svg`: the pictures; `index.html`: the page.
- Reuses `../general-relativity/grtools.py`.
