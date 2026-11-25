# General relativity through the circle operation (Nov 25, 2026)
Page: https://claude.ai/artifact/SHrK88RpTBLLfTaeAxgi37 (source index.html)

Question (Brennan, Nov 25): special relativity's velocity addition is the circle operation ∘₁ (called mult-add before 1.0); can we get from the
circle-operation geometry to Einstein's field equations G = 8πG T, or find the exact step where it fails?

## Answer

Not from the circle operation alone. The circle operation gives the local rules of motion (the dial q sets the signature and
the speed of light, and twist = curvature × area). Einstein's equations need a rule for how matter
bends spacetime, and every route tested has to supply that rule from outside. What we found is the
exact place each route needs it, and one route that comes close in a picture that is ours to use.

## Route (a): let the dial vary from point to point

- A dial field k(x) = 1/c(x)² changes the local speed of light. On its own it only touches clocks (or
  clocks and rulers equally, the conformal case). That is a scalar theory of gravity (Einstein 1911–12,
  Nordström 1913).
- Exact geodesic integration at the Sun and for Mercury (checks.py, section B):

| Metric | γ | β | Light bending | Mercury perihelion |
|---|---|---|---|---|
| GR (Schwarzschild) | 1 | 1 | 1.751″ | 42.98″/century |
| dial in time only, A = 1 − 2U | 0 | 0 | 0.876″ | 28.65″ |
| dial in time only, A = e^(−2U) | 0 | 1 | 0.876″ | 14.33″ |
| conformal dial (Nordström) | −1 | 1 | 0 | −14.33″ |

- Every dial-only theory bends light by exactly half, whatever its β. Cassini (γ − 1 = (2.1 ± 2.3)·10⁻⁵)
  puts γ = 0 about 43,000 σ away.
- **Where it fails:** space has to curve as much as time slows. The dial alone can't do that.

## Route (c): the river (weak field and beyond)

Two postulates, both stated in circle-operation language:
1. At each event there is a free-fall frame (the river) in which space is flat and special relativity
   holds, so speeds compose by the circle operation on the Doppler-squared ruler.
2. The river flows past static observers at Newton's escape speed, v² = 2m/r (r = areal radius).

Results (checks.py, sections C to C3; all exact, sympy):
- The metric −dt² + (dr + v dt)² + r²dΩ² with v² = 2M/r has zero Einstein tensor and is Schwarzschild
  after the clock change t_S = t − ∫ v/(1 − v²) dr. River observers fall freely and are at rest at infinity.
- Static clocks and static radial rulers are boosted through the river, and one boost gives both: clock
  factor 1 − v² = −g_tt, ruler factor 1/(1 − v²) = g_rr. So the time deficit and the radial-space deficit
  are circle inverses, (−g_tt − 1) ∘₁ (g_rr − 1) = 0. This is why γ = 1 and why light bends twice the
  Newtonian amount: half from clocks, half from rulers, because a boost dilates time and contracts length
  by the same γ.
- Weak field: −g_tt = 1 − 2U + 2U² (β = 1), g_ij = (1 + 2U)δ_ij (γ = 1). So light bending 1.751″,
  perihelion 42.98″/century, Shapiro delay and redshift all follow.
- The horizon r = 2M is where the river reaches light speed: the wall of the velocity law.
- With dust, the spherical Einstein equations are exactly Newton's: G^t_t gives the shell theorem
  m′ = 4πr²ρ with v²/2 = m/r, G^r_r gives "m is carried along by the flow", and G^r_t and G^θ_θ add
  nothing new. The Hubble river v = −Hr gives Friedmann's equations (H² = 8πρ/3 is Newton's escape speed
  for a uniform ball).
- Gravity Probe B's exact factor √(1 − 3M/r) is 1/γ of the gyroscope's speed through the river (river
  radial, orbit tangential, composed at right angles): −w² = (−v_r²) ∘₁ (−v_o²). This replaces the earlier
  twisted-boosts bookkeeping, which took Schwarzschild's space curvature as an input; the river supplies it.
- Known: Painlevé 1921, Gullstrand 1922, Milne–McCrea 1934, Hamilton & Lisle 2008 ("river model").

**Where it stops:** gravitational waves. A plane wave has real tidal curvature (R^x_txt = −h″/2), and its
transverse-traceless part can't be removed by any change of coordinates, so no slicing makes space flat.
A river in flat space can't carry the two polarizations, and the data need them (Hulse–Taylor orbital
decay 0.9983 ± 0.0016 of GR's prediction; LIGO). Pressure was not tested here.

## Route (b): holonomy, the full curvature and the Bianchi identity

- A small loop in the t–r plane outside a mass returns a vector boosted, with rapidity (2M/r³) × area
  (numerical parallel transport, 0.5% at the loop size used). Curvature in a plane that contains time
  gives a boost; that is the twisted-boosts law with q < 0.
- The Bianchi identity ∇·G = 0 holds automatically (checked on every river v(t, r)). The contracted
  Bianchi ∇·Ric = ½ dR forces the left side to be R − ½Rg (b = −a/2). Lovelock (1971): in four dimensions
  that, plus Λg, is the only choice for a second-order local law. A metric field and "second order,
  local" are inputs put in by hand.
- Section A: the dial's boosts close on −q × a rotation and commute at q = 0. Hojman, Kuchař and
  Teitelboim (1976) showed that if these boosts are made local (a lapse that varies from place to place)
  and must compose the same way whichever path is taken, the only second-order dynamics for a spatial
  metric is Einstein's, with G and Λ left free. The coefficient in their bracket is ε = sign(q). This is
  the closest thing to a derivation from the composition law: the twist rule, made local, plus a metric.

## Leads, judged on their own merits

- **A derivation chain worth writing up, not new physics.** Dial boosts → twist rule (ours) → made local
  → HKT → Einstein's equations. Every link is proved, ours or known; the result is GR exactly.
- **No new prediction.** Every route either reproduces GR or is ruled out by data. GW170817 also says the
  speed of gravity equals the speed of light to about 10⁻¹⁵, so light and gravity share one dial.
- **The one place a varying dial goes beyond GR:** letting q vary inside the HKT bracket ("deformed
  general covariance"). This is active research: in loop quantum cosmology the bracket's coefficient turns
  negative above half the critical density, which reads as a switch from Lorentzian to Euclidean
  signature (Mielczarek, arXiv:1304.2208; arXiv:1606.03271; disputed with self-dual variables,
  arXiv:1610.07467). Bojowald and Duque's emergent modified gravity builds consistent theories where that
  coefficient is a free function (arXiv:2310.06798, arXiv:2404.06375, arXiv:2407.13583). It is a lead
  worth following as theory if the dial-as-field idea matters to the project; its effects appear only near
  Planck density or through free functions, so there is no data test today.

## Files
- checks.py → checks_output.txt, results.json: 44 of 44 checks pass (sympy, mpmath, scipy, numpy).
- grtools.py: small sympy helper (Christoffels, Riemann, Ricci, Einstein, divergence).
- lean/GeneralRelativity.lean: 28 theorems, core Lean 4.37.0-rc1, no sorry; compiled here and merged into
  the shared library as FusedAlgebra/Topics/GeneralRelativity (library builds clean, 2,715 theorems).
