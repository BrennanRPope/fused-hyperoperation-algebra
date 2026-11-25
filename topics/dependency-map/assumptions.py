"""Assumptions of each write-up on the quantum-gravity path, drawn from its own text, its checks and
its Lean hypotheses. kind: K = rests on a known published result, O = a step of ours,
in = an input put in by hand. 'stated' = the write-up states it; 'implicit' = it uses it without saying."""

CATS = [
 ("sec_a", "Has an Assumptions section"),
 ("sec_k", "Has a Known vs new section"),
 ("dial", "Dial value set by hand"),
 ("dyn", "Dynamics (H, action, V) as input"),
 ("hbar", "Size of ħ as input"),
 ("grav", "Metric, 2nd order or source law as input"),
 ("ruler", "Leans on a 1D match (ruler change)"),
 ("limit", "Infinite limit cited, not proved"),
 ("toy", "Tested on a restricted case"),
 ("key", "Key step: published or ours"),
 ("lean", "Lean backs its identities"),
]

W = []
def add(name, folder, checks, lean, assumptions, known, ours, cats):
    W.append(dict(name=name, folder=folder, checks=checks, lean=lean, assumptions=assumptions,
                  known=known, ours=ours, cats=cats))

add("order dial", "order-dial", "71/71", "OrderDial, 33 theorems (Int); centre_fixed assumes r·x = -d",
 [("Maps are x ↦ (1+a)x + c with integer or real coefficients", "K", "implicit"),
  ("The circle operation is composition of stretches sharing the centre −1/k (T_comp)", "K", "stated"),
  ("The route to quantum mechanics uses the Heisenberg group on functions; i and ħ are put in by hand", "in", "stated"),
  ("Category-theory framing: a monoid is a one-object category", "K", "stated")],
 "affine group, Heisenberg group, Weyl relation, monoids as categories",
 "the order gap read as the cost of order; the circle operation as the single-centre slice",
 dict(sec_a="no", sec_k="no", dial="input", key="known", lean="yes"))

add("why i (+ path integrals, weak)", "why-i", "44/44", "WhyI, 21 theorems (Int); stands alone, imports nothing",
 [("Real coefficients; half-plane uses orientation-keeping maps (1 + ka > 0)", "K", "stated"),
  ("'Length' = left-invariant metric; 'keeps size' = unitary representation; continuous characters", "K", "stated"),
  ("i needs the twist plus one more ingredient (length, size-keeping, or a loop)", "K", "stated"),
  ("Path integrals: finite matrices and lattices; unbounded H rests on Chernoff's theorem (cited)", "K", "stated"),
  ("Hamiltonian, size of ħ and a real-time path measure are not supplied", "in", "stated"),
  ("Wall-flips-i: 4096-point FFT, 'keeps sign' = under 1e-3 of weight crosses zero", "O", "stated")],
 "Kähler triples, Stone–von Neumann, Kostant–Souriau, Hestenes, Stueckelberg, Volterra, Chen, Feynman–Kac, Wick",
 "half-plane picture of the dials; wall flips i; Wick rotation read as a turn of k",
 dict(sec_a="yes", sec_k="no", dial="input", dyn="input", hbar="input", limit="known", toy="known", key="known", lean="yes"))

add("Cayley–Klein", "why-i/cayley-klein", "13/13", "none",
 [("Real Cayley–Klein geometries; imaginary k only as a complex coordinate for the elliptic case", "K", "stated"),
  ("Sign convention of Herranz–Ortega–Santander (κ₁ > 0 is anti-de Sitter)", "K", "stated"),
  ("Cayley map normalized so a ≈ T for small T; only the sign of κ is invariant", "O", "stated"),
  ("'No free label per pair' uses motions acting on one space, not abstract graded contractions", "K", "stated")],
 "Yaglom, Herranz et al., Gromov, Lévy-Leblond, Bacry–Lévy-Leblond, Coleman–McGuinness",
 "κ = −k²/4 via the Cayley map; (k, d) collapse to κ = d − k²/4",
 dict(sec_a="yes", sec_k="no", ruler="ours", toy="known", key="ours", lean="no"))

add("dial map (Bronstein cube)", "why-i/dial-map", "15/15", "none",
 [("Stability test = H²(g, g) = 0 for finite-dimensional Lie algebras", "K", "stated"),
  ("Quantum case: a 3D stand-in; Moyal quoted, not computed", "K", "stated"),
  ("3+1 algebras built with labels (κ₁, κ₂, 1, 1): de Sitter family only", "O", "stated"),
  ("Newton binding (G1) is a first-order reading at fixed separation", "O", "stated"),
  ("Bootstrap toy (G3) is scalar, not Deser's argument", "O", "stated")],
 "Inönü–Wigner, Lévy-Nahas, Bacry–Lévy-Leblond, Vey, BFFLS, Mendes, Lovelock, Deser, MacDowell–Mansouri",
 "the dial lattice; 'stable dials need non-commuting elements'; the Newton-binding reading",
 dict(sec_a="yes", sec_k="no", ruler="ours", toy="ours", key="known", lean="no"))

add("κ-Poincaré", "why-i/kappa-poincare", "10/10", "none",
 [("Majid–Ruegg bicrossproduct basis; other bases give other formulas", "K", "stated"),
  ("Massless, collinear, positive branch", "K", "stated"),
  ("Magueijo–Smolin via the Judes–Visser U-map law (debated at the time)", "K", "stated"),
  ("Torsion and curvature read numerically at the origin (relative-locality dictionary)", "K", "stated"),
  ("κ set to the Planck energy", "in", "implicit")],
 "κ-Poincaré, DSR, de Sitter momentum space, relative locality",
 "collinear law written as a + b + kab with k = −1/κ; Magueijo–Smolin on the merged-wall line",
 dict(sec_a="yes", sec_k="yes", dial="input", ruler="known", key="known", lean="no"))

add("e, i, π", "why-i/e-i-pi", "42/42", "EIPi, 17 theorems (Int, Gaussian integers); stands alone",
 [("Dial and lift taken from order-dial and Cayley–Klein, not re-derived", "O", "stated"),
  ("'Size-keeping' = lift of modulus 1; its physical reading (probability kept) is interpretation", "O", "stated"),
  ("Parameter counts K(n) = n², n(n+1)/2, n(2n−1)", "K", "stated"),
  ("Picking complex numbers needs Hardy / Masanes–Müller axioms", "in", "stated"),
  ("Physics values from CODATA as floats: formulas checked, not experiments", "K", "stated")],
 "Euler, product integrals, local tomography, KMS, Unruh, Hawking, Bisognano–Wichmann",
 "e, i, π as dial directions; K_AB = K_A K_B as the circle operation; the dial ledger; caustic/conformal-factor parallel (untested)",
 dict(sec_a="yes", sec_k="yes", dial="input", hbar="input", key="ours", lean="yes"))

add("chain resolves", "why-i/e-i-pi/chain-resolves", "18/18", "none",
 [("Gaussian actions only (free particle, oscillator, two straight routes)", "O", "stated"),
  ("'Resolving' = ħ → 0 picks a classical path; 'reversible' = inverse of the same kind", "O", "stated"),
  ("Dial applied as e^{κS/ħ} with κ on the unit circle; ħ's size not tested", "in", "stated"),
  ("Checkerboard: turn rate 1, step 0.05, 401 sites", "O", "stated")],
 "stationary phase, Laplace, Morse index, Maslov, Hopf–Cole, Crandall–Lions, Kac",
 "'allows vs forces' framing; dial-plane picture of allowed κ",
 dict(sec_a="yes", sec_k="yes", dyn="input", hbar="input", toy="ours", key="known", lean="no"))

add("checkerboard", "why-i/checkerboard", "24/24", "Checkerboard, 4 theorems (Int); stands alone",
 [("c = ħ = 1, equal space and time steps", "K", "stated"),
  ("Tick = turn then move; start facing right", "O", "stated"),
  ("Dial multiplies only the turns: a continuation in mass, not a Wick rotation", "K", "stated"),
  ("Finite ring larger than the light cone", "O", "stated")],
 "Feynman–Hibbs, Kac telegraph, GJKS 1984, Ichinose, Meyer, Bialynicki-Birula",
 "writing both models as one weighted-circle chain (exact identity B1)",
 dict(sec_a="yes", sec_k="no", dial="input", toy="known", limit="known", key="known", lean="yes"))

add("checkerboard + potential", "why-i/checkerboard-potential", "17/17", "CheckerboardPotential, 3 theorems (Int); stands alone",
 [("1+1D, c = ħ = 1", "K", "stated"),
  ("Real, time-independent scalar potential", "in", "stated"),
  ("Limit proofs for bounded smooth V; unbounded V checked numerically in a box", "K", "stated"),
  ("Random setting: V ≥ 0 keeps a sub-probability", "K", "stated")],
 "Feynman–Hibbs 2-6, Ichinose–Tamura, Cedzich et al., Feynman–Kac, Trotter",
 "the potential as one chain element; the chain identity with potential",
 dict(sec_a="yes", sec_k="no", dyn="input", toy="known", limit="known", key="known", lean="yes"))

add("path integrals (strong)", "why-i/path-integral-strong", "28/28", "none (uses WhyI, uncompiled when written; now compiled)",
 [("'Strong sense' = each of H, ħ, measure from a structure with a published theorem (own definition)", "O", "stated"),
  ("Curvature trade assumes a closed genus-g quotient, which needs a lattice in SL(2,R)", "O", "stated"),
  ("Klauder's curvature read as scalar curvature", "K", "stated"),
  ("Numerics: 1200-point grid, Monte Carlo at 4 standard errors", "O", "stated")],
 "Bargmann, Lévy-Leblond, Groenewold–van Hove, Kostant–Souriau, Cameron, Nelson, Daubechies–Klauder–Paul",
 "the curvature trade (ħ cancels, only β fixed); wedge on exponents",
 dict(sec_a="yes", sec_k="no", dyn="input", hbar="input", limit="known", key="known", lean="no"))

add("special relativity", "relativity", "37/37", "Relativity, 30 theorems (Rat); hypotheses: 1 − u ≠ 0, 1 + uv ≠ 0",
 [("Velocities in units of c; the ruler d = 2u/(1 − u) is chosen", "O", "implicit"),
  ("Collinear motion only", "K", "implicit"),
  ("Hypotheses keep away from the walls (u ≠ 1, uv ≠ −1)", "K", "implicit")],
 "Einstein velocity addition, rapidity, Doppler",
 "velocity addition as the circle operation on the Doppler-squared ruler; landscape placement",
 dict(sec_a="no", sec_k="no", ruler="ours", key="known", lean="yes"))

add("twisted boosts (+ GP-B, pulsar)", "twisted-boosts", "36/36 + 21/21 + 13/13", "TwistedBoosts 36, GravityProbeB 16, BinaryPulsar 10 (Rat)",
 [("Half-velocities t = tanh(φ/2); dial boost B_q(u)", "O", "implicit"),
  ("GP-B uses exact Schwarzschild transport: GR's space curvature taken as input in the first pass", "in", "stated"),
  ("Exact velocity factor 1/γ vs Thomas's γ − 1: equal to first order only", "K", "stated"),
  ("Pulsar: masses and orbit from Breton et al. 2008", "K", "implicit")],
 "Thomas 1926, Wigner 1939, Borel, Ungar gyrogroups",
 "landscape q = curvature; twist = phase of the split corner; twist bookkeeping of precession",
 dict(sec_a="no", sec_k="no", grav="input", key="known", lean="yes"))

add("general relativity (river)", "general-relativity", "44/44", "GeneralRelativity, 28 theorems (Rat); hypotheses: 1 − 2x ≠ 0",
 [("Postulate 1: a free-fall frame at every event (equivalence principle)", "in", "stated"),
  ("Postulate 2: the river flows at Newton's escape speed v² = 2GM/r", "in", "stated"),
  ("Spherical symmetry, dust; pressure not tested", "O", "stated"),
  ("Metric field and 'second order, local' put in by hand (route b)", "in", "stated")],
 "Painlevé–Gullstrand, Hamilton–Lisle river, Lovelock, Hojman–Kuchař–Teitelboim",
 "clock and ruler deficits as circle-operation inverses; GP-B factor as river composition",
 dict(sec_a="no", sec_k="no", grav="input", toy="ours", key="known", lean="yes"))

add("GR bridges", "gr-bridges", "65/65", "GRBridges, 22 theorems (Int, Rat)",
 [("From 3D on, Schur forces a direction-dependent dial = a metric", "K", "stated"),
  ("Jacobson: entropy ∝ area with η = 1/(4Għ) is the outside input", "in", "stated"),
  ("Lovelock: second order is an input; river as gauge, not field", "in", "stated"),
  ("Deser: massless spin-2 needed; traded for the observed fact (LIGO)", "in", "stated")],
 "Jacobson 1995, TEGR, Lovelock, Deser, Brans–Dicke, Einstein–aether, Lasenby–Doran–Gull",
 "Raychaudhuri focusing as the circle-operation cross term; ADM constraints as e₂ of strain rates",
 dict(sec_a="no", sec_k="no", grav="input", key="known", lean="yes"))

add("electromagnetism", "electromagnetism", "49/49", "Electromagnetism, 41 theorems (Rat)",
 [("Wave equation kept under a dial boost only at d = −1/c²", "K", "implicit"),
  ("Source law ∇F = J taken from outside", "in", "stated"),
  ("Vacuum superposition is plain addition (k = 0)", "K", "stated")],
 "Maxwell, Riemann–Silberstein, Hestenes 1966, Wilson loops",
 "gauge phase and boost both obey a·conj(a) = 0; EM vs gravity = commuting vs non-commuting",
 dict(sec_a="no", sec_k="no", dyn="input", key="known", lean="yes"))

WRITEUPS = W
