# DRAFT · Tunnelling, worked up with our own tools

Nov 26, 2026. Brennan asked to work on the tunnelling lead that `closed-leads-fresh-look/` (F4) reopened, and to
explore it freely: "just trying out the model and maybe seeing if there's anything that comes of it".
Checks: `checks.py` → `checks_output.txt`, **37/37 pass** (numpy, scipy, mpmath, sympy; key laws polished to 30
digits). Pictures: `figures.py` → `figures/`. Lean: `lean/Tunnelling.lean`, 7 polynomial identities, compiled in the shared library as
`FusedAlgebra.Topics.Tunnelling` (as written).

## What came out, in one paragraph

A barrier is a **boost of our dial**: its reflection amplitude is twisted-boosts' half-velocity, with curvature −1.
Its lossless cousin, a two-mode mixer (beam splitter, Rabi pulse, polarizer), is the same thing on the sphere,
curvature +1. Once you see that, one interference of two paths gives **three arithmetics at once**: at a fixed phase
the amplitudes compose by our (k = 0, d = K) law and pick up a twist equal to half the triangle's area; averaged
over the phase, probabilities compose by ∘ₖ with **k = −2K**; and the log-average picks the stronger path exactly,
which is the tropical (max) arithmetic. On the sphere that last law has a **wall** at s₁ + s₂ = 1; in the hyperbolic
plane it has none, which is why "blocked" glues by ∘₋₁ exactly along any chain of barriers. Two more things fell
out: the mean transmission of a chain falls per barrier by **Legendre's elliptic integral K of the reflection
amplitude**, and **rewinds** work for single numbers but not for the full spread, which loses a comb of
components whose positions are all real.

---

## 1. One law, three curvatures

![Fig 1](figures/fig1-one-law-three-curvatures.png)

**Fig 1.** Two steps of lengths D₁, D₂ meeting at angle φ, on a sphere, a flat plane and the hyperbolic plane.
Average the third side over φ. The squared half-chord t = (1 − cos √K·D)/(2K) composes by ∘ₖ with k = −2K: ∘₋₂ on
the sphere, plain addition on the flat plane, ∘₂ in the hyperbolic plane.

- **A barrier is a dial boost (1.1, 1.2).** Every delta barrier's transfer matrix is B₋₁(u) times a phase, where
  B_q(u) = [[1, u], [−q ū, 1]] is the boost from `twisted-boosts/` and u is the reflection amplitude. So
  |u| = tanh(D/2) is a half-velocity, and the odds t = R/T are its flip, |u|²/(1 − |u|²).
- **The mean law, every curvature (4.1–4.4).** Averaging the law of cosines over the angle makes the cos φ term
  vanish, and what is left is ∘ₖ with k = −2K, exactly, for every K tested (−2, −1, −0.3, 0, 0.5, 1).
  - K = −1: t is the barrier's odds R/T and the law is ∘₂, the mean-resistance rule.
  - K = 0: t = D²/4 and the law is plain addition, the rule that variances of independent steps add.
  - K = +1: t is a mixer's switch probability s and the law is ∘₋₂. In bias form, 1 − 2s multiplies. That is the
    rule for independent bit flips, and the piling-up lemma of cryptanalysis. It is also the freely rotating chain
    of polymer physics, where bond correlations multiply.
- **Only the sign is physics (4.5).** Rescaling lengths changes K and k together, so at d = 0 only the sign of k
  carries meaning. That is F1's "k alone: only the sign survives", now with a reason: the sign is the curvature's.
- **Known versus ours.** The averages are classical (each field has its own: Anderson, Thouless, Abrahams and
  Fisher 1980 for barriers; Matsui 1993 for the piling-up lemma; Flory for the freely rotating chain). The barrier
  as an SU(1,1) element and a Lorentz-type boost is known (Monzón, Sánchez-Soto and coauthors, multilayer optics).
  **Ours:** one formula, k = −2K, that puts all four on the homothety dial, and the identification of the barrier
  with twisted-boosts' B_q at q = −1.

## 2. One interference, three arithmetics

![Fig 2](figures/fig2-three-arithmetics.png)

**Fig 2.** (a) Two paths with a random relative phase. The mean of |A + Be^(iφ)|² is |A|² + |B|² (green). The mean
of its log is exactly ln max(|A|², |B|²), with a corner (blue on the dashed line). (b) Two mixers: whichever of
stay–stay or switch–switch is likelier sets the typical stay probability; the two tie on the red wall s₁ + s₂ = 1.
(c) Two barriers over the loop phase: the pair's odds (black) and the three averaged laws.

- **The log-average is the max, exactly (5.1).** ⟨ln|A + Be^(iφ)|²⟩ = ln max(|A|², |B|²). This is Jensen's formula
  (1899). Read in our terms: averaging over the phase sends the quantum arithmetic (amplitudes) to the classical one
  (probabilities add, the mean), and the log of it to the tropical one (max). The semiring audit had the tropical
  arithmetic only as the cold end of a limit; here a random phase produces it with no limit at all.
- **Mixers: the typical law has a wall (5.2, 5.3).** ⟨ln(1 − s)⟩ = ln max((1 − s₁)(1 − s₂), s₁s₂). Below the wall
  s₁ + s₂ = 1 the stay probability multiplies, which is ∘₋₁ on s. Above it, the typical stay probability is the
  product of the two switch probabilities: two switches make a stay. The wall is where plain addition reaches 1,
  because (1 − s₁)(1 − s₂) − s₁s₂ = 1 − s₁ − s₂ (Lean `sphere_paths`).
- **Barriers: no wall, so ∘₁ is exact forever (5.5, 5.6).** The direct path always outweighs the double reflection,
  since (1 + t₁)(1 + t₂) > t₁t₂ (Lean `hyperbolic_direct_wins`). So ⟨ln(1 + t)⟩ adds at every barrier, for any
  earlier history, and "blocked" glues by ∘₋₁ exactly along any chain (checked on six different barriers). This
  is the cleanest answer yet to why the old T1 test found ∘₋₁ only approximately: it never averaged the phase.
- **The pair's mean transmission is ∘₀ in both worlds (6.1–6.3).** Barriers: 1/⟨T⟩ − 1 = t₁ + t₂. Mixers:
  ⟨1/(1 − s)⟩ = 1/|1 − (s₁ + s₂)|, which blows up on the same wall. One identity gives both (Lean `pair_identity`).
- **Known versus ours.** Jensen's formula and the phase averages are known. **Ours:** reading the three averages as
  the three arithmetics of one interference, and the wall on the sphere as the reason the hyperbolic law is exact.
  We did not find the mixer wall stated as such, but it is a two-line consequence of Jensen, so it may well be in
  print somewhere.

## 3. Fixed phase: the edges read d, the twist reads the holonomy

![Fig 3](figures/fig3-edges-and-twist.png)

**Fig 3.** (a) The pair's reflection amplitude as the spacing changes. Its top edge is (r₁ + r₂)/(1 + r₁r₂), our
(k = 0, d = −1) law, relativity's velocity addition; the bottom edge is the same law run in reverse. (b) The pair's
transmission phase minus each barrier's phase and the propagation phase (black) sits on ± half the area of the
barrier triangle (red).

- **Edges (2.1–2.4).** Barriers: the extreme amplitudes are the (k = 0, d = −1) law on |r| = tanh(D/2). Mixers: the
  extreme tan(θ/2) is the (k = 0, d = +1) law, tangent addition. In both, the curvature reading −Δ/4 = d equals the
  curvature of the geometry, sign included.
- **Twist (3.1, 3.2).** The extra transmission phase is ± half the hyperbolic triangle's area for barriers and ± half
  the spherical triangle's area for mixers, to 1e-14. That is twisted-boosts' holonomy ω = q × area, with the phase
  of X + iY equal to ω/2.
- **Known versus ours.** Known: Wigner's angle as a geometric phase in multilayer stacks (Monzón, Sánchez-Soto and
  coauthors, around 2001), Pancharatnam's phase on the sphere, and Ungar's gyrogroups. **Ours:** that the same dial
  boost B_q covers relativity, multilayers and beam splitters, with q read as the curvature.
- **One honest wrinkle for the design rule.** The edge laws read the curvature through d (K = d). The averaged laws
  sit at d = 0 and read it through k (K = −k/2). The reading −Δ/4 of an averaged law is −k²/4, which is −1 for both
  ∘₂ and ∘₋₂. So for the averaged laws the reading has the right size but loses the sign; the sign lives in k. This
  is not a contradiction (the reading belongs to the law's own conic, not to the geometry it averages over), but it
  is worth a sentence wherever the curvature reading is defined.

## 4. The Legendre ladder, and Legendre's K

![Fig 4](figures/fig4-legendre-ladder.png)

**Fig 4.** (a) Each Legendre function P_ν(1 + 2t) multiplies on average and grows at rate ν(ν + 1) per unit of
typical length. The mean odds are ν = 1 (rate 2), the typical value is the slope at ν = 0 (rate 1), and the mean
transmission sits at ν = −½ (rate ¼, the bottom of the hyperbolic plane's spectrum). (b) The per-barrier factor of the
mean transmission, P₋½(1 + 2t), equals (2/π)·√T·K(|r|) (dots). (c) Mean transmission from the Legendre spectrum
(orange) against simulated chains (dots), with the typical and mean-odds curves for comparison.

- **One family (8.1, 8.2).** ν = 1 gives ∘₂; the slope in ν at ν = 0 is ln(1 + t), so ∘₁ is the same family's
  first derivative.
- **Mean transmission (8.3–8.6).** T = 2/(1 + z) has Mehler–Fock coefficient 2π/cosh(πλ), so after N identical
  barriers ⟨T⟩ = ∫ λ tanh(πλ)·2π/cosh(πλ)·P₋½₊ᵢλ(1 + 2t₁)ᴺ dλ; it matches the simulations at N = 2, 4, 8. At large
  N the decay per barrier is P₋½(1 + 2t₁) = (2/π)·√T₁·K(|r₁|), Legendre's complete elliptic integral of the first
  kind with modulus equal to the reflection amplitude, with the usual N^(−3/2) in front.
- **Rates 2, 1, ¼ (8.7).** For weak barriers the mean odds run twice as fast as the typical value (the ratio of their
  k's), and the mean transmission a quarter as fast. The ¼ is the bottom of the spectrum of the curvature −1 plane.
- **Known versus ours.** Known: the product formula, Mel'nikov's equation (1981), the continuum mean transmission
  with its e^(−s/4), and Ng and Englert's Legendre-function treatment of phase disorder (2013, arXiv:1212.1951; we
  could not confirm whether they write the per-barrier factor as an elliptic integral). P₋½(cosh η) as an elliptic
  integral is a standard identity. **Ours:** reading the three averages as three points of one ladder, and the
  link to the elliptic lead (`closed-leads-fresh-look/` F3): tunnelling reaches Legendre's K through the
  mean transmission, where F3 reached it through Jacobi's addition law.

## 5. Rewinds

![Fig 5](figures/fig5-rewinds.png)

**Fig 5.** (a) Adding one barrier multiplies each spectral component of the spread by P₋½₊ᵢλ(1 + 2t₁); at its zeros
(red) the component is erased. (b) The same function over complex λ: every zero is on the real line. (c) Rewinding the
mean divides by 1 + 2t₁ for barriers (never zero) and by 1 − 2s₁ for mixers (zero at s₁ = ½: the wall).

- **Single numbers rewind (9.1, 9.2).** Removing a barrier from the mean odds is the ∘₂ inverse (t − t₁)/(1 + 2t₁),
  always defined. On the sphere the same rewind has a wall: a 50/50 mixer erases the bias and cannot be undone.
- **The full spread does not (9.3, 9.4).** Each barrier's character P₋½₊ᵢλ(1 + 2t₁) has only real zeros in λ (box
  count equals real count), because Mehler's kernel is increasing and Pólya's 1918 theorem applies. So every barrier
  erases a comb of components exactly, and no rewind of the spread can bring them back.
- **Fit with the project.** This matches the Riemann work's finding that our exact rewinds are single-number flows:
  here the mean and the typical value rewind; the distribution does not. In the continuum the flow is backward heat
  flow on the hyperbolic plane (Section 6), the same obstruction in another geometry.
- **Known versus ours.** Pólya's theorem and the conical functions are known; applying them to "which parts of a
  barrier chain's statistics can be rewound" is ours, and modest.

## 6. The continuum flow picks k

- **The flow (7.1–7.4).** Many weak elements give the flow ∂ₜ t(1 + εt) ∂ₜ with ε = −K. At ε = 1 it is Mel'nikov's
  one-channel equation (DMPK for one channel); at ε = −1, heat flow on the sphere.
- **It keeps exactly two settings.** Among all powers (1 + kt)^α and logs ln(1 + kt), the only exact laws are the
  mean at k = 2ε and the log at k = ε (sympy, Lean `continuum_mean_unique`, `continuum_log_unique`). So the physics
  does not just allow k = 2 and k = 1 for barriers; it forces them, and forces −2 and −1 for mixers.

## 7. Where tunnelling puts the dial

![Fig 6](figures/fig6-dial-line.png)

**Fig 6.** The homothety line. Means sit at k = −2K, typical values at k = −K, flat geometry and the pair's mean
transmission at k = 0.

## Leads, judged on their own merits

1. **The dial boost covers barriers and mixers (open, strong).** It joins tunnelling to twisted boosts, relativity and
   the holonomy with no new parts, and it makes predictions anyone can check in a multilayer or an interferometer.
   Next: the three-barrier case, where order matters at fixed phase, read through the order dial's oriented area.
2. **Phase averaging as exact dequantization (open, worth following).** A random phase turns amplitudes into
   probabilities (+) and log-probabilities into max, with no limit taken. This could give the semiring audit a
   physical source, and the wall rule (stronger path flips at a tie) is the kind of thing that carries to other
   areas: interferometers, random networks, and possibly decoding, where the same max appears.
3. **Legendre's K in the mean transmission (open, small).** A clean bridge between tunnelling and the elliptic lead.
   Next: see whether the other conical points of the ladder have similar closed forms.
4. **Rewinds (closed for now).** Single numbers rewind, the spread does not; consistent with the Riemann pages and
   not a new road.
5. **The curvature reading (a note for the design rule).** The averaged laws carry their curvature's sign in k, not
   in −Δ/4. Worth one sentence on the dial pages.

## Inputs and readings

Inputs (base dials and settings): each element's odds t (or switch probability s), the loop phase φ, and the
curvature K of the element's geometry. Readings: the pair's odds, the lengths D, the twist ω, the growth rates, and
the curvature reading −Δ/4 of each law.
