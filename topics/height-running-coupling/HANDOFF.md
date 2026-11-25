# Handoff: the height flow read as a running coupling

*Nov 25, 2026. Stopped partway because Brennan merged the height-flow work into the entropy stream. This file is everything found so far. Script: `running_coupling.py` → `running_coupling_output.txt` (14 of 14 checks pass), numbers in `running_coupling.json`. No figures made yet; `figures/` is empty.*

## Verdict so far

The running-coupling reading is **exactly right but generic**. Any map with a double fixed point (like E(x) = e^x − 1 at 0) has a flow that looks like a marginal coupling at one and two loops; nothing here is special to E. The renormalon resemblance is **real in mechanism** (same object, same exponent law) but **differs in where the singularities sit** (imaginary axis for E, real axis for QFT renormalons). I found **no physics test**. As a route toward physics or the North Star it is weak; as a correct dictionary entry between iteration theory and RG it is worth keeping. Most of it is known on both sides; the side-by-side reading is ours.

## 1. The height flow as a beta function (checked)

- **β(x) = x²/2 − x³/12 + x⁴/48 − x⁵/180 + 11x⁶/8640 − …** This is Julia's iterative logarithm of E, from v(E(x)) = E′(x)·v(x). Exact coefficients to x²⁴ in the JSON, 150 terms at 260 digits. Known (Julia, Erdős–Jabotinsky 1960).
- **Scheme-free numbers.** In QFT, b0 and b1 are the same in every scheme and b2 onward are not. The same holds here, and for the same reason: under any change of coordinate x = y + a·y² + …, b0 stays put and so does **ρ = −b1/b0² = Res(dx/β) = 1/3**, while b2 changes. Checked with x = y + 0.7y² − 1.3y³: b2 went from 0.0208 to 0.858, while b0 and ρ stayed the same. Known on both sides: the formal normal form of a double-zero vector field is the 't Hooft scheme, β = b0x² + b1x³ exactly.
- **Two-loop running is visible.** With z = −2/x the Fatou coordinate (Abel function) is Φ = z − (1/3)·log z + …. The log term is the two-loop term. For x_n = E^n(−1), the one-loop formula −2/n is off by 1.75e-4 at n = 200; the two-loop formula is off by 1.3e-8.
- **Sign pattern.** E's b0 and b1 have opposite signs, like φ⁴ in 4d (β = 3λ²/16π² − (17/3)λ³/(16π²)²). QED and QCD have b0 and b1 with the same sign (table in the output). The size of ρ is **not** comparable, because the time unit is free (one height = ln μ or ln μ²?): rescaling time multiplies ρ by any positive number, so a numerical "match" with any theory can always be manufactured. Only the sign carries meaning.
- **Key dictionary line (known in principle, ours to state):** E's Fatou coordinate Φ plays exactly the role of ln Λ in QFT (the RG-invariant scale, ln(Λ/μ) = −∫dg/β).
- **Known toy-model literature:** Curtright & Zachos, "Renormalization Group Functional Equations", Phys. Rev. D 83, 065019 (2011), arXiv:1010.5174, treats RG flows as Schröder/Abel continuous iteration. This is the closest known prior work; read it before claiming anything new.

## 2. The renormalon test (checked numerically)

- **Divergent, with singularities spaced along a lattice.** Borel transforms in z of both the generator V(z) = 1 + (1/3)/z + … and the formal half-iterate E^(1/2) have their nearest singularities at **±2πi**. Padé poles land at ±6.273i; the expected 4πi, 6πi, … follow. Known: Écalle; Dudko & Sauzin, arXiv:1307.8093 (2πiℤ-resurgent Fatou coordinates).
- **The singularity exponent is set by ρ, as in QFT.** The large-order coefficients oscillate with a phase that drifts like j^(−iσ). Fitted σ = 2.0932 for the generator and 2.0940 for the half-iterate; the prediction is **2πρ = 2π/3 = 2.0944**. In QFT, the exponent of a renormalon singularity is also set by b1/b0² (Beneke, Phys. Rept. 317, 1999). Mechanism, the same in both: corrections are exponentials of the Fatou coordinate, e^(cΦ), and the log ρ·log z inside Φ turns into a power.
- **What differs:**
  - In QFT the factor c is real (power corrections Λ^p, from operator dimensions), so renormalons sit on the **real** Borel axis and the series is genuinely ambiguous there.
  - For E, c = 2πik is imaginary (from the step Φ ↦ Φ + 1), so the singularities sit on the **imaginary** axis. The series is Borel-summable along the real line, and for real x the regular iterate is unique, with no ambiguity.
  - E's "ambiguity" is the classical freedom Φ → Φ + p(Φ), with p 1-periodic, in choosing a fractional iterate, and the Écalle–Voronin invariants measure how the two sides fail to glue. In short: renormalon structure rotated by 90°.

## 3. Testable consequence

None found for physics. The only setting where this kind of mathematics predicts something physical is a **discrete** RG step: log-periodic critical amplitudes in hierarchical models (Derrida, Itzykson & Luck, Commun. Math. Phys. 94, 115 (1984)). That is not a continuous running coupling, and I didn't find a way to turn it into a test for us. I rate it low.

## Open questions for the entropy stream

1. Read Curtright–Zachos 2011 and say exactly what is new, if anything, in our reading.
2. The log-periodic door: would a "height = discrete RG step" model give log-periodic wiggles of measurable size? This needs the Écalle–Voronin invariants of e^x − 1 (Dudko–Sauzin 1307.8095 gives them as convergent series). Probably tiny.
3. How this sits with the earlier asymptotic-safety finding (G's one-loop flow matched the circle operation and broke past one loop). Agreement at one loop is generic for any double-zero flow, so it's expected that it breaks at two loops unless ρ matches. Worth one line in the merged write-up.
4. Figures still to draw: the β phase line with the flows; x_n against the one-loop and two-loop curves; coefficient growth against (k−1)!/(2π)^k; Borel-plane Padé poles for E next to a cartoon of QFT renormalons on the real axis; the b0/b1 sign table.

## Known vs ours

- **Known:** iterative logarithm and its divergence; regular iteration (Szekeres, Écalle); resurgence of Fatou coordinates at 2πiℤ (Écalle; Dudko–Sauzin); b0 and b1 scheme-independent, 't Hooft scheme; renormalon exponent from b1/b0² (Beneke); RG as functional conjugation (Curtright–Zachos); log-periodic amplitudes (Derrida–Itzykson–Luck).
- **Ours:** the side-by-side dictionary for E (ρ = 1/3 as its two-loop invariant; Φ ↔ ln Λ; the φ⁴-like sign pattern; the time-unit caveat; the numerical checks σ = 2πρ for this map; renormalons "rotated 90°").
