# The repetition dial (FHA 1.0 checklist, item 1) — Nov 25, 2026

Question (Brennan): the framework isn't complete "until we figure out how repetition can fall out of a greater framework". A separate dial is fine "if that's what makes sense". Repetition here means how many times one element may be used, the leftover found on the size-dial page (circle-operation counts each element at most once). The rung kind of repetition (adding repeated is multiplying) is a different thing and already lives in rung flow.

## Result in plain words
- **Repetition inside one element falls out of k. No new dial is needed for it.** Give an element the weight Z_k(x) = (1 + kx)^(1/k). Then 1/k is the number of cells the element has. Positive k: cells that hold one copy each, so k = 1/m allows at most m copies (k = 1 is the circle operation's "at most once"). Negative k: bins that hold any number, so repetition is unlimited (k = −1 is multisets). k = 0 is infinitely many tiny cells, the classical count e^x.
- **k sweeps it smoothly.** Every coefficient of Z_k is a polynomial in k, ∏_{j<n}(1 − jk)/n!, so nothing jumps, including at k = 0. Every value of k composes by the circle operation a + b + kab: read each element through the k-log (Z^k − 1)/k and its value is its own weight x.
- **One honest catch: on the "at most" side, only k = 1, 1/2, 1/3, … give real counts.** In between (say k = 2/3) some counts come out negative (the coefficient of x³ is −1/54). On the "any number" side every k gives positive counts. So the dial is smooth, but its counting meaning is quantised on one side (Figure 4).
- **k can be read off one measurement.** Average copies n = x/(1 + kx) = 1/(e^cost + k) (Fermi-Dirac at k = 1, Boltzmann at 0, Bose-Einstein at −1), so 1/n − 1/x = k, and the spread is n(1 − kn): k is how much copies avoid (k > 0) or bunch with (k < 0) each other.
- **Two other kinds of "repetition" are not a single k** (Figure 3: a true k is a flat line at every crowding):
  - A hard cap of m copies with equal weight (1 + x + … + x^m, Gentile 1940) slides from k = −1 to 1/m as it fills. But it is built from pieces already in the framework: an "any number" element at size 1 divided by one at size m + 1.
  - Haldane's exclusion statistics (g between 0 and 1: anyons, fractional quantum Hall, Calogero-Sutherland), where copies of one element use up room for others, slides from k = 2g − 1 when sparse to k = g when crowded. That is crowding between elements. The circle operation ∘₁ (called mult-add before 1.0) composes independent elements by design, so this is genuinely outside it: it needs its own dial g, or a k that changes with crowding.

## Verdict for 1.0
- Settled: repetition within an element is k (1/k cells), proved and pictured.
- Open, outside the composition law: crowding between elements (Haldane's g). Adding it later as a new dial would not change any existing definition, notation or theorem name, so under the versioning rule it can come in a 1.x minor release. Whether it belongs before 1.0 is Brennan's call.

## Known vs new (honest)
- All the pieces are known: the interpolating formula 1/(e^x + a) is textbook; (1 + x/m)^m as m exclusive states is the binomial theorem; (Z^k − 1)/k is the Box-Cox transform / Tsallis q-logarithm (Box & Cox 1964; Tsallis 1988, with Tsallis q_T = 1 − k); the fluctuation n(1 ∓ n) is textbook; Haldane 1991 (PRL 67, 937) and Wu 1994 (PRL 73, 922) give the g statistics, and the sparse and crowded limits used here follow from Wu's equation.
- Already in this project, found during the names review (item 4): topics/dial-survey item 4 and lean-library DialSurvey.lean had the per-element copies count first, with the opposite sign (factor (1 − kz)^(−1/k), so k = −1 is subsets and k = −1/m is at most m copies), and proved those cases in Lean (kRise_fermi, kRise_bose, kRise_parafermi, fall_vanish). The pages agree once the sign is read as weight versus probability (names-review/lean/Conventions.lean, weight_prob_flip). This page duplicated that result; it should cite dial-survey.
- What this page adds on top: atoms compose by the k-rule (atoms_compose), the occupation and fluctuation in the k form, the 1/n − 1/x = k gap, the negative-count check for k = 2/3, and the one-dial test (Haldane's g as a k that slides with crowding). Packaging only, not found stated, not a discovery.
- Naming hazard for checklist item 4: Tsallis's index is q_T = 1 − k for the same family.

## Figures
- fig1-repetition-dial.png: the dial as cells, and the chance of 0 to 5 copies at the same average for k = 1, 1/2, 0, −1/2, −1.
- fig2-average-copies.png: n = 1/(e^cost + k) for five k; Fermi, classical, Bose on one fan.
- fig3-one-dial-test.png: measured k = 1/n − 1/x against crowding; the k family is flat, Haldane's g and the hard cap drift.
- fig4-honest-counts.png: sign of each count across k; all positive for k ≤ 0, honest on the exclusive side only at k = 1/m.

## Files (for the repo push)
- README.md, checks.py, checks_output.txt (29/29; needs sympy, numpy, mpmath), figures.py, fig1–fig4 PNGs.
- lean/RepetitionDial.lean: 22 theorems, core Lean 4.37.0-rc1, no sorry, standard axioms only. Imports Base only. To add to the library: put it at FusedAlgebra/Topics/RepetitionDial.lean, add `FusedAlgebra.Topics.RepetitionDial` to build.sh after DialPotential and an import line in FusedAlgebra.lean. Key theorems: atoms_compose, gbinom_cells, cap_at_m, cells_exclusive, unlimited_when_nonpos, bose_pos, two_thirds_negative, fluctuation, inverse_gap, semion_no_constant_k.
