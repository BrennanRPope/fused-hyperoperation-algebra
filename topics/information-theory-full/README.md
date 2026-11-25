# Information theory through mult-add: the full pass

Page: index.html (self-contained). Checks: checks.py → checks_output.txt (28/28 pass). Lean: lean/InformationFull.lean (31 theorems, compiles clean), lean/README.md.

Builds on math-theory/information-theory/ (page SBsWmduwSNADpPVbeRDq2v, InfoFused.lean), entropy-noise (1To24w7SazEuaE5Pjsj4YM), quantum-noise (X6u4iNGch19w41ioqJ5VoW), thermodynamics (SKKLwKiEg4nyYGfrw32sp5), q (9AtexUZrzvsPfzSA7g91dL), and the measure-probability page of this pass (its E[X ✳ₖ Y] = EX ✳ₖ EY engine).

## Findings
- Grafting a D-ary decision tree onto every leaf of another composes question counts by the dial with k = D − 1; leaves are the lift; the depth-n complete tree asks [n]_D questions (q = D). Counts classical; reading ours. Binary case proved on real trees in Lean.
- Tsallis q-log is the dial's logarithm: ln_q(xy) = ln_q x ✳ₖ ln_q y, k = 1 − q; e_q is its inverse (Tsallis 1988/1994).
- Tsallis' pseudo-additivity = "expectation passes through the dial" applied to q-surprise; one Lean lemma covers Shannon (k = 0) and Tsallis. Correlation adds exactly k·Cov. Proof route ours.
- Conventions: project q = k + 1 = 2 − q_Tsallis, which is Tsallis' known additive duality q ↔ 2 − q, i.e. the dial's sign flip.
- Sharma–Mittal H_{α,β} (1975) = Rényi_α pushed through the dial's exponential, composes by k = 1 − β; Tsallis is β = α, Landsberg–Vedral β = 2 − α.
- Any binary channel: the crossover sum a + b chains by the union rule (det multiplies); BSC and Z-channel are special cases.
- Dropped: the q-exponential as a rung flow (half-steps match e^x − 1 through x², fail at x³: 1/8 vs 1/6).

## Verdict
No new lead on its own merits; the pass explains where each weight comes from and makes information theory a derived branch of the measure picture. Worth revisiting only if a bent landscape law describes a real coding cost.
