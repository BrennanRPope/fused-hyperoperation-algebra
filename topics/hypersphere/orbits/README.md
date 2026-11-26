# DRAFT · Newton's orbits and hydrogen's levels on one 3-sphere: does the swap rule's area law give Mercury and fine structure?

**DRAFT, Nov 26 2026.** Part of Affine Composition Theory (repo: fused-hyperoperation-algebra), a section of
`hypersphere/`. Brennan asked to follow the lead "about Newton's orbits and hydrogen's energy levels": do GR's
perihelion shift (Mercury's 43″ a century) and hydrogen's fine structure both come out as the swap rule's area law
(turn = d × area) on the 3-sphere they share?

Checks: `checks_orbits.py` → `checks_orbits_output.txt`, **23 of 23 pass**, and `six/checks_six.py` →
`six/checks_six_output.txt`, **11 of 11 pass**, and `twov/checks_twov.py` → `twov/checks_twov_output.txt`, **8 of 8 pass**, and
`dcopy/checks_dcopy.py` → `dcopy/checks_dcopy_output.txt`, **11 of 11 pass** (numpy, scipy; real orbit simulations, not
formulas fed back). Pictures: `figures_orbits.py` → `figures/`. Lean: `../lean/Hypersphere.lean` (§5–§9, 33 theorems), **not
compiled here**. Colour key as before: **orange = known, blue = ours**.

## The result in three sentences

1. **The 3-sphere they share is real and known** (Fock 1935): every Kepler orbit and every hydrogen level is a great
   circle on one 3-sphere, which is why the orbit closes and why level n holds n² states.
2. **For the electric force, what breaks that sphere is exactly one area, on a different sphere**: the relativistic
   orbit's turn per lap, and hydrogen's fine-structure splitting, are both the area of one disc in relativity's
   velocity space, and that is the same area that turns a gyroscope (Thomas precession).
3. **For gravity, Einstein's 1 − 6M/r splits exactly into two twists** (section 5): (1 − 2M/r) for curved space, times
   (1 − (2v)²), a velocity disc at **twice** the orbital speed. A gyroscope on the same orbit gets the same space
   factor times a disc at the orbital speed itself. The "twice the speed" is really 1 + 3 (section 6): v² from
   special relativity, as in the electric orbit, plus 3v² because gravity pulls on energy, and a falling body's
   energy grows. So the 6 is 2 (space) + 1 (special relativity) + 3 (energy). Mercury's 42.98″ follows. Through
   string theory's double copy (section 7), 1 + 3 is one square: gravity's coupling is the electric one squared,
   (E + p)², and its cross term 2E·p doubles Sommerfeld's electric term, so the 2 in 2v is the 2 of the square.

**Does it open a lead?** Yes. For the electric force it turns the Coulomb orbit and Dirac's fine structure into one
swap-rule area. For gravity the 6 is no longer put in by hand: it is an exact product of the same two twists the
Gravity Probe B page found for the gyroscope, with the velocity disc's radius doubled. That joins orbits and
gyroscopes into one swap-rule picture, and section 6 says where each piece of the 6 comes from.

## 1. Two spheres, not one

![Two spheres](figures/figO1-two-spheres.png)

*Figure O1. A: a relativistic orbit; the ellipse turns a little each lap. B: on Fock's momentum sphere each lap is a
great circle, and the circle tilts lap by lap. C: the amount it tilts is the area of a disc in relativity's velocity
space (a sphere of imaginary radius c, the curvature dial at d = −1/c²).* (The curvature dial was formerly called the plane dial.)

**Known.** Fock (1935) mapped momentum p to the unit 3-sphere by stereographic projection with radius p0 = √(−2E)
(Lean `fock_on_sphere`). Each Kepler orbit's momenta then lie on one great circle (check K3), tilted by the
eccentricity: (h·p0)² + e² = 1 (Lean `kepler_tilt`). Rotating the 3-sphere is the hidden SO(4) symmetry carried by
the Runge–Lenz vector (Pauli 1926; Bargmann 1936; check K1); regularizing collisions on it is Moser (1970). Level n
holds 1 + 3 + … + (2n − 1) = n² orbital states (K2; Lean `deg_eq`).

**What we tested.** If the precession were an area on Fock's sphere, different orbits would enclose different
areas. They don't: on the unit Fock sphere every orbit is a great circle and bounds half the sphere, 2π, for all
orbits (P7). So **the precession's area is not on Fock's sphere**. Fock's sphere only shows the effect: the circle
turns round the energy axis by exactly the advance each lap (P8: 0.31000 vs 0.31001).

## 2. The electric orbit: the turn per lap is one disc's area

![Advance is the disc area](figures/figO2-advance-is-the-disc-area.png)

*Figure O2. A: eight simulated orbits with different eccentricities (dots) sit on one curve, the area of the velocity
disc of speed-radius k/L. B: that curve is identical to the Thomas turn of a gyroscope on the circular orbit with the
same angular momentum.*

**Known.** Sommerfeld (1916) found the relativistic Coulomb orbit's advance per lap,
2π(1/√(1 − (k/cL)²) − 1). Thomas (1926) found a gyroscope on a circular orbit turns 2π(γ − 1) per lap, and that this
is the curvature of velocity space times the area its velocity encloses (the swap rule, already on our
`precession/` page).

**Ours (as stated; check the literature).** Two exact facts tie them together:
- The advance is **exactly** the area of a disc of speed-radius k/L in velocity space (curvature −1/c²), at any
  eccentricity (P1, P2: simulated 0.31000577, area 0.31000577 at c = 3).
- The circular orbit with that L moves at speed exactly k/L (Lean `circular_L_v`), so **an orbit precesses exactly
  as much as a gyroscope on the circular orbit with the same angular momentum** (R1, exact at β = 0.1, 0.5, 0.9;
  Lean `sommerfeld_eq_thomas`).

**Literature check (Nov 26, at Brennan's request).** We found this equality **implied but not stated**:
- **Biedenharn (1983)**, "The Sommerfeld puzzle revisited and resolved", shows that Sommerfeld's rotating frame for
  the precessing orbit is exactly paralleled by the Dirac transformation that carries the electron's Thomas spin
  precession, and turns into it in the classical limit. That is the same correspondence in qualitative form. We read
  it through Vickers (2018), who quotes it; we could not open Biedenharn's paper here, and no source we read gives the
  equality as numbers.
- **Belloni and Reina (1986)** is about Sommerfeld's own geometric derivation of Thomas precession. Its abstract says
  nothing about the orbit, and its full text wouldn't open here.
- **Vankov (arXiv:1008.1811)** links Thomas precession to perihelion advance only to first order, in a non-mainstream
  argument that rejects general relativity. It does not state this equality.
- **Gibbons (2015)**, on Hamilton's hodograph in relativity, treats the same orbits without mentioning Thomas.
- Other searches (velocity-space area, same angular momentum, exact equality) found nothing.

So the exact statement (advance per lap = velocity-disc area = the Thomas turn of the circular orbit with the same
L) looks unstated in what we could read. It is a short calculation, though, so it should be called a small new
observation in a long-studied area, not a discovery. Reading Biedenharn's paper in full is the one check left.

## 3. Hydrogen: the same area gives Dirac's fine structure, exactly

![Fine structure from the area](figures/figO4-fine-structure-from-the-area.png)

*Figure O4. A: the 3-sphere makes each Bohr level n²-fold (grey); one area per κ splits it into Dirac's levels
(blue). B: the area for each κ = j + ½, the disc of speed-radius αc/κ.*

**Ours, as an exact rewrite of Dirac.** Take x = (the disc area for κ)/2π and the quantum defect δ = κx/(1 + x). Then
(κ − δ)² = κ² − α², which is exactly Dirac's defect (Lean `dirac_defect`), so E = mc²/√(1 + α²/(n − δ)²) is Dirac's
formula. Rebuilt this way, every level up to n = 7 matches Dirac to the last digit (F1), and 2P3/2 − 2P1/2 comes out
10.950 GHz, Dirac's value (F2).

**Limits, said plainly.** The measured splitting is 10.969 GHz; the 0.2% gap is QED (the electron's anomalous
moment), which this law doesn't contain. 2S1/2 and 2P1/2 share κ = 1, so the area law leaves them equal; the Lamb
shift (1.058 GHz) is outside it (F5). And this is a rewrite of Sommerfeld–Dirac, not a derivation from the sphere.

## 4. Gravity: one area, but a coefficient of 6, and only at first order

![One area, three coefficients](figures/figO3-one-area-three-coefficients.png)

*Figure O3. A: advance ÷ (hodograph area/c²) from simulated orbits: scalar gravity −1, the electric force +1, general
relativity +6. B: GR's coefficient drifts away from 6 as the field gets stronger.*

![Mercury's hodograph](figures/figO5-mercury-hodograph.png)

*Figure O5. Mercury's velocities trace a circle of radius 48.9 km/s (Hamilton's hodograph). 6 × its area/c² per orbit
gives 42.98″ a century.*

**What holds.** Every theory's first-order perihelion shift is a whole number times the same area, the hodograph
disc over c² (P4): Nordström's scalar gravity −1, the electric force +1, Schwarzschild +6, at any eccentricity (P5).
For Mercury, 6 × area/c² = **42.98″ per century** (P6), the general-relativity value.

**What doesn't.** That number is GR's own first-order formula, rewritten; the 6 is put in, not derived. And GR is not
one exact area law: at stronger fields the coefficient drifts (8.89 at c = 5, P5b). So the honest answer to "does
Mercury come out as the swap rule's area law on the 3-sphere?" is **no, not exactly, and not on Fock's sphere**: it
comes out as 6 times an area on the velocity sphere, at first order.

**The sharpened question.** Compare the orbit with a gyroscope on it, both in units of that area (known PPN
formulas, R3; Lean §7):

| | orbit turns | gyroscope turns | ratio |
|---|---|---|---|
| electric force (exact) | 1 | 1 | 1 |
| Nordström scalar gravity (first order) | −1 | −1 | 1 |
| general relativity (first order) | 6 | 3 | 2 |

Brennan said to look for an exact split of 1 − 6M/r. Section 5 has it.

## 5. Gravity's 6, split exactly

![The six splits](figures/figO6-the-six-splits.png)

*Figure O6. A: the orbit's factor 1 − 6M/r (blue dashed) is exactly the space factor 1 − 2M/r (grey) times
1 − (2v)² (blue), v being the local orbital speed. The gyroscope's 1 − 3M/r (orange dashed) is the same space factor
times 1 − v² (orange). Dots are simulated orbits. B: in velocity space the gyroscope turns by the disc of radius v,
and the orbit by the disc of radius 2v. In O6 and O7 the colours mark the gyroscope and the orbit, not known vs ours.*

**The result.** For a near-circular orbit, write T = (radial frequency ÷ angular frequency)²; the advance per lap is
2π(1/√T − 1). Measure v as a local observer at rest would. Then, exactly:

- **orbit:** T = 1 − 6M/r = (1 − 2M/r) · (1 − (2v)²) (S1; Lean `orbit_split`)
- **gyroscope:** 1 − 3M/r = (1 − 2M/r) · (1 − v²) (S2; Lean `gyro_split`, as on the Gravity Probe B page)

Both are products of factors, so in the theory's language the fractions fuse by the circle operation (called mult-add before 1.0): −6M/r = (−2M/r) ✳ (−4v²) (Lean
`orbit_fractions_fuse`). In area terms, the 6 = **2** (curved space, the same for orbit and gyroscope) + **4** (a
velocity disc of radius 2v has four times Thomas's area) (S8).

**Two landmarks fall out.** The orbit's disc reaches light speed (2v = c) exactly at r = 6M, the innermost stable
orbit; the gyroscope's disc reaches light speed (v = c) exactly at r = 3M, the photon sphere (S3, S4; Lean
`isco_two_v`, `photon_sphere_v`). This also explains Figure O3B: the coefficient drifts from 6 in strong fields
because the true law is a product, not a sum.

![One rule, three forces](figures/figO7-one-rule-three-forces.png)

*Figure O7. A: one rule for circular orbits of three kinds. N is Newton's apsidal number, 3 + r·(force)′/(force),
which is 1 for the inverse square, so Newton's orbits close. B: simulated orbits for five electric-type potentials
and three metrics land on the line.*

**The general rule.** The split is not special to Schwarzschild. For **any** static round spacetime
(−A dt² + B dr² + r² dφ²), T = (N − 4v²)/B, where N is Newton's closure number built from the clock rate A (S5,
symbolic). For a relativistic particle in **any** potential, T = N − v² (S6). Simulated geodesics in Schwarzschild,
in a clock-warp-only spacetime (B = 1) and in a charged black hole (Reissner–Nordström) all match to a few parts in a
million (S7). The clock-warp-only case shows where each piece lives: 1/B (space) supplies the 2, and the clock rate
supplies N − 4v².

**Limits, said plainly.**
- This is the near-circular law, small wobbles about a circle. The electric result in section 2 holds at every
  eccentricity; for gravity, eccentric orbits have no closed form, and the 6 is only first order there.
- Mercury (S9) tests only the first order, which was already known to give 42.98″.
- The frequencies themselves are textbook; what we add is the grouping into twist factors and the reading.
  We have not searched the literature for this grouping yet.
- **Why 2v?** Section 6.

## 6. Where the 2v comes from, and what string theory says

![Where 2v comes from](figures/figO8-where-2v-comes-from.png)

*Figure O8. A: the first-order budget. The electric orbit has only special relativity (1). The gyroscope adds curved
space (2 + 1 = 3). The orbit adds a third piece, gravity pulling on energy (2 + 1 + 3 = 6). B: the potential the
orbit actually feels is steeper than Newton's 1/r. C: the exact factors at every field strength.*

**The answer (ours, exact).** It is not a doubling. Inside the clock warp, a geodesic is exactly a
special-relativity orbit in the potential E(1 − 1/√A) (D1, and simulated orbits agree to 3 parts in 100 million, D7),
and that picture's speed is the local speed (D2). Two things then spoil Newton's closed orbit:

- **Special relativity, −v²:** the same velocity disc as the electric orbit and the Thomas turn (section 2).
- **Gravity pulls on energy, −3v²:** a falling body's energy grows as 1/√(1 − 2M/r), and gravity pulls on that energy.
  So the pull is steeper than 1/r: 1 − 1/√(1 − 2M/r) = −M/r − (3/2)M²/r² − … (D6). That steepness lowers Newton's
  closure number from 1 to exactly 1 − 3v², for any clock rate (D3, D4; Lean `energy_pull_closure`).

v² + 3v² = 4v² = (2v)², which is why it looked like a disc of radius 2v. Exactly:
**1 − 6M/r = (1 − 2M/r) · [(1 − 3v²) − v²]** (D5; Lean `six_is_two_one_three`). The 6 is 2 (space) + 1 (special
relativity) + 3 (energy) (D8). The gyroscope gets only the first two, 2 + 1 = 3, because a gyroscope's turn doesn't
care how steep the pull is.

**String theory, as Brennan suggested.** We looked in three places:
- **T-duality and open/closed strings:** no factor of 2 on a speed there that bears on orbits. Nothing found.
- **The double copy** (gravity = gauge theory squared, from string theory's KLT relations, 1986): this is the right
  place. It builds Schwarzschild from the Coulomb field (Monteiro, O'Connell and White 2014), maps a charge's orbits
  onto geodesics (Gonzo and Shi 2021), and its central rule is that gravity's charge is energy and momentum rather
  than a fixed charge. That is the same physics as our 3v². No paper we found derives the 6 from the double copy;
  section 7 now does.
- **String theory's own gravity can remove the 6 entirely.** In N = 8 supergravity, the low-energy limit of string
  theory, orbits of extremal black holes don't precess at all: the extra scalar and vector forces cancel it
  (Caron-Huot and Zahraee 2018–19). That fits our table: scalar forces push the coefficient negative (Nordström −1),
  and the budget can sum to zero. So the 6 is specific to pure general relativity, not to string theory's gravity.

**Is the 1 − 6M/r split new? Not really.** Both ingredients are textbook: Ω_r² = Ω_φ²(1 − 6M/r) and a local orbital
speed of c/2 at the innermost stable orbit. The split is a one-line rewrite of them. We did not find it written this
way, nor the general rule T = (N − 4v²)/B, but published formulas for epicyclic frequencies in general static
spacetimes contain it. Call it a new way of reading known results. The 2 + 1 + 3 budget is the useful part.

## 7. The 6 through the double copy: the 2 in 2v is the 2 of a square

![Gravity's 6 through the double copy](figures/figO9-double-copy-six.png)

*Figure O9. A: gravity's coupling is the electric coupling squared, and the square (E + p)² has three pieces. B: the
same 6 read two ways. Both give curved space 2; the double copy turns special relativity's 1 and the energy pull's 3
into one square, 2² = 4. C: switching the pieces on one at a time, exactly.*

Brennan asked to try deriving the 2 + 1 + 3 through string theory's double copy. **It works, and it explains why
1 + 3 is a perfect square.**

**The setup (known).** Written in Kerr–Schild form, the double copy's own coordinates, Schwarzschild is flat space plus
(2M/r) k k, with k a light-like direction, and the inverse is exact (K1). So a falling body obeys
η(p, p) − (2M/r)(k·p)² = −1, with k·p = E + p (K2). The single copy, Coulomb, has the same k but one power:
η(p, p) − (2Q/r)(k·p) = −1. Gonzo and Shi (2021) use exactly this map, which "replaces the color charge of the test
particle by its momentum", to rebuild Schwarzschild geodesics from charged orbits.

**What we added (ours).** Expand the square: (k·p)² = E² + 2E·p + p² (K3). Each piece does one job, and we switched
them on one at a time, exactly and in simulated orbits (K4–K10):

- **E², Newton's pull:** on its own, the orbit closes exactly, at every field strength: no turn at all (K4, K8).
- **p², space:** with Newton's pull it gives exactly 1 − 2M/r, the same space factor as the gyroscope (K5). It adds 2.
- **2E·p, the cross term:** in the single copy this term is what makes Sommerfeld's electric turn, the 1 (K9). In
  gravity the square doubles its strength, from M to 2M, so it gives 2² × 1 = 4 (K10).

So **6 = 2 (space) + 2² (the cross term)**, and all three together give exactly 1 − 6M/r = (1 − 2M/r)(1 − (2v)²)
(K6; Lean `ks_all`, `six_is_space_plus_square`). The "2v" is the 2 of the binomial square. The clock-warp reading of
section 6 found the same 6 as 2 + 1 + 3; the two readings agree on space's 2, and the double copy shows that special
relativity's 1 and the energy pull's 3 are one square, 1 + 3 = 2² (Lean `one_plus_three_is_square`).

**One caution.** How the 6 splits into pieces depends on the coordinates; only the total is the same in every frame.
The Kerr–Schild frame is the natural one here because it is where the double copy is exact and linear.

**Is it published?** The Kerr–Schild form, its exact inverse and the charge-to-momentum map are known (Monteiro,
O'Connell and White 2014; Gonzo and Shi 2021). Gonzo and Shi do not compute a precession, and we found no paper that
splits the 6 into the three pieces of (k·p)², or ties gravity's 2v to the binomial 2. So the derivation is ours, built
on known parts.

## Known vs ours

- **Known:** Fock's 3-sphere and SO(4) (Fock 1935, Pauli 1926, Bargmann 1936), Moser regularization (1970),
  Runge–Lenz, Hamilton's hodograph (1847), Sommerfeld's advance (1916), Thomas precession (1926), Dirac's levels
  (1928), the PPN perihelion and geodetic formulas, Mercury's 43″.
- **Ours, exact, implied by Biedenharn (1983) but not stated as an equality in anything we read:** Sommerfeld's
  advance per lap = the velocity-disc area = the Thomas turn of the same-L circular orbit (R1, P2).
- **Ours, exact rewrite:** Dirac's fine structure from one disc area per κ (F1, F4).
- **Ours, first order, a reading:** every theory's perihelion shift is an integer times the hodograph area/c²
  (−1, +1, +6), and GR's orbit turns twice its gyroscope (R2, R3).
- **Known ingredients, our grouping (not found written this way):** 1 − 6M/r = (1 − 2M/r)(1 − (2v)²), beside the
  gyroscope's (1 − 2M/r)(1 − v²); the general rules T = (N − 4v²)/B and T = N − v² (S1–S7).
- **Ours, exact:** the 6 = 2 (space) + 1 (special relativity) + 3 (gravity pulling on energy), via the
  special-relativity picture of a clock-warp geodesic (D1–D8).
- **Known parts, our derivation:** the 6 = 2 (space) + 2² (the cross term 2E·p of the double copy's (k·p)²), with
  Newton's pull alone giving no turn at all; 1 + 3 = 2² (K1–K10). The Kerr–Schild form and the charge-to-momentum
  map are known; the split is not, in anything we found.
- **Known, string theory:** the double copy (KLT 1986; Monteiro–O'Connell–White 2014; Gonzo–Shi 2021); no precession
  in N = 8 supergravity (Caron-Huot–Zahraee 2019).
- **Known, now read as discs:** v = c/2 at the innermost stable orbit, v = c at the photon sphere.
- **Does not work:** the precession as an area on Fock's sphere (P7); GR's perihelion as one exact sum of areas
  (P5b); the split for eccentric orbits in closed form.

## Next (offered, not started)

1. Try the same three-piece split for a spinning black hole (Kerr), whose Kerr–Schild form is also exact.
2. Read Biedenharn (1983) in full for the R1 equality, and search for the T = (N − 4v²)/B grouping.
3. Compile `../lean/Hypersphere.lean` in the shared library (33 theorems).

## Files

`checks_orbits.py`, `checks_orbits_output.txt`, `results.json`, `figures_orbits.py`, `figures/figO1…figO9`,
`six/` (`checks_six.py` 11/11, `checks_six_output.txt`, `results_six.json`, `figures_six.py`, `explore.py`,
`general.py`, `srnum.py`, `nordnum.py`), `twov/` (`checks_twov.py` 8/8, `checks_twov_output.txt`, `derive.py`,
`figure_twov.py`), `dcopy/` (`checks_dcopy.py` 11/11, `checks_dcopy_output.txt`, `results_dcopy.json`, `derive.py`,
`numcheck.py`, `figure_dcopy.py`), and `../lean/Hypersphere.lean`.

## References

V. Fock, Z. Phys. 98, 145 (1935). W. Pauli, Z. Phys. 36, 336 (1926). V. Bargmann, Z. Phys. 99, 576 (1936).
J. Moser, Comm. Pure Appl. Math. 23, 609 (1970). W. R. Hamilton, Proc. R. Irish Acad. 3, 344 (1847).
A. Sommerfeld, Ann. Phys. 51, 1 (1916). L. H. Thomas, Nature 117, 514 (1926). P. A. M. Dirac, Proc. R. Soc. A 117,
610 (1928). C. M. Will, Theory and Experiment in Gravitational Physics (PPN perihelion and geodetic formulas).
J. Bertrand, C. R. Acad. Sci. 77, 849 (1873) (closed orbits; the apsidal number). L. Belloni and C. Reina, Eur. J. Phys. 7, 55 (1986). A. A. Vankov, arXiv:1008.1811 (2010). L. C. Biedenharn, Found.
Phys. 13, 13 (1983). P. Vickers, Brit. J. Phil. Sci. (2018), arXiv:1912.04698. G. W. Gibbons, arXiv:1509.01286 (2015).
H. Kawai, D. C. Lewellen and S.-H. H. Tye, Nucl. Phys. B 269, 1 (1986) (KLT). R. Monteiro, D. O'Connell and C. D. White,
JHEP 12, 056 (2014). R. Gonzo and C. Shi, Phys. Rev. D 104, 105012 (2021), arXiv:2109.01072. S. Caron-Huot and
Z. Zahraee, JHEP 07, 179 (2019), arXiv:1810.04694.
