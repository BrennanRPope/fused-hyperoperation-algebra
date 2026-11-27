# DRAFT: The spinning river, and every idea of a frame (Nov 26, 2026)

Asked by Brennan (voice call). He said "Yeah, that is a, a good idea" to giving each river cell a rotation as well as a
velocity. He added: "And I'd also like to figure out if there's any other concept of frame. I feel like procession is
where we started talking about frame dragging at first." He then said "Yes" to widening the search to frame-dragging
look-alikes outside gravity: Larmor's theorem, spinning fluids and Foucault's pendulum.

Follows `../jet-river/` (the river as a jet field).

## The result in plain words

**One move gives the whole spinning black hole.** Take the source of a non-spinning black hole and turn its height on
the spin axis by i, the dial's rotation number, from 0 to i·a. That one move sets everything the river needs: the flat
lattice's squashed spheres, each cell's speed, and each cell's turn. Every field equation of empty space then holds,
to 50 digits, and the curvature is Kerr's own. What is still put in is the mass, the river's rules, and the reading
that spin is this shift. Frame dragging far away sets its size: a = J/M, spin divided by mass. The shift is known as
the Newman–Janis shift, and nobody has settled why it works. In this theory i is a setting of the dial, not a trick,
which makes it a lead.

How it got here:

1. **Tested (Part 1):** a velocity plus a turn for each cell gives Kerr, and only the right turn works.
2. **Turn derived (Part 3):** from shells that turn rigidly and a clock gap that each falling ring keeps, plus
   Sagnac's area law and frame dragging.
3. **Speed and lattice derived (Part 4):** from the complex distance to the shifted source.
4. **Turn from the same shift (Part 5):** the shifted source's light rays lean sideways by exactly the turn.
5. **Charge (Part 6):** the same move with a charge gives the charged spinning black hole, and a magnet with the
   electron's ratio g = 2.
6. **Two bodies (Part 7):** one shear per cell can't hold two bodies, and two shears don't combine. What adds for
   bodies at rest is the log of the clock rate, by velocity addition, with a strut whose push is Newton's force.
7. **Two spinning bodies in orbit (Part 8):** three rules (fields sent at light speed, clock rates multiply, masses
   seen through each other's clocks) give gravity's slow-orbit field exactly. The orbits turn and decay as binary
   pulsars show, the spins don't line up but sort into two families, and the wave's run-up and ring-down come out.
   The merger is still a gap.
8. **The next level of detail (Part 9):** the three rules fail at 2PN. A clock riding on one body runs at the wrong
   rate at that order, by a formula in the mass ratio. The bodies only feel each other's own fields, and the clock
   rule turns out to shape only the field between them. What's missing is likely the field's own energy between the
   bodies, or the field of a body that is speeding up.
9. **A fourth rule (Part 10):** adding gravity's own field energy as a source closes that gap exactly, at every mass
   ratio. The last piece is a speeding-up body's own field multiplied by itself. With it, both body 1's clock rate and
   the orbit's size match general relativity's 2PN answers.

The turn shows up in a measurable place. The cell's clock leg doesn't close around the spin axis: after one lap at
the equator it is off by 2πa. That is exactly the shift in Kerr's orbital periods: prograde orbits take 2πa longer,
retrograde ones 2πa shorter, at every radius. This is the known "gravitomagnetic clock effect". Reading it as the
cell's clock gap is ours.

**Frames.** Our pages already hold 17 ideas of a frame. Foucault's pendulum, Larmor precession and the Berry phase
were on the precession page all along. They line up along the dial with one law, turn per loop = d × area. Of the
three outside-gravity look-alikes:

- **Larmor's theorem (new here):** a magnetic field acts on charges exactly like a turning frame. Its ½ is the same
  ½ as in the gyroscope rule of frame dragging.
- **Spinning fluids:** sound in a flowing fluid sees exactly our river. A draining, spinning bathtub drags waves the
  way a spinning star drags light. It does so with negative energy in our Q layer, though, and a fluid can't make
  the turn, so no fluid can make full Kerr.
- **Foucault:** this one was already ours. Earth's own frame dragging adds only 4.6 × 10⁻¹⁰ to the pendulum's turn,
  at every latitude.

![cell legs](figures/fig1_cell_legs.png)

*Fig 1. Left: a non-spinning hole. Each cell's radial leg is pushed by v along its time leg, which is one shear.
Right: a spinning hole. The radial leg is also pushed by −v·α along the azimuth leg: the turn, a second shear. The
two shears commute. Both sit at the dial's flat middle, where boosts and turns both become shears.*

All numbers are from `checks.py` (**73/73 pass**, about 17 minutes). Curvature is computed with exact symbolic
derivatives and 50-digit arithmetic at sample points (`numgr.py`).

---

## Part 1. The spinning river

### 1.1 The lattice: flat space in the shape of a spinning body (fundamental)

A spinning body is squashed at the poles, so the natural lattice is made of squashed spheres (oblate spheroids).
Their semi-axes are √(r² + a²) across and r along the axis, where a = J/M is the spin per unit mass. At r = 0 they
shrink to a ring of radius a. This is still ordinary flat space, only labeled differently:

> x = √(r²+a²) sin θ cos φ, y = √(r²+a²) sin θ sin φ, z = r cos θ  gives dx² + dy² + dz² exactly (K1b)

Every curvature component of this lattice is zero (K1a).

![spheroid lattice](figures/fig2_spheroid_lattice.png)

*Fig 2. The lattice in a slice through the spin axis. Black ellipses are the squashed spheres; dotted curves run
outward. The red curve is the horizon, the orange dots are the ring, and the shading shows how much each cell turns.*

### 1.2 Each cell: a velocity and a turn

The cell's legs are the lattice's legs, changed in one place. The radial leg picks up

> e_r = (lattice radial leg) + v·(time leg) − v·α·(azimuth leg),  v = √(2Mr)/ρ,  α = a sin θ / √(r² + a²)

with ρ² = r² + a² cos²θ. The first push is the river's boost, unchanged from the non-spinning case. The second is the
turn: the radial leg leans toward the azimuth by α. The map has the form I + (shear) with (shear)² = 0, and the two
shears commute (K5a; Lean `shear_sq`, `cellMap_mul`). That is the dial's flat middle, q = 0, the dual numbers, where
a boost and a turn are both shears and simply add.

### 1.3 Kerr comes out, exactly

| check | result |
|---|---|
| every field equation of empty space (Ricci = 0), a = 0.7 M | max error 4.5 × 10⁻⁵⁰ (K2a) |
| the same near the fastest spin, a = 0.99 M | 3.3 × 10⁻⁵⁰ (K2b) |
| curvature invariant R_abcd R^abcd vs Kerr's formula | agree to every digit shown (K2c) |
| horizon: g^rr ρ² = r² − 2Mr + a² | Kerr's horizon, r₊ = M + √(M² − a²) (K6b) |
| slow spin: the turn's first-order part | the swirl 2J/r³ of `jet-river` (K6a) |

This metric is known: it is Doran's form of Kerr (2000). Hamilton and Lisle (2008) read it as a river whose cells
carry a Lorentz transformation, a boost plus a rotation. What is ours is the reading at the dial's flat middle, two
commuting shears on a flat spheroidal lattice, and the forcing tests below.

### 1.4 The turn is forced

![turn is forced](figures/fig3_turn_is_forced.png)

*Fig 3. Largest field-equation error for each variant (log scale). Only the right turn with the right speed law is
zero (green). Velocity alone fails, any other size or shape of turn fails, and any other speed law fails (red).*

| variant | max error |
|---|---|
| velocity only, no turn | 0.53 (K3a) |
| half the turn / 1.5× the turn | 0.40 / 0.66 (K3b) |
| turn ∝ sin θ instead of sin²θ | 0.26 (K3b) |
| speed √(2Mr) or √(2Mr)/ρ² instead of √(2Mr)/ρ | 22.8 / 2.4 (K3c) |

So the turn isn't decoration. On this lattice, a spinning hole needs it, at exactly this size.

### 1.5 The velocity rule is unchanged

The cell's 4-velocity is u = −dt, and g^tt = −1 (K4a). So u has unit length and is a gradient, which makes it a free
fall. Its energy is E = 1, starting from rest far away, and it carries no angular momentum (K4b). This is the same
rule as in `velocity-lines` and `jet-river`: every cell is a zero-energy free faller. The speed law v² = 2Mr/ρ² is
that rule with the squashed-sphere "radius" ρ²/r in place of r.

### 1.6 The turn is a clock gap: Kerr's clock effect

The cell's time leg is dt − a sin²θ dφ. It does not close: going once around the axis, the clocks disagree by
2πa sin²θ.

Circular orbits in the equatorial plane have periods (K7a; Lean `kerr_period`)

> T± = 2π (r^(3/2)/√M ± a)  (prograde +, retrograde −)

So each orbit is shifted by exactly 2πa from Kepler's period, at every radius. That is the clock gap of one lap at
the equator (K7b). Their difference T₊ − T₋ = 4πa (4πJ/Mc² in ordinary units) is the known gravitomagnetic clock
effect (Cohen & Mashhoon, 1993). It doesn't depend on the radius, or even on G.

![clock gap](figures/fig4_clock_gap.png)

*Fig 4. Left: clocks carried around the spin axis fail to agree after one lap. Right: prograde and retrograde orbital
periods differ from Kepler's by +2πa and −2πa at every radius.*

### 1.7 What this changes from the last page

`jet-river` said full spin needed "cells that carry a rotation and a lattice that curves". The first half holds. The
second half needs a correction: the lattice the cells move through stays flat. What curves is the slice the cells
themselves see, the curvature 3R = 2Ma²(3 sin²θ − 2)/r⁵ found there. The input list is the same as for the non-spinning
case, a flat lattice and cells that are free fallers. The new structure is that each cell carries a turn as well as a
velocity.

Still not covered: gravitational waves. A wave is a ripple that travels, and a fixed lattice with fixed cell maps
can't carry one. It would need the cell maps themselves to change in time.

---

## Part 2. Every idea of a frame

### 2.1 Pass one: our own pages

The precession page is where frame dragging first came up. Its swap rule says b ∘ a = X⁻¹(a ∘ b)X: swapping the
inputs turns the answer and keeps its size. Over a loop that becomes **turn = d × area**. The precession page then
named frame dragging as "the only precession here that the theory can't produce on its own". Part F of the gravity
page closed that gap. This page closes it for full spin.

Here are all 17 frame ideas on our pages, sorted by the dial (survey of the repository):

![frame family](figures/fig5_frame_family.png)

*Fig 5. Every idea of a frame in our work, by dial setting. Left, d = +1, turning: Foucault, Coriolis, the tumbling
body, Larmor and Berry, Fuller's coil. Middle, d = 0, the flat middle: Galilean boosts (nothing precesses), the river
cells, the acoustic metric. Right, d = −1, relativity: Thomas–Wigner, geodetic, Lense–Thirring, rotating BTZ, Kerr's
clock effect. Purple: added on this page.*

| idea | where on our pages | dial | tie to frame dragging |
|---|---|---|---|
| swap rule, turn = d × area | precession §1 (A1–A7, Lean `swap_is_turn`) | any | the parent rule |
| holonomy ω (an output reading) | twisted-boosts; affine-names §6 | any | frame dragging's turn is a holonomy |
| Thomas–Wigner rotation | precession §2; twisted-boosts (36/36) | −1 | added to the drag in the double pulsar |
| geodetic (de Sitter) precession | twisted-boosts GP-B | −1 | frame dragging's sister effect |
| Foucault's pendulum, Coriolis | precession §4 (C5–C7) | +1 | turning frame, the everyday version |
| tumbling body's phase | precession §5 (D1–D3) | +1 | none |
| torqued gyroscope | precession §6 (E1–E6) | +1 | the 90° of the swap rule |
| Larmor precession, Berry phase | precession §7 (F1–F4) | +1 | half the d = +1 twist |
| Fuller's coil twist | fusion-magnets §3.5 | +1 | none |
| rotating BTZ drag (Fermi–Walker cross-check) | precession §8 (41/41) | −1 | a 2+1 frame dragging, exact |
| Lense–Thirring from boosting Newton | gravity-from-theory Part F (26/26) | −1 | derives frame dragging |
| river swirl 2J/r³ | jet-river J7 | 0 and −1 | slow-spin frame dragging |
| teleparallel tetrad, torsion | gr-bridges Part T | 0 | the river's frame and its defects |
| gauge connection, Aharonov–Bohm | electromagnetism | commuting case | same loop law, phases commute |
| Purcell: magnetism from a moving frame | electromagnetism | −1 | B as a frame effect |
| Bargmann phase m·area/ħ | bridge-survey (prediction) | 0 | the d = 0 version of the loop law |
| perihelion turn ≈ 6 × area/c² | river-rebuild §6 | −1 | shaped like the swap rule (observation) |

Missing before this page: Larmor's theorem, analogue gravity (spinning fluids), the Sagnac effect, and the Hannay
angle (cited once, never used).

### 2.2 Pass two: the literature, and the three look-alikes outside gravity

**Larmor's theorem (1897).** A charge in a weak magnetic field moves exactly as it would with no field, watched from
a frame turning at −qB/2m, up to terms of order B². Our check: the leftover error grows 4.18 times when B doubles, as a
B² error should, while without the turn it grows only 2 times (F1a).

- Our tools express it: the magnetic force qv × B has the Coriolis form 2m v × Ω, and Coriolis is the swap rule at
  d = +1 (precession page).
- The ½ in −qB/2m is the same ½ as in Part F's gyroscope rule, Ω = −½ curl h. In gravity, the field B_g = curl h plays
  the magnetic role. This is the known gravitational Larmor theorem (Mashhoon). Ours: the line through the swap rule.
- **Lead?** Modest. It explains a factor we already had, rather than opening something new.

![larmor and foucault](figures/fig6_larmor_foucault.png)

*Fig 6. Left: a charge orbiting in a magnetic field (blue) and the field-free orbit watched from a frame turning at
Larmor's rate (dashed orange) lie on top of each other. Right: Foucault's turn against latitude (green). Earth's own
frame dragging (dashed purple) has the same shape, at 4.6 × 10⁻¹⁰ of the size.*

**Spinning fluids (analogue gravity).** Sound in a moving fluid travels on the metric −c²dt² + (dx − u dt)². That is
exactly our river, with the fluid's velocity as the river's (Unruh 1981; Visser 1998). A draining, spinning bathtub
has u_r = −A/r and u_φ = B/r. Sound is trapped inside r = A, and between r = A and r = √(A² + B²) waves must move
with the spin, as light must inside a spinning hole's ergoregion. Spinning water tanks have shown waves being
amplified this way (Torres et al., 2017).

- **Ours, the field-law side:** for any 2D flow, the energy density the river's cells measure is the Q layer of the
  strain, G_nn = det(sym ∇u) (F2a). For the bathtub that is −(A² + B²)/r⁴, which is negative (F2b). So the vortex
  drags waves like a spinning star, but in gravity's terms it would need negative energy. The analogy is in the
  velocities, not in the field law.
- **Ours, the Kerr side:** a fluid has only a velocity per point, no turn. By K3a, velocity alone can't make Kerr.
  This matches the known result that fluid vortices reproduce Kerr only in the equatorial plane (Visser & Weinfurtner,
  2005). The turn is the reason.
- **Lead?** Yes, a small one. "Fluids lack the turn" is a sharp, testable statement in our language. A system whose
  cells carry an internal turn, such as a fluid with a spin field (micropolar fluids, superfluid helium-3 textures),
  is where a full Kerr analogue might live. That is untested.

![bathtub](figures/fig7_bathtub.png)

*Fig 7. Left: streamlines of a draining, spinning bathtub, with its sound horizon (black) and ergoregion edge (red
dashed). Right: the vortex's swirl rate B/r² and a spinning mass's swirl rate 2J/r³ both drag waves. The vortex's Q
layer is negative, though, while gravity's empty space has Q = 0.*

**Foucault's pendulum as a geometric phase.** Foucault's turn, 360° × sin(latitude) per sidereal day, is the
holonomy of carrying a direction around a circle of latitude on the sphere: d = +1 times the enclosed area. This was
already checked on our precession page (C5–C7), with the Berry and Hannay phases as its family (Berry 1984; Hannay
1985).

- New here: Earth's own frame dragging adds a vertical turn 2 sin(lat)·GJ/(c²R³). That is the fraction
  2(I/MR²)(GM/Rc²) = 4.6 × 10⁻¹⁰ of Foucault's turn, the same at every latitude, about 155 mas/yr at 45° (F3a).
- **Lead?** No. It is a nice number, but known physics with nothing new to test.

**Sagnac effect (from the survey's gap).** Light sent both ways around a turning loop returns with a time difference
4ΩA/c², a frame effect proportional to area: the same loop law again. The Kerr clock effect in §1.6 is gravity's
Sagnac effect, and here it is the cell clock's gap.

- **Lead?** Yes. The clock gap is a direct, quantitative reading of "the turn each cell carries". If the theory can say
  why the gap is a sin²θ per radian, it would derive the turn instead of finding it. (Done in Part 3.)

## Part 3. Deriving the turn (round 7)

Asked by Brennan, who answered "Yes" to: "Should I try to derive the turn's size from the dial, so Kerr's clock gap
becomes an output instead of something we found by testing?"

**Result: the turn's shape and size come out of four short steps, and the field equations then only confirm.** It
takes two assumptions, one known law, one far-field match and one reading. Kerr's a = J/M is an output, and so is its
clock gap, 2πJ/M per lap. Nothing is fitted. What is still put in by hand: the flat spheroidal lattice, and the speed
law v² = 2Mr/ρ² from Part 1, which is still matched by the field equations rather than derived.

![derivation chain](figures/fig8_turn_derivation.png)

*Fig 8. The chain. Blue: the area law gives the shape. Purple: the circulation law gives the rate in each shell.
Green: frame dragging far away gives the size. Orange: the rim at light speed ties the lattice's ring to the turn.
Red: the output, Kerr, with every empty-space field equation satisfied.*

### 3.1 Shape: the area law (D1)

Suppose each lattice spheroid turns as a whole, at its own rate Ω(r) (assumption 1: rigid shells). At the dial's flat
middle, the turning shell's time leg is the lattice time leg pushed sideways by the shell's speed: one more shear.
Around a latitude circle, that leg's clock gap is exactly 2 × Ω × (the circle's area). This is Sagnac's law, known,
and at q = 0 it holds exactly, with no γ factor (D1a). The circle's area is π(r² + b²) sin²θ, so

  turn f = Ω(r) · (r² + b²) · sin²θ.

**The sin²θ shape is just the area of a latitude circle.** And Part 1's turn amount, α = a sinθ/√(r² + a²), turns out
to be exactly the speed of the shell at that point (D1b). So the "turn" is the cell being pushed along the turning
shell's own time leg.

### 3.2 Rate per shell: the circulation law (D2)

For any turn function, the cells fall along cones of fixed θ (D2a). A latitude ring of cells therefore stays a
latitude ring as it falls. The ring's clock gap changes at the rate −2π u^r ∂f/∂r (D2b). Assumption 2 is that a
falling ring keeps its gap. This is a circulation law, like Kelvin's for a swirling fluid. It holds exactly when f
doesn't depend on r. With 3.1 that forces

  Ω(r) = κ / (r² + b²),  f = κ sin²θ.

### 3.3 Size: frame dragging far away (D3)

Far out, this cell gives g_tφ → −2Mκ sin²θ / r. Part F derived frame dragging's far field, −2J sin²θ / r, and checked
it against Gravity Probe B and LAGEOS. Matching the two gives **κ = J/M**.

### 3.4 The ring: the rim at light speed (D4)

The innermost shell, r = 0, is a flat disk of radius b. Its rim moves at κ/b. Reading the dial's edge as the fastest
a shell can go puts the rim at light speed, which gives **b = κ**. The lattice's ring and the turn are then the same
length, and every other shell moves slower than light: 1 − α² = ρ²/(r² + κ²) ≥ 0, which is zero only on the ring.
This matches a known picture of Kerr's source: Israel (1970) read it as a flat disk of radius a spanning the ring,
with its material reaching light speed at the rim (and with negative energy density, so it isn't a realistic star).
Our r = 0 shell is a lattice shell in flat space, not Israel's material disk. The match is in the radius and the rim
speed.

### 3.5 The field equations confirm, and every step matters (D5, D6)

The derived cell is identical to Part 1's tested one, with a = J/M (D5a). It solves every empty-space field equation,
max |Ricci| = 5 × 10⁻⁵⁰ (D5b). Change any step and empty space breaks (D5c):

![each step matters](figures/fig10_each_step_matters.png)

*Fig 10. The largest field-equation error for the derived chain (green) and for five changes to it (red). Those are
a rim faster or slower than light, two ways of letting the clock gap change between shells, and a gap that drifts
with r. Log scale; the dotted line is the 10⁻³⁰ "zero" threshold.*

There is also a bonus, which is a check we didn't aim for. At the horizon, Kerr's dragging rate −g_tφ/g_φφ equals the
derived shell rate a/(r₊² + a²) (D6b). That is Kerr's known horizon spin, Ω_H. It is also the rate of Carter's frame.
The t–φ block closes exactly there: g_tφ² − g_tt g_φφ = sin²θ (r² − 2Mr + a²) (D6a).

![turning shells](figures/fig9_turning_shells.png)

*Fig 9. Left: the flat lattice, each spheroid turning rigidly at Ω(r) = a/(r² + a²). The arrows show each shell's
equator speed, and only the r = 0 disk's rim reaches c. Right: one curve carries the light-speed rim at r = 0, Kerr's
horizon spin at r₊, and frame dragging's 2J/r³ far away.*

### 3.6 What is still assumed, plainly

| step | status |
|---|---|
| flat spheroidal lattice | put in (Part 1.1); derived in Part 4 |
| speed v² = 2Mr/ρ² | put in, matched by the field equations (K3c); derived in Part 4 |
| shells turn rigidly | **assumption 1** |
| clock gap = 2 Ω × area | known law (Sagnac), exact at q = 0 |
| a falling ring keeps its gap | **assumption 2** (circulation law) |
| κ = J/M | derived, from Part F's far field |
| rim at light speed, b = κ | a **reading** of the dial's edge; the field equations confirm only b = κ works |
| Kerr, a = J/M, clock gap 2πJ/M, Ω_H | **output**, checked |

The next gap is the speed law. If it can be derived too, perhaps as a Gauss law on the spheroids the way rv² = 2Gm
came out of the jet river, then Kerr would follow from the lattice shape plus the two assumptions.

## Part 4. Deriving the speed law, and the lattice's shape (round 8)

Asked by Brennan, who answered "Yes" to: "Should I try to derive the speed law next, so Kerr follows from just the
lattice's shape and those two assumptions?"

**Result: one move gives both remaining inputs.** Take the non-spinning source and turn its height on the axis by
the dial's rotation number i, moving it from 0 to i·a. Its distance from each point then becomes the complex number
D = r − i a cos θ. Three things fall out at once. The real part of D is the lattice label r, so the squashed spheres
are the level sets of Re D. D = 0 is a ring of radius a, the same ring the turn needed in Part 3. And Newton's
potential of the moved source, Re(M/D) = Mr/ρ², is exactly the speed law v² = 2Mr/ρ² through the jet river's rule
v²/2 = −Φ. This is a reading, not a proof: it assumes the spin is that one shift. Within the reading, nothing else
is put in except the mass M.

![complex shift](figures/fig11_complex_shift.png)

*Fig 11. Left: the one move. The source's height is turned a quarter turn in the complex plane, from 0 to i·a, and
not slid along the real axis by a (grey, which fails). Right: what that does in real space. Purple shows the level
sets of Re D, which are the lattice's squashed spheres. Red marks D = 0, the ring of radius a, which is the edge of
the flat disk r = 0. Blue shows the cell speed² = 2 Re(M/D).*

### 4.1 The complex distance (S1)

Move a point source up the spin axis by i·a. Its squared distance x² + y² + (z − ia)² is exactly (r − i a cos θ)² in
the lattice's coordinates (S1a, also `complex_distance` in Lean). So the lattice's spheroids are not a choice: they
are where the shifted source's distance has a fixed real part. The distance vanishes on the ring r = 0, θ = π/2,
radius a (S1b).

**Two routes meet at one ring.** Part 3 found the turn's ring at radius b = κ = J/M, from frame dragging. Here the
ring sits at radius a, from the shift. They agree only when a = J/M, and that is Kerr's own relation.

### 4.2 The speed law (S2, S3)

Newton's potential of the moved source is Re(M/D) = M r/(r² + a² cos²θ) (S2a). Far away it is M/r, the ordinary
potential of mass M. The jet river's rule says the cell's speed sets the potential, Φ = −v²/2 (jet-river J6b). So
v² = 2Mr/ρ², exactly the law Part 1 found by testing (S2b). Re(M/D) is an exact flat-space potential: Laplace's
equation holds everywhere off the disk, not just far away (S3). That is the empty-space form of the jet river's
trace law.

### 4.3 Being a potential of mass M is not enough (S4)

Many flat-space potentials have mass M. A point at the centre, a point slid up the axis by a real a, and a point
moved by i·a/2 are all harmonic, and all tend to M/r far away (S4a). With the lattice and turn from Part 3, every one
of them breaks empty space. Only the i·a shift works (S4b).

![only the i shift](figures/fig12_only_the_i_shift.png)

*Fig 12. The largest field-equation error for four speed laws, each a flat-space potential of mass M. Only the
source moved by i·a gives empty space. The i·a/2 shift fails because its ring (radius a/2) is not the turn's ring.*

### 4.4 What is assumed now, plainly

| piece | status |
|---|---|
| spin = the source's height turned by i, from 0 to i·a | **reading** (the dial's rotation number); known as the Newman–Janis shift |
| lattice's squashed spheres, and their ring | **derived** from the shift (S1) |
| speed law v² = 2Mr/ρ² | **derived**: the shifted source's potential, through the jet river's v²/2 = −Φ (S2, S3) |
| shells turn rigidly; a falling ring keeps its clock gap | the two assumptions of Part 3 |
| a = J/M | derived twice (frame dragging far away, and the two rings meeting) |
| Kerr | output, checked to 50 digits |

**Lead, on its own merits.** The shift by i·a is known: Newman and Janis (1965) used it to get Kerr from
Schwarzschild, and Appell (1887) and Newman (1973) read it as a source at a complex position. Kerr's Kerr–Schild
function 2Mr/ρ² is also known to be this harmonic function. Why the shift works has never had a settled explanation.
In the theory, i is not a trick: it is the dial's rotation setting. Spin enters the river twice, as a shear at the
dial's flat middle (the turn, q = 0) and as a quarter turn of the source (q = +1). If the dial can say why those two
settings go together, that would explain the Newman–Janis shift. That lead is worth following. Part 5 derives the turn from
the shift as well, so one move gives everything.

## Part 5. The turn from the same shift (round 9)

Asked by Brennan: "Yes fix the sentences and try to get the turn from the same shift".

**Result: the shifted source's own light rays carry the turn.** Through each point, the source moved to i·a sends two
real light rays, one out and one in. Both lean sideways by the same amount, and that sideways part is −a sin²θ dφ:
exactly the turn round 7 derived from rigid shells and the circulation law. So the one move of Part 4 gives the
lattice, the speed and the turn. Together these give the whole spinning black hole, with every field equation of empty
space holding to 50 digits.

![shifted rays](figures/fig13_shifted_rays.png)

*Fig 13. Rays in the plane of the equator. Left: no spin, so the rays run straight out from the centre. Right: the
source is moved by i·a, and every ray is a straight line tangent to the ring of radius a. The lean, a/√(r² + a²), is
the speed of round 7's turning shell.*

### 5.1 The rays of a shifted source (T1–T3)

For a real source, the ray through a point is the unit direction to the source. For the shifted source that direction,
E = grad D, is a complex unit vector, E = p + i q. A real ray k must be a unit vector with k·E = 1, which means k·p = 1
and k·q = 0. Those are a circle and a plane, so there are exactly two solutions: k = (p ± q × p)/|p|² (T1). With no
spin, q = 0 and both are the ordinary radial ray.

On the lattice, the outgoing ray is k = (ρ²/(r²+a²)) dr − a sin²θ dφ. That is a radial part along the shell's normal,
plus a sideways part (T2). The ingoing ray has the same sideways part and the opposite radial part (T3). The lean is
exactly α = a sin θ/√(r² + a²), the speed of round 7's turning shell (T2b). In the plane of the equator, every ray is
a straight line tangent to the ring of radius a.

### 5.2 The turn (T4)

The river rule reads: push each cell along its ray, without the ray's radial part. With no spin the ray is purely
radial, and this is the plain river, w = dt. With the shift it gives w = dt − a sin²θ dφ. That is round 7's turn
f = κ sin²θ, with κ = a. This rule is a reading. The calculation behind it is T2.

### 5.3 The whole black hole from one move (T5)

The same rays and speed give Kerr in its ray form, flat + v²(dt + k)(dt + k). This is the known Kerr–Schild form, and
it solves every empty-space field equation (T5a). Other light rays with the same speed law fail: a purely radial ray,
and a ray with half the sideways part (T5b). Built from the shift alone, the river cell is Part 1's Kerr cell (T5c).

![one move](figures/fig14_one_move.png)

*Fig 14. One move gives everything. Turning the source's height by i gives the complex distance D. From D come the
lattice (Re D), the speed (Re M/D) and the turn (the light rays with k·E = 1). Together they give Kerr. Frame dragging
far away sets a = J/M.*

### 5.4 What is assumed now, plainly

| piece | status |
|---|---|
| spin = the source's height turned by i, from 0 to i·a | **reading** (the dial's rotation number); known as the Newman–Janis shift |
| lattice, speed law, turn | **derived** from the shift (Parts 4 and 5) |
| push the cell along its ray without the radial part; v²/2 = −Φ | river rules (the second from the jet river) |
| a = J/M | from frame dragging far away (Part 3) |
| Kerr | output, checked to 50 digits |

The two assumptions of Part 3 (rigid shells, a kept clock gap) are no longer needed. They now come out as properties
of the shift's rays: the ray's lean is the rigid shell's speed, and the turn has no r in it.

**Lead, on its own merits.** In this theory, spin is one setting of the dial applied to the source's position. Rotation
numbers (q = +1) act on where the source is. Shears at the flat middle (q = 0) act on each cell. The first produces
the second. That two-setting structure is a candidate explanation of why the Newman–Janis shift works, and it is worth
following. The next tests: charge (Kerr–Newman), where the shift is known to work and the river should follow (done
in Part 6); then a case where the shift isn't known to work, such as two bodies.

## Part 6. Charge on the same shifted source (round 10)

Asked by Brennan, who answered "Yes" to: "Should I test the same move with charge next, the charged spinning black
hole?"

**Result: the same move works with charge.** Put a charge Q on the shifted source and change nothing else. The charged
spinning black hole (Kerr–Newman) comes out. Einstein's equations, now including the electric field's own energy, hold
to 50 digits, and so do Maxwell's equations. The shifted charge also turns out to be a magnet, with the electron's
magnetic ratio, g = 2.

![charge](figures/fig15_charge.png)

*Fig 15. Left: the shifted charge's one complex potential, Q/D. Its real part (blue) is the electric potential. Its
imaginary part (orange) is a magnetic potential, whose far field is a magnet of moment Qa. Right: the largest
Einstein–Maxwell error for the shifted charge (green) and for three rivals (red). The rivals fail at about 10⁻³ against
10⁻⁵⁰. Their errors are small only because the charge is small, Q = 0.5.*

### 6.1 No spin first: the field's energy in the jet river's Gauss law (C1)

The jet river's Gauss law says r v² = 2 m(r), where m(r) is the mass inside radius r. A charge's electric field holds
energy Q²/(2r) outside radius r, so the mass inside is m(r) = M − Q²/(2r). That gives v² = 2M/r − Q²/r² (C1a). With
the field A = (Q/r) dt, this river is the charged black hole (Reissner–Nordström). Both sets of equations hold (C1b).

### 6.2 Shift the charge: electric and magnetic together (C2)

Moved by i·a, the charge's potential Q/D is one complex number. Its real part is the electric potential Qr/ρ². Its
imaginary part, Qa cos θ/ρ², is a magnetic potential, which far away is a magnet of moment Qa (C2). The black hole's
spin is J = Ma, so the ratio of magnetism to spin, in the units used for particles, is g = 2. That is the electron's
value. For Kerr–Newman this has been known since Carter (1968). Here it comes from the same quarter turn that gives
the spin.

### 6.3 The speed law with charge (C3)

Shift both terms of 6.1. The mass term becomes 2 Re(M/D), as in Part 4. The field-energy term becomes |Q/D|², the size
of the complex potential, not its real part. That is because the field's energy counts the electric and magnetic
parts together: (E² + B²)/8π = |E + iB|²/8π. Together, v² = 2 Re(M/D) − |Q/D|² = (2Mr − Q²)/ρ² (C3).

### 6.4 The charged spinning black hole (C4, C5)

Take that speed, the lattice and turn from the shift, and the field A = Re(Q/D) (dt − a sin²θ dφ). That field is the
charge's potential carried on the cell's turned time leg. Einstein's equations with the field's energy hold, and so do
Maxwell's (C4). The horizon sits where r² − 2Mr + a² + Q² = 0, Kerr–Newman's (C4b). Three rivals fail (C5): taking
the field energy from Re(Q/D)² instead of |Q/D|², dropping the turn from the field, or leaving the charge at the
centre.

### 6.5 What this adds, on its own merits

- **A real test passed.** The shift was built for spin, and charge was not used to set it up. It then gives the
  charged case with no new input except Q. The one new choice, |Q/D|² for the field's energy, is forced: the real
  part fails.
- **Known vs ours.** The charged spinning black hole was found in 1965 with this same complex trick (Newman and
  co-workers), and its g = 2 is known (Carter 1968). Ours: the river form built from the shift, the field-energy term
  as |Q/D|² through the jet river's Gauss law, the field carried on the turned leg, and the rivals failing.
- **Lead.** The electron's g = 2 is also what Dirac's equation gives, and people have long noticed that Kerr–Newman
  shares it. In this theory both the spin and the magnetism come from one quarter turn. That is worth following,
  carefully. With an electron's mass, charge and spin there would be no horizon at all, since the spin and charge are
  far too large for the mass. That would leave a bare ring about the electron's Compton size. So this is a shared
  number, not a model of the electron.
- **Next test, where the shift is not known to work:** two bodies, or a moving one.

## Part 7. Two bodies (round 11)

Asked by Brennan, who answered "Yes!" to: "Should I try a case where the shift isn't known to work, like two bodies?
It's a much harder test, since even standard gravity has no exact answer for two orbiting bodies." He then asked how
we could verify it, and whether the gravitational waves detected from colliding black holes could be used.

**Result: the river as it stands can't hold two bodies, and the reason is plain.** One shear per cell only ever makes
curvature of a special kind, and two bodies' field isn't of that kind. Giving each cell two shears, one per body, also
fails, however the two are combined. What does add exactly, for bodies at rest, is a different quantity: the log of
each cell's clock rate. Along the axis between the bodies, those combine by Einstein's velocity-addition rule. The
exact field then needs a strut to hold the bodies apart, and the strut's push is Newton's inverse-square law once they
are far apart.

![two bodies](figures/fig16_two_bodies.png)

*Fig 16. Left: two black holes held at rest on one axis, masses 1 and ½, in Weyl's coordinates, where each is a rod
of length 2M. Blue: ψ = ½ ln(1 − v²), the log of the clock rate, which is the two bodies' ψ added. Orange: the strut
between them. Right: the strut's push from the exact field (solid) against Newton's M₁M₂/ℓ² (dashed). The red dot is
measured from the metric in check B4.*

![one shear cannot](figures/fig17_one_shear_cannot.png)

*Fig 17. Left: how far the curvature is from special, where 0 means two of its three stretch rates are equal. One body
gives 0. Two bodies give values up to 1, except close to one of them. Right: the leftover in the empty-space equations
when each body brings its own shear, for three ways of combining the two shears. Either body alone leaves none.*

### 7.1 Why one shear per cell can't (B1)

Curvature in empty space stretches and squeezes in three directions, at three rates. For one body, at rest or
spinning, two of the three rates are equal: the curvature is "algebraically special" (B1a). That is not an accident of
these examples. Every empty space made by one shear per cell is special: the Goldberg–Sachs theorem, which applies to
every Kerr–Schild metric, the form the river takes. Two bodies at rest, in gravity's exact solution, are not special:
the three rates all differ, by an order-one amount away from the bodies (B1b). So no choice of one shear per cell can
give two bodies.

### 7.2 Two shears, one per body, don't combine (B2)

Each body alone is exact. Put both in, each with its own shear, and combine them three ways: add the two metrics, add
the two shears, or compose the two cell maps one after the other. All three leave a leftover in the empty-space
equations (B2). It grows as M₁·M₂: make both masses ten times smaller and the leftover is a hundred times smaller.
That product is exactly what the interaction between two bodies looks like, so the leftover is the part of gravity
that comes from the bodies acting on each other, which none of the combinations supplies.

### 7.3 What does add: clock rates, by velocity addition (B3)

For bodies at rest on one axis, Weyl found in 1917 that one quantity does add exactly: ψ = ½ ln(1 − v²), the log of
the river cell's clock rate. In Weyl's own flat coordinates, ψ is a sum of Newtonian potentials, one per body, each
spread along a rod of length 2M (B1b). On the axis, each body's ψ is −artanh(u), with u = M/(z − c) its potential.
Adding artanh's is how speeds combine in relativity, so the two bodies' u's combine by Einstein's rule
(u₁ + u₂)/(1 + u₁u₂) (B3; Lean `velocity_addition`). So for two bodies the quantity that adds is a clock-rate angle, not
the cell's speed. In the dial's terms, one body's cell sits at the flat middle, where shears add, while two bodies
combine at the boost setting, where angles add. That last sentence is a reading, ours.

### 7.4 The strut, and Newton's law (B4)

Two bodies can't stay at rest without something holding them apart, and the exact field shows it. On the axis between
them, the field has the same cone-shaped defect at every height: a strut. Its push is (e^−γ − 1)/4, which equals
M₁M₂/(ℓ² − (M₁ + M₂)²) exactly, with ℓ the distance between the centres (B4; Lean `strut_push`). Far apart that is
Newton's M₁M₂/ℓ², so the inverse-square law comes back out of the exact two-body field. This is known: Bach and Weyl
(1922) and Israel and Khan (1964). It is the cleanest place where this page meets the inverse-square law of the thread's
title.

### 7.5 What this means, on its own merits

- **Fails, cleanly:** the river with one shear per cell is a one-body picture. That is a structural fact, not a
  tuning problem, so two bodies need a new ingredient. Adding a second shear isn't it.
- **What works:** for bodies at rest, adding clock-rate angles (Weyl) gives the exact field, with a strut whose push is
  Newton's force. All of that is known. Ours is checking it in this page's terms, the velocity-addition reading on the
  axis, and the failure of each way of combining shears.
- **Lead:** the shift's potential M/D also turns up on the axis of the exact spinning black hole, in the form gravity
  uses for fields that don't change in time. On the axis, in Weyl's coordinate z, Ernst's potential for Kerr is
  (1 − u)/(1 + u) with u = M/(z + ia): a point source shifted by i·a, the same u as in B3 with c = −ia (B5). A known method
  (Sibgatullin's) turns values on the axis into a whole field, and the exact two-spinning-body family is known
  (Kramer and Neugebauer 1980). So the natural next test with spin: combine two shifted sources by velocity addition
  on the axis, build the field, and see whether it is two spinning black holes with a strut. That is checkable, and
  a pass would make the shift a rule for combining bodies, not only for making one.
- **How to verify against nature.** The detected gravitational waves are the real test, and they need moving bodies:
  1. two bodies at rest: done here, by the known route;
  2. two spinning bodies at rest on one axis: the lead above;
  3. slow orbits, the corrections physicists call post-Newtonian, checked against binary pulsars, whose orbits shrink
     as they give off waves (done in Part 8);
  4. the waveform itself, against the LIGO–Virgo–KAGRA signals, hundreds of them since GW150914 in 2015. The
     ringdown at the end of each comes from one spinning black hole, which this page already matches.

## Part 8. Two spinning bodies in orbit (round 12)

Asked by Brennan, who said "Yes" to two spinning bodies on one axis, then asked instead: "Please try two spinning
bodies in a spinning orbit around each other". He added: "I'd be interested to see if the spins start lining up to
each other if they start different", "And then we need to simulate the orbits decaying", and "And at some point the
merge will create a graviton wave".

**Result: three river rules give gravity's field for two bodies in a slow orbit, exactly.** The rules: each body
sends out its own far field at light speed, clock rates multiply, and each mass is seen through the other body's
clock. Leave any one out and the field is wrong. From that field the oval orbit turns and the spins wobble at the
known rates. Adding the waves' known back-push (put in, not yet derived from the rules), the orbit shrinks at the rate
measured in real binary pulsars. All of this downstream is standard physics. Spins that start pointing different ways **do not line up on average**.
They sort into two families instead. The wave's run-up and the merged hole's ring-down come out. The merger in
between needs a supercomputer simulation and is not done here.

![three rules](figures/fig18_three_rules.png)

*Fig 18. Left: the three rules for two moving bodies. Right: the largest mismatch with gravity's published slow-orbit
field (Blanchet, Faye and Ponsot 1998) at random points, to 40 digits. With all three rules it is zero to rounding.
Dropping any one leaves an order-one error.*

### 8.1 The field of two moving bodies, from three rules (O1)

Gravity has no exact two-body solution, but when the bodies move slowly compared with light it has a field good to
the first correction beyond Newton. That is the "1PN" field that fixes the Einstein–Infeld–Hoffmann motion, and it is
what binary pulsars test. Check O1 builds it from three rules and compares it term by term, at random configurations
with 40-digit arithmetic:

1. **A, light-speed lag:** each body's own far field is the field of one moving mass, sent out at light speed. The
   field at a point comes from where the body was a moment ago, not where it is now.
2. **B, clock rates multiply:** the rate from body 1 times the rate from body 2. This is the same rule Part 7 found for
   bodies at rest (B3). In the slow-orbit form it supplies the cross term 4U₁U₂ exactly (Lean `clock_multiply`).
3. **C, dressed masses:** each mass is seen through the other's clock, m₁ → m₁(1 − m₂/r₁₂).

The frame-dragging part of the field is the one from Part F: −4mv/r for each moving body, plus each spin's own drag.
All three rules together match to 2 × 10⁻⁴¹. Without A the mismatch is 0.06, without B 0.25, without C 0.25.

### 8.2 The spins wobble (O2)

Each spin turns in this field by two pieces already on our pages: the geodetic turn, (3/2) v × ∇U (1 from curved space
plus ½ from Thomas), and half the curl of the frame-dragging field. Together they give the known spin-orbit and
spin-spin precession exactly (Barker and O'Connell; Kidder 1995), symbolically. On a circular orbit, spin 1 turns
about the orbit's axis at (2 + 3m₂/(2m₁)) L/r³.

### 8.3 The orbit turns, shrinks and rounds off (O3, O4)

![orbit decays](figures/fig19_orbit_decays.png)

*Fig 19. Left: the first six laps of an oval orbit of two equal masses. Red: the closest points, turning about 24° per
lap. Middle: the separation over the whole run. The waves carry energy off, so the swing between near and far dies out
and the bodies spiral together. Right: real binary pulsars, observed over predicted.*

The motion follows from the field, plus the waves' back-push on the bodies (the 2.5PN term), integrated numerically.
The back-push is not derived from the rules here: it is the known term (harmonic gauge), put in, with its sign checked
by energy balance.

- **Turning:** the closest point moves on by 6πM/p per lap, p the orbit's width parameter, whatever the mass ratio
  (O3: 1.0014 and 1.0005 of that, for mass ratios 0 and 1, at p ≈ 840 M). This is Mercury's turning, for two bodies.
- **Shrinking:** the energy the orbit loses matches the quadrupole formula (32/5) ν²M⁵/r⁵ to 0.7% at r = 80 M and 2.4%
  at 40 M, closer the wider the orbit, as it should be (O3).
- **Against real pulsars (O4):** for the Hulse–Taylor pulsar the formula gives dP/dt = −2.402 × 10⁻¹², and the
  measured value is 0.9983 ± 0.0016 of that (Weisberg and Huang 2016). For pulsar B of the double pulsar the spin
  wobble comes out at 5.07° per year, and the measured value is 4.77 ± 0.66° (Breton et al. 2008).

### 8.4 Do two spins that start different line up? (O5)

![spins](figures/fig20_spins.png)

*Fig 20. Left: one pair, mass ratio 0.7, both spins 0.8 of the maximum, from 200 M down to 10 M. Both spins wobble
around the orbit's axis while the orbit's axis itself wobbles. Middle: 80 pairs started with spins pointing every which
way at r = 1000 M, and the angle between their spins at the start and at r = 10 M. Right: the angle between the
spins' shadows on the orbit plane at 10 M, split by which spin started nearer the orbit's axis.*

No, not on average. The mean cosine of the angle between the two spins was −0.04 at the start and −0.08 at the end:
no drift toward lined up. But the spins do sort. When the heavier body's spin starts nearer the orbit's axis, the
spins' shadows on the orbit plane end up close together (median 41°). When it starts farther, they end up nearly
opposite (median 138°). These are the spin-orbit resonances found by Schnittman (2004), and they come out of the
precession in 8.2 with no extra input. A quick version runs in O5 (40 pairs, mass ratio 0.9); `orbits.py` makes the
80-pair run in Fig 20.

### 8.5 The wave (O6)

![wave](figures/fig21_wave.png)

*Fig 21. The wave from a collision like GW150914. Blue: the run-up as the two bodies spiral in, at leading order.
Grey: the merger, which needs a full supercomputer simulation and is not done here. Red: the merged spinning hole's
own tone. The run-up uses the known quadrupole formula, not a derivation from the rules. The heights of the two pieces
are matched by hand; only their shapes and frequencies are computed.*

The waves that shrink the orbit are what the detectors see. A "graviton wave" in the request is this wave: the graviton
is its quantum particle, the way a photon is light's. The classical wave is what LIGO measures, and that is what is
computed here.

- **Run-up:** from 35 to 150 Hz the leading-order formula gives 0.17 s and 8.7 cycles, against about 0.2 s and 8 cycles
  published for GW150914.
- **Ring-down:** the merged hole (62 solar masses, spin 0.67) rings at 250 Hz and fades in about 4 ms. That is Kerr's
  fundamental tone of a spinning black hole, the one this page builds (Part 5). The number is from the known fit of
  Berti, Cardoso and Will (2006), not computed here.
- **The merger:** the last few orbits and the collision itself are fully nonlinear. Only numerical relativity can do
  them, so that piece is a gap here.
- **The real data:** the detector data could not be downloaded from this session, so no direct comparison with the
  recorded signal was made.

### 8.6 What this means, on its own merits

- **Ours:** the decomposition into three rules, and especially that the clock-multiplication rule found for bodies at
  rest (Part 7) carries over to moving bodies and gives the interaction term with nothing added. Reading the light-speed
  lag as "each body's field is sent out" is the obvious picture, but writing the 1PN field as exactly these three rules
  is not something we found stated.
- **Known:** everything downstream: the 1PN field itself, the motion, the precession, the decay, the pulsar numbers, the
  resonances and the ringdown. These are textbook general relativity, and they come out because the field is exact at
  this order.
- **Put in, not derived:** the waves' back-push on the orbit and the quadrupole formula for the wave. Both should
  follow from rule A carried one step further, and that is a separate test.
- **Tested since (Part 9):** whether the three rules keep working at the next order, 2PN. They don't: body 1's clock
  rate comes out wrong at x³. Rule B turns out to shape only the field between the bodies, not their motion.
- **Still out of reach:** the merger. Matching the detected waveforms end to end would need the rules at higher order
  plus numerical relativity.

## Part 9. The next level of detail (round 13)

Asked by Brennan, who tapped "Test it" on: "Should it test the next level of detail? That's the first place the three
rules could fail."

**Result: the three rules fail at the next level of detail (2PN).** They give body 1's clock rate exactly at the first
correction beyond Newton, and exactly at the next one when body 1 is too light to pull back. With two real masses the
next correction comes out wrong, by ν(15/2 + 9Δ/2 − 23ν/3) in the x³ coefficient (ν = m₁m₂/m², Δ = (m₂ − m₁)/m, with
m = m₁ + m₂). The test also shows why. Each body only feels the other body's own field, lagged and dressed. The clock
rule B, which Part 8 needed for the field, drops out of what the bodies feel. What general relativity adds at this
order must come from the field between the bodies, or from the field of a body that is speeding up, and the rules
have neither.

![next order](figures/fig22_next_order.png)

*Fig 22. Left: body 1's clock rate z₁ on a circular orbit, term by term, for body 1 at ⅓ of the mass. The Newton and
1PN terms match general relativity. The 2PN term is off by 1.62, with the wrong sign. Middle: the 2PN term across the
mass ratio. Grey is general relativity, the red dots are the rules, the dashed curve is the fitted miss, and the open
dot was computed after the fit and lands on it. At m₁ → 0 the two agree. Right: what each body feels, and the two
candidates for what is missing.*

### 9.1 The test: a clock on the orbit (O7)

Fields depend on the coordinates you write them in, so comparing two fields term by term needs the same coordinates on
both sides. Part 8 had that, because the 1PN field is published in harmonic coordinates (Blanchet, Faye and Ponsot
1998). The 2PN field is in the same paper, but its source couldn't be fetched from this session. So this round compares
something that doesn't depend on coordinates: the rate of a clock riding on body 1 around a circular orbit, as a
function of the orbit's frequency. Write it z₁ = 1/u^t and expand it in x = (mω)^(2/3), so x is about v².

- **General relativity's answer:** the first law of binary mechanics (Le Tiec, Blanchet and Whiting 2012) gives
  z₁ = 1 + ∂E/∂m₁ − ω ∂J/∂m₁ at fixed ω. The orbit's energy E and angular momentum J are known to 2PN. Out comes
  z₁ = 1 − (m₂(2m₁ + 3m₂)/2) x − … through x³ (checked symbolically in O7).
- **The rules' answer:** build the field of two bodies on a circle from the rules (`next_order.py`). Each body's own
  field is the moving mass's field, lagged at light speed, with its own nonlinear part (a moving black hole's).
  Clock factors multiply, the drag and space parts add, and each mass is dressed by the field it sits in. Then let
  each body move on a free-fall path of that field, with its own infinite part removed by Hadamard's finite part, the
  method Blanchet, Faye and Ponsot use. Solve for the circular orbit order by order, and read off body 1's clock rate.
  Everything is exact fractions.

The method passes both checks it should pass. Through x² it gives general relativity's clock rate for every mass ratio,
as it must, since Part 8's field is exact there. Through x³ it gives the exact answer for a test mass (m₁ → 0), as it
must, since then body 1 only rides in body 2's own field, which is a black hole's.

### 9.2 The miss (O8)

With two real masses the x³ term is wrong. For body 1 at ⅓ of the mass the rules give +1205/2916 and general
relativity gives −3523/2916. Over seven mass ratios the miss is exactly

  rules − general relativity = m₁m₂(23m₁² − 50m₁ + 36)/3 = ν(15/2 + 9Δ/2 − 23ν/3)   (total mass 1).

That formula was fitted to six mass ratios (it is fixed by five, and the sixth agrees), then confirmed on a seventh,
m₁ = 2/5, computed afterwards (984/625 both ways). The two ways of writing it are one identity (Lean `two_pn_miss`).
The miss is already there at first order in the lighter mass: about 12m₁ when m₁ is small. So it is not a small
correction that only shows up for equal masses. It is the light body's own pull on the heavy one, at 2PN.

### 9.3 Why: the bodies don't feel rule B (O7)

Swapping rule B for its rival, adding the clock rates instead of multiplying them, leaves body 1's clock rate unchanged
through 2PN, to every digit. The reason is how the infinite part of a body's own field is removed. The leading product
of the two fields, the cross term 4U₁U₂, goes as 1/r₁ near body 1 times something smooth, and the finite part of that
averages to zero over directions. The higher products drop out of the clock rate too, as the run shows. So:

- **Rule B shapes the field between the bodies** (that is where Part 8's O1 needed it, at random points of space),
  **but not their motion.** "Each rule is needed" in Part 8 holds for the field, not for the orbits.
- **What a body feels** is the other body's own field, lagged (rule A) and with its mass dressed (rule C), plus its
  own field's finite part. At 1PN that is enough. At 2PN it isn't.
- **The space rule doesn't matter either:** multiplying or adding the two bodies' space stretch (the 4U₁U₂ term in the
  space part) gives the same clock rate (next_order.json, kappa0).

So whatever general relativity has at 2PN that the rules lack must survive the finite part at a body. Two pieces do,
and the rules have neither:

1. **Gravity's own field energy between the bodies, acting as a source.** In general relativity the field's energy
   also pulls. Between two bodies that gives a term that isn't a product of the two fields at the same point but is
   spread out between them, built from ln(r₁ + r₂ + r₁₂) (r₁₂ the distance between the bodies). A spread-out term
   doesn't drop out at a body the way a product does.
2. **The extra field of a body that is speeding up.** Each body here carries the field of a body moving at constant
   velocity. On a circle the bodies accelerate, and at 2PN the nonlinear part of the field feels that.

Which of the two closes the gap, or whether both are needed, can't be told apart without general relativity's 2PN field
itself, which this session couldn't fetch.

### 9.4 What this means, on its own merits

- **The three rules as they stand stop at 1PN.** That is an honest limit. The slow-orbit field and everything Part 8
  built on it (orbit turning, spin wobble, pulsar decay with the back-push put in) are untouched. Effects that need 2PN,
  like the last orbits before a merger, are not given by the rules.
- **Ours:** the test itself (the rules' field run through the clock-rate comparison), the miss formula, and the
  finding that rule B drops out of the motion. The extension of the rules to 2PN is also ours, and it involves choices:
  each body's own nonlinear field is a moving black hole's with no acceleration terms, and the parts of the field that
  aren't clock rates add.
- **Known:** the first law and the clock rate (Le Tiec, Blanchet and Whiting 2012; Detweiler 2008), the 2PN energy and
  angular momentum (Blanchet's review), Hadamard's finite part for point bodies (Blanchet, Faye and Ponsot 1998), and
  that 2PN gravity has field-energy terms spread between the bodies.
- **The next test:** add gravity's own field energy as a fourth rule (the spread-out term between the bodies) and see
  whether the miss closes. If it does, the rules become "fields are sent out, rates multiply, masses are dressed, and
  the field's own energy is a source", which is general relativity's own structure, written in the river's terms.

## Part 10. A fourth rule: the field's own energy (rounds 14 and 15)

Asked by Brennan, who said "Yes" to: "Should I add gravity's own field energy as a fourth rule and see whether it
closes the gap?" (round 14), and then "Yes" to "Should I try to derive the leftover quarter point term?" and "Also see
if you can find the paper" (round 15).

**Result: adding gravity's own field energy closes the 2PN gap exactly, at every mass ratio.** The rules plus the
field's energy rebuild general relativity's 2PN field piece by piece. For body 1 at ⅓ of the mass the miss in the x³
coefficient of its clock rate goes 1.62 → 0.44 → 0.13 → 0. Round 14 left the last 0.13 as one point-like term at each
body. Round 15 derived it: each body speeds up toward the other, its field shifts with it, and the −2V² term of
general relativity's clock field multiplies that shift by the body's own field. Round 14's difference D had left that
product out. With it the miss is 0 at all six mass ratios, a second, separate quantity, the orbit's size at a given
orbital rate, matches the known 2PN law at all six, and the clock field at a body equals the paper's own formula
term by term.

![fourth rule](figures/fig23_fourth_rule.png)

*Fig 23. Left: the miss in body 1's 2PN clock-rate term, step by step, for body 1 at ⅓ of the mass. Red: the three
rules. Orange: plus the closed-form pieces of the field's energy. Blue: plus the four integrals where the field acts on
the field. Green: plus each speeding-up body's own field squared (round 15), which brings the miss to 0. Middle: the same steps across the mass ratio; dots are computed,
curves are the formulas, and the open dots (m₁ = 2/5) were never used in a fit. Right: the picture. The field's energy
between the bodies is itself a source; the last piece sits at each body, from its own field shifting as it speeds
up toward the other.*

### 10.1 What general relativity has that the rules lack

General relativity's 2PN field of two point bodies is published in harmonic coordinates (Blanchet, Faye and Ponsot
1998, arXiv gr-qc/9804079). Round 14 wrote its formulas from the standard form; round 15 read the paper and confirmed
them against its Eqs. (2.4)–(2.5) and its values of Ŵᵢᵢ and the dressed mass at a body, Eqs. (4.7)–(4.8):
g₀₀ = −1 + 2V − 2V² + 8(X̂ + V_iV_i + V³/6), g₀ᵢ = −4Vᵢ − 8R̂ᵢ, gᵢⱼ = δᵢⱼ(1 + 2V + 2V²) + 4Ŵᵢⱼ, with X̂, R̂ᵢ and
Ŵᵢⱼ sourced partly by the field itself. The paper's equations are also on main, copied from its arXiv source
([`topics/gravity-from-theory/sources/bfp1998-2pn-metric.md`](../sources/bfp1998-2pn-metric.md)). Checked during the work on one moving body: these formulas give
the moving black hole that Part 9's rules use (that check is not in `checks.py` yet).

Subtracting the rules' field from that one leaves the difference D = general relativity − rules. Every piece of D is
the field's energy acting as a source, and they fall into two kinds:

- **Closed-form pieces** (`fourth_rule_closed.py`): the cross part of Ŵᵢⱼ, of R̂ᵢ and of X̂ (moving and speeding-up
  bodies), plus 16V₁ₖV₂ₖ. Each is a derivative, with respect to a body's position, of Fock's function
  g = ln(r₁ + r₂ + r₁₂), whose Laplacian is 1/(r₁r₂): the field of two point sources' energy spread between them.
  The speeding-up body's own extra field (Part 9's second candidate) is one of these pieces.
- **The field acting on the field** (`fourth_rule_integrals.py`): the term Ŵᵢⱼ ∂ᵢⱼV in X̂'s source, the 1PN
  space-space field (itself made of field energy) acting through the Newtonian potential's curvature. With both parts
  static it gives four products, and these have no closed form here, so they are integrated directly.

### 10.2 The closed-form pieces: 1.62 → 0.44 (O10)

Each closed-form piece is valued at each body with Hadamard's finite part, exactly as Part 9 values the rules, and
added to the rules' field before the same clock-rate pipeline runs. For body 1 at ⅓ of the mass the miss drops from
394/243 to 4/9. Over six mass ratios the new miss is exactly

  m₁m₂(5 − 10m₁ + 3m₁²)   (total mass 1),

fitted on three ratios and holding on the other three, including m₁ = 2/5 (222/625), which was never used in a fit.
So the closed-form pieces take away m₁m₂(14m₁² − 20m₁ + 21)/3 of Part 9's miss.

### 10.3 The field acting on the field (O9)

The four products are T1 = Ŵself₁ : ∂∂V₂, T2 = Ŵself₂ : ∂∂V₁, T3 = Ŵcross : ∂∂V₁, T4 = Ŵcross : ∂∂V₂. Each gives X̂
and its gradient at body 1 as a Poisson integral over all space, done in two-centre coordinates: the inner integral
exactly, the outer one numerically after taking out its infinite parts, whose finite part is added exactly. Two
details matter:

- **The shape of the cut-out.** Ordering the integral around body 1 cuts a thin slab, not a ball, out around body 2.
  For a leading term A(u)/ρ³ there, ball minus slab is a₂/3 + b/4, where A = c₀ + b u + c₂u² and a₂ = 2c₂/3. The
  rule was checked on integrands with known answers, and T4 ordered around body 2 instead (a ball there, no
  correction) gives the same numbers.
- **The point part.** ∂ᵢⱼ(1/r) has a point part, −(4π/3)δᵢⱼ δ³. At body 2 it picks up the trace of Ŵ there, taken as
  its average over directions (Hadamard), which adds (1/3)m₂Ŵₖₖ(y₂) times the kernel.

| product | masses | X̂ at body 1 | its gradient | point part (value, gradient) |
|---|---|---|---|---|
| T1 | m₁²m₂ | −1/6 | 1/10 | −1/6, 1/6 |
| T2 | m₁m₂² | 0 | 0 | 0, 0 |
| T3 | m₁²m₂ | −1/3 | −1/10 | 0, 0 |
| T4 | m₁m₂² | 1 | −2 | 1/3, −1/3 |

(totals, point part included). So X̂ at body 1 gains (αm₁²m₂ + βm₁m₂²)/r₁₂³ with **α = −1/2, β = 1**, and its gradient
(γm₁²m₂ + δm₁m₂²)/r₁₂⁴ with **γ = 0, δ = −2**. Every finite part lands on a simple fraction to better than 10⁻⁶. T3's
value was also checked by direct 3D integration with a small ball cut out (−0.3331 at a cut of 10⁻³). T1's gradient
from the slab-corrected integral (−1/15 before its point part) agrees with a separate analytic route through the
inverse Laplacian of Ŵself's gradient (a scratch check, not in the package).

### 10.4 What round 14 left: one point-like term (O10)

With α..δ added, the miss at m₁ = ⅓ is 32/243, and at every mass ratio it is

  m₁m₂²(3 − m₁)/3.

This also follows by pure algebra: a hand derivation at Newtonian order of how the four numbers move body 1's clock
rate (`fourth_rule.newtonian_map`), added to 10.2's formula, gives the same polynomial (Lean `fourth_rule_residual`),
and the map matches the pipeline at all six ratios.

The leftover has exactly the shape of one extra point-like source at each body: X̂ ⊃ m_a m_b²/(4 r₁₂² r_b), which is
β → 5/4 and δ → −9/4. With it the miss is 0 at all six mass ratios (Lean `fourth_rule_point_term`). Equivalently,
T4's point part would need Ŵcross's trace at body 2 counted 7/4 times instead of once. Round 14 guessed that the gap
came from how that trace is valued at a point. Round 15 showed the guess was wrong: T4's value holds, and the
1/4 is a missing piece of D (10.5).

### 10.5 The leftover derived: a speeding-up body's own field squared (O11, round 15)

Each body's Newtonian potential V is sent out at light speed. Expanded to the order needed, the delay adds a piece set
by the body's own acceleration a: V_b ⊃ M_b/r_b − (M_b/2) a_b·n_b, with n_b pointing from the body to the field point.
(This is the same term that dresses each body's mass in the paper's Eq. (4.8) with −½ m₁m₂/r₁₂.) Round 14 put the
X̂ part of this into D, as Xacc = −(M_b²/2) a_b·(x − y_b)/r_b². But g₀₀ also has −2V², and squaring V gives the cross
term of the body's two pieces:

  −2V² ⊃ +2 M_b² a_b·(x − y_b)/r_b² = −4 Xacc.

So in g₀₀ the speeding-up body's self part is 8 Xacc − 4 Xacc = 4 Xacc, not 8 Xacc. The rules' clock field has no
−2V² cross of this kind, so the term belongs in D, and round 14 had left it out.

At body a, with the other body pulling it at a = m_b/r₁₂², this piece has value 2 m_a m_b²/r₁₂³ and x-gradient
−2 m_a m_b²/r₁₂⁴. That is exactly 8 times the point term of 10.4, β → 5/4 and δ → −9/4, at all six ratios. So the
"quarter point term" is this product, and nothing about point values was needed.

Three checks, the first two exact at m₁ = ⅓, ½, ¼, ⅔, ⅕ and ⅖:

- **Body 1's clock rate:** with the piece added, the x³ miss is 0 against the first law of binary mechanics.
- **The orbit's size:** at orbital rate ω = 1 (total mass 1) the distance between the bodies comes out as
  r₁₂ = 1 + (−1 + ν/3)x + B x², and B matches the known harmonic-coordinate 2PN law
  ω² = (m/r³)[1 + (ν − 3)γ + (6 + 41ν/4 + ν²)γ²], γ = m/r, at all six ratios (e.g. B = 1547/1458 at ⅓, 43/36 at ½).
  Round 14's version missed it by the same point term.
- **Term by term against the paper (O12):** the paper gives general relativity's g₀₀ at body 1 explicitly, its
  Eq. (7.6a). Along the solved circular orbit our g₀₀ there (rules + field energy + own field squared) equals it
  exactly through 2PN at m₁ = ½ and ⅖ (and ⅓ in a scratch run). Without the new piece it misses by −2m₁m₂²x³,
  which is the point term.

### 10.6 What this means, on its own merits

- **What it contributes:** general relativity's 2PN field of two bodies on a circle splits exactly as "the three
  rules + the field's own energy", and each piece can be added and measured on its own: Fock's closed-form function
  does most of the work, the field acting on the field most of the rest, and a body's own field squared the last
  piece. Two separate quantities, a clock rate and an orbit size, close at six mass ratios. So the four rules, "fields
  are sent out at light speed, clock rates multiply, masses are dressed, and the field's own energy is a source",
  reproduce gravity's field through 2PN for circular orbits. That is a lead that held.
- **What is still put in:** the shape of the field-energy sources is read from general relativity's 2PN equations, not
  derived from the dial. Deriving that shape, rather than checking it, is the next step for the fourth rule. Orbits
  that aren't circles, spins, and 2.5PN (the back-push of the waves) are untested at this order.
- **Ours:** the split D = general relativity − rules, the step-by-step closing (1.62 → 0.44 → 0.13 → 0 at m₁ = ⅓),
  the formulas for each stage, the four field-on-field numbers as computed here, the slab-to-ball rule, the finding
  that round 14's leftover was one point-like term, and identifying it as a speeding-up body's own field squared.
- **Known:** the 2PN harmonic field and its point-body regularization (Blanchet, Faye and Ponsot 1998), Fock's function
  ln(r₁ + r₂ + r₁₂) (Fock 1939), Hadamard's finite part, and the idea that the field's energy gravitates and rebuilds
  general relativity (Feynman's lectures; Deser 1970).
- **Checked against the paper:** the field's formulas match Blanchet, Faye and Ponsot's Eqs. (2.4)–(2.5), the
  dressed mass matches their Eq. (4.8), and g₀₀ at a body matches their Eq. (7.6a) term by term (O12). Only the
  value at the body is compared that way; its gradient, which the orbit also uses, is checked through the orbit size
  and clock rate.

## Verdict, on its own merits

- **Holds:** full Kerr comes out of a flat lattice plus cells that each carry a velocity and a turn. The equations
  hold to 50 digits, the curvature is Kerr's own, and the turn's size, shape and speed law are all forced. The cells
  obey the same zero-energy rule as before. Brennan's "more structure than velocity" is exactly one extra shear per
  cell.
- **Opens a lead:** the turn is visible as a clock gap, 2πa per lap, equal to Kerr's clock effect. Deriving that gap
  from the dial, rather than reading it off Doran's metric, would make the turn an output. A second, smaller lead: no
  fluid can make full Kerr because fluids lack the turn.
- **Derived since (Part 3):** the turn's form a sin²θ and its size a = J/M now come from rigid shells, the circulation
  law, Sagnac's area law and frame dragging's far field. That makes the clock gap 2πJ/M an output. Kerr's horizon spin
  comes out as a bonus.
- **Derived since (Part 4):** the speed law v² = 2Mr/ρ² and the lattice's squashed spheres both come from one reading:
  the source's height turned by i, from 0 to i·a. Other potentials of mass M fail, and the shift's ring is the turn's
  ring only when a = J/M.
- **Derived since (Part 5):** the turn comes from the same shift, through its light rays. So one move gives the
  lattice, the speed and the turn, and Part 3's two assumptions become properties of the rays.
- **Charge (Part 6):** the same shift with a charge gives Kerr–Newman, with Einstein's and Maxwell's equations
  holding. The field's energy enters as |Q/D|², and the real part fails. The shifted charge is a magnet with g = 2.
- **Two bodies (Part 7):** one shear per cell can't hold two bodies (Goldberg–Sachs), and none of three ways of
  combining two shears works: the leftover grows as M₁·M₂. For bodies at rest, the log of the clock rate adds
  exactly (Weyl), by velocity addition on the axis, and the strut between them pushes with Newton's force.
- **Two bodies in orbit (Part 8):** three rules give gravity's slow-orbit field of two moving bodies exactly, and
  each is needed. Clock-rate multiplication, found for bodies at rest, carries over. The known orbits, spin wobbles,
  pulsar decay rates and spin resonances follow. Spins don't line up on average.
- **The next level (Part 9):** the three rules fail at 2PN. Body 1's clock rate on a circular orbit is right through
  1PN and for a test mass, and wrong at 2PN by ν(15/2 + 9Δ/2 − 23ν/3). Rule B drops out of the bodies' motion
  entirely, so "each is needed" holds for the field, not the orbits. The likely missing piece is gravity's own field
  energy between the bodies.
- **A fourth rule (Part 10):** adding the field's own energy closes the 2PN miss exactly at every mass ratio
  (1.62 → 0.44 → 0.13 → 0 at m₁ = ⅓). The last piece is a speeding-up body's own field squared, round 14's
  "quarter point term". The orbit's size then also matches the known 2PN law. The lead held for circular orbits; the
  sources' shape is still read from general relativity, not derived.
- **Not derived:** the shift itself (a reading, known as Newman–Janis) and the river's push-along-the-ray rule. The
  waves' back-push and the wave formula are put in from known results. 2PN is reached for circular orbits with the
  fourth rule (Part 10); 2.5PN and beyond, and the merger, are still out of reach.
- **Frames:** one area law runs through all of them, turn = d × area. Thomas, geodetic, Lense–Thirring and Kerr sit at
  d = −1; Foucault, Coriolis, Larmor and Berry at d = +1; the river and the acoustic metric at d = 0. Larmor's theorem
  explains the ½ in frame dragging's gyroscope rule.

## Files

- `checks.py`: 73 checks (K1–K7 the spinning river, F1–F3 frames, D1–D6 deriving the turn, S1–S4 deriving the
  speed law, T1–T5 the turn from the shift, C1–C5 charge, B1–B5 two bodies, O1–O6 two bodies in orbit, O7–O8 the
  next level of detail, O9–O12 the fourth rule), writes `checks_output.txt`, `results.json`.
- `next_order.py`: the 2PN test of Part 9, body 1's clock rate from the rules and from the first law, exact fractions;
  `python3 next_order.py` writes `next_order.json` (seven mass ratios, the fit, the held-out ratio, the two rule
  swaps; about 5 minutes).
- `fourth_rule.py`: Part 10, the rules plus the field's own energy at six mass ratios; writes `fourth_rule.json`
  (stages, formulas, the point-term closure, and round 15's own-field-squared stage with the orbit-size check against
  the known 2PN law). Uses `fourth_rule_cache/` (the closed-form pieces at each body); with
  `--fresh` it recomputes them (about 6 minutes per mass ratio).
- `fourth_rule_closed.py`: the closed-form field-energy pieces (Fock's function) and their finite parts at each body.
- `fourth_rule_integrals.py`: the four field-on-field integrals; writes `fourth_rule_integrals.json`; `--cross` adds
  the cross-checks (about 2 minutes, a few more with `--cross`).
- `orbits.py`: orbits, spins and the spin ensemble for Part 8; `python3 orbits.py` writes `orbits.json` (80 pairs,
  about 5 minutes on 4 cores).
- `numgr.py`: curvature at a point, exact derivatives and 50-digit arithmetic.
- `grtools.py`: symbolic helper shared with the earlier pages.
- `figures.py`: the twenty-three figures in `figures/` (Figs 19–21 need `orbits.json`, Fig 22 `next_order.json`,
  Fig 23 `fourth_rule.json`).
- `lean/SpinningRiver.lean`: `shear_sq`, `cellMap_mul`, `shear_layers`, `kerr_period`; round 7 adds `area_law`,
  `rim_at_light`, `tphi_block`, `horizon_rate`; round 8 adds `complex_distance`, `conj_product`; round 9 adds `ray_norm`; round 11 adds `velocity_addition`, `strut_push`; round 12 adds `clock_multiply`; round 13 adds `two_pn_miss`; round 14 adds `fourth_rule_residual`, `fourth_rule_point_term`; round 15 adds `own_field_squared`. Written without a Lean
  toolchain; the repo copy is the shared library's `FusedAlgebra.Topics.SpinningRiver` (namespace
  `CircleOp.SpinningRiver`, `ring` written as `grind`), which compiles clean, rounds 13 to 15 included from Nov 27. It
  builds on `FusedAlgebra.Topics.JetFields`.

## Known vs ours

| piece | status |
|---|---|
| Kerr metric (1963); Doran's river form (2000); Hamilton–Lisle river with boost + rotation (2008) | known |
| oblate spheroidal coordinates are flat space | known |
| reading the cell map as two commuting shears at the dial's flat middle | **ours** |
| forcing tests: no turn, other sizes, other shapes, other speed laws all fail | **ours** (checks) |
| cells are zero-energy, zero-angular-momentum fallers | known (Doran); **ours** as "the velocity rule, unchanged" |
| clock effect T₊ − T₋ = 4πa | known (Cohen & Mashhoon 1993); **ours** as the cell clock's gap per lap |
| Larmor's theorem; gravitational Larmor theorem | known (Larmor 1897; Mashhoon); **ours**: the ½ tied to the swap rule and Part F |
| acoustic metric = river; bathtub horizon and ergoregion | known (Unruh 1981; Visser 1998) |
| G_nn = det(sym ∇u) for any 2D acoustic river; the bathtub's Q < 0 | **ours** as the Q-layer reading (the energy constraint itself is known) |
| fluids can't make full Kerr | known (Visser & Weinfurtner 2005); **ours**: the reason is the missing turn |
| Foucault as holonomy; Earth's drag at 4.6 × 10⁻¹⁰ of it | known physics; the ratio's form 2(I/MR²)(GM/Rc²) worked out here |
| sin²θ shape as the area of a latitude circle; α as the turning shell's speed | **ours** (D1); Sagnac's law known |
| circulation law for the cell clock gap, forcing Ω(r) = κ/(r²+b²) | **ours** (D2), an assumption, modelled on Kelvin |
| κ = J/M from the far field | **ours** as a derivation step; a = J/M itself known |
| light-speed rim of a flat disk of radius a | known picture (Israel 1970); **ours** as the lattice's r = 0 shell and as the reason b = κ |
| horizon spin Ω_H = a/(r₊²+a²) = the shell rate; Carter's frame rate | known values; **ours** as the derived shell rate |
| complex shift by i·a giving Kerr | known (Newman & Janis 1965); complex-position source known (Appell 1887; Newman 1973) |
| Kerr–Schild function 2Mr/ρ² = 2 Re(M/D), harmonic | known (Kerr & Schild 1965) |
| shift as the dial's rotation number; deriving the lattice shape and the speed law from it through v²/2 = −Φ | **ours** (S1–S3) |
| rival potentials of mass M fail; the shift's ring = the turn's ring ⇔ a = J/M | **ours** (S4, S1b with D3–D4) |
| Kerr congruence as the rays of a complex-point source; Kerr–Schild form | known (Kerr & Schild 1965; Newman 1973) |
| rays as the real unit k with k·E = 1, E = grad D; their lean = the turning shell's speed; the turn read off as the ray's sideways part | **ours** as written (T1–T4); likely equivalent to known constructions |
| Kerr–Newman by the complex shift; its g = 2; Reissner–Nordström mass function M − Q²/2r | known (Newman et al. 1965; Carter 1968; standard) |
| Kerr–Newman river from the shift: field energy as \|Q/D\|² in the jet river's Gauss law, field on the turned leg, rivals failing | **ours** (C1–C5); the river form itself may be in the Doran-coordinate literature |
| vacuum Kerr–Schild metrics are algebraically special; two static bodies are not | known (Goldberg–Sachs 1962; the Bach–Weyl field is type I) |
| two-body static field, strut, push M₁M₂/(ℓ² − (M₁+M₂)²) | known (Weyl 1917; Bach & Weyl 1922; Israel & Khan 1964) |
| three ways of combining two cell shears all fail, leftover ∝ M₁M₂; the velocity-addition reading on the axis; two bodies at the boost setting | **ours** (B2, B3); the first likely known for Kerr–Schild sums |
| Ernst potential of Kerr on the axis = M/(z − ia); axis data → field (Sibgatullin); double Kerr (Kramer & Neugebauer 1980) | known; **ours** only as the proposed next test |
| two-body 1PN field (harmonic gauge); EIH motion; 6πM/p; quadrupole decay; 2.5PN back-push | known (Einstein, Infeld & Hoffmann 1938; Blanchet, Faye & Ponsot 1998; Peters & Mathews 1963; Iyer & Will 1993) |
| that field as exactly three rules: light-speed lag, clock rates multiply, dressed masses; each needed | **ours** (O1) |
| spin precession as geodetic turn + half the curl of the drag field | known (Barker & O'Connell 1975; Kidder 1995); **ours** only as assembled from our pages' pieces (O2) |
| pulsar decay and double-pulsar wobble; spin-orbit resonances; GW150914 run-up and Kerr ringdown | known (Weisberg & Huang 2016; Breton et al. 2008; Schnittman 2004; Abbott et al. 2016; Berti, Cardoso & Will 2006) |
| first law of binary mechanics; clock rate z₁ = 1/u^t; 2PN circular-orbit energy and angular momentum | known (Le Tiec, Blanchet & Whiting 2012; Detweiler 2008; Blanchet's Living Review) |
| Hadamard's finite part for point bodies; 2PN field-energy terms spread between the bodies | known (Blanchet, Faye & Ponsot 1998) |
| the rules' 2PN extension, run through the clock-rate test; the miss ν(15/2 + 9Δ/2 − 23ν/3); exact at 1PN and for a test mass | **ours** (O7, O8) |
| rule B (and the space rule) drop out of the bodies' motion through 2PN; B shapes only the field between them | **ours** (O7) |
| 2PN harmonic field of two point bodies, its field-energy sources and point-body regularization; Fock's function ln(r₁ + r₂ + r₁₂); the field's energy gravitating (Feynman–Deser) | known (Blanchet, Faye & Ponsot 1998; Fock 1939; Deser 1970) |
| D = general relativity − rules at 2PN, added piece by piece; the stage formulas m₁m₂(5 − 10m₁ + 3m₁²) and m₁m₂²(3 − m₁)/3; the four field-on-field numbers α..δ as computed here; the slab-to-ball rule | **ours** (O9, O10) |
| round 14's leftover as one point-like term m_a m_b²/(4 r₁₂² r_b) per body, closing the miss at every ratio | **ours** (O10) |
| the −2V² cross of a body's M/r with its own acceleration term −(M/2) a·n (the retarded potential) | known (part of the 2PN field; Blanchet, Faye & Ponsot 1998) |
| identifying the leftover as that cross term, so the rules + field energy close the 2PN clock rate and orbit size exactly, and g₀₀ at a body matches the paper's Eq. (7.6a) | **ours** (O11, O12) |
