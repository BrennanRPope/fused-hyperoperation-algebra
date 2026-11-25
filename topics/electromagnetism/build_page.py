"""Build index.html (single file, SVGs inlined). Run after figures.py and checks.py."""
import os, json, re

HERE = os.path.dirname(os.path.abspath(__file__))
F = lambda n: open(os.path.join(HERE, "figures", n + ".svg")).read()
res = json.load(open(os.path.join(HERE, "results.json")))
npass, nfail = res["passed"], res["failed"]
nlean = 41


def fig(name, cap):
    return f'<figure>{F(name)}<figcaption>{cap}</figcaption></figure>'


html = f'''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Electromagnetism through mult-add</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Hanken+Grotesk:wght@300;400;700&display=swap" rel="stylesheet">
<style>
/* VOICE: a careful lab notebook for one person with good intuition; pictures first, claims plain.
   IDEA: the next plate in the physics atlas: electromagnetism first, then special relativity, then gravity.
   PICTURE: the history strip (what each theory hands the next), then one picture per fit.
   COLOUR and TYPE: deep sea ground; sky = same structure, amber = renaming, pink = fails or comes from outside; Hanken Grotesk.
   LAYOUT: one column of strips up to 1100px, figures full width, tables stack on phones. */
:root{{color-scheme:dark;--bg:#0a2540;--ink:#f4f8fb;--soft:#b9cbd9;--sky:#8fd3ff;--amb:#ffc66b;--pink:#ff5c8a;--line:rgba(255,255,255,.16);--panel:rgba(255,255,255,.04)}}
:root[data-theme=light],:root[data-theme=dark]{{color-scheme:dark}}
*{{box-sizing:border-box;min-width:0}}
body{{margin:0;background:var(--bg);color:var(--ink);font:17px/1.6 "Hanken Grotesk",Helvetica,Arial,"DejaVu Sans",sans-serif;overflow-x:clip}}
.band{{padding:clamp(28px,4vw,60px) max(5vw,16px)}}.wrap{{max-width:1100px;margin:0 auto}}
.band+.band{{border-top:1px solid var(--line)}}
h1{{font-weight:300;font-size:clamp(32px,5vw,50px);line-height:1.05;text-wrap:balance;margin:0 0 18px}}
h2{{font-weight:700;font-size:26px;line-height:1.15;margin:0 0 12px;text-wrap:balance}}
h3{{font-size:18px;margin:22px 0 6px}}
p,li{{max-width:72ch}}.lede{{font-size:19px;color:var(--soft)}}
.sky{{color:var(--sky)}}.amb{{color:var(--amb)}}.pink{{color:var(--pink)}}
figure{{margin:22px 0;padding:14px;border:1px solid var(--line);border-radius:16px;background:var(--panel)}}
figure svg{{display:block;width:100%;height:auto;color:var(--ink)}}
figcaption{{font-size:15px;color:var(--soft);margin-top:10px}}
.answer{{display:grid;gap:12px;padding:0;list-style:none;margin:20px 0}}
.answer li{{padding:14px 16px;border-left:3px solid var(--sky);background:var(--panel);border-radius:0 12px 12px 0;max-width:none}}
.answer li.ren{{border-left-color:var(--amb)}}.answer li.no{{border-left-color:var(--pink)}}
table{{border-collapse:collapse;width:100%;font-size:15px;margin:12px 0}}
th,td{{text-align:left;vertical-align:top;padding:9px 10px;border-bottom:1px solid var(--line)}}
th{{color:var(--soft);font-weight:700}}
.tag{{white-space:nowrap;font-weight:700}}
.scroll{{overflow-x:auto;display:flex;gap:18px;padding:12px 0;border-bottom:1px solid var(--line)}}
.scroll a{{flex:none;color:var(--sky);text-decoration:none;font-size:15px}}
a{{color:var(--sky)}}a:focus-visible{{outline:2px solid var(--sky);outline-offset:3px}}
@media(min-width:900px){{nav.scroll{{position:sticky;top:0;background:var(--bg);z-index:2}}}}
code{{overflow-wrap:anywhere;font-size:.92em;color:var(--sky)}}
.small{{font-size:14px;color:var(--soft)}}
html{{scroll-padding-top:72px}}
@media(max-width:600px){{table{{font-size:13px}}.stack thead{{display:none}}.stack tr{{display:grid;padding-block:10px}}.stack td{{display:flow-root;border:0;padding:3px 0;padding-inline-start:8em;text-align:start}}.stack td::before{{content:attr(data-l);float:inline-start;width:7.5em;margin-inline-start:-8em;white-space:normal;color:var(--soft)}}.stack td:empty{{display:none}}}}
</style><script>new MutationObserver(()=>{{for(const c of document.querySelectorAll('.stack td:not([data-l])')){{let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}}}}).observe(document,{{childList:true,subtree:true}})</script></head><body>
<header class="band"><div class="wrap">
<h1>Electromagnetism through the mult-add lens</h1>
<p class="lede">Electromagnetism fits under the same lens as gravity, in three exact places, and it comes first in the story: Maxwell's waves are what set the dial that special relativity and then gravity run on. Nothing here is new physics; it is known physics seen through one shared lens, with the renamings marked as renamings.</p>
{fig("history", "History in order. Maxwell's equations fix the speed c = 1/√(ε₀μ₀), which is the dial of the velocity law. Special relativity composes boosts by mult-add. General relativity makes the twist rule local. Electromagnetism and gravity share the loop rule in the bottom band.")}
</div></header>
<nav class="scroll band" style="padding-block:0"><a href="#answer">The answer</a><a href="#dial">1. The dial</a><a href="#boosts">2. Boosting E and B</a><a href="#gauge">3. The gauge phase</a><a href="#loops">4. Loops and lattices</a><a href="#corner">5. F as the corner</a><a href="#one">6. One equation</a><a href="#walls">7. Where it stops</a><a href="#gravity">8. On to gravity</a><a href="#table">Scorecard</a><a href="#files">Checks and proofs</a></nav>

<section class="band" id="answer"><div class="wrap">
<h2>The answer, judged on its own merits</h2>
<ul class="answer">
<li><b class="sky">Same structure: electromagnetism sets the dial.</b> The wave equation that comes out of Maxwell's equations survives a dial boost only when q = −1/c², and c² = 1/(ε₀μ₀) (CODATA: c²ε₀μ₀ = 1.000000000000). A Galilean boost (q = 0) breaks it. This is the historical reason relativity grew out of electromagnetism.</li>
<li><b class="sky">Same structure: the field rides the boost.</b> Split E and B along the two light lines and a boost multiplies one part by the Doppler factor D and the other by 1/D. Their ratio is multiplied by D² = 1 + d(v), mult-add's lift, and two boosts multiply the D's, so the d's mult-add. A light wave sits on one light line (the wall) and its energy density falls by exactly that lift.</li>
<li><b class="sky">Same structure: loop turn = curvature × area.</b> The phase a charge picks up around a loop is charge × enclosed flux, the same rule as gravity's "a loop returns boosted by curvature × area". On a lattice each plaquette is the mult-add of its four links, and neighbouring plaquettes mult-add into the bigger loop (discrete Stokes).</li>
<li class="ren"><b class="amb">Renamings, said plainly.</b> Writing a phase as a = e^(iθ) − 1 turns multiplication into mult-add, but that works for any group of numbers under multiplication, so it adds nothing by itself. Calling F "the ab corner" of ∇A is notation; the content underneath is Clifford's split of a product into a dot part and a wedge part, and ∇F = J is Hestenes's 1966 spacetime algebra.</li>
<li><b class="sky">A small frame of our own:</b> the gauge phase and the boost are the same rule, a ✳ ā = 0, at mirror dials: q = +1 (a circle that wraps and never meets the wall) and q = −1 (a hyperbola that runs up to the light walls). Exact and tidy; a way of seeing, not a prediction.</li>
<li><b class="sky">Why electromagnetism is linear and gravity isn't:</b> phases commute, so plain mult-add (which is commutative) is exactly electromagnetism's group law, a constant potential makes no field, and the field carries no charge. Boosts and Yang–Mills phases don't commute; the gap between the two orders is the commutator, which is the field's self-coupling. Known in substance (abelian versus non-abelian gauge theory); the mult-add reading is ours.</li>
<li class="no"><b class="pink">Where it stops.</b> Mult-add does not give the field equation itself; the source law ∇F = J comes from outside, as Einstein's equations did on the gravity page. Fields in empty space add by plain addition (the k = 0 rung), not mult-add. Born–Infeld's field cap has the exact shape of the light-speed wall but does not compose by Einstein's rule.</li>
<li class="no"><b class="pink">No new prediction and no lead toward new physics here.</b> What is worth keeping is the bridge: electromagnetism is the commuting half of the same loop rule gravity uses, which is why the river picture's swirl turning out to be teleparallel torsion (a gauge field) fits the family.</li>
</ul>
<p class="small">{npass} of {npass + nfail} exact checks pass (sympy, numpy, scipy); {nlean} Lean theorems, core Lean 4.37.0-rc1, no sorry, standard axioms only.</p>
</div></section>

<section class="band" id="dial"><div class="wrap">
<h2>1. Electromagnetism sets the dial</h2>
<p>Take a boost with dial q: x′ = γ(x − vt), t′ = γ(t + q v x). Put it into the wave operator ∂²/∂t² − c²∂²/∂x² that Maxwell's equations give for light. A cross term ∂²/∂x′∂t′ appears with coefficient −2γ²v(1 + c²q). It vanishes only at <b>q = −1/c²</b>, and then the operator comes back with the same wave speed. At q = 0 (Galileo) the cross term is −2γ²v and the wave equation changes form.</p>
<p>Maxwell's constants fix c: c² ε₀ μ₀ = 1, so in SI units the dial is q = −ε₀μ₀ = −1.1127 × 10⁻¹⁷ s²/m². From there the relativity page takes over: velocities compose by (u + v)/(1 − q u v), and d(u) = 2u/(1 − u) composes by mult-add.</p>
{fig("dial", "Left: a Galilean boost shears the light lines, so light would move at c − v one way and c + v the other. Right: a Lorentz boost squeezes along the light lines by D and 1/D; the light lines stay put and areas are kept.")}
<p class="small">Known: Lorentz 1904, Einstein 1905. The dial reading is the project's landscape law (relativity page). Checks A1 to A3; Lean <code>wave_cross_zero_iff</code>, <code>galilean_cross</code>, <code>wave_speed_kept</code>.</p>
</div></section>

<section class="band" id="boosts"><div class="wrap">
<h2>2. Boosting E and B: the field on the light cone</h2>
<p>Write the boost with the Doppler factor D instead of v. Then γ = (D² + 1)/2D and γv = (D² − 1)/2D, everything is rational, and the textbook boost of the fields along x splits cleanly:</p>
<ul>
<li>E_y + B_z and E_z − B_y are multiplied by <b>1/D</b>; E_y − B_z and E_z + B_y by <b>D</b>.</li>
<li>So the ratio of the two parts is multiplied by D² = (1 + v)/(1 − v) = 1 + d(v): mult-add's lift. Their product does not change (a squeeze, the same as T² − X² for clocks and rulers).</li>
<li>Two boosts make one with D = D₁D₂. The speed of D₁D₂ is Einstein's sum, and D² − 1 composes by mult-add: (D₁D₂)² − 1 = (D₁² − 1) ✳ (D₂² − 1).</li>
<li>E² − B² and E·B are the same in every frame. A light wave has both equal to 0 (a null field): it sits on a light line, the wall. Boosted along its direction its amplitude is divided by D and its energy density by D², the inverse lift.</li>
<li>A pure electric field, boosted, gains B_z = −γv E_y: magnetism is electricity seen from a moving frame (Purcell's picture).</li>
</ul>
{fig("field", "The plane of E_y and B_z. The light lines (dashed) are where light waves live. A boost slides a field along its hyperbola E² − B² = constant; a light wave on the light line slides toward zero by 1/D per boost.")}
<p>The same boost is the Clifford sandwich F′ = R F R̃ with rotor R = cosh(φ/2) − sinh(φ/2) σₓ, where D = e^φ: the half-velocities of the Clifford page, now acting on the field.</p>
<p class="small">Known: the transformation of fields (Lorentz, Einstein 1905), light-cone components (standard in radiation theory), rotor form (Hestenes 1966). Checks B0 to B7; Lean <code>lc_plus</code>, <code>lc_minus</code>, <code>lc_pair2</code>, <code>lc_ratio</code>, <code>doppler_multAdd</code>, <code>dop_compose</code>, <code>vD_mul</code>, <code>inv1_boost</code>, <code>inv2_boost</code>, <code>null_stays_null</code>, <code>wave_boost</code>, <code>purcell</code>.</p>
</div></section>

<section class="band" id="gauge"><div class="wrap">
<h2>3. The gauge phase as a dial</h2>
<p>A charged particle's wavefunction carries a phase e^(iθ), and electromagnetism is the freedom to change that phase from place to place. Write the phase as 1 + a. Then:</p>
<ul>
<li>Phases compose by mult-add: (e^(iθ₁) − 1) ✳ (e^(iθ₂) − 1) = e^(i(θ₁+θ₂)) − 1. <span class="amb">This is a renaming of multiplication.</span></li>
<li>A number is a unit phase exactly when <b>a ✳ ā = 0</b>: its complex conjugate is its mult-add inverse, because a ✳ ā = |1 + a|² − 1. Unit phases are closed under ✳, and the wall a = −1 is never one.</li>
<li>A boost, written as a split-complex number 1 + a = cosh φ + j sinh φ with j² = +1, obeys the very same rule a ✳ ā = 0.</li>
<li>In half-angle charts the two compose by the landscape law at mirror dials: tan(θ/2) at <b>q = +1</b>, (t₁ + t₂)/(1 − t₁t₂); tanh(φ/2) at <b>q = −1</b>, (s₁ + s₂)/(1 + s₁s₂).</li>
</ul>
{fig("circle", "Left: gauge phases on the unit circle; a and its conjugate are mult-add inverses; the circle wraps around and never reaches the wall 1 + a = 0. Right: boosts on the unit hyperbola obey the same rule with j² = +1; its arms run up to the light walls.")}
<p>Gauge invariance itself checks out exactly: shifting the potential by a gradient, A → A + ∇χ and φ → φ − ∂χ/∂t, leaves E and B unchanged.</p>
<p class="small">The circle-and-hyperbola pairing is the project's landscape picture (q page, Clifford page item 02); the a ✳ ā = 0 statement for both is this page's. Checks C1 to C5; Lean <code>maC_conj</code>, <code>maC_lift</code>, <code>phase_closed</code>, <code>wall_not_phase</code>, <code>maS_conj</code>, <code>boost_closed</code>, <code>half_charts</code>.</p>
</div></section>

<section class="band" id="loops"><div class="wrap">
<h2>4. Loops and lattices: loop turn = curvature × area</h2>
<p>Carry a charge around a closed loop and its phase turns by charge × (magnetic flux through the loop). With a uniform field this is exact for any polygon: the loop integral of A equals B × area (20 random polygons, exact rationals). Around a thin solenoid the field on the path is zero, yet every loop that goes round it picks up the full flux and every loop that misses it picks up nothing (Aharonov and Bohm 1959; seen by Tonomura and others in 1986).</p>
<p>On a lattice (Wilson 1974) each link carries a phase 1 + a. A plaquette's value is the mult-add of its four links, the two backward links entering as mult-add inverses −a/(1 + a). For a potential that varies linearly it equals e^(iε²F) − 1 exactly. Two plaquettes side by side mult-add into the big loop because the shared link runs both ways and cancels: Stokes's theorem in mult-add form.</p>
{fig("holonomy", "Left: two plaquettes, each the mult-add of four links; the shared link cancels, so the big loop is a₁ ✳ a₂. Middle: the Aharonov–Bohm phase. Right: the gravity page's loop in the t–r plane returns boosted by curvature × area.")}
<p>This is the place where electromagnetism and gravity share the most: both are connections, and in both the turn around a small loop is the curvature times the area. The difference is what turns. In electromagnetism it is a phase, and phases commute.</p>
{fig("commute", "Left: with commuting phases, a loop of constant potentials closes, so a constant potential makes no field. Right: with non-commuting ones (boosts, or Yang–Mills), the same loop leaves a gap equal to the commutator, and the field carries its own charge.")}
<p class="small">Known: Weyl 1929 (gauge phase), Yang–Mills 1954, Aharonov–Bohm 1959, Wilson 1974, Wu–Yang 1975 (loop phases as the full content of a gauge field). The mult-add form of plaquettes and Stokes is this page's rewriting. Checks C6 to C8; Lean <code>plaquette</code>, <code>stokes</code>, <code>stokes_multAdd</code>.</p>
</div></section>

<section class="band" id="corner"><div class="wrap">
<h2>5. F as the corner of a Clifford product</h2>
<p>The Clifford thread found that with multivector inputs mult-add stays associative and the gap between the two orders is exactly the wedge: x ✳ y − y ✳ x = 2 x∧y. In spacetime algebra a vector product splits the same way, ab = a·b + a∧b, so the corner of the lift square (1 + a)(1 + b) is a dot part plus a wedge part.</p>
<p>Put the derivative ∇ in for a and the potential A in for b. The corner ∇A splits into ∇·A, which is a gauge choice (the Lorenz gauge sets it to 0), and ∇∧A, which is <b>exactly the field F</b>, with E = −∇φ − ∂A/∂t and B = ∇×A. The field strength is the wedge, the order-gap part of the corner.</p>
{fig("corner", "The lift square's corner splits into a·b (the same in both orders) and a∧b (flips with order). For the derivative and the potential, the wedge part is the electromagnetic field.")}
<p class="amb">Honest status: this is notation that lines up. ∇ is an operator, not a number, so the "lift" (1 + ∇)(1 + A) has no meaning of its own; the content is the dot-plus-wedge split, which is Clifford's.</p>
<p class="small">Checks D1, D2; Clifford page item 08 (Lean <code>G4.vector_comm_wedge</code>).</p>
</div></section>

<section class="band" id="one"><div class="wrap">
<h2>6. Maxwell as one equation</h2>
<p>In spacetime algebra F = E + IB is one object, and Maxwell's four equations are one: <b>∇F = J</b>, with J = ρ plus the current. Its vector part is Gauss's law and the Ampère–Maxwell law; its trivector part is "no magnetic charge" and Faraday's law. Checked symbolically with general E(t, x, y, z) and B(t, x, y, z).</p>
{fig("one", "One equation, two parts, four laws. Below: conservation comes free in both theories.")}
<ul>
<li><b>Charge conservation is automatic.</b> ∇∇ is the wave operator, a scalar, so ∇(∇F) is again a field-type object with no scalar part; hence ∇·J = 0 for every F. In momentum space: k_μ k_ν F^μν = 0 for any antisymmetric F (Lean <code>charge_conserved</code>).</li>
<li>The same geometry gives gravity's ∇·G = 0 (the Bianchi identity, checked on every river on the GR page), which forces conservation of energy.</li>
<li>FF = (E² − B²) + 2(E·B) I: the two invariants are the scalar and the pseudoscalar parts of the field's square, and light is the case FF = 0.</li>
<li>In the Lorenz gauge, ∇∇A = J is the wave equation, and ∇∧(∇∧A) = 0 always (no magnetic charge is an identity, not a law).</li>
</ul>
<p class="small">Known: Hestenes 1966 (Space-Time Algebra); Doran and Lasenby 2003. Checks D3 to D8, F1.</p>
</div></section>

<section class="band" id="walls"><div class="wrap">
<h2>7. Where the lens stops</h2>
<ul>
<li><b class="pink">The source law comes from outside.</b> Mult-add explains how fields transform and how phases compose, not why the field obeys ∇F = J (or, as a Lagrangian, why the action is F²). That is the same gap the gravity page found for Einstein's equations.</li>
<li><b class="pink">Superposition is plain addition.</b> Maxwell in empty space is linear: two solutions add by ordinary addition, the k = 0 rung, with no corner term. Mult-add enters through the phases and the boosts, not through how fields pile up.</li>
<li><b class="pink">Born–Infeld shares the wall's shape, not the law.</b> Born and Infeld (1934) capped the field the way relativity caps speed: E = D/√(1 + D²/b²), the same curve as v = p/√(1 + p²). But Born–Infeld fields add by adding D (Gauss's law), and E(D₁ + D₂) ≠ E₁ ⊕ E₂ (0.9029 against 0.9479 in the example). ATLAS light-by-light data rule out a Born–Infeld scale below about 100 GeV (Ellis, Mavromatos and You 2017).</li>
<li>Real electrodynamics turns nonlinear only near the Schwinger field m²c³/(eħ) = 1.32 × 10¹⁸ V/m, through quantum effects, not through a mult-add rule.</li>
</ul>
{fig("walls", "One curve, two readings: speed against momentum and Born–Infeld's field against D. The wall is shared; the rule for adding is not.")}
<p class="small">Checks E1 to E3; Lean <code>bi_cap</code>, <code>bi_den_pos</code>, <code>bi_inverse</code>. Source for the bound: <a href="https://arxiv.org/abs/1703.08450">arXiv:1703.08450</a>.</p>
</div></section>

<section class="band" id="gravity"><div class="wrap">
<h2>8. On to gravity</h2>
<p>Einstein built on electromagnetism, and the lens shows the steps. Electromagnetism hands special relativity its dial, and special relativity hands gravity the twist rule. Gravity then keeps electromagnetism's shape and changes one thing:</p>
<table class="stack"><thead><tr><th>Piece</th><th>Electromagnetism</th><th>Gravity (GR and bridges pages)</th></tr></thead><tbody>
<tr><td>What turns around a loop</td><td>a phase (U(1), commutes)</td><td>a boost or rotation (Lorentz group, doesn't commute)</td></tr>
<tr><td>Loop rule</td><td>phase = charge × flux</td><td>rapidity = curvature × area</td></tr>
<tr><td>Group law</td><td>plain mult-add (commutative)</td><td>the sandwich a + b + a·k·b; order gap = twist</td></tr>
<tr><td>Field from the potential</td><td>F = ∇∧A</td><td>curvature = dω + ω∧ω (extra term from the order gap)</td></tr>
<tr><td>Conservation for free</td><td>∇·J = 0</td><td>∇·G = 0 (Bianchi)</td></tr>
<tr><td>Linear?</td><td>yes: the field carries no charge</td><td>no: gravity gravitates</td></tr>
<tr><td>Gauge form</td><td>Maxwell, ∇F = J</td><td>teleparallel: the river's swirl is the torsion; the momentum rule is a gravitomagnetic Ampère law</td></tr>
</tbody></table>
<p>So the river-equals-teleparallel finding is the natural next plate: gravity written as a gauge field in the same family as electromagnetism, with the one difference that its "phases" do not commute.</p>
<p class="small">Known: Utiyama 1956 and Kibble 1961 (gravity as a gauge theory). Cited from the <a href="https://claude.ai/artifact/SHrK88RpTBLLfTaeAxgi37">general relativity page</a>, the <a href="https://claude.ai/artifact/LhAbnuPMJnNoq9cRYBz3xc">GR bridges page</a> and the <a href="https://claude.ai/artifact/NU8rLX3jBBeaojcq8g3hMW">Clifford page</a>.</p>
</div></section>

<section class="band" id="table"><div class="wrap">
<h2>Scorecard</h2>
<table class="stack"><thead><tr><th>Piece of electromagnetism</th><th>Mult-add reading</th><th>Status</th><th>Known or ours</th></tr></thead><tbody>
<tr><td>Speed of light from ε₀, μ₀</td><td>sets the dial q = −1/c²; only that q keeps the wave equation</td><td class="tag sky">same structure</td><td>known; dial reading from the relativity page</td></tr>
<tr><td>Velocity and rapidity composition</td><td>d = D² − 1 composes by mult-add</td><td class="tag sky">same structure</td><td>known; relativity page</td></tr>
<tr><td>Boosts of E and B</td><td>light-cone parts × D and × 1/D; ratio × lift; D's multiply</td><td class="tag sky">same structure</td><td>known; mult-add phrasing ours</td></tr>
<tr><td>Light as a null field</td><td>on the wall; energy density × 1/lift</td><td class="tag sky">same structure</td><td>known; phrasing ours</td></tr>
<tr><td>U(1) phase composition</td><td>(e^(iθ) − 1) mult-adds</td><td class="tag amb">renaming</td><td>true for any multiplicative group</td></tr>
<tr><td>Phase and boost together</td><td>a ✳ ā = 0 at q = +1 and q = −1</td><td class="tag sky">exact framing</td><td>ours (small)</td></tr>
<tr><td>Aharonov–Bohm, holonomy</td><td>loop turn = curvature × area, as in gravity</td><td class="tag sky">same structure</td><td>known (Wu–Yang)</td></tr>
<tr><td>Lattice plaquettes</td><td>mult-add of four links; neighbours mult-add (Stokes)</td><td class="tag sky">same structure</td><td>Wilson; mult-add form ours</td></tr>
<tr><td>F = ∇∧A</td><td>the wedge (order-gap) part of the corner ∇A</td><td class="tag amb">notation</td><td>content is Clifford's split</td></tr>
<tr><td>Maxwell as ∇F = J</td><td>one equation, grade by grade</td><td class="tag amb">known rewrite</td><td>Hestenes 1966</td></tr>
<tr><td>Charge conservation</td><td>automatic from the geometry, like the Bianchi identity</td><td class="tag sky">same structure</td><td>known</td></tr>
<tr><td>Linearity of Maxwell</td><td>U(1) commutes, so no self-coupling</td><td class="tag sky">same structure</td><td>known; mult-add reading ours</td></tr>
<tr><td>The source law ∇F = J itself</td><td>not derived</td><td class="tag pink">from outside</td><td>as with Einstein's equations</td></tr>
<tr><td>Field superposition</td><td>plain addition (k = 0), not mult-add</td><td class="tag pink">doesn't fit</td><td></td></tr>
<tr><td>Born–Infeld field cap</td><td>wall shape matches; composition doesn't</td><td class="tag pink">partial</td><td>known (1934); bound from ATLAS data</td></tr>
</tbody></table>
</div></section>

<section class="band" id="files"><div class="wrap">
<h2>Checks and proofs</h2>
<p>Exact checks: <code>checks.py</code> ({npass} of {npass + nfail} pass), in six groups: A the dial, B boosts of the fields, C the gauge phase, loops and lattices, D Clifford form and Maxwell as one equation, E where it stops, F the bridge to gravity. Outputs in <code>checks_output.txt</code> and <code>results.json</code>.</p>
<p>Lean: <code>lean/Electromagnetism.lean</code>, {nlean} theorems, core Lean 4.37.0-rc1 (no mathlib), no sorry, only the standard axioms. It imports the shared library's Base, Thermodynamics, Relativity and TwistedBoosts. What Lean does not check: derivatives, the Clifford products and the numerics; those are in the checks.</p>
<p class="small">Files in math-theory/electromagnetism in the project folder. Pictures: figures.py → figures/*.svg.</p>
</div></section>
</body></html>
'''
open(os.path.join(HERE, "index.html"), "w").write(html)
print("index.html", len(html))
