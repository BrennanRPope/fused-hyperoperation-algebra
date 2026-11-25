"""Assemble index.html with the SVG plates inline. Run: python3 build_page.py"""
import os, re
D = os.path.dirname(os.path.abspath(__file__))
def plate(name):
    s = open(os.path.join(D, "figures", name + ".svg")).read()
    return s.replace('<svg xmlns="http://www.w3.org/2000/svg" ', '<svg ')
css = open(os.path.join(D, "..", "fractional-powerset-dials", "index.html")).read()
style = css[css.index("<style>")+7: css.index("</style>")]
style = re.sub(r"/\* VOICE.*?\*/", """/* VOICE a careful reading document in the infinities series, honest about known vs open
   IDEA a plate from a geometry atlas: five fractals of one dimension, each with its hidden dials read off
   PICTURE five half boxes, wobble curves, weighted bars, dust against gasket, the binary tree, the dial map
   COLOUR white ground, near-black ink, ultramarine for what the circle operation has, amber partial, vermilion missing
   TYPE Newsreader for titles, Inter for text; LAYOUT verdict and first plate first, one plate per section */""", style, flags=re.S)
style += "\n.plate svg{display:block;width:100%;height:auto}\n.verdict li{margin:6px 0}\n"
names = ["p1_five_half_boxes","p2_wobble","p3_weights_t","p4_topology","p5_joining_tree","p6_dial_map","p7_growth_dial"]
P = {n: plate(n) for n in names}
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Hidden dials of a half box</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Inter:wght@400;600&display=swap" rel="stylesheet">
<style>{style}</style>
</head>
<body>
<header class="band top"><div class="wrap">
<p class="mute" style="margin:0 0 8px">Fused Hyperoperation Algebra · fractional power sets · Nov 25, 2026</p>
<h1>Hidden dials of a half box</h1>
<p class="lede">Many fractals share dimension ½. What else tells them apart? Six known dials do: zoom, placement, wobble, weights, topology and branching growth. The circle operation already has the two that matter for counting (dimension and weights), has a cousin of a third (wobble), and is missing one clean one: a zoom dial for size. No mix of shape dials gives the natural iterating half of the power set; that still needs a dial on how fast branching grows.</p>
<div class="box verdict"><strong>Verdict, on its own merits: a small, clean addition; no new lead on the iterating half.</strong>
<ol class="key">
<li><strong>All of these dials are known math.</strong> Zoom ratios and Moran's equation (1946), lacunarity (Mandelbrot; Allain and Cloitre 1991), Minkowski content and complex dimensions (Lapidus and van Frankenhuijsen), multifractal spectra (Halsey and others 1986), and topology (Brouwer 1910). None is new here.</li>
<li><strong>The circle operation already has two.</strong> Dimension is the fractional circle power (1 + a)^D. Unequal weights p : 1 − p are the per-element dial t = p/(1 − p), and the multifractal spectrum is what t does to the counting.</li>
<li><strong>One is a cousin.</strong> Lattice fractals wobble forever as you zoom, because their dimension has imaginary partners ½ + i·2πn/ln(1/r). That is a phase turning with zoom, like the circle operation's twist or a Clifford rotor. It is an analogy, not the same thing.</li>
<li><strong>One is missing: zoom.</strong> The circle operation multiplies sides 1 + a and has no notion of size. A zoom dial would make dimension an output, D = log(pieces) / log(1/zoom), instead of an input. This is the one dial worth adding, and it is modest.</li>
<li><strong>The rest are geometry.</strong> Placement (lacunarity, Minkowski content) moves no counts. Topology is silent below dimension 1, since every dust is the same Cantor set.</li>
<li><strong>Which combination pins one half box?</strong> Binary branching, equal zooms and fair weights pin exactly one: two pieces at ¼. Joined with itself it is the binary full box as a tree, and gluing ends gives the stick. That pins the joining half (already the dial √2 − 1), not an iterating half.</li>
<li><strong>Why no shape dial iterates.</strong> Every shape dial keeps level counts of the form c·bⁿ, and joining keeps that form. The power set's counts 2^(2ⁿ) fit no such form (Lean: <code>powerset_not_selfsimilar</code>). A half step needs a dial on branching growth, which is the circle operation's rung flow at t = ½; it exists but is unique only under extra conditions.</li>
</ol></div>
<div class="plate">{P['p1_five_half_boxes']}</div>
<p class="cap">Plate 1. Five sets of dimension ½. Rows A to E differ in number of pieces, zoom, placement and evenness, and the dials on the right read those differences off. "Content" is the blurred size rescaled by ε^½, averaged over many blur sizes; lacunarity is the gliding-box value at three box sizes.</p>
</div></header>
<nav class="nav" aria-label="Sections"><a href="#table">The dials</a><a href="#wobble">Wobble</a><a href="#weights">Weights</a><a href="#topology">Topology</a><a href="#pin">Pinning one</a><a href="#map">circle operation's map</a><a href="#iterate">Iterating</a><a href="#known">Known vs new</a></nav>

<section class="band" id="table"><div class="wrap">
<h2>Each dial, and what it does to the corner count</h2>
<p>A box of dimension D has 2^D corner pieces in the counting sense: each level of the fractal multiplies the count by a fixed factor. The question for each hidden dial is whether it changes or pins down how those pieces are counted.</p>
<table class="stack"><thead><tr><th>Dial</th><th>Tells apart</th><th>Fixes the count?</th><th>In the circle operation</th></tr></thead><tbody>
<tr><td>Dimension D</td><td>nothing at fixed D</td><td>Only the number 2^D</td><td>Has it: (1 + a)^D</td></tr>
<tr><td>Zoom ratios r</td><td>A, B, D, E (2, 3, 2, 4 pieces)</td><td>Fixes piece sizes; with piece count, gives D</td><td>Missing: one candidate</td></tr>
<tr><td>Placement (lacunarity, content)</td><td>B from C (same pieces, middle moved)</td><td>No: counts unchanged</td><td>Outside the algebra</td></tr>
<tr><td>Wobble (complex dimensions)</td><td>lattice A, B, C, E from non-lattice D</td><td>No: a periodic correction on top</td><td>Partial: like twist or a rotor phase</td></tr>
<tr><td>Weights p : 1 − p</td><td>the same set carrying different mass</td><td>Weights the corners; fair weights give plain 2^D</td><td>Has it: t = p/(1 − p)</td></tr>
<tr><td>Topology (touching, gluing)</td><td>nothing below D = 1; gasket from dust above</td><td>No</td><td>Outside; only matters at D ≥ 1</td></tr>
<tr><td>Branching growth</td><td>shapes from the power set</td><td>Yes, for iterating</td><td>Has it, separately: rung flow at t = ½</td></tr>
</tbody></table>
</div></section>

<section class="band" id="wobble" style="background:var(--panel)"><div class="wrap">
<h2>Wobble: the dial that turns as you zoom</h2>
<p>Blur a fractal by ε, measure the blurred set, and divide by ε^½. For a non-lattice set like D, whose two zooms 9/100 and 49/100 share no common power, this settles to a number (about 3.00). For lattice sets it never settles: it repeats every time ε shrinks by the zoom.</p>
<div class="plate">{P['p2_wobble']}</div>
<p class="cap">Plate 2. The wobble is the trace of complex dimensions ½ + i·2πn/ln(1/r). They are evenly spaced, 4.53 apart for A, 2.86 for B and C, 2.27 for E. Lean proves D is non-lattice by showing 9^p·100^q is never 49^q·100^p.</p>
<p>In circle-operation terms this is a phase on the zoom, which is why it sits next to the twist dial and the Clifford rotors in the map below. It is a resemblance worth noting, not a proved identity.</p>
</div></section>

<section class="band" id="weights"><div class="wrap">
<h2>Weights are the per-element dial t</h2>
<p>Keep the set A fixed and share its mass unevenly between the two pieces at every level. The set, and so its dimension, does not move. What moves is how much each corner piece counts. That is exactly the per-element dial t: a piece's weight is t/(1 + t).</p>
<div class="plate">{P['p3_weights_t']}</div>
<p class="cap">Plate 3. At t = 1 every level-5 piece carries the same mass and all Rényi dimensions are ½. As t moves off 1 the spectrum of local dimensions α spreads. The collision sum p² + (1 − p)² is smallest, ½, at fair weights (Lean: <code>weights_odds</code>).</p>
</div></section>

<section class="band" id="topology" style="background:var(--panel)"><div class="wrap">
<h2>Topology is silent below dimension 1</h2>
<div class="plate">{P['p4_topology']}</div>
<p class="cap">Plate 4. Left: B and C look different, but every compact dust with no isolated points is the same Cantor set up to stretching. Right: at dimension log 3 / log 2 the gasket's pieces share corners at every level, while the 9-at-¼ dust's pieces never touch. Lean counts the gasket's (3^(m+1) + 3)/2 corners, fewer than its 3·3ᵐ triangles have between them.</p>
<p>So topology is a real dial only once a set is big enough for pieces to touch, which needs dimension at least 1. For the half box it adds nothing.</p>
</div></section>

<section class="band" id="pin"><div class="wrap">
<h2>Pinning one definite half box</h2>
<p>Ask for the plainest half box: two pieces per level (binary), both shrunk the same, mass split fairly. Moran's equation then forces the zoom to ¼, and stretching it to fill [0, 1] forces the pieces to the two ends. That is set A, and it is the only one.</p>
<div class="plate">{P['p5_joining_tree']}</div>
<p class="cap">Plate 5. A × A is the four-corner set: four pieces at ¼ each level, which is two binary choices per level, so it is the binary tree (distances agree within a factor of about 2 either way at every depth, checked to depth 4). Gluing the two codes of each dyadic point, 0.0111… = 0.1000…, turns the tree into the stick.</p>
<p>This is the joining half, P½ × P½ = P, which the dial √2 − 1 already gave. It answers "which half box", not "what is half of taking a power set".</p>
</div></section>

<section class="band" id="map" style="background:var(--panel)"><div class="wrap">
<h2>Which dials the circle operation already has</h2>
<div class="plate">{P['p6_dial_map']}</div>
<p class="cap">Plate 6. Blue: covered. Amber: a resemblance. Red: the one gap. Grey: geometry that does not change counts. The circle operation's numeric dials from the dial survey are k and the per-element t; dimension is the fractional circle power; the rung flow handles growth between rungs.</p>
<p>Adding a zoom dial would be small: give each element a size r alongside its value a, so that a fractal's dimension comes out of its piece count and zoom. Whether that earns a place depends on whether other parts of the theory want sizes too; on this page it is only what lets circle-operation talk about fractals directly.</p>
</div></section>

<section class="band" id="iterate"><div class="wrap">
<h2>Why no combination iterates</h2>
<div class="plate">{P['p7_growth_dial']}</div>
<p class="cap">Plate 7. On a log scale, self-similar shapes are straight lines (counts c·bⁿ), and joining two lines gives a line (Lean: <code>pieces_join_closed</code>). The power set's counts 2^(2ⁿ) curve away; each level squares the last (Lean: <code>powerset_branch</code>). The dashed curve is one exact half step, h(h(n)) = 2ⁿ, built to be linear on [0, ½] and forced from there.</p>
<p>Any half step must bend between the lines and the power set, so it needs a dial on how fast branching grows. That is the half-exponential, which in the circle operation is the rung flow at t = ½. Kneser (1950) built a smooth one; it becomes unique only once extra conditions are added, as Paulsen and Cowgill did for the related tetration. Shape dials cannot supply that choice.</p>
</div></section>

<section class="band" id="known" style="background:var(--panel)"><div class="wrap">
<h2>Known vs new</h2>
<p>Every dial on this page is known fractal geometry. What is new is only the placement in this project's picture: which of them the circle operation already has, the one it lacks, and the reason none of them can give an iterating half. Nothing here goes into the README list of things for mathematicians. The lead that remains is not in the fractal shapes; it is the growth dial, which the project already has as the rung flow.</p>
<h3>Files</h3>
<ul>
<li><code>fractals.py</code>: the five half-dimension sets and the measuring tools.</li>
<li><code>checks.py</code>: 25 of 25 pass (dimension certificates, lattice tests, Minkowski wobble, lacunarity, multifractal spectra, tree distances, topology, closure under joining).</li>
<li><code>lean/FractalDials.lean</code>: 14 theorems, compiled in the shared library as <code>Topics/FractalDials</code> (the whole shared library builds cleanly, no sorry).</li>
<li><code>figures.py</code> and <code>figures/*.svg</code>: the seven plates.</li>
</ul>
<h3>Sources</h3>
<p class="mute">Moran, Proc. Cambridge Phil. Soc. 42 (1946), similarity dimension. Lapidus and van Frankenhuijsen, Fractal Geometry, Complex Dimensions and Zeta Functions (Springer, 2nd ed. 2013): lattice vs non-lattice, Minkowski measurability. Allain and Cloitre, Phys. Rev. A 44 (1991), gliding-box lacunarity. Halsey, Jensen, Kadanoff, Procaccia and Shraiman, Phys. Rev. A 33 (1986), f(α). Falconer, Fractal Geometry (Wiley). Brouwer (1910), uniqueness of the Cantor set. Kneser, J. reine angew. Math. 187 (1950). Paulsen and Cowgill, Adv. Comput. Math. 43 (2017). Earlier pages: fractal boxes, power set dials, dial survey.</p>
</div></section>
<script>{{for(const c of document.querySelectorAll('.stack td')){{let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}}}}</script>
</body></html>
"""
open(os.path.join(D, "index.html"), "w").write(html)
print(len(html))
