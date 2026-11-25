"""Assemble index.html with the SVG plates inline. Run: python3 build_page.py"""
import os, re
D = os.path.dirname(os.path.abspath(__file__))
def plate(name):
    s = open(os.path.join(D, "figures", name + ".svg")).read()
    return s.replace('<svg xmlns="http://www.w3.org/2000/svg" ', '<svg ')
css = open(os.path.join(D, "..", "fractional-powerset-dials", "index.html")).read()
style = css[css.index("<style>")+7: css.index("</style>")]
style = re.sub(r"/\* VOICE.*?\*/", """/* VOICE a careful reading document in the infinities series, honest about known vs open
   IDEA a plate from a geometry atlas: the corner pieces of a box, then the same box at half a dimension
   PICTURE signed corner layers of the half box, the Cantor half box, the pieces of a shape, the gap at infinity
   COLOUR white ground, near-black ink, ultramarine for pieces that count +, vermilion for pieces that count -
   TYPE Newsreader for titles, Inter for text; LAYOUT verdict and first plate first, one plate per section */""", style, flags=re.S)
style += "\n.plate svg{display:block;width:100%;height:auto}\n.verdict li{margin:6px 0}\n"
P = {n: plate(n) for n in ["p1_half_box_layers","p2_two_halves_join","p3_half_box_fractal","p4_many_halves","p5_pieces_of_a_shape","p6_infinite_sizes"]}
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>A box of half a dimension</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Inter:wght@400;600&display=swap" rel="stylesheet">
<style>{style}</style>
</head>
<body>
<header class="band top"><div class="wrap">
<p class="mute" style="margin:0 0 8px">Fused Hyperoperation Algebra · fractional power sets · Nov 25, 2026</p>
<h1>A box of half a dimension</h1>
<p class="lede">The power set of n things is the corner pieces of an n-dimensional box. A box of dimension ½ does have √2 corner pieces in an exact sense, and two of them join into a stick. But that is the joining half we already had (the dial √2 − 1), and a fractal box cannot make the iterating half: the shape is not unique, its own power set needs a half-exponential, and at infinity every fractal has the same number of points.</p>
<div class="box verdict"><strong>Verdict, on its own merits: no new lead.</strong> The fractal route is consistent for joining and closed for iterating, for clear reasons.
<ol class="key">
<li><strong>Joining (P½ × P½ = P): works, known.</strong> Corner layers of a half box are Newton's binomial C(½,k): 1, ½, −⅛, 1/16, …; they add to √2, and two half boxes join into a stick exactly. A half-dimension box with sides 2 counts the same as one stick of side √2, so it is the per-element dial √2 − 1 again.</li>
<li><strong>Not canonical as a shape.</strong> Many fractals have dimension ½; each has its own wobbling corner count. Only the number 2^½ is forced.</li>
<li><strong>Iterating (P½ ∘ P½ = P): fractals can't do it.</strong> A shape's own power set is the space of all its pieces, and at each scale it counts 2^N − 1 when the shape counts N. A half step must be a half-exponential on counts: such functions exist but none is canonical and none is built from exp and log. Products, unions and stretching change size only linearly, and twice-linear never reaches 2^N (Lean: <code>no_linear_half_pow2</code>).</li>
<li><strong>The choice-free open lead is out of reach of geometry.</strong> Every Cantor set, of any dimension, is P(ℕ) as a set, and no closed set has a size strictly between ℕ and P(ℕ). The open half lives only among sets that are not geometric in this sense.</li>
</ol></div>
<div class="plate">{P['p1_half_box_layers']}</div>
<p class="cap">Plate 1. Left: the stick, square and cube from the call, each corner piece counted once. Right: a box of dimension ½. Its layers are C(½,k); from layer 2 on they alternate in sign, and the running total settles on √2.</p>
</div></header>
<nav class="nav" aria-label="Sections"><a href="#join">Two halves join</a><a href="#shape">As a real shape</a><a href="#iterate">Iterating</a><a href="#infinite">At infinity</a><a href="#tests">Consistency</a><a href="#known">Known vs new</a></nav>

<section class="band" id="join"><div class="wrap">
<h2>Two half boxes make a stick</h2>
<p>Joining two boxes multiplies them: a box of dimension D next to one of dimension e is a box of dimension D + e, and its layer k is the sum of (layer i of the first) × (layer j of the second) over i + j = k. Do that with two half boxes and every layer past the first two cancels to zero. What is left is 1 empty piece and 1 side: a stick.</p>
<div class="plate">{P['p2_two_halves_join']}</div>
<p class="cap">Plate 2. Each coloured diagonal adds to one layer of the joined box. The negative pieces of the half box are exactly what make the cancellation work.</p>
<div class="eq">box of dimension D, sides 1 + a: &nbsp;(1 + a)^D = 1 + (a ∘₁ a ∘₁ … ∘₁ a, D times) = 1 + a^(∘₁D)</div>
<p>So "fractional dimension" is the fractional circle power, and at a = 1 it is 2^D. A box of dimension ½ with sides 2 has the same count as one stick of side √2: the half box is the dial √2 − 1 from the dials page, in a different drawing. Three third-boxes make a stick the same way. The numbers are forced: any continuous rule with f(D + e) = f(D)·f(e) and f(1) = 2 is 2^D.</p>
</div></section>

<section class="band" id="shape" style="background:var(--panel)"><div class="wrap">
<h2>The half box as a real shape</h2>
<p>A fractal makes the half box visible. Keep the two outer quarters of a stick, then of each piece, forever. Each level has twice the pieces at a quarter the size, so doubling the size gives √2 times as much: dimension ½. Square it and you get the four-corner set, which has dimension 1 like a stick but is dust.</p>
<div class="plate">{P['p3_half_box_fractal']}</div>
<p class="cap">Plate 3. The half box and its square. The square has the stick's dimension and the stick's count per scale, but it is not a stick: no two pieces touch, and it casts almost no shadow (its Favard length is zero, a known fact about this set).</p>
<div class="plate">{P['p4_many_halves']}</div>
<p class="cap">Plate 4. Three different halves, all dimension ½. A whole box has an exact corner count; a half box only has a growth rate, and the leftover wobble depends on the shape. The shape does not pin down one half.</p>
</div></section>

<section class="band" id="iterate"><div class="wrap">
<h2>Iterating needs something geometry doesn't supply</h2>
<p>The power set of a shape is the space of all its closed pieces, measured by how far apart two pieces are. At a scale where the shape looks like N separate points, its pieces look like 2^N − 1 separate points (subsets of a net are a net; subsets of a spread-out set stay spread out). So geometry's own power set really is the power set, scale by scale. Counting by growth rate also removes the finite wall at n = 0 from the earlier page.</p>
<div class="plate">{P['p5_pieces_of_a_shape']}</div>
<p class="cap">Plate 5. A half step between a shape and its pieces would turn a count N into h(N) with h(h(N)) = 2^N. Two different whole-number answers are drawn; there are many.</p>
<p>That lands exactly on the half-exponential: functions with h(h(x)) = 2^x. They exist (Kneser built a smooth one for base e in 1950), but none is canonical without extra conditions, and none can be built from arithmetic, exp and log. Products, unions and stretching change the log of the count at most linearly, and doing a linear thing twice is still linear. So fractal geometry turns the iterating half into a known problem with a known answer: possible, not unique.</p>
</div></section>

<section class="band" id="infinite" style="background:var(--panel)"><div class="wrap">
<h2>At infinity, dimension can't see size</h2>
<p>A point of a Cantor set is an infinite path of left/right choices, which is a subset of ℕ. So every Cantor set, of every dimension, is P(ℕ) as a set; the Lean file proves the left/right code is one to one at each level in any base. Dimension measures how spread out the points are, not how many.</p>
<div class="plate">{P['p6_infinite_sizes']}</div>
<p class="cap">Plate 6. The open lead from the earlier page needs a size strictly between ℕ and P(ℕ), which can exist only without choice. A closed set is either countable or contains a whole Cantor set (Cantor–Bendixson, no choice needed), so no fractal can sit in the gap.</p>
</div></section>

<section class="band" id="tests"><div class="wrap">
<h2>The consistency tests</h2>
<table class="stack"><thead><tr><th>Test</th><th>Fractal box</th><th>Why</th></tr></thead><tbody>
<tr><td>Joining: P½ × P½ = P</td><td>Yes</td><td>Layers C(½,k) convolve to 1, 1, 0, 0, … (exact through layer 10 in Lean)</td></tr>
<tr><td>Iterating: P½ ∘ P½ = P</td><td>No</td><td>Needs a half-exponential on counts; linear-size operations can't</td></tr>
<tr><td>Canonical number</td><td>Yes</td><td>2^D is forced by continuity</td></tr>
<tr><td>Canonical shape</td><td>No</td><td>Many dimension-½ fractals, each with its own count wobble</td></tr>
<tr><td>Fits the circle operation</td><td>Yes</td><td>Dimension D is the circle power a^(∘₁D) = (1 + a)^D − 1</td></tr>
<tr><td>Merging and complement</td><td>Same as the dials</td><td>Union half √2 − 1, intersection half √2 + 1, as before</td></tr>
<tr><td>Reaches the choice-free lead</td><td>No</td><td>All fractals have exactly the size of P(ℕ)</td></tr>
</tbody></table>
<h3>Other ways to make "half an element"</h3>
<p><strong>Magnitude</strong> (Leinster's size of a metric space) grows at the rate set by the dimension for compact sets in space (Meckes 2015), so it adds nothing beyond dimension here. <strong>Groupoid size</strong>: one point with a two-fold symmetry counts ½, but its maps into a two-element set number 2, not √2, so it gives no half power set. <strong>Species</strong>: split a set into a part and the rest and the power set is sets × sets, so the half under that kind of joining is "one structure per set" (count 1), not √2 per element. A different way of joining gives a different half, which is the same lesson as Plate 4.</p>
</div></section>

<section class="band" id="known" style="background:var(--panel)"><div class="wrap">
<h2>Known vs new</h2>
<p>Everything here is known mathematics put in this project's picture: Newton's binomial series, self-similar dimension, the four-corner set, hyperspaces with the Hausdorff distance, Kneser's half-exponential, Cantor–Bendixson. Nothing here goes into the README list of things for mathematicians. What it does is close one route clearly: a fractional-dimension box gives the joining half and cannot give the iterating half, so the open question stays where the earlier page left it, among permutation models without choice.</p>
<h3>Files</h3>
<ul>
<li><code>checks.py</code>: 27 of 27 pass (exact fractions, box counting, Hausdorff distances, Cantor codes).</li>
<li><code>lean/FractalBoxes.lean</code>: 15 theorems, compiled in the shared library as <code>Topics/FractalBoxes</code> (library now 2,915 theorems, all build cleanly, no sorry).</li>
<li><code>figures.py</code> and <code>figures/*.svg</code>: the six plates.</li>
</ul>
<h3>Sources</h3>
<p class="mute">Kneser, J. reine angew. Math. 187 (1950), half iterate of e^x. Wikipedia, "Half-exponential function" (no exp-log half-exponential; non-analytic constructions). Nazarov–Peres–Volberg and Łaba–Zhai on the four-corner set's Favard length (arXiv 0902.0964). McClure, Real Analysis Exchange 21 (1995/96), entropy dimensions of the hyperspace of compact sets. Meckes, Potential Analysis 42 (2015), magnitude and dimension. Earlier pages: fractional power set, power set dials.</p>
</div></section>
<script>{{for(const c of document.querySelectorAll('.stack td')){{let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}}}}</script>
</body></html>
"""
open(os.path.join(D, "index.html"), "w").write(html)
print(len(html))
