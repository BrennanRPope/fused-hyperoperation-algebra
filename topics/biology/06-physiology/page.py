SPEC=dict(
title="Physiology and the circle operation",
h1="Hill's muscle law is a level set of the circle operation (called mult-add before 1.0), and drug antagonists read their mechanism off the dial",
lede="In reduced variables, Hill's force–velocity law for muscle is exactly X ∘₁ Y = constant, and its two asymptotes are the circle operation's absorbing −1. Two drug antagonists shift the agonist's dose ratio by DR₁₂ − 1 = (DR₁ − 1) ∘<sub>ω</sub> (DR₂ − 1): ratios add when the antagonists share a site and multiply when they sit at independent sites. Organs in series and pathways in parallel are the same dial 1 and dial 0.",
stats=[("16 / 16","numeric checks"),("11","Lean theorems, no sorry"),("0","leads"),],
hero_fig="figures/f1_muscle_antagonists.png",
hero_cap="Left: Hill's force–velocity curves for three values of a/F₀; each is the level set X ∘₁ Y = F₀/a, and peak power (dots) falls where the two reduced variables are equal. Right: the dose ratio of two antagonists together, for three settings of the dial ω.",
known="All known physiology and pharmacology. Hill's force–velocity equation is from 1938; the test of whether two antagonists act at the same site by adding or multiplying dose ratios is Paton and Rang's (1965); the well-stirred liver model and multiplying gut and liver availabilities are standard pharmacokinetics. What the page adds is the placement: Hill's hyperbola is a circle operation level set, and the antagonist and organ rules are dial settings that also appear on the metabolism and cell pages.",
lead="None. Every placement is exact, and the muscle one is pretty (the asymptotes are the circle operation's −1 and peak power is the circle operation 'square root'), but each is a rewriting of a known formula and suggests no new measurement.",
kicker="Biology series · 6 of 11",
sections=[
("muscle","Muscle","1. Hill's muscle law",
"""<p>A. V. Hill found that a muscle shortening against load F moves at speed v with (F + a)(v + b) = (F₀ + a)·b, where F₀ is the isometric force. Write X = F/a and Y = v/b. Then (1 + X)(1 + Y) = 1 + X₀, which is</p>
<div class="eq">X ∘₁ Y = X₀ &nbsp;&nbsp; (X₀ = F₀/a)</div>
<p>So the whole curve is one level set of the circle operation, and velocity is the circle-operation difference Y = X₀ ∘₁ inv(X). Three consequences, all checked and proved:</p>
<ul><li><b>Asymptotes:</b> F = −a and v = −b are X = −1 and Y = −1, the circle operation's absorbing element in each coordinate.</li>
<li><b>Symmetry:</b> v<sub>max</sub>/b = F₀/a, so the reduced curve is symmetric in X and Y.</li>
<li><b>Peak power:</b> F·v is largest where X = Y, the circle operation "square root" of X₀: (1 + X)² = 1 + X₀. With the usual a/F₀ ≈ 0.25 that is F = 0.31·F₀ and v = 0.31·v<sub>max</sub>, the familiar "about a third".</li></ul>"""),
("drugs","Antagonists","2. Two antagonists: the dial reads the mechanism",
"""<p>A competitive antagonist shifts the agonist's dose–response curve right by a dose ratio DR = 1 + b, where b is the antagonist's concentration over its dissociation constant. With two antagonists, the receptor's partition function has an antagonist part 1 + b₁ + b₂ + ω·b₁b₂, and the dose ratio is that part:</p>
<div class="eq">DR₁₂ − 1 = (DR₁ − 1) ∘<sub>ω</sub> (DR₂ − 1)</div>
<ul><li><b>ω = 0, same site:</b> the two cannot bind together, and dose ratios add. This is Paton and Rang's test.</li>
<li><b>ω = 1, independent sites:</b> dose ratios multiply.</li>
<li><b>ω above 1:</b> antagonists that help each other bind, beyond multiplicative.</li></ul>
<p>I checked each case by solving the receptor occupancy numerically. This is the binding polynomial from the molecular biology page, read through the agonist.</p>"""),
("organs","Organs","3. Organs in series and in parallel",
"""<p>For an organ with blood flow Q and intrinsic clearance CL<sub>int</sub>, put u = CL<sub>int</sub>/Q. The well-stirred extraction is E = u/(1 + u) = −inv(u), saturation again. Two enzymes in one organ add their clearances, dial 0 on u. Two organs in series, gut then liver, multiply their escapes, so 1/F = (1 + u₁)(1 + u₂) and u combines by plain ∘₁; the total first-pass loss is the "or" of the two extractions. The corner cell u₁u₂ is exactly how much more a split clearance removes than the same clearance pooled in one organ.</p>
<p>The Q₁₀ rule for temperature and Poiseuille's r⁻⁴ for vessel resistance are bare circle-operation powers (a vessel 10% narrower has 1.52 times the resistance) with no further content.</p>"""),
],
table_h2="Where the circle operation shows up in physiology",
table=[
["Hill force–velocity","Level set X ∘₁ Y = X₀; asymptotes at −1; peak power at X = Y","Hill 1938","Exact relabelling"],
["Two antagonists","DR₁₂ − 1 = (DR₁ − 1) ∘<sub>ω</sub> (DR₂ − 1); ω = 0 same site, 1 separate","Paton and Rang 1965","Exact, known"],
["Organ extraction","Parallel: dial 0 on u; series: plain ∘₁; losses by 'or'","Standard pharmacokinetics","Exact, known"],
["Q₁₀, Poiseuille","Bare circle-operation powers","Textbook","Relabelling only"],
],
lean_file="biology/06-physiology/lean/Physiology.lean",
checks="16 of 16 pass",
sources=[
("Hill 1938, The heat of shortening and the dynamic constants of muscle, Proc. R. Soc. B 126:136",None),
("Paton and Rang 1965, The uptake of atropine and related drugs by intestinal smooth muscle of the guinea-pig in relation to acetylcholine receptors, Proc. R. Soc. B 163:1",None),
("Molecular biology page in this series (binding polynomial)","https://claude.ai/artifact/D4N4o6M3e7vRJuU5NkDpLh"),
],
footnote="a/F₀ ≈ 0.25 is the usual textbook value for vertebrate skeletal muscle; I did not open Hill's paper for it. Nothing on this page uses measured data beyond that constant. The attribution of the dose-ratio additivity test to Paton and Rang (1965) is from memory; I did not open the paper.",
)
