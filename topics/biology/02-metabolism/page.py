SPEC=dict(
title="Metabolism and mult-add",
h1="Two-substrate enzymes sit on the dial, and their plotted lines cross at its wall",
lede="Write each substrate as u = K/S and the textbook rate law for an ordered two-substrate enzyme is exactly 1 + (u_A ✳_k u_B), with dial k = K_ia/K_a. Ping-pong enzymes are dial 0 and independent random binding is dial 1. The point where every double-reciprocal line crosses, which biochemists use to read off K_ia, is the wall of the dial. Two-nutrient growth laws sit on the same dial; Liebig's minimum does not.",
stats=[("19 / 19","numeric checks"),("19","Lean theorems, no sorry"),("0","new leads; one useful reading")],
hero_fig="figures/f1_lines_meet_at_wall.png",
hero_cap="Reciprocal rate against reciprocal first substrate, one line per level of the second substrate. At dial 0 the lines are parallel; at any other dial they all pass through u_A = −1/k, the wall. The middle panel uses measured constants for rat liver mitochondrial malate dehydrogenase (K_a 72, K_ia 12, K_b 86 µM), which put it at k = 1/6.",
known="Everything here is known enzymology and microbiology. Cleland's rate laws for ordered and ping-pong mechanisms date from 1963, and the crossing point of double-reciprocal lines at 1/A = −1/K_ia is how K_ia is measured. The multiplicative and additive two-nutrient forms and Liebig's law of the minimum are standard growth models. What this page adds is the placement: the mechanism type is a dial setting, and the crossing point is the dial's wall.",
lead="No new lead. The reading is clean and exact: one dial number, k = K_ia/K_a, places every two-substrate mechanism between ping-pong (0) and independent random binding (1), and the wall theorem from the library is the crossing-point rule enzymologists already use. That makes it good teaching glue for the dial and a nice confirmation that the wall shows up in a real lab plot. It does not predict anything a biochemist does not already compute. The nutrient forms are the same algebra one level up, and Liebig's minimum, the law ecologists often prefer, is provably off the dial, which is worth knowing but not a lead.",
kicker="Biology series · 2 of 11",
sections=[
("enzymes","Enzymes","1. Two substrates: the mechanism is a dial setting",
"""<p>An enzyme that joins two substrates A and B has the textbook (Cleland) rate law</p>
<div class="eq">V / v = 1 + K<sub>a</sub>/A + K<sub>b</sub>/B + K<sub>ia</sub>·K<sub>b</sub> / (A·B)</div>
<p>Put u<sub>A</sub> = K<sub>a</sub>/A and u<sub>B</sub> = K<sub>b</sub>/B, how far each substrate is below its own half-saturation level. The last term is then k·u<sub>A</sub>·u<sub>B</sub> with k = K<sub>ia</sub>/K<sub>a</sub>, so</p>
<div class="eq">V / v = 1 + (u<sub>A</sub> ✳<sub>k</sub> u<sub>B</sub>)</div>
<p>I did not take the textbook formula on trust. I wrote down the elementary steps of each mechanism, solved the steady state exactly for 300 random sets of rate constants, and compared with the dial form: the ordered mechanism matches to 2×10⁻¹² and the ping-pong mechanism to 9×10⁻¹⁴. Ping-pong has no product term at all, dial 0. Random binding to two independent sites factors as (1 + u<sub>A</sub>)(1 + u<sub>B</sub>), dial 1, plain ✳. Lean proves the identity over any field.</p>
<p>A real enzyme: rat liver mitochondrial malate dehydrogenase has K<sub>a</sub> = 72 µM, K<sub>ia</sub> = 12 µM, K<sub>b</sub> = 86 µM, so k = 1/6. It is a sequential enzyme, but nearer the ping-pong end of the dial than the independent end.</p>"""),
("wall","The wall","2. The crossing point is the wall",
"""<p>Enzymologists plot 1/v against 1/A, one line for each level of B. Read the rate law as a line in u<sub>A</sub>: its height at zero is 1 + u<sub>B</sub> and its slope is 1 + k·u<sub>B</sub>. Changing B moves the whole line by (change in u<sub>B</sub>) × (1 + k·u<sub>A</sub>). That factor is zero at one place only:</p>
<div class="eq">u<sub>A</sub> = −1/k &nbsp; ⇔ &nbsp; 1/A = −1/K<sub>ia</sub></div>
<p>So every line passes through the same point, at height 1 − 1/k. That is the dial's wall: in the library, wAdd(k, a, −1/k) = −1/k whatever a is, and here it means that at the wall the second substrate stops mattering. Lean proves both that the lines meet there and that no other point works. At dial 0 the wall runs off to infinity and the lines are parallel, which is exactly the textbook signature of a ping-pong enzyme.</p>
<p>The crossing point sits at negative 1/A, so it is never a real concentration. Like the helix–coil wall on the previous page, it is a point of the algebra that shows up only by extending lines past the data, and that is how biochemists use it.</p>"""),
("nutrients","Nutrients","3. Two scarce nutrients, and saturation",
"""<p>A microbe short of two nutrients grows at some fraction of its maximum. Each nutrient alone gives θ = S/(K + S) = 1/(1 + u). The two common ways to combine them are the same dial:</p>
<ul><li><b>Multiplicative:</b> growth ∝ θ₁θ₂, so 1/growth = 1 + (u₁ ✳ u₂), dial 1.</li>
<li><b>Additive resistance:</b> 1/growth = 1/θ₁ + 1/θ₂ − 1 = 1 + u₁ + u₂, dial 0.</li></ul>
<p>Seen from the other side, the shortfall from full growth, 1 − θ₁θ₂, is the "or" of the two shortfalls, dial −1, the same law as two proofreading checks.</p>
<p>Liebig's law of the minimum, growth = the smaller θ, is not on the dial. With both nutrients equally scarce, the best single dial is k = 0 and it still misses Liebig by 17% at the worst point; negative dials send growth through a pole. Lean proves the dial forms; the miss is a numerical check.</p>
<p>Saturation itself is the mult-add inverse lifted: S/(K + S) = 1 + inv(K/S), and the free fraction is −inv(K/S). Kleiber's ¾ law is also a mult-add power (doubling body mass raises basal rate by 2<sup>¾</sup> − 1 = 68%), but that is a relabelling with no content, so it gets no more than this line.</p>
<figure><img src="figures/f2_nutrient_dial.png" alt="Left: growth rate against scarcity for multiplicative, additive and Liebig laws. Right: the worst gap between the dial forms and Liebig, smallest at k = 0 and still 17%."><figcaption>Left: the dial-1 and dial-0 laws fall away from Liebig's minimum as both nutrients get scarce. Right: the worst gap to Liebig for each dial setting from 0 to 1.4; no setting closes it.</figcaption></figure>
<p>Single-substrate binding, inhibition and the binding polynomial belong to the chemistry series' biochemistry page, so they are left out here.</p>"""),
],
table_h2="Where mult-add shows up in metabolism",
table=[
["Ordered two-substrate enzyme","V/v = 1 + (u<sub>A</sub> ✳<sub>k</sub> u<sub>B</sub>), dial k = K<sub>ia</sub>/K<sub>a</sub>; MDH sits at k = 1/6","Cleland 1963","Exact, known"],
["Ping-pong / random independent","Dial 0 (parallel lines) and dial 1 (plain ✳)","Cleland 1963","Exact, known"],
["Crossing of double-reciprocal lines","The wall u<sub>A</sub> = −1/k, i.e. 1/A = −1/K<sub>ia</sub>; unique","Standard way to measure K<sub>ia</sub>","Exact, known; wall seen in a lab plot"],
["Two-nutrient growth","Multiplicative is dial 1, additive resistance dial 0; shortfalls combine by 'or'","Standard growth models","Exact, known"],
["Liebig's minimum","Not on the dial; best k = 0 misses by 17%","Liebig 1840","Negative result"],
["Saturation / Kleiber","θ = 1 + inv(K/S); Kleiber is a bare mult-add power","Michaelis–Menten 1913; Kleiber 1932","Relabelling only"],
],
lean_file="biology/02-metabolism/lean/Metabolism.lean",
checks="19 of 19 pass",
sources=[
("Cleland 1963, The kinetics of enzyme-catalyzed reactions with two or more substrates or products, Biochim. Biophys. Acta 67:104",None),
("Rat liver mitochondrial malate dehydrogenase constants (K_a 72, K_ia 12, K_b 86 µM), as recorded in my working notes from Massey's 1988 thesis",None),
("Molecular biology page in this series","https://claude.ai/artifact/D4N4o6M3e7vRJuU5NkDpLh"),
],
footnote="The malate dehydrogenase constants come from a 1988 thesis I read earlier in this work; I did not re-open it for this page. The rate-law identities do not depend on them: they were checked by solving each mechanism's steady state directly.",
)
