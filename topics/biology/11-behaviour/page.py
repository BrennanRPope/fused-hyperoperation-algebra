SPEC=dict(
title="Behaviour and the circle operation",
h1="Learning curves are the 'or' power, and discount functions are sorted by how they compose over split delays",
lede="Rescorla–Wagner learning is the 'or' law applied once per trial, so the learning curve is wpow(−1, αβ, n) and training phases compose by the same law. In discounting, write the 'interest' a delay demands as u = 1/d − 1. Hyperbolic interest adds over split delays (dial 0), exponential interest composes by ∘₁ (dial 1), and a family in between composes by a fixed dial κ. That family is exponential measured from the wall −1/κ. Many-eyes vigilance and Herrnstein's hyperbola are the 'or' power and saturation again.",
stats=[("13 / 13","numeric checks"),("12","Lean theorems, no sorry"),("0","leads"),],
hero_fig="figures/f1_discount_learning.png",
hero_cap="Left: discount factors that compose over split delays by the dial κ, all with the same short-delay slope; κ → 0 is hyperbolic and κ = 1 exponential. Right: a Rescorla–Wagner learning curve and the 'or' power that it is.",
known="All known. Rescorla and Wagner's rule (1972), its exponential approach to the asymptote, and blocking are textbook learning theory. Hyperbolic discounting (Mazur 1987), the generalised hyperbola (Loewenstein and Prelec 1992), subadditive discounting (Read 2001), Pulliam's many-eyes model (1973) and Herrnstein's hyperbola (1970) are standard. The dial-stationary discount family is a one-line consequence of the dial algebra; I have not checked whether it appears in the discounting literature.",
lead="None. The placements are exact and tidy: subadditivity is exactly the dial's corner term k²D₁D₂, and the κ family joins hyperbolic to exponential through the dial. But there are already many two-parameter discount functions, and nothing here predicts a choice pattern that they do not.",
kicker="Biology series · 11 of 11",
sections=[
("rw","Learning","1. Rescorla–Wagner is the 'or' law per trial",
"""<p>On each trial associative strength moves a fraction αβ of the way to the asymptote λ. In units of λ that step is V ∘₋₁ αβ, so after n trials</p>
<div class="eq">V<sub>n</sub> = λ · wpow(−1, αβ, n) = λ(1 − (1 − αβ)ⁿ)</div>
<p>Two training phases of n₁ and n₂ trials compose by the 'or' law. A compound of two cues learns its summed strength as a single cue of rate α₁ + α₂, because they share one prediction error. That shared error is also why blocking happens: in the simulation, a cue pretrained near λ left 1% of λ for a new cue added later.</p>"""),
("disc","Discounting","2. Discount functions, sorted by the dial",
"""<p>Write the interest a delay D demands as u = 1/d(D) − 1, so that a delayed reward is worth A/(1 + u). How u behaves when the delay is split into D₁ then D₂ sorts the classic forms:</p>
<ul><li><b>Hyperbolic, d = 1/(1 + kD):</b> u = kD adds, dial 0. Discounting each piece separately and multiplying raises 1/d by exactly k²D₁D₂, the dial's corner term. That is Read's subadditivity: people discount a delay more when it is presented in pieces.</li>
<li><b>Exponential:</b> u composes by plain ∘₁, dial 1. This is the only form with no preference reversals.</li>
<li><b>Dial κ:</b> u(D) = (e<sup>cD</sup> − 1)/κ composes by ∘<sub>κ</sub>, because 1 + κu is exponential: this is exponential discounting measured from the wall −1/κ. κ = 1 is exponential, and κ → 0 with c/κ fixed is hyperbolic. If u is continuous and composes by a fixed κ, this is the only form, by Cauchy's equation. I did not prove that uniqueness in Lean.</li></ul>
<p>Loewenstein and Prelec's generalised hyperbola (1 + αD)<sup>−β/α</sup> is 1/(1 + mpow(αD, β/α)): a circle operation (called mult-add before 1.0) power of the hyperbolic interest, with α = β hyperbolic and α → 0 exponential.</p>"""),
("other","Groups","3. Vigilance and matching",
"""<p>A group of n animals, each spotting an approaching predator with chance p, spots it with wpow(−1, p, n): 74% for six animals at p = 0.2. Herrnstein's response-rate hyperbola R = k·r/(r + r<sub>e</sub>) is k times −inv(r/r<sub>e</sub>), the saturation that runs through the whole series.</p>"""),
],
table_h2="Where the circle operation shows up in behaviour",
table=[
["Rescorla–Wagner","Learning curve wpow(−1, αβ, n); phases compose by 'or'","Rescorla and Wagner 1972","Exact, known"],
["Hyperbolic discounting","Interest adds (dial 0); subadditivity = corner k²D₁D₂","Mazur 1987; Read 2001","Exact, known"],
["Exponential discounting","Interest composes by ∘₁ (dial 1)","Samuelson 1937","Exact, known"],
["Dial-κ discounting","Exponential from the wall −1/κ; joins the two","Follows from the dial","Exact; adds no prediction"],
["Generalised hyperbola","1/(1 + mpow(αD, β/α))","Loewenstein and Prelec 1992","Relabelling"],
["Many eyes, Herrnstein","'Or' power; saturation","Pulliam 1973; Herrnstein 1970","Relabelling"],
],
lean_file="biology/11-behaviour/lean/Behaviour.lean",
checks="13 of 13 pass",
sources=[
("Read 2001, Is time-discounting hyperbolic or subadditive? J. Risk Uncertainty 23:5",None),
("Loewenstein and Prelec 1992, Anomalies in intertemporal choice, Q. J. Econ. 107:573",None),
("Ecology page in this series","https://claude.ai/artifact/Xw4bKkWgC9u7G2EsK7vnwW"),
],
footnote="All numbers are simulated or algebraic; no behavioural data. Citations and dates are from memory and I did not open them.",
)
