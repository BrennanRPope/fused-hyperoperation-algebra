import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
iso = r["isoboles"]
cols = {"1000000000.0": ('var(--mute)', "α → ∞, one site (Loewe)"), "5": ('var(--green)', "α = 5"), "1": ('var(--blue)', "α = 1, independent (Bliss)"), "0.2": ('var(--red)', "α = 0.2, synergy")}
ser = [dict(pts=iso[k], c=c, w=3, name='', short='') for k, (c, _) in cols.items()]
picture = "plot(document.getElementById('pic'),{x:[0,1],y:[0,1],ml:52,mr:20,xticks:[0,0.25,0.5,0.75,1],yticks:[0,0.25,0.5,0.75,1],xlab:'first inhibitor, [I]/K_I',ylab:'second inhibitor, [J]/K_J, for 50%% inhibition',series:%s,extra:({X,Y})=>`%s`});" % (
 json.dumps(ser), "".join(f'<text x="${{X(0.62)}}" y="${{Y({0.95 - 0.07*i})}}" fill="{c}">{lab}</text>' for i, (c, lab) in enumerate(cols.values())))
ns = [0, 1, 2, 3]
pser = []
for D, c in ((10, 'var(--green)'), (100, 'var(--blue)'), (1000, 'var(--red)')):
    pser.append(dict(pts=[[n, 1/(1 + D**(n+1))] for n in ns], c=c, dots=True, r=5, w=2.5, name='', short=''))
pr_js = """
plot(document.getElementById('proof-fig'),{x:[-0.2,3.2],y:[1e-13,1],ylog:true,ml:56,mr:20,xticks:[0,1,2,3],yticks:[1,1e-3,1e-6,1e-9,1e-12],xlab:'proofreading stages added',ylab:'error rate',series:%s,
 extra:({X,Y})=>`<text x="${X(-0.1)}" y="${Y(1e-8)}" fill="var(--green)">10 : 1 per stage</text><text x="${X(-0.1)}" y="${Y(1e-9)*0.5+Y(1e-10)*0.5}" fill="var(--blue)">100 : 1 per stage</text><text x="${X(-0.1)}" y="${Y(1e-11)}" fill="var(--red)">1000 : 1 per stage</text>`});
""" % json.dumps(pser)
sections = [
 dict(id="inhibitors", nav="Two inhibitors", h="Two inhibitors: Yonetani and Theorell's α is the dial, between Loewe and Bliss",
  body="""<p>When two inhibitors I and J act on one enzyme, the four enzyme states (free, with I, with J, with both) give</p>
<div class="eq">v₀/vᵢ − 1 = a + b + ab/α = a ⊞ₖ b,  with a = [I]/K<sub>I</sub>, b = [J]/K<sub>J</sub>, k = 1/α</div>
<p>This is the Yonetani–Theorell analysis, and its interaction constant α is the dial (k = 1/α), measured routinely in enzyme kinetics. The two ends are the two standard ways pharmacologists judge drug combinations:</p>
<ul><li><b>α → ∞, one shared site (k = 0):</b> plain addition, so the 50% isobole is Loewe's straight line.</li>
<li><b>α = 1, independent sites (k = 1):</b> the circle operation (called mult-add before 1.0), so the remaining activities multiply, which is Bliss independence.</li>
<li><b>α &lt; 1:</b> the bound inhibitors help each other, which is synergy.</li></ul>
<p>Two side facts follow. By Loewe's yardstick independent inhibitors already look synergistic, which is part of the long Loewe versus Bliss debate. And no binding α can push the isobole beyond Loewe's line, because k = 1/α is never negative; antagonism of that kind needs another mechanism.</p>"""),
 dict(id="mm", nav="Michaelis–Menten", h="Michaelis–Menten is Langmuir, and steps in series add their times",
  body="""<p>Michaelis–Menten saturation has odds v/(V<sub>max</sub> − v) = S/K<sub>m</sub>, the Langmuir isotherm of subject 6. Irreversible steps in a catalytic cycle add their mean times, so the overall turnover rate is the reciprocal of the summed times (the k = 0 law, on 1/rate), checked by simulation.</p>"""),
 dict(id="proof", nav="Proofreading", h="Kinetic proofreading multiplies discrimination odds",
  body="""<p>Each proofreading stage multiplies the odds of right over wrong, so the error falls as f₀<sup>n+1</sup>: a 100 : 1 discriminator gives errors of 1%, 0.01% and 0.0001% with 0, 1 and 2 extra stages (Hopfield's scheme). In the language of subject 12 the fidelity ee composes by the velocity-addition law. The odds multiplying is the whole content; nothing further comes out.</p>""",
  fig=True, figcls="", cap="Error rate against the number of proofreading stages added, for discriminators of 10 : 1, 100 : 1 and 1000 : 1 per stage. On a log scale each stage takes off a fixed number of decades.", js=pr_js),
 dict(id="n2", nav="Nitrogenase", h="Nitrogenase spends 16 ATP and a quarter of its electrons to do what Haber does with heat",
  body=f"""<p>Nitrogenase runs N₂ + 8H⁺ + 8e⁻ + 16ATP → 2NH₃ + H₂ + 16ADP + 16Pᵢ. Six of the eight electrons reach ammonia and two go to the obligatory hydrogen, so 75% of the electrons are used. The 16 ATP supply about {r['nitrogenase_ATP_kJ']:.0f} kJ at the standard −30.5 kJ/mol, although the overall conversion of N₂ and H₂ to ammonia is downhill by 32.9 kJ. As subjects 4, 7, 8 and 9 showed, the difficulty is the barrier of the first N≡N step. Haber pays for it with heat and pressure; the enzyme pays with ATP.</p>"""),
 dict(id="photo", nav="Photosynthesis", h="Photosynthesis: the Calvin cycle stores 34% of the photons' energy at the 8-photon limit",
  body=f"""<p>Fixing one CO₂ takes 3 ATP and 2 NADPH, and at least 8 photons, a maximum quantum yield of 0.125. Eight 680 nm photons carry {8*r['E680']:.0f} kJ per mole; one CO₂ fixed into glucose stores 2870/6 = 478 kJ, so {100*r['co2_eff']:.0f}% at that limit. This sits inside the stage chain of subject 10, where losses compose at k = −1 down to the 4.6% and 6.0% ceilings for whole plants.</p>"""),
]
claims = [
 ("Yonetani–Theorell: v₀/vᵢ − 1 = a ⊞ₖ b with k = 1/α", "known (as α); dial reading here", "four-state equilibrium; Lean <code>yonetani_theorell</code>"),
 ("k = 0 is Loewe additivity, k = 1 is Bliss independence", "known", "Python; Lean <code>exclusive_loewe</code>, <code>independent_lift</code>, <code>bliss</code>, <code>bliss_fraction</code>"),
 ("50% isoboles for every α; binding cannot reach beyond Loewe's line", "tidy reading", "Python; Lean <code>isobole</code>"),
 ("Michaelis–Menten is Langmuir; series steps add times", "known", "simulation; Lean <code>series_rate</code>"),
 ("Proofreading multiplies discrimination odds", "known", "Python; Lean <code>proofreading</code>"),
 ("Nitrogenase: 75% of electrons to NH₃, 8 ATP per NH₃", "known", "Python; Lean <code>nitrogenase_electrons</code>, <code>nitrogenase_atp</code>"),
 ("Calvin cycle stores 34% of 8 photons at 680 nm", "known", "Python; Lean <code>calvin_efficiency</code>"),
]
lean = [("yonetani_theorell","the two-inhibitor law is the dial"),("exclusive_loewe, independent_lift","k = 0 and k = 1 ends"),
 ("bliss, bliss_fraction","Bliss independence"),("isobole","the 50% isobole for each k"),("series_rate","steps in series"),
 ("proofreading","odds multiply"),("nitrogenase_electrons, nitrogenase_atp, calvin_efficiency","bookkeeping")]
n = build("index.html", title="Biochemistry", tile=("13","E·S","Biochem"),
 h1="In enzymes the dial is a measured constant: Yonetani and Theorell's α runs from Loewe additivity to Bliss independence",
 lead="Two inhibitors on one enzyme give v₀/vᵢ − 1 = a ⊞ₖ b with k = 1/α, so mutually exclusive inhibitors follow Loewe's addition and independent ones follow Bliss's multiplication. Michaelis–Menten is Langmuir, proofreading multiplies odds, nitrogenase spends 16 ATP and a quarter of its electrons per N₂, and the Calvin cycle stores 34% of its photons at the 8-photon limit.",
 figs=[("13 / 13","checks pass"),("11","Lean theorems, compiled"),("k = 1/α","Yonetani–Theorell"),("75%","nitrogenase electrons to NH₃")],
 picture=picture, caption="Fifty-percent isoboles for two inhibitors: the doses of each that together halve the enzyme's rate. The straight grey line is mutual exclusion (Loewe), blue is independence (Bliss), and red is synergy.",
 verdict_h="Lead check: not a lead",
 verdict="This is the clearest place the dial appears in chemistry: it is a constant enzymologists measure, and its two ends are the two standard rules for drug combinations. But all of that is established, and the dial reading predicts nothing beyond the Yonetani–Theorell equation. The rest is bookkeeping, so I didn't follow it further.",
 sections=sections, claims=claims, lean_file="chemistry/13-biochemistry/lean/ChemBiochem.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/13-biochemistry/: <code>checks.py</code> (13 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemBiochem.lean</code>, <code>make_page.py</code>.",
 sources='Sources: Yonetani and Theorell, Arch. Biochem. Biophys. 106, 243 (1964); Loewe (1926) and Bliss (1939) on drug combinations; Greco, Bravo and Parsons, Pharmacol. Rev. 47, 331 (1995); Hopfield, PNAS 71, 4135 (1974); standard biochemistry texts for the nitrogenase and Calvin-cycle stoichiometry.',
 series_nav="Chemistry pass, part 13 of 13: biochemistry. Before: organic. This is the last subject.")
print("page bytes", n)
