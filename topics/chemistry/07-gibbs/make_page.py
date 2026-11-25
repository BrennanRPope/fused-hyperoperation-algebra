import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
ks = [10**(-3 + 0.02*i) for i in range(0, 275)]
def cov(k, x=1.0, y=1.0):
    Z = 1 + x + y + k*x*y
    return x*y*(k-1)/Z**2
picture = """
plot(document.getElementById('pic'),{x:[0.001,300],y:[-0.13,0.08],xlog:true,ml:62,mr:%MR%,xticks:[0.001,0.01,0.1,1,10,100],yticks:[-0.1,-0.05,0,0.05],xlab:'dial k = interaction Boltzmann factor e^(−ε/RT)',ylab:'covariance of the two occupancies (x = y = 1)',
 bands:[{x0:0.001,x1:1,fill:'var(--red)',label:'repulsion'},{x0:1,x1:300,fill:'var(--blue)',label:'attraction'}],
 hlines:[{y:0,c:'var(--mute)'}],
 series:[{pts:%C%,c:'var(--ink)'},{pts:[[0.001,%C0%]],c:'var(--red)',line:false,dots:true,r:6,labels:[{x:0.001,y:%C0%,t:'k → 0: exclusion, Z = 1 + x + y (addition)',dx:10,dy:18}]},
 {pts:[[1,0]],c:'var(--yel)',line:false,dots:true,r:8,labels:[{x:1,y:0,t:'k = 1: independent, Z = (1+x)(1+y)',dx:-10,dy:-14,a:'end'}]}]});
""".replace("%C%", json.dumps([[k, cov(k)] for k in ks])).replace("%C0%", str(cov(0.001))).replace("%MR%", "40")
dH, dS = r["dH"], r["dS"]
Ts = [250 + 5*i for i in range(0, 131)]
haber_js = """
plot(document.getElementById('haber-fig'),{x:[250,900],y:[-50,100],ml:52,mr:40,xticks:[300,400,500,600,700,800,900],xlab:'temperature, K',ylab:'ΔG° for N₂ + 3H₂ → 2NH₃, kJ',
 bands:[{x0:650,x1:770,fill:'var(--yel)',op:.25,label:'Haber reactors'}],hlines:[{y:0,c:'var(--mute)'}],
 series:[{pts:%G%,c:'var(--blue)'},{pts:[[%TS%,0]],c:'var(--red)',line:false,dots:true,r:7,labels:[{x:%TS%,y:0,t:'K = 1 at %TSS% K',dx:-10,dy:-14,a:'end'}]},
 {pts:[[298.15,%G298%]],c:'var(--blue)',line:false,dots:true,r:6,labels:[{x:298.15,y:%G298%,t:'−32.7 kJ at 298 K',dx:10,dy:16}]}]});
""".replace("%G%", json.dumps([[T, (dH - T*dS)/1000] for T in Ts])).replace("%TS%", f"{r['T_K1']:.1f}").replace("%TSS%", f"{r['T_K1']:.0f}").replace("%G298%", f"{r['dG298']/1000:.2f}")
K = r["K"]
sections = [
 dict(id="lift", nav="ΔG and K", h="Free energies add, so equilibrium constants multiply: K is the lift side",
  body=f"""<p>When two reactions are added, their free energies add and their equilibrium constants multiply, because K = e<sup>−ΔG/RT</sup>. In the circle operation terms the excess K − 1 combines by the circle operation, (K₁ − 1) ∘₁ (K₂ − 1) = K₁K₂ − 1, exactly as for any quantity whose logarithm adds.</p>
<p>This is how cells run uphill reactions. Making glutamine from glutamate and ammonium costs +14.2 kJ/mol (K′ = {r['K_gln']:.1e}). Coupled to ATP hydrolysis at −30.5 kJ/mol (K′ = {r['K_atp']:.1e}), the combined step has K′ = {r['K_coupled']:.0f}. The energies add and the constants multiply.</p>"""),
 dict(id="dial", nav="The dial as interaction", h="The dial is the Boltzmann factor of an interaction",
  body="""<p>This is where the dial shows up as a physical quantity. Take two binding sites with weights x and y (concentration times affinity), and let a pair of occupants interact with energy ε, so the doubly-occupied state carries an extra factor ω = e<sup>−ε/RT</sup>. The partition function is</p>
<div class="eq">Z = 1 + x + y + ω·xy,  so  Z − 1 = x ⊞<sub>ω</sub> y</div>
<p>which is the weighted circle operation with the dial set to ω. The three familiar cases are the three special settings:</p>
<ul><li><b>k = 0:</b> hard exclusion. Only one site can be filled, Z = 1 + x + y, and the occupations simply add.</li>
<li><b>k = 1:</b> no interaction. Z = (1 + x)(1 + y), the sites are independent, and the law is the circle operation.</li>
<li><b>k &gt; 1:</b> attraction, so the occupants help each other (positive cooperativity).</li></ul>
<p>The dial's distance from 1 is the correlation between the sites: the covariance of the two occupancies is exactly xy(k − 1)/Z², which Lean proves. At room temperature a 5 kJ/mol repulsion sets k = 0.13, and an equal attraction sets k = 7.5. Equilibrium (next subject) takes this to many sites.</p>"""),
 dict(id="haber", nav="Nitrogen fixation", h="Making ammonia is downhill at room temperature and uphill where it is made",
  body=f"""<p>From standard tables, N₂ + 3H₂ → 2NH₃ has ΔH° = {dH/1000:.1f} kJ and ΔS° = {dS:.1f} J/K. At 298 K, ΔG° = {r['dG298']/1000:.1f} kJ, matching twice the tabulated −16.4 kJ/mol for ammonia, and K ≈ {K['298.15']:.1e}. The reaction gives off heat but loses entropy (four gas molecules become two), so ΔG° crosses zero near {r['T_K1']:.0f} K in the constant-ΔH approximation. At the 700 K of a Haber reactor, K has fallen to about {K['700']:.0e} per bar².</p>
<p>So thermodynamics wants nitrogen fixed cold, and the reason for running hot is the rate (subject 9): the N≡N bond's hard first step (subject 4) is too slow at room temperature even on iron. The high pressure pushes the equilibrium back (subject 8).</p>""",
  fig=True, figcls="", cap="ΔG° for ammonia synthesis against temperature, from standard ΔH° and ΔS°. It is negative (downhill) below about 463 K, where industrial reactors cannot run fast enough.", js=haber_js),
]
claims = [
 ("ΔG adds, K multiplies, K − 1 composes by the circle operation", "tidy reading", "Python; Lean <code>K_excess</code>, <code>K_three</code>, <code>gibbs_adds</code>"),
 ("Glutamine synthesis coupled to ATP: K′ 3.2e−3 × 2.2e5 = 7.2e2", "known", "Python (textbook ΔG°′)"),
 ("Two interacting sites: Z − 1 = x ⊞_ω y; ω = 0 adds, ω = 1 composes by the circle operation", "tidy reading", "Python; Lean <code>Z2_dial</code>, <code>exclusion_adds</code>, <code>independent_factors</code>"),
 ("Occupancy covariance = xy(k − 1)/Z², so the sign of k − 1 is the sign of the correlation", "tidy reading", "Python; Lean <code>covariance</code>, <code>covariance_pos</code>"),
 ("Haber: ΔH° −91.8 kJ, ΔS° −198.1 J/K, ΔG°(298) −32.7 kJ, K = 1 near 463 K", "known", "standard tables; Lean <code>haber_enthalpy</code>, <code>haber_entropy</code>"),
]
lean = [("K_excess, K_three","equilibrium constants' excesses the circle operation"),("Z2, Z2_dial","the two-site partition function is the dial"),
 ("exclusion_adds, independent_factors, independent_circ","k = 0 adds, k = 1 factors"),("covariance, covariance_pos","the covariance is xy(k − 1)/Z²"),
 ("haber_entropy, haber_enthalpy, gibbs_adds","Haber's ΔS° and ΔH°; ΔG adds")]
n = build("index.html", title="Chemical thermodynamics", tile=("7","ΔG","Gibbs"),
 h1="Free energies add so equilibrium constants multiply, and the dial is the Boltzmann factor of an interaction",
 lead="Adding reactions adds ΔG and multiplies K, so K − 1 composes by the circle operation (called mult-add before 1.0). Two binding sites with an interaction have Z − 1 = x ⊞ₖ y with k = e^(−ε/RT): k = 0 is exclusion (addition), k = 1 independence (the circle operation), and k − 1 is the correlation between the sites. Ammonia synthesis is downhill only below about 463 K.",
 figs=[("11 / 11","checks pass"),("11","Lean theorems, compiled"),("k = e^(−ε/RT)","what the dial is"),("463 K","where K = 1 for ammonia")],
 picture=picture, caption="How correlated two sites are, against the dial k. Below 1 the occupants repel (k → 0 is full exclusion, ordinary addition); above 1 they attract. At k = 1 they are independent and the law is the circle operation.",
 verdict_h="Lead check: a useful reading, not yet a lead",
 verdict="The dial gets a plain physical meaning here, the Boltzmann factor of an interaction, with addition and the circle operation as its exclusion and independence settings. That is exact and it organizes the next subjects, but on two sites it is ordinary statistical mechanics. Whether it says anything new depends on many sites, which subject 8 tests.",
 sections=sections, claims=claims, lean_file="chemistry/07-gibbs/lean/ChemGibbs.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/07-gibbs/: <code>checks.py</code> (11 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemGibbs.lean</code>, <code>make_page.py</code>.",
 sources='Sources: NIST Chemistry WebBook and CRC Handbook standard values for NH₃, N₂, H₂; Lehninger, Principles of Biochemistry, for the glutamine synthetase and ATP ΔG°′.',
 series_nav="Chemistry pass, part 7 of 13: chemical thermodynamics. Before: gas laws. Next: equilibrium and mass action.")
print("page bytes", n)
