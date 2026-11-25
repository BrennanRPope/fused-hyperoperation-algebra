import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
alpha = r["alpha"]
pts = [("HCl bond", 9.856e-05), ("H₂ bond", r["table"][0][1]), ("hydrogen atom", r["hydrogen_theta"]),
       ("positronium", alpha*1.000004), ("He⁺ ion", 2*alpha*1.000031), ("deuteron", r["table"][3][1])]
rows = [[p[1], i+1] for i, p in enumerate(pts)]
labels = [dict(x=p[1], y=i+1, t=p[0], dx=12) for i, p in enumerate(pts)]
picture = f"""
plot(document.getElementById('pic'),{{xlog:true,x:[3e-5,1],y:[0,7],mr:40,ml:20,xticks:[1e-4,1e-3,1e-2,1e-1,1],yticks:[],
 xfmt:v=>({{'0.0001':'10⁻⁴','0.001':'10⁻³','0.01':'10⁻²','0.1':'10⁻¹','1':'1 rad'}})[String(v)],
 xlab:'binding angle θ, where the whole weighs (m₁² + m₂² + 2m₁m₂ cos θ)^½',
 bands:[{{x0:3e-5,x1:3e-4,fill:'var(--blue)',label:'chemical bonds'}},{{x0:5e-3,x1:2e-2,fill:'var(--yel)',op:.35,label:'Coulomb atoms: θ = Zα'}},{{x0:5e-2,x1:0.3,fill:'var(--red)',label:'nuclei'}}],
 series:[{{pts:{json.dumps(rows)},line:false,dots:true,r:7,c:'var(--ink)',labels:{json.dumps(labels)}}}],
 extra:({{X,Y,H,m}})=>`<line x1="${{X({alpha})}}" x2="${{X({alpha})}}" y1="${{m.t+24}}" y2="${{H-m.b}}" stroke="var(--red)" stroke-dasharray="4 4"/><text x="${{X({alpha})+6}}" y="${{H-m.b-8}}" fill="var(--red)">α = 1/137.036</text>`}});
"""
curve = []
for i in range(0, 301):
    s = -math.pi + i*(math.pi+2.0)/300
    c = math.cos(-s) if s < 0 else math.cosh(s)
    curve.append([round(s, 4), round(math.sqrt((1+c)/2), 5)])
countries_js = f"""
plot(document.getElementById('countries-fig'),{{x:[-3.1416,2],y:[0,2],mr:20,ml:44,xticks:[-3,-2,-1,0,1,2],xlab:'angle θ (bound, left)   ·   rapidity η (free, right)',ylab:'M ÷ (m₁ + m₂), equal masses',
 bands:[{{x0:-3.1416,x1:0,fill:'var(--blue)',label:'turning: bound'}},{{x0:0,x1:2,fill:'var(--red)',label:'two walls: free'}}],
 series:[{{pts:{json.dumps(curve)},c:'var(--ink)'}}],
 extra:({{X,Y,H,m}})=>`<line x1="${{X(0)}}" x2="${{X(0)}}" y1="${{m.t+20}}" y2="${{H-m.b}}" stroke="var(--yel)" stroke-width="4"/><text x="${{X(0)-8}}" y="${{Y(1.55)}}" text-anchor="end">border: Lavoisier,</text><text x="${{X(0)-8}}" y="${{Y(1.55)+18}}" text-anchor="end">M = m₁ + m₂</text>`}});
"""
sections = [
 dict(id="countries", nav="Three countries", h="Invariant mass is a law of cosines, and Lavoisier sits on the border",
  body="""<p>Two parts with masses m₁ and m₂ make a whole whose mass obeys</p>
<div class="eq">M² = m₁² + m₂² + 2 m₁ m₂ c</div>
<p>This is the law of cosines, the line q = 1 on the landscape. The value of c decides the country:</p>
<ul><li><b>c = 1, the border.</b> M = m₁ + m₂ exactly: Lavoisier's conservation of mass. The triangle is flat.</li>
<li><b>c = cos θ &lt; 1, turning.</b> A bound state: the whole weighs less than its parts, and θ is a real angle.</li>
<li><b>c = cosh η &gt; 1, two walls.</b> Two free particles flying apart: the pair weighs more, and η is their relative rapidity, the ruler from the relativity page.</li></ul>
<p>So binding is a rotation and free motion is a boost, the same split the relativity and early-laws pages found between Euclid's and Einstein's countries. Lean proves the bookkeeping: bound pairs weigh less, free pairs weigh more, and the binding energy B satisfies B(2(m₁ + m₂) − B) = 2m₁m₂(1 − c) exactly.</p>""",
  fig=True, figcls="", cap="The mass of a pair of equal parts, as a fraction of the parts' sum, across the three countries. Left of the border the pair is bound by angle θ; right of it the parts fly apart with rapidity η.", js=countries_js),
 dict(id="alpha", nav="The angle α", h="Hydrogen's binding angle is the fine-structure constant",
  body=f"""<p>Measure the angle for real bound states and the Coulomb ones line up on one value. Hydrogen's 1s state gives θ = {r['hydrogen_theta']:.6e}, against α = {alpha:.6e}: they agree to 8 parts per million. Positronium and muonium, with very different masses, give the same angle. He⁺ gives 2α. To leading order θ = √(2B/μ), and a Coulomb binding energy is ½μ(Zα)², so θ = Zα for every hydrogen-like system, whatever its masses.</p>
<p>Taking sin θ = α exactly, the compact two-body Dirac form, predicts hydrogen's ionisation energy as {r['B_pred_eV']:.6f} eV against the measured 13.598434 eV. The difference, 3.4 × 10⁻⁵ eV, is the 1s Lamb shift (a quantum-electrodynamics correction) to within 2%.</p>
<p>Chemistry sits much nearer the border: an H₂ bond has θ ≈ 1.4 × 10⁻⁴, so burning hydrogen loses only 1.8 × 10⁻¹⁰ of the water's mass. Lavoisier's balances, good to about one part in ten thousand, could never see it. That is why conservation of mass is exact for chemists and approximate for physicists: chemistry lives within a hair of the border.</p>"""),
 dict(id="charge", nav="Charge", h="Charge is the k = 0 law, and its phase is mult-add on the circle",
  body="""<p>Charges add with no corner: the k = 0 line, mult-add's own identity line. Their exponentials multiply. The phase e<sup>iQφ</sup> of a charge Q, or the fugacity z<sup>Q</sup> in a grand ensemble, multiplies as charges add, so its excess mult-adds: e<sup>i(Q₁+Q₂)φ</sup> − 1 = (e<sup>iQ₁φ</sup> − 1) ✳ (e<sup>iQ₂φ</sup> − 1). On the unit circle that is rotation, the turning country again.</p>
<p>Every additive conservation law in chemistry has this shape: the conserved quantity is the log of a lift whose excess mult-adds. A balanced redox equation such as MnO₄⁻ + 5Fe²⁺ + 8H⁺ → Mn²⁺ + 5Fe³⁺ + 4H₂O balances atoms and charge (+17 on each side) by plain addition.</p>"""),
]
claims = [
 ("Invariant mass is a law of cosines, M² = m₁² + m₂² + 2m₁m₂c", "known", "500 random four-vector pairs"),
 ("Lavoisier = the border c = 1; bound = turning; free = two walls", "tidy reading", "Lean <code>lavoisier</code>, <code>bound_lighter</code>, <code>free_heavier</code>"),
 ("Exact binding relation B(2(m₁+m₂) − B) = 2m₁m₂(1 − c)", "known", "Lean <code>binding_exact</code>"),
 ("Hydrogen-like binding angle θ = Zα, whatever the masses", "tidy reading", "H, positronium, muonium, He⁺ to ≤ 3e-4"),
 ("sin θ = α gives hydrogen's ionisation energy to 2.5 ppm; the rest is the Lamb shift", "known", "mpmath, CODATA 2018"),
 ("Chemical bonds sit at θ ~ 10⁻⁴, so chemistry is within a hair of the border", "tidy reading", "H₂ and HCl bond energies"),
 ("Burning hydrogen loses 1.8 × 10⁻¹⁰ of the mass", "known", "ΔH = −571.66 kJ per 2 mol"),
 ("Charge phases mult-add on the circle; charges add at k = 0", "tidy reading", "300 cases; Lean <code>fugacity_multAdd</code>, <code>circle_closed</code>"),
]
lean = [("lavoisier","c = 1 gives (m₁ + m₂)² exactly"),("deficit","(m₁ + m₂)² − M² = 2m₁m₂(1 − c)"),("binding_exact","B(2(m₁+m₂) − B) = 2m₁m₂(1 − c)"),
 ("bound_lighter","c &lt; 1 and positive masses: the whole is lighter"),("free_heavier","c &gt; 1: the pair is heavier"),("cosine_extra","the extra piece beyond m₁² + m₂² is 2c·m₁m₂"),
 ("small_angle","c = 1 − θ²/2 gives deficit m₁m₂θ²"),("fugacity_mul, fugacity_multAdd","z<sup>Q</sup> multiplies, its excess mult-adds"),("neutral_identity","charge 0 is mult-add's identity"),
 ("circle_closed","products of points on the unit circle stay on it")]
n = build("index.html", title="Conservation of mass and charge", tile=("2","M","Conservation"),
 h1="Conservation of mass is the border of a law of cosines, and binding is a turn by Zα",
 lead="The mass of two parts is M² = m₁² + m₂² + 2m₁m₂c. Lavoisier's exact conservation is c = 1, the border; bound atoms and molecules turn by a real angle, and free particles boost. Hydrogen's binding angle is the fine-structure constant.",
 figs=[("15 / 15","numeric checks pass"),("11","Lean theorems, compiled"),("θ = α","hydrogen's angle, to 8 ppm"),("1.8 × 10⁻¹⁰","mass lost burning hydrogen")],
 picture=picture, caption="Binding angles of real bound states on one log scale: chemical bonds near 10⁻⁴, every hydrogen-like atom at Zα, the deuteron near 0.1 rad. Computed in checks.py from measured binding energies and CODATA masses.",
 verdict_h="Lead check: a sturdy picture, not a lead to new results",
 verdict="The angle picture is exact and geometric, and it ties Lavoisier, atomic binding and free collisions into one law of cosines; it predicts hydrogen's ionisation energy to 2.5 parts per million from α alone. But it is known physics reparametrised: θ = Zα is the Coulomb binding energy written as an angle. I found no prediction it makes that the usual formulas don't. Its value for the project is the map it gives the rest of the chemistry pass: bound means turning, free means two walls.",
 sections=sections, claims=claims, lean_file="chemistry/02-conservation/lean/ChemConservation.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/02-conservation/: <code>checks.py</code> (15 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemConservation.lean</code>, <code>make_page.py</code>.",
 sources="Sources: CODATA 2018 masses and α; NIST hydrogen ionisation energy 13.598434005 eV and 1s Lamb shift 8172.874 MHz; NIST-JANAF and CRC bond energies for H₂ (4.478 eV) and HCl (4.434 eV); deuteron binding 2.224566 MeV.",
 series_nav="Chemistry pass, part 2 of 13: conservation of mass and charge. Before: quantum chemistry. Next: the periodic law.")
print("page bytes", n)
