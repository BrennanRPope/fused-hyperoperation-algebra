import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
from scipy.integrate import quad
r = json.load(open("results.json"))
eps, sig, NA = 119.8, 3.405e-8, 6.02214076e23
def B2(T):
    f = lambda x: (math.exp(-4*eps/T*((sig/x)**12 - (sig/x)**6)) - 1)*x*x
    return -2*math.pi*NA*(quad(f, 1e-10, sig, limit=200)[0] + quad(f, sig, 20*sig, limit=200)[0])
Ts = [100 + 10*i for i in range(0, 91)]
lj = [[T, B2(T)] for T in Ts]
a, b = 1.355, 0.0320; Rg = 0.083144626
vdw = [[T, 1000*(b - a/(Rg*T))] for T in Ts]
TB = r["Boyle"]
picture = """
plot(document.getElementById('pic'),{x:[100,1000],y:[-160,40],ml:56,mr:%MR%,xticks:[100,200,400,600,800,1000],xlab:'temperature, K',ylab:'argon second virial coefficient B₂, cm³/mol',
 hlines:[{y:0,c:'var(--yel)'}],
 series:[{pts:%LJ%,c:'var(--blue)',name:'Lennard-Jones, from Π(1 + f)',short:'LJ'},{pts:%VDW%,c:'var(--red)',dash:'6 5',w:2,name:'van der Waals b − a/RT',short:'vdW',dy:16},
 {pts:[[%TB%,0]],c:'var(--blue)',line:false,dots:true,r:7,labels:[{x:%TB%,y:0,t:'Boyle point %TBS% K (ab initio 408 K)',dx:10,dy:-12}]}]});
""".replace("%LJ%", json.dumps(lj)).replace("%VDW%", json.dumps(vdw)).replace("%TB%", f"{TB:.1f}").replace("%TBS%", f"{TB:.0f}").replace("%MR%", "200")

# CO2 at 350 K: the four pieces vs molar volume
aC, bC, T = 3.640, 0.04267, 350
Vs = [0.10 + 0.005*i for i in range(0, 181)]
pieces = {"inv": [], "al": [], "mb": [], "co": []}
for V in Vs:
    P = Rg*T/(V-bC) - aC/V**2
    al = aC/(P*V*V); be = bC/V
    for k_, v_ in (("inv", (1/(P*V/(Rg*T))) - 1), ("al", al), ("mb", -be), ("co", -al*be)):
        if -0.6 <= v_ <= 1.2: pieces[k_].append([V, v_])
vdw_js = """
plot(document.getElementById('vdw-fig'),{x:[0.1,1.0],y:[-0.6,1.2],ml:48,mr:16,xlab:'CO₂ molar volume at 350 K, L/mol',ylab:'pieces of 1/Z − 1',
 hlines:[{y:0,c:'var(--mute)'}],
 series:[{pts:%A%,c:'var(--blue)'},{pts:%B%,c:'var(--red)'},{pts:%C%,c:'var(--yel)'},{pts:%D%,c:'var(--ink)'}],
 extra:({X,Y,W,m})=>[['var(--blue)','α = a/PV²'],['var(--red)','−β = −b/V'],['var(--yel)','corner −αβ'],['var(--ink)','1/Z − 1 = α ∘₁ (−β)']].map((c,i)=>`<line x1="${W-m.r-170}" x2="${W-m.r-146}" y1="${m.t+14+i*20}" y2="${m.t+14+i*20}" stroke="${c[0]}" stroke-width="4"/><text x="${W-m.r-140}" y="${m.t+19+i*20}">${c[1]}</text>`).join('')});
""".replace("%A%", json.dumps(pieces["al"])).replace("%B%", json.dumps(pieces["mb"])).replace("%C%", json.dumps(pieces["co"])).replace("%D%", json.dumps(pieces["inv"]))

zc = r["Zc"]
cubic_js = """
plot(document.getElementById('cubic-fig'),{x:[-1,10],y:[0.25,0.40],ml:52,mr:120,xticks:[0,1,8],xlab:'discriminant u² − 4w of the attraction denominator',ylab:'critical compressibility Zc',
 hlines:[{y:0.29,c:'var(--green)',label:'real gases ≈ 0.29',short:'real'}],
 series:[{pts:%P%,c:'var(--blue)',line:false,dots:true,r:7,labels:[{x:0,y:%Z0%,t:'van der Waals',dx:10,dy:-8},{x:1,y:%Z1%,t:'Redlich–Kwong',dx:10,dy:-8},{x:8,y:%Z8%,t:'Peng–Robinson',dx:-6,dy:-12,a:'end'}]}]});
""".replace("%P%", json.dumps([[0, zc["van der Waals"]], [1, zc["Redlich–Kwong / SRK"]], [8, zc["Peng–Robinson"]]])).replace("%Z0%", str(zc["van der Waals"])).replace("%Z1%", str(zc["Redlich–Kwong / SRK"])).replace("%Z8%", str(zc["Peng–Robinson"]))
c = r["CO2_terms"]; h = r["haber"]
sections = [
 dict(id="vdw", nav="van der Waals", h="van der Waals' own equation is a circle operation (called mult-add before 1.0) of attraction and excluded volume",
  body=f"""<p>Van der Waals wrote his law as a product, (P + a/V²)(V − b) = RT. Divide both sides by RT·V and it says, with no approximation,</p>
<div class="eq">1/Z = (1 + α)(1 − β),  α = a/PV²,  β = b/V</div>
<p>so 1/Z − 1 = α ∘₁ (−β): the attraction correction composed by the circle operation with minus the excluded-volume correction. Checked for five gases in four states each, and proved in Lean from the equation itself.</p>
<p>The corner −αβ is not a rounding term. For CO₂ at 300 K and 0.2 L/mol, α = {c['alpha']:.2f}, β = {c['beta']:.2f} and the corner is {c['corner']:.2f}. To first order the two corrections simply add, which is why the second virial coefficient is b − a/RT. The corner first enters at the third virial coefficient.</p>""",
  fig=True, figcls="", cap="The three pieces of 1/Z − 1 for CO₂ at 350 K. At large volume α and −β nearly cancel; squeezed, the corner grows to a real share.", js=vdw_js),
 dict(id="mayer", nav="Mayer clusters", h="The Boltzmann factor is the circle-operation lift, and Mayer's clusters are its corners",
  body=f"""<p>For particles with pair energies u, the Boltzmann factor is a product over pairs of e<sup>−u/kT</sup> = 1 + f. A product of (1 + f)'s is the lift of a circle operation fold, so</p>
<div class="eq">Π (1 + f<sub>ij</sub>) = 1 + (f₁ ∘₁ f₂ ∘₁ … ∘₁ f<sub>m</sub>)</div>
<p>and expanding the fold gives one term for each nonempty set of pairs, 2<sup>m</sup> − 1 in all. Those sets are Mayer's cluster graphs, so the cluster expansion is the circle operation's corner expansion. Lean proves the three-pair case term by term and the general list version.</p>
<p>The first term gives the second virial coefficient. For argon with a Lennard-Jones potential it changes sign at {r['Boyle']:.0f} K, within half a percent of the 408 K from the 2023 ab initio argon potential (the picture at the top).</p>"""),
 dict(id="cubic", nav="Cubic equations", h="The cubic equations of state differ only in a quadratic denominator, placed here for the record",
  body=f"""<p>Redlich–Kwong and Peng–Robinson keep van der Waals' repulsion and change the attraction's denominator from V² to V² + u·bV + w·b²: (0, 0) for van der Waals, (1, 0) for Redlich–Kwong and Soave, and (2, −1) for Peng–Robinson. The discriminants are 0, 1 and 8. Their critical compressibilities come out at {zc['van der Waals']:.3f}, {zc['Redlich–Kwong / SRK']:.3f} and {zc['Peng–Robinson']:.3f}, moving toward the real value near 0.29.</p>
<p>A quadratic form's discriminant is what sorts the landscape's three countries, so van der Waals sits on the border and the other two in the two-walls country. That placement is formal: I found no physics in it beyond the known fact that a third parameter (Patel–Teja and others) is needed to match Zc.</p>""",
  fig=True, figcls="", cap="Critical compressibility of the three classic cubic equations against the discriminant of their attraction denominator, with real gases near 0.29.", js=cubic_js),
 dict(id="haber", nav="Haber conditions", h="At Haber conditions the gases are nearly ideal, and the corrections multiply",
  body=f"""<p>For nitrogen fixation the gas laws matter at the Haber process's 700 K and 200 bar. Van der Waals gives fugacity coefficients of {h['N₂']['phi']:.3f} for N₂, {h['H₂']['phi']:.3f} for H₂ and {h['NH₃']['phi']:.3f} for NH₃. They enter the equilibrium as a product, φ<sub>NH₃</sub>²/(φ<sub>N₂</sub>φ<sub>H₂</sub>³) = {r['K_phi']:.2f}, so their logarithms add: the lift once more. On this van der Waals estimate, non-ideality raises the ammonia yield at a given K by about {100*(1/r['K_phi']-1):.0f}%. Equilibrium itself comes in subject 8.</p>"""),
]
claims = [
 ("vdW's (P + a/V²)(V − b) = RT is exactly 1/Z = (1 + α)(1 − β), so 1/Z − 1 = α ∘₁ (−β)", "tidy reading", "5 gases × 4 states; Lean <code>vdw_inverse_Z</code>, <code>two_corrections</code>"),
 ("The corner −αβ first enters at the third virial coefficient; B₂ = b − a/RT", "known", "Python"),
 ("Π(1 + f) = 1 + fold of f under the circle operation; Mayer clusters are the corners", "tidy reading", "Python; Lean <code>mayer_three</code>, <code>mayer_clusters</code>, <code>mayer_list</code>"),
 ("Lennard-Jones argon: Boyle temperature 409 K vs 408 K ab initio", "known", "numerical integration"),
 ("vdW Zc = 3/8; RK 1/3; PR 0.307", "known", "Python; Lean <code>critical_Z</code>"),
 ("Cubic denominators have discriminants 0, 1, 8 (placement only)", "tidy reading", "Lean <code>cubic_discriminants</code>"),
 ("Haber at 700 K, 200 bar: vdW K_φ ≈ 0.59", "known", "Python (vdW estimate)"),
]
lean = [("two_corrections, corner","(1 + α)(1 − β) = 1 + α ∘₁ (−β), with corner −αβ"),("vdw_inverse_Z","van der Waals' equation gives the factorization"),
 ("critical_Z","Zc = 3/8"),("mayer_three, mayer_clusters, mayer_list","the Boltzmann product is the lift"),
 ("vdw_denominator, rk_denominator, pr_denominator, cubic_discriminants","the three cubic denominators"),("dalton, fugacity_product","partial pressures add, fugacity coefficients multiply")]
n = build("index.html", title="Gas laws", tile=("6","Ar","Gas laws"),
 h1="Van der Waals' own equation says 1/Z − 1 = attraction ∘₁ (−excluded volume), and Mayer's clusters are the circle operation's corners",
 lead="Written the way van der Waals wrote it, his equation factors exactly as 1/Z = (1 + α)(1 − β). The Boltzmann factor of a gas is a product of (1 + f) over pairs, so the cluster expansion is the circle-operation corner expansion. The cubic equations of state sit on the landscape only formally.",
 figs=[("15 / 15","checks pass"),("13","Lean theorems, compiled"),(f"{r['Boyle']:.0f} K","argon Boyle point, LJ (408 ab initio)"),("2ᵐ − 1","Mayer clusters = corners")],
 picture=picture, caption="Argon's second virial coefficient from the Mayer product with a Lennard-Jones potential (blue), against van der Waals' b − a/RT with textbook constants (red, dashed). Below zero attraction wins; above it the particles' size wins.",
 verdict_h="Lead check: not a lead",
 verdict="The van der Waals factorization is exact, but it is van der Waals' own product form divided through, so it predicts nothing new. The Mayer reading is exact and is what the cluster expansion already is. The cubic-equation placement is formal.",
 sections=sections, claims=claims, lean_file="chemistry/06-gas-laws/lean/ChemGas.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/06-gas-laws/: <code>checks.py</code> (15 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemGas.lean</code>, <code>make_page.py</code>.",
 sources='Sources: <a href="https://arxiv.org/pdf/2304.14719">Thermophysical properties of argon gas from an improved two-body potential (2023)</a> for the 408 K sign change; van der Waals (1873); Mayer and Mayer, Statistical Mechanics (1940); Redlich and Kwong (1949); Peng and Robinson (1976); textbook van der Waals constants.',
 series_nav="Chemistry pass, part 6 of 13: gas laws. Before: stoichiometry. Next: chemical thermodynamics.")
print("page bytes", n)
