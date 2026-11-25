import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
F, Rg, T = 96485.33212, 8.314462618, 298.15
f = F/(Rg*T)
y = r["kok_yields"]
damp = r["kok_damp"]; mean = (1 - 0.10 + 0.05)/4
env_hi = [[n+1, mean + (y[2]-mean)*damp**(n-2)] for n in range(2, 16)]
picture = "plot(document.getElementById('pic'),{x:[0.5,16.5],y:[0,0.62],ml:52,mr:120,xticks:[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],yticks:[0,0.2,0.4,0.6],xlab:'flash number',ylab:'O₂ per flash (fraction of centres)',series:%s,hlines:[{y:%f,c:'var(--mute)',dash:'4 4',label:'long run'}],extra:({X,Y})=>`<text x=\"${X(3.3)}\" y=\"${Y(0.57)}\">peak on flash 3, then about every 4th</text><text x=\"${X(8)}\" y=\"${Y(0.5)}\">the swing decays by a factor 0.851 per flash</text>`});" % (
 json.dumps([dict(pts=[[i+1, v] for i, v in enumerate(y)], c='var(--blue)', dots=True, r=5, name='', short='')]), mean)
Es = [-150 + 2*i for i in range(151)]
def frac(E_mV, k):
    x = math.exp(-f*E_mV/1000)/math.sqrt(k); Z = 1 + 2*x + k*x*x
    return (x + k*x*x)/Z
curves = []
for k, c, lab in ((0.02, 'var(--green)', 'k = 0.02 (ΔE = 136 mV)'), (1, 'var(--blue)', 'k = 1 (ΔE = 35.6 mV)'), (25, 'var(--red)', 'k = 25 (inverted)')):
    curves.append(dict(pts=[[E, frac(E, k)] for E in Es], c=c, w=2.5, name='', short=''))
two_js = """
plot(document.getElementById('two-fig'),{x:[-150,150],y:[0,1],ml:48,mr:20,xticks:[-150,-100,-50,0,50,100,150],xlab:'E − E_mid (mV), midpoint where half the electrons are on',ylab:'fraction of electrons taken up',series:%s,
 extra:({X,Y})=>`<text x="${X(8)}" y="${Y(0.93)}" fill="var(--red)">k = 25: inverted, n = 1.67</text><text x="${X(8)}" y="${Y(0.84)}" fill="var(--blue)">k = 1: independent, n = 1</text><text x="${X(8)}" y="${Y(0.75)}" fill="var(--green)">k = 0.02: two waves</text>`});
""" % json.dumps(curves)
sections = [
 dict(id="nernst", nav="Nernst and Frost", h="ΔG adds and K multiplies, but potentials are averaged, not composed by the circle operation",
  body=f"""<p>Coupling two half-reactions adds their ΔG = −nFE, so their equilibrium constants multiply and K − 1 composes by the circle operation, exactly as in subject 7. The same holds for the forward/backward rate ratio of each step, e<sup>A/RT</sup>: in series the ratios multiply, so ratio − 1 composes by the circle operation (De Donder's relation).</p>
<p>The potentials themselves do something else. Combining two couples gives the electron-weighted mean, Luther's rule E = (n₁E₁ + n₂E₂)/(n₁ + n₂). From standard tables, Cu²⁺/Cu⁺ (0.153 V) and Cu⁺/Cu (0.521 V) give {r['cu2_cu']:.3f} V for Cu²⁺/Cu (tables list 0.342), and iron gives {r['fe3_fe']:.3f} V (tables list −0.037). A mean is idempotent, a couple combined with itself returns its own potential, so this is averaging, not compose by the circle operation. On a Frost diagram it is the slope rule, and a middle state above the chord disproportionates: Cu⁺ does, with K ≈ {r['cu_disp_K']:.1e}.</p>
<p>Butler–Volmer current is the circle operation in the overpotential only in the one-sided limit α = 1, where i/i₀ = e<sup>Fη/RT</sup> − 1. At α = ½ the gap is large ({r['bv_gap_half']:.2f} in the test case); only the forward/backward ratio composes. The Tafel slope there is {r['tafel_mV']:.0f} mV per decade.</p>"""),
 dict(id="two", nav="Two redox sites", h="Two redox sites: Z − 1 = x ⊞ₖ x, and 35.6 mV marks independence",
  body=f"""<p>A molecule with two one-electron sites has the electron binding polynomial Z = 1 + 2x + kx², with x the electron activity and k the interaction factor. So Z − 1 = x ⊞ₖ x, the dial again. Matching it to the two measured potentials gives E₁ − E₂ = (RT/F) ln(4/k):</p>
<ul><li><b>k = 1, independent sites:</b> the potentials sit {r['dE_indep_mV']:.1f} mV apart, purely from counting.</li>
<li><b>k &gt; 1, cooperative:</b> ΔE &lt; 35.6 mV, the Nernst curve is steeper than one electron, and Z has complex zeros.</li>
<li><b>k &gt; 4, inverted:</b> the second electron goes on more easily than the first.</li></ul>
<p>The midpoint steepness is n = 2√k/(1 + √k). This is subject 8's diacid with electrons in place of protons (there k = 4·10<sup>−ΔpKa</sup>; here k = 4/K<sub>c</sub>, the comproportionation constant). The 35.6 mV statistical spacing is well known in electrochemistry.</p>""",
  fig=True, figcls="", cap="Fraction of electrons taken up by a two-site molecule against potential, for three interaction factors. Small k gives two separate one-electron waves; k = 1 is the independent case; large k merges them into one steep two-electron wave.", js=two_js),
 dict(id="photo", nav="Photosynthesis", h="Photosynthesis: losses compose at k = −1, the Z-scheme spends two photons per electron, and Kok's clock is a circulant",
  body=f"""<p>Any chain of stages loses energy by failure-union: total loss = L₁ ⊞₋₁ L₂ ⊞₋₁ …, the same law as yields in subject 5. The published ceilings for converting sunlight into biomass are 4.6% for C3 and 6.0% for C4 plants, or 9.4% and 12.3% of photosynthetically active light; both pairs imply that about {100*r['par_fraction']['C3']:.0f}% of sunlight is usable, so they are consistent.</p>
<p>The Z-scheme uses 8 photons per O₂ (4 at 680 nm and 4 at 700 nm, 14.4 eV) to lift 4 electrons from water (+0.815 V) to NADP⁺ (−0.320 V), storing {100*r['zscheme_eff']:.1f}%. Thermodynamics alone would allow one 680 nm photon per electron ({100*r['one_photon_eff']:.0f}%), so the scheme, not the energy, sets the cost.</p>
<p>Kok's water-splitting clock advances 0, 1 or 2 steps per flash (misses α, double hits β). Its transition matrix is circulant, so the oscillation decays by |α − β + i(1 − α − β)| per flash, {damp:.3f} for α = 10%, β = 5% (picture at the top), and settles at (1 − α + β)/4 O₂ per flash. Composing flashes multiplies the one-flash generating function, which is the lift, but that is all it is.</p>"""),
]
claims = [
 ("Coupled steps: K and forward/backward ratios multiply, so ratio − 1 composes by the circle operation", "tidy reading", "Python; Lean <code>ratio_series</code>"),
 ("Luther's rule averages potentials; it is idempotent, so not compose by the circle operation", "known", "copper and iron tables; Lean <code>luther_balance</code>, <code>luther_idem</code>"),
 ("Butler–Volmer current circle operations only at α = 1", "tidy reading", "Python; Lean <code>bv_one_sided</code>"),
 ("Two redox sites: Z − 1 = x ⊞ₖ x, ΔE = (RT/F) ln(4/k), 35.6 mV at k = 1", "known", "Python; Lean <code>two_site</code>, <code>two_site_disc</code>, <code>nernst_midpoint</code>"),
 ("Stage losses compose at k = −1; published C3/C4 ceilings are mutually consistent", "tidy reading", "Python; Lean <code>stage_losses</code>, <code>stage_losses3</code>"),
 ("Z-scheme stores 31.6% of 8 photons; one photon per electron would allow 62%", "known", "Python; Lean <code>zscheme_bounds</code>"),
 ("Kok clock damps by |α − β + i(1 − α − β)| per flash, at most 1", "known", "Python; Lean <code>kok_mean</code>, <code>kok_damped</code>"),
]
lean = [("ratio_series, bv_one_sided","ratios multiply, so ratio − 1 composes by the circle operation"),("luther_balance, luther_idem","Luther's mean, and why it is not compose by the circle operation"),
 ("two_site, two_site_disc","two redox sites are the dial, complex zeros for k > 1"),("nernst_midpoint","midpoint at x = 1/√k"),
 ("stage_losses, stage_losses3","losses compose at k = −1"),("kok_mean, kok_damped","Kok clock mean advance and damping"),("zscheme_bounds","Z-scheme efficiency between 30% and a third")]
n = build("index.html", title="Electrochemistry and Nernst", tile=("10","E","Electro"),
 h1="Electrochemistry averages potentials and multiplies constants, two redox sites are the dial, and photosynthesis loses energy by failure-union",
 lead="Coupled half-reactions multiply K, so K − 1 composes by the circle operation (called mult-add before 1.0), but the potentials combine by an electron-weighted mean. Two one-electron sites give Z − 1 = x ⊞ₖ x, with 35.6 mV marking independent sites. In photosynthesis stage losses compose at k = −1, the Z-scheme stores 31.6% of its photons, and Kok's four-flash clock is a damped circulant.",
 figs=[("22 / 22","checks pass"),("12","Lean theorems, compiled"),("35.6 mV","independent two-site spacing"),("31.6%","Z-scheme storage")],
 picture=picture, caption="Oxygen released on each flash by dark-adapted photosystem II in Kok's model (25% S₀, 75% S₁, 10% misses, 5% double hits). Yield peaks on flash 3 and about every fourth flash after, decaying toward the long-run value.",
 verdict_h="Lead check: not a lead",
 verdict="Everything here checks out, but each piece is either a restatement (K − 1 composes by the circle operation, losses fold at k = −1) or a known result (Luther's rule, the 35.6 mV statistical spacing, Kok's model). The two-site dial is the same structure as subject 8 and adds nothing beyond it. I didn't follow any of it further.",
 sections=sections, claims=claims, lean_file="chemistry/10-electrochemistry/lean/ChemElectro.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/10-electrochemistry/: <code>checks.py</code> (22 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemElectro.lean</code>, <code>make_page.py</code>.",
 sources='Sources: standard reduction potential tables (CRC Handbook); Zhu, Long and Ort, Curr. Opin. Biotechnol. 19, 153 (2008) for the 4.6% and 6.0% ceilings, with the PAR-based 9.4% and 12.3% as listed by BioNumbers; Kok, Forbush and McGloin, Photochem. Photobiol. 11, 457 (1970) for the S-state clock; standard texts for Butler–Volmer and the two-site 35.6 mV spacing.',
 series_nav="Chemistry pass, part 10 of 13: electrochemistry. Before: rates. Next: acids and bases.")
print("page bytes", n)
