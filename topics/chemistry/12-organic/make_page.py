import sys, json, math
sys.path.insert(0, "..")
from page_builder import build
r = json.load(open("results.json"))
def levels(N, mob):
    ev = sorted(round(2*math.cos(((2*j+1) if mob else 2*j)*math.pi/N), 9) for j in range(N))
    out = []; i = 0
    while i < len(ev):
        j = i
        while j < len(ev) and abs(ev[j] - ev[i]) < 1e-6: j += 1
        m = j - i
        for t in range(m): out.append((-0.07*(m-1)/1 + 0.14*t if m > 1 else 0, ev[i]))
        i = j
    return out
hu, mo = [], []
for N in range(3, 9):
    hu += [[N - 0.18 + dx, e] for dx, e in levels(N, False)]
    mo += [[N + 0.18 + dx, e] for dx, e in levels(N, True)]
picture = "plot(document.getElementById('pic'),{x:[2.5,8.5],y:[-2.3,2.5],ml:52,mr:20,xticks:[3,4,5,6,7,8],yticks:[-2,-1,0,1,2],xlab:'ring size N (left: Hückel, right: Möbius)',ylab:'orbital energy, units of β (bonding up)',series:%s,hlines:[{y:0,c:'var(--mute)'}],extra:({X,Y})=>`<text x=\"${X(2.6)}\" y=\"${Y(2.35)}\" fill=\"var(--blue)\">● Hückel ring: 2cos(2πj/N)</text><text x=\"${X(5.6)}\" y=\"${Y(2.35)}\" fill=\"var(--red)\">● Möbius ring: 2cos((2j+1)π/N)</text>`});" % json.dumps(
 [dict(pts=hu, c='var(--blue)', line=False, dots=True, r=5, name='', short=''), dict(pts=mo, c='var(--red)', line=False, dots=True, r=5, name='', short='')])
def land(a, b): return (a + b)/(1 + a*b)
es = [i/100 for i in range(0, 101)]
cols = [(0.0, 'var(--mute)', 'first step 0'), (0.5, 'var(--green)', 'first step 50%'), (0.8, 'var(--blue)', 'first step 80%'), (0.95, 'var(--red)', 'first step 95%')]
ser = [dict(pts=[[e, land(e1, e)] for e in es], c=c, w=2.5, name='', short='') for e1, c, _ in cols]
ee_js = """
plot(document.getElementById('ee-fig'),{x:[0,1],ml:44,mr:20,xticks:[0,0.2,0.4,0.6,0.8,1],yticks:[0,0.2,0.4,0.6,0.8,1],xlab:'ee of the second step',ylab:'ee of the matched pair',y:[0,1.05],series:%s,
 extra:({X,Y})=>`<text x="${X(0.55)}" y="${Y(0.36)}" fill="var(--red)">first step 95%%</text><text x="${X(0.55)}" y="${Y(0.28)}" fill="var(--blue)">first step 80%%</text><text x="${X(0.55)}" y="${Y(0.20)}" fill="var(--green)">first step 50%%</text><text x="${X(0.55)}" y="${Y(0.12)}" fill="var(--mute)">first step 0 (diagonal)</text>`});
""" % json.dumps(ser)
hor = r["horeau"]; kj = r["kJ_for_ee"]
sections = [
 dict(id="rings", nav="Rings", h="Closing the chain: Hückel and Möbius rings, and the 4n + 2 rule",
  body="""<p>Subject 4 found open Hückel chains on the landscape. Closing the chain into a ring gives levels 2cos(2πj/N) in units of β; putting a half-twist in the ring (a Möbius ring) shifts them to 2cos((2j + 1)π/N). The characteristic polynomials are 2T<sub>N</sub>(x/2) − 2 and 2T<sub>N</sub>(x/2) + 2, with T<sub>N</sub> the Chebyshev polynomial, checked for N = 3 to 10.</p>
<p>Filling the levels gives the 4n + 2 rule. Benzene has π energy 8β, 2β more than three ethylenes. Cyclobutadiene gains nothing and leaves two electrons in a degenerate pair. A Möbius four-ring reverses this: levels ±√2, π energy 5.66β and a full gap. Rings with 6 or 10 electrons are closed-shell as Hückel rings; rings with 4, 8 or 12 are closed-shell as Möbius rings (picture at the top). This is standard Hückel and Heilbronner theory, now checked and proved for the small rings.</p>"""),
 dict(id="hammett", nav="Hammett", h="Hammett: σ adds, so rate ratios multiply",
  body="""<p>Hammett's rule log(k/k₀) = ρσ makes substituent constants add (the k = 0 law), so the rate ratios of two substituents multiply and k/k₀ − 1 composes by the circle operation (called mult-add before 1.0). This is the same log-of-the-lift fact as Gibbs energies and pKas, and it adds nothing beyond them.</p>"""),
 dict(id="ee", nav="Selectivity", h="Selectivity: ee composes by the velocity-addition law, the landscape point (0, −1)",
  body=f"""<p>An enantiomer ratio r gives ee = (r − 1)/(r + 1), which is tanh of half of ln r, so ee is a hyperbolic tangent of the free-energy difference. When two selective steps act together in the matched sense their ratios multiply, so</p>
<div class="eq">ee = (e₁ + e₂)/(1 + e₁e₂)</div>
<p>which is the landscape law at (k, q) = (0, −1), and also Einstein's addition of velocities with ee playing v/c. It sits in the two-walls country, with walls at ee = ±1: a perfectly selective step absorbs anything. A 10:1 step (ee 0.818) matched with a 5:1 step (ee 0.667) gives 0.961; mismatched they give 0.333.</p>
<p>Horeau's duplication is the same law applied to itself, e ⊕ e = 2e/(1 + e²): 20% ee becomes {100*hor[1]:.0f}%, {100*hor[2]:.0f}%, {100*hor[3]:.0f}% and {100*hor[4]:.1f}% in four rounds. Free energy plays the role of rapidity, so each extra nine of ee costs close to RT ln 10 = 5.7 kJ/mol at room temperature ({kj['0.9']:.1f}, {kj['0.99']:.1f} and {kj['0.999']:.1f} kJ/mol for 90%, 99% and 99.9%).</p>""",
  fig=True, figcls="", cap="ee of a matched pair of selective steps against the second step's ee, for four values of the first step's ee. The curves follow the velocity-addition law and all bend toward the wall at ee = 1.", js=ee_js),
 dict(id="yields", nav="Yields", h="Yields: losses fold at k = −1, and only the longest branch counts",
  body="""<p>A linear route's losses compose by failure-union, as in subject 5: six steps at 90, 85, 95, 80, 90 and 70% lose 64% overall. Eight 90% steps in a line keep 43%, while two parallel four-step branches joined by a 90% coupling keep 59%, because only the longest branch counts. This is the familiar case for convergent synthesis.</p>"""),
]
claims = [
 ("Hückel and Möbius ring levels; characteristic polynomials 2T_N(x/2) ∓ 2", "known", "N = 3 to 10; Lean <code>benzene_charpoly</code>, <code>cyclobutadiene_charpoly</code>, <code>mobius4_charpoly</code>"),
 ("4n + 2 rule for Hückel rings, 4n for Möbius rings; benzene resonance 2β", "known", "Python; Lean <code>benzene_resonance</code>"),
 ("Hammett ratios multiply, so ratio − 1 composes by the circle operation", "tidy reading", "Python; Lean <code>hammett_compose</code>"),
 ("ee composes by (e₁ + e₂)/(1 + e₁e₂), the landscape point (0, −1) and the velocity-addition law", "tidy reading", "Python; Lean <code>ee_compose</code>, <code>ee_is_landscape</code>, <code>ee_two_walls</code>, <code>ee_wall</code>"),
 ("Horeau duplication is e ⊕ e", "known", "Python; Lean <code>horeau</code>"),
 ("Route losses fold at k = −1; convergent routes keep only the longest branch", "known", "Python; Lean <code>route_losses</code>"),
]
lean = [("benzene_charpoly, cyclobutadiene_charpoly, mobius4_charpoly","ring characteristic polynomials"),("benzene_resonance","8β against 6β"),
 ("hammett_compose","Hammett ratios the circle operation"),("ee_compose, ee_is_landscape","ee composes by the (0, −1) landscape law"),
 ("ee_two_walls, ee_wall","two-walls country, ee = 1 absorbs"),("horeau","Horeau duplication is e ⊕ e"),("route_losses","losses fold at k = −1")]
n = build("index.html", title="Organic chemistry", tile=("12","C","Organic"),
 h1="Rings close Hückel chains into the 4n + 2 rule, Hammett ratios multiply, and enantiomeric excess composes by the velocity-addition law",
 lead="Hückel rings and their Möbius twins give the 4n + 2 and 4n rules, with characteristic polynomials from Chebyshev polynomials. Hammett constants add, so rate ratios multiply. Matched selective steps compose ee by (e₁ + e₂)/(1 + e₁e₂), the landscape point (0, −1), the same law as relativistic velocity addition, and Horeau's duplication is that law applied to itself.",
 figs=[("16 / 16","checks pass"),("11","Lean theorems, compiled"),("(0, −1)","ee law on the landscape"),("0.961","10:1 matched with 5:1")],
 picture=picture, caption="π orbital energies for rings of 3 to 8 atoms: Hückel rings in blue, Möbius rings in red, with degenerate levels side by side. Every Hückel ring has one lowest-energy level at 2β; Möbius levels come in degenerate pairs, apart from a single level at −2β in odd rings.",
 verdict_h="Lead check: not a lead",
 verdict="The ee law is exact and neat, and it puts selectivity in the two-walls country beside the relativistic velocities, but it follows directly from ee being a tanh and ratios multiplying, which chemists already use. The rings and Hammett parts are textbook. Nothing here predicts something new, so I didn't follow it further.",
 sections=sections, claims=claims, lean_file="chemistry/12-organic/lean/ChemOrganic.lean", lean_items=lean,
 files="In the project folder under math-theory/chemistry/12-organic/: <code>checks.py</code> (16 checks), <code>checks_output.txt</code>, <code>results.json</code>, <code>lean/ChemOrganic.lean</code>, <code>make_page.py</code>.",
 sources='Sources: Hückel (1931); Heilbronner, Tetrahedron Lett. 5, 1923 (1964) for Möbius rings; Hammett, J. Am. Chem. Soc. 59, 96 (1937); Masamune et al. (1985) on double asymmetric induction; Horeau, Tetrahedron 31, 1307 (1975) on duplication.',
 series_nav="Chemistry pass, part 12 of 13: organic. Before: acids and bases. Next: biochemistry, with photosynthesis and nitrogenase.")
print("page bytes", n)
