"""Draw one labeled picture pair per Clifford piece: left = Clifford, right = the mult-add analog.

Run: python3 figures.py  -> figures/NN-a.svg (Clifford), figures/NN-b.svg (mult-add), figures/hero.svg
Colours come from CSS variables with fallbacks, so the same SVG reads in the page (light or dark)
and as a standalone file.
"""
import math, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

STYLE = """<style>
.k{stroke:var(--ink,#141414);fill:none;stroke-width:2}
.kt{stroke:var(--ink,#141414);fill:none;stroke-width:1.2}
.m{stroke:var(--ma,#1f5fa8);fill:none;stroke-width:2.4}
.mf{fill:var(--ma,#1f5fa8)}
.kf{fill:var(--ink,#141414)}
.hl{fill:var(--hl,#ffe600)}
.hls{stroke:var(--hl,#ffe600);fill:none;stroke-width:7;stroke-linecap:round}
.pl{fill:var(--plane,#e4e4e0)}
.dash{stroke:var(--mute,#8a8a86);fill:none;stroke-width:1.4;stroke-dasharray:5 4}
.mute{fill:var(--mute,#8a8a86)}
text{font:14px 'Hanken Grotesk','DejaVu Sans',system-ui,sans-serif;fill:var(--ink,#141414)}
.tt{font-weight:600;font-size:15px}
.sm{font-size:13px}
.mt{fill:var(--ma,#1f5fa8)}
</style>"""

W, H = 320, 230


def svg(body, title, w=W, h=H):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{title}">'
            f'{STYLE}<title>{title}</title>{body}</svg>')


def t(x, y, s, cls="", anchor="start"):
    c = f' class="{cls}"' if cls else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}"{c}>{s}</text>'


def ln(x1, y1, x2, y2, cls="k"):
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="{cls}"/>'


def circ(cx, cy, r, cls="k"):
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" class="{cls}"/>'


def path(d, cls="k"):
    return f'<path d="{d}" class="{cls}"/>'


def poly(pts, cls="k"):
    return '<polyline points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f'" class="{cls}"/>'


def arrow(x1, y1, x2, y2, cls="k", head=8):
    a = math.atan2(y2 - y1, x2 - x1)
    hx1, hy1 = x2 - head * math.cos(a - 0.4), y2 - head * math.sin(a - 0.4)
    hx2, hy2 = x2 - head * math.cos(a + 0.4), y2 - head * math.sin(a + 0.4)
    fill = "mf" if cls.startswith("m") else "kf"
    return ln(x1, y1, x2, y2, cls) + f'<polygon points="{x2:.1f},{y2:.1f} {hx1:.1f},{hy1:.1f} {hx2:.1f},{hy2:.1f}" class="{fill}"/>'


def dot(x, y, r=4.5, cls="kf"):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" class="{cls}"/>'


FIGS = []  # (id, piece, a_title, a_body, b_title, b_body)


def add(n, piece, at, ab, bt, bb):
    FIGS.append((n, piece, at, ab, bt, bb))


# ------------------------------------------------------------------ 01 geometric product
a = (t(14, 24, "Two vectors make a bivector", "tt")
     + f'<polygon points="70,170 190,170 240,110 120,110" class="pl"/>'
     + arrow(70, 170, 190, 170) + arrow(70, 170, 120, 110)
     + t(196, 186, "u") + t(98, 112, "v")
     + t(130, 148, "u∧v", "sm")
     + t(14, 212, "uv = u·v + u∧v,  u² = 1,  v² = q", "sm"))
b = (t(14, 24, "Two lifts multiply to the law", "tt")
     + f'<rect x="20" y="58" width="120" height="44" rx="6" class="pl"/>' + t(80, 86, "1 + aθ", "", "middle")
     + f'<rect x="180" y="58" width="120" height="44" rx="6" class="pl"/>' + t(240, 86, "1 + bθ", "", "middle")
     + t(160, 86, "×", "tt", "middle")
     + arrow(160, 108, 160, 140, "m")
     + f'<rect x="20" y="146" width="280" height="44" rx="6" class="hl"/>'
     + t(160, 174, "(1 − qab)(1 + F(a,b)·θ)", "", "middle")
     + t(14, 214, "θ = uv,  θ² = pθ − q", "sm mt"))
add("01", "Geometric product, Cl(p,q)", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 02 three 2D algebras
def three(cx, cy, s, which):
    out = ""
    if which == "c":
        out += circ(cx, cy, s, "k")
    elif which == "h":
        pts1 = [(cx + s * math.cosh(u), cy - s * math.sinh(u)) for u in [i / 20 for i in range(-24, 25)]]
        pts2 = [(cx - s * math.cosh(u), cy - s * math.sinh(u)) for u in [i / 20 for i in range(-24, 25)]]
        out += poly(pts1) + poly(pts2) + ln(cx - 1.4 * s, cy - 1.4 * s, cx + 1.4 * s, cy + 1.4 * s, "dash") + ln(cx - 1.4 * s, cy + 1.4 * s, cx + 1.4 * s, cy - 1.4 * s, "dash")
    else:
        out += ln(cx - s * 0.9, cy - 1.3 * s, cx - s * 0.9, cy + 1.3 * s) + ln(cx + s * 0.9, cy - 1.3 * s, cx + s * 0.9, cy + 1.3 * s)
    return out
a = (t(14, 24, "The three 2D Clifford algebras", "tt")
     + three(60, 112, 30, "c") + three(160, 112, 22, "h") + three(262, 112, 30, "d")
     + t(60, 182, "e² = −1", "sm", "middle") + t(160, 182, "e² = +1", "sm", "middle") + t(262, 182, "e² = 0", "sm", "middle")
     + t(60, 204, "complex", "sm", "middle") + t(160, 204, "split", "sm", "middle") + t(262, 204, "dual", "sm", "middle"))
b = (t(14, 24, "The landscape's three countries", "tt")
     + three(60, 112, 30, "c").replace('class="k"', 'class="m"') + three(160, 112, 22, "h").replace('class="k"', 'class="m"')
     + three(262, 112, 30, "d").replace('class="k"', 'class="m"')
     + dot(160 + 22 * math.cosh(0.6), 112 - 22 * math.sinh(0.6), 6, "hl")
     + t(60, 182, "circle", "sm", "middle") + t(160, 182, "line", "sm", "middle") + t(262, 182, "border", "sm", "middle")
     + t(60, 204, "tan law", "sm", "middle") + t(160, 204, "mult-add here", "sm mt", "middle") + t(262, 204, "addition", "sm", "middle"))
add("02", "The three 2D algebras", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 03 k as a multivector
a = (t(14, 24, "A multivector K = k₀ + k₁e", "tt")
     + arrow(40, 190, 290, 190) + arrow(60, 205, 60, 40)
     + t(282, 210, "k₀", "sm") + t(66, 48, "k₁", "sm")
     + dot(200, 90) + ln(60, 190, 200, 90, "dash")
     + t(208, 88, "K", "") + t(14, 222, "every K obeys K² = tr(K)·K − N(K)", "sm"))
b = (t(14, 24, "Its trace and norm are (p, q)", "tt")
     + arrow(40, 190, 290, 190, "kt") + arrow(60, 205, 60, 40, "kt")
     + ln(60, 190, 210, 40, "hls") + ln(60, 190, 210, 40, "m")
     + t(214, 50, "light cone, N = 0", "sm mt")
     + dot(150, 100, 6, "mf") + t(160, 112, "mult-add's k", "sm mt")
     + dot(250, 150, 5) + t(234, 140, "q > 0 off it", "sm", "end")
     + t(14, 222, "p = tr K,   q = N K,   e² = +1", "sm"))
add("03", "k as a multivector (your question)", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 04 grades / blades
def cube(cx, cy, s, cls, labels, hl=None):
    # subset lattice of {1,2,3} drawn as a cube in oblique projection
    P = {}
    for m in range(8):
        x = (m & 1) * s + ((m >> 2) & 1) * s * 0.55
        y = -((m >> 1) & 1) * s - ((m >> 2) & 1) * s * 0.45
        P[m] = (cx + x, cy + y)
    out = ""
    for m in range(8):
        for bit in (1, 2, 4):
            if not m & bit:
                out += ln(*P[m], *P[m | bit], cls)
    for m in range(8):
        out += dot(*P[m], 4, "hl" if hl and m in hl else "kf")
        out += t(P[m][0] + 7, P[m][1] + 15 if m in (1, 5) else P[m][1] - 8, labels[m], "sm")
    return out
lab = ["1", "e₁", "e₂", "e₁₂", "e₃", "e₁₃", "e₂₃", "e₁₂₃"]
a = (t(14, 24, "Grades: 2ⁿ blades", "tt") + cube(70, 180, 90, "kt", lab)
     + t(14, 222, "each subset of directions is its own part", "sm"))
lab2 = ["", "a₁", "a₂", "a₁a₂", "a₃", "a₁a₃", "a₂a₃", "a₁a₂a₃"]
b = (t(14, 24, "a₁e₁ ✳ a₂e₂ ✳ a₃e₃ keeps all 7", "tt") + cube(70, 180, 90, "m", lab2, hl=set(range(1, 8)))
     + t(14, 222, "add the parts up: ordinary a₁ ✳ a₂ ✳ a₃", "sm mt"))
add("04", "Grades and multivectors", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 05 involutions
a = (t(14, 24, "Conjugation and reversion", "tt")
     + ln(160, 40, 160, 200, "dash")
     + dot(110, 100) + t(80, 90, "x + yθ", "sm") + dot(210, 140) + t(212, 160, "x + py − yθ", "sm")
     + arrow(118, 104, 200, 136, "kt")
     + t(14, 222, "X·X̄ = norm,  (XY)~ = Ỹ X̃", "sm"))
b = (t(14, 24, "In the chart: the inverse", "tt")
     + arrow(30, 130, 300, 130, "kt") + dot(80, 130, 6, "hl") + t(80, 155, "−1 wall", "sm", "middle")
     + dot(230, 130, 5, "mf") + t(230, 118, "a", "sm mt", "middle")
     + dot(140, 130, 5, "mf") + t(140, 118, "−a/(1+a)", "sm mt", "middle")
     + path("M230,140 Q185,185 140,140", "m")
     + t(14, 214, "conjugate of a lift = (1+a)·(lift of a⁻¹)", "sm"))
add("05", "Reversion and involutions", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 06 rotors
cx, cy, r = 110, 120, 70
a = (t(14, 24, "A rotor holds half the angle", "tt") + circ(cx, cy, r, "kt")
     + arrow(cx, cy, cx + r, cy) + arrow(cx, cy, cx + r * math.cos(-1.2), cy + r * math.sin(-1.2))
     + path(f"M{cx + 30},{cy} A30,30 0 0 0 {cx + 30 * math.cos(-0.6):.1f},{cy + 30 * math.sin(-0.6):.1f}", "hls")
     + t(cx + 36, cy - 8, "φ/2", "sm") + t(cx + r * math.cos(-1.2) + 4, cy + r * math.sin(-1.2) - 6, "turn φ", "sm")
     + t(200, 100, "R v R̃", "") + t(200, 124, "R and −R:", "sm") + t(200, 142, "same turn", "sm")
     + t(14, 222, "R = cos(φ/2) + sin(φ/2) e₁₂", "sm"))
b = (t(14, 24, "Half-velocities are rotor charts", "tt")
     + arrow(20, 120, 300, 120, "kt")
     + dot(110, 120, 5, "mf") + t(110, 108, "t₁", "sm mt", "middle") + dot(170, 120, 5, "mf") + t(170, 108, "t₂", "sm mt", "middle")
     + dot(250, 120, 6, "hl") + t(250, 108, "t₁ ⊕ t₂", "sm", "middle")
     + t(160, 160, "t = tan(φ/2):  circle law", "sm", "middle") + t(160, 182, "t = tanh(φ/2):  velocity law", "sm", "middle")
     + t(14, 216, "the chart can't see R vs −R; the lift can", "sm mt"))
add("06", "Rotors and spinors", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 07 quaternions
a = (t(14, 24, "Hamilton's quaternions", "tt") + circ(110, 122, 72, "kt")
     + path("M38,122 A72,24 0 0 0 182,122", "kt") + path("M38,122 A72,24 0 0 1 182,122", "dash")
     + t(200, 80, "i² = j² = k²", "sm") + t(200, 100, "  = ijk = −1", "sm")
     + t(14, 222, "unit quaternions turn 3D space", "sm"))
b = (t(14, 24, "The q = +1 dial boost", "tt")
     + f'<rect x="30" y="54" width="120" height="70" rx="6" class="pl"/>'
     + t(90, 82, "[ 1    u ]", "", "middle") + t(90, 106, "[ −ū   1 ]", "", "middle")
     + t(170, 92, "= 1 + u j", "mt")
     + t(20, 160, "product's corner:", "sm") + t(20, 182, "(1 − u·v)  +  i (u×v)", "")
     + t(14, 214, "scalar = dot, i-part = twist", "sm mt"))
add("07", "Hamilton's quaternions", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 08 Grassmann
a = (t(14, 24, "Grassmann's wedge: oriented area", "tt")
     + f'<polygon points="60,180 190,180 250,100 120,100" class="pl"/>'
     + arrow(60, 180, 190, 180) + arrow(60, 180, 120, 100) + t(196, 196, "x") + t(96, 100, "y")
     + path("M140,150 A20,20 0 1 1 160,140", "kt") + t(170, 146, "x∧y", "sm")
     + t(14, 222, "x∧y = −y∧x", "sm"))
b = (t(14, 24, "Order matters by twice the area", "tt")
     + t(20, 70, "x ✳ y = x + y + x·y + x∧y", "")
     + t(20, 100, "y ✳ x = x + y + x·y − x∧y", "")
     + f'<rect x="16" y="120" width="290" height="40" rx="6" class="hl"/>'
     + t(26, 146, "x ✳ y − y ✳ x = 2 x∧y", "tt")
     + t(14, 200, "the same corner split as twisted boosts", "sm mt"))
add("08", "Grassmann's exterior product", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 09 biquaternions
def space(cx, cy, kind, cls):
    if kind == "e":
        return circ(cx, cy, 34, cls) + path(f"M{cx - 34},{cy} A34,11 0 0 0 {cx + 34},{cy}", "kt")
    if kind == "p":
        return f'<rect x="{cx - 34}" y="{cy - 30}" width="68" height="60" class="{cls}"/>' + ln(cx - 34, cy, cx + 34, cy, "kt")
    pts = [(cx + 34 * math.sin(u), cy - 30 + 60 * (u + 1.2) / 2.4 - 18 * (1 - math.cos(u))) for u in [i / 10 for i in range(-12, 13)]]
    return poly(pts, cls) + path(f"M{cx - 30},{cy + 26} Q{cx},{cy - 6} {cx + 30},{cy + 26}", "kt")
a = (t(14, 24, "Biquaternions q + ωr (1873)", "tt")
     + space(55, 110, "e", "k") + space(160, 110, "p", "k") + space(265, 110, "h", "k")
     + t(55, 168, "ω² = +1", "sm", "middle") + t(160, 168, "ω² = 0", "sm", "middle") + t(265, 168, "ω² = −1", "sm", "middle")
     + t(55, 188, "elliptic", "sm", "middle") + t(160, 188, "Euclidean", "sm", "middle") + t(265, 188, "hyperbolic", "sm", "middle")
     + t(14, 216, "ω² is the sign of space's curvature", "sm"))
b = (t(14, 24, "θ is Clifford's (1 + ω)/2", "tt")
     + circ(80, 112, 40, "m") + circ(80, 112, 40, "kt")
     + path("M40,112 L120,112", "dash")
     + t(80, 102, "left", "sm mt", "middle") + t(80, 132, "right", "sm", "middle")
     + t(150, 90, "elliptic motion", "sm") + t(150, 110, "= left half × right half", "sm")
     + t(150, 140, "θ² = θ: the left", "sm mt") + t(150, 160, "projector", "sm mt")
     + t(14, 210, "the coefficient ring is our theta-lift, q = −ω²", "sm"))
add("09", "Biquaternions and Clifford's ω", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 10 dual numbers
xs = [20 + i * 2.8 for i in range(101)]
curve = [(x, 180 - 0.0025 * (x - 40) ** 2) for x in xs]
a = (t(14, 24, "Dual numbers: ε² = 0", "tt") + poly(curve, "kt")
     + ln(80, 176, 250, 108) + dot(160, 144)
     + t(176, 140, "a + a′ε", "sm")
     + t(14, 222, "the ε part rides along as a tangent", "sm"))
b = (t(14, 24, "Mult-add carries its derivative", "tt")
     + t(20, 76, "(a₀ + a₁ε) ✳ (b₀ + b₁ε) =", "")
     + t(20, 106, "a₀ ✳ b₀", "mt")
     + t(20, 136, "+ (a₁(1 + b₀) + b₁(1 + a₀)) ε", "")
     + t(14, 210, "border country = re-charted dual numbers", "sm mt"))
add("10", "Dual numbers", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 11 Clifford parallels / Hopf
a = (t(14, 24, "Clifford parallels in S³", "tt")
     + "".join(path(f"M{60 + 20 * i},{60 + 6 * i} C{120 + 20 * i},{30 + 6 * i} {200 + 20 * i},{90 + 6 * i} {150 + 20 * i},{170 + 6 * i}"
                    + f" C{110 + 20 * i},{200 + 6 * i} {40 + 20 * i},{150 + 6 * i} {60 + 20 * i},{60 + 6 * i}", "kt") for i in range(3))
     + t(14, 222, "great circles, equally spaced everywhere", "sm"))
cx, cy = 110, 125
tri = [(cx - 40, cy + 30), (cx + 50, cy + 20), (cx + 5, cy - 45)]
b = (t(14, 24, "The twist's phase is a Hopf phase", "tt") + circ(cx, cy, 78, "kt")
     + f'<polygon points="{" ".join(f"{x},{y}" for x, y in tri)}" class="hl"/>'
     + poly(tri + [tri[0]], "m")
     + t(cx - 44, cy + 48, "rest", "sm") + t(cx + 52, cy + 38, "u", "sm") + t(cx + 10, cy - 50, "u⊕v", "sm")
     + t(200, 100, "arg(X + iY)", "sm") + t(200, 120, "= area / 2", "sm mt") + t(200, 140, "= Hopf holonomy", "sm")
     + t(14, 222, "checked on 300 triangles, error 3e−16", "sm"))
add("11", "Clifford parallels and the Hopf fibration", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 12 Clifford torus
a = (t(14, 24, "The flat Clifford torus", "tt")
     + f'<ellipse cx="160" cy="120" rx="110" ry="50" class="k"/>' + f'<ellipse cx="160" cy="116" rx="42" ry="14" class="k"/>'
     + t(14, 212, "|z₁| = |z₂| in S³: two circles' product", "sm"))
b = (t(14, 24, "One circle, two meanings", "tt")
     + circ(80, 120, 55, "m") + t(80, 194, "q = −1: light wall", "sm mt", "middle")
     + circ(240, 120, 55, "kt") + f'<ellipse cx="240" cy="120" rx="55" ry="15" class="m"/>'
     + t(240, 194, "q = +1: torus's base", "sm", "middle")
     + t(160, 124, "|t| = 1", "sm", "middle"))
add("12", "The Clifford torus", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 13 Clifford–Klein
a = (t(14, 24, "Clifford–Klein: fold a space", "tt")
     + circ(100, 120, 62, "kt") + dot(100 + 62 * math.cos(0.8), 120 - 62 * math.sin(0.8)) + dot(100 - 62 * math.cos(0.8), 120 + 62 * math.sin(0.8))
     + ln(100 + 62 * math.cos(0.8), 120 - 62 * math.sin(0.8), 100 - 62 * math.cos(0.8), 120 + 62 * math.sin(0.8), "dash")
     + t(180, 110, "S³ / ±1", "") + t(180, 132, "= elliptic space", "sm")
     + t(14, 222, "glue opposite points: R ~ −R", "sm"))
pts = [(110 + 70 * math.cos(2 * math.pi * j / 12), 120 + 70 * math.sin(2 * math.pi * j / 12)) for j in range(12)]
b = (t(14, 24, "The chart is that fold", "tt") + circ(110, 120, 70, "m")
     + "".join(dot(x, y, 4.5, "kf") for x, y in pts) + dot(*pts[0], 6, "hl")
     + t(196, 110, "tan(πj/12)", "sm") + t(196, 130, "a clock of 12", "sm mt")
     + t(14, 222, "partial: 1-D forms only so far", "sm"))
add("13", "Clifford–Klein forms, elliptic space", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 14 space theory of matter
hills = "M20,170 C60,170 70,120 100,120 C130,120 140,170 180,170 C210,170 220,140 240,140 C260,140 270,170 300,170"
a = (t(14, 24, "Space theory of matter (1870)", "tt") + path(hills, "k")
     + t(100, 108, "curved", "sm", "middle") + arrow(140, 70, 220, 70) + t(150, 60, "the hill moves = motion", "sm")
     + t(14, 212, "matter is space bending, passing on", "sm"))
b = (t(14, 24, "Light bending, from the GR page", "tt")
     + f'<rect x="20" y="70" width="140" height="30" class="pl"/>' + t(26, 90, "time dial: half", "sm")
     + f'<rect x="160" y="70" width="140" height="30" class="hl"/>' + t(166, 90, "space curve: half", "sm")
     + t(20, 128, "Clifford's space curvature is the", "sm") + t(20, 146, "half a dial alone misses", "sm mt")
     + t(20, 180, "waves: where the river picture stops", "sm"))
add("14", "The space theory of matter", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 15 Bott
vals = [1, 1.414, 1, 0, -1, -1.414, -1, 0]
names = ["ℝ", "ℝ⊕ℝ", "ℝ", "ℂ", "ℍ", "ℍ⊕ℍ", "ℍ", "ℂ"]
cx, cy = 110, 112
a = t(14, 24, "Bott: Cl(p,q) repeats every 8", "tt") + circ(cx, cy, 58, "kt")
for r in range(8):
    ang = -math.pi / 2 + 2 * math.pi * r / 8
    x, y = cx + 58 * math.cos(ang), cy + 58 * math.sin(ang)
    a += dot(x, y, 5, "kf") + t(cx + 80 * math.cos(ang), cy + 80 * math.sin(ang) + 5, names[r], "sm", "middle")
a += t(cx, cy + 5, "p − q mod 8", "sm", "middle") + t(14, 222, "matrix type by (p − q) mod 8", "sm")
b = (t(14, 24, "The wall stays a point only in", "tt")
     + "".join(f'<rect x="{30 + 92 * i}" y="70" width="76" height="60" rx="8" class="hl"/>' + t(68 + 92 * i, 108, s, "tt", "middle")
               for i, s in enumerate(["ℝ", "ℂ", "ℍ"]))
     + t(20, 160, "Cl(0,0), Cl(0,1), Cl(0,2)", "sm mt")
     + t(20, 184, "everywhere else a blade squares to +1", "sm")
     + t(20, 202, "and the wall grows into a cone", "sm"))
add("15", "Bott periodicity", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 16 chain of circles
a = (t(14, 24, "Chain of circles (1871)", "tt")
     + "".join(circ(140 + 34 * math.cos(2 * math.pi * i / 4 + 0.3), 118 + 34 * math.sin(2 * math.pi * i / 4 + 0.3), 52, "kt") for i in range(4))
     + dot(140, 118, 6, "hl") + dot(140, 118, 2.5)
     + t(14, 216, "4 lines: 4 circles through one point", "sm"))
b = (t(14, 24, "Survives every landscape move", "tt")
     + "".join(circ(130 + 30 * math.cos(2 * math.pi * i / 4 + 0.9), 112 + 26 * math.sin(2 * math.pi * i / 4 + 0.9), 40 + 6 * i, "m") for i in range(4))
     + t(14, 200, "Möbius maps keep it, so it holds in every chart;", "sm")
     + t(14, 218, "no ✳ statement of it yet", "sm mt"))
add("16", "Clifford's chain of circles", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ 17 special divisors
a = (t(14, 24, "Special divisors (1878)", "tt")
     + f'<ellipse cx="160" cy="118" rx="120" ry="56" class="k"/>'
     + path("M90,118 Q110,104 130,118 Q110,128 90,118", "k") + path("M190,118 Q210,104 230,118 Q210,128 190,118", "k")
     + t(14, 210, "on a curve: dim |D| ≤ deg D / 2", "sm"))
b = (t(14, 24, "No analog found", "tt")
     + f'<rect x="20" y="60" width="280" height="90" rx="8" class="pl"/>'
     + t(30, 92, "would need ✳ acting on divisors", "sm") + t(30, 116, "(sums of points on a curve)", "sm")
     + t(30, 140, "no route seen", "sm mt"))
add("17", "Clifford's theorem on special divisors", "Clifford", a, "mult-add", b)

# ------------------------------------------------------------------ hero: K plane with three regions and the cone
HW, HH = 640, 300
hero = t(16, 26, "Where the dial K lives", "tt")
# axes: k0 horizontal, kappa*k1^2 sign as regions: draw three strips
for i, (lbl, sub, cls) in enumerate([("e² = −1", "circle: rotations", "pl"), ("e² = 0", "border: addition", "pl"), ("e² = +1", "line: boosts", "pl")]):
    x0 = 20 + i * 205
    hero += f'<rect x="{x0}" y="44" width="190" height="220" rx="10" class="{cls}"/>'
    hero += t(x0 + 95, 70, lbl, "tt", "middle") + t(x0 + 95, 252, sub, "sm", "middle")
    cx, cy = x0 + 95, 160
    hero += arrow(cx - 80, cy, cx + 80, cy, "kt") + arrow(cx, cy + 70, cx, cy - 70, "kt")
    hero += t(cx + 72, cy + 18, "k₀", "sm") + t(cx + 6, cy - 62, "k₁", "sm")
    if i == 0:
        hero += circ(cx, cy, 45, "k")
    elif i == 1:
        hero += ln(cx - 45, cy - 60, cx - 45, cy + 60, "k") + ln(cx + 45, cy - 60, cx + 45, cy + 60, "k")
    else:
        hero += ln(cx - 62, cy + 62, cx + 62, cy - 62, "hls") + ln(cx - 62, cy + 62, cx + 62, cy - 62, "m") + ln(cx - 62, cy - 62, cx + 62, cy + 62, "dash")
        hero += dot(cx + 34, cy - 34, 7, "mf") + t(cx + 22, cy - 50, "mult-add", "sm mt", "end")
open(os.path.join(OUT, "hero.svg"), "w").write(svg(hero, "Where the dial K lives: three kinds of direction, mult-add on the light cone", HW, HH))

index = []
for n, piece, at, ab, bt, bb in FIGS:
    open(os.path.join(OUT, f"{n}-a.svg"), "w").write(svg(ab, f"{piece}: Clifford"))
    open(os.path.join(OUT, f"{n}-b.svg"), "w").write(svg(bb, f"{piece}: mult-add analog"))
    index.append({"id": n, "piece": piece})
json.dump(index, open(os.path.join(OUT, "index.json"), "w"), indent=1)
print(f"wrote {2 * len(FIGS) + 1} SVG files to {OUT}")
