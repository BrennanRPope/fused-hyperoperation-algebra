"""Labeled SVG pictures for the electromagnetism page. Run: python3 figures.py -> figures/*.svg"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
SKY, AMB, PINK, SOFT = "#8fd3ff", "#ffc66b", "#ff5c8a", "#b9cbd9"
FONT = 'font-family="Hanken Grotesk, Helvetica, Arial, DejaVu Sans, Noto Sans Math, STIX Two Math, Cambria Math, sans-serif"'


def T(x, y, s, size=13, anchor="middle", fill="currentColor", weight=400, italic=False):
    st = ' font-style="italic"' if italic else ""
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}" font-weight="{weight}" fill="{fill}" stroke="none"{st}>{s}</text>'


def L(x1, y1, x2, y2, col="currentColor", w=1.5, dash=None, arrow=False, op=1):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = ' marker-end="url(#ar)"' if arrow else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{col}" stroke-width="{w}"{d}{m} stroke-opacity="{op}"/>'


def P(pts, col="currentColor", w=1.8, dash=None, fill="none", arrow=False, close=False, op=1):
    d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y in pts) + (" Z" if close else "")
    da = f' stroke-dasharray="{dash}"' if dash else ""
    m = ' marker-end="url(#ar)"' if arrow else ""
    return f'<path d="{d}" stroke="{col}" stroke-width="{w}" fill="{fill}"{da}{m} stroke-opacity="{op}"/>'


def C(x, y, r, col="currentColor", fill="none", w=1.5, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" stroke="{col}" fill="{fill}" stroke-width="{w}"{da}/>'


def R(x, y, w, h, col="currentColor", fill="none", rx=10, sw=1.5, fop=1, dash=None):
    da = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" stroke="{col}" fill="{fill}" fill-opacity="{fop}" stroke-width="{sw}"{da}/>'


def save(name, w, h, label, body):
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{label}" {FONT} '
           f'font-size="13" fill="none" stroke="currentColor"><defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" '
           f'markerUnits="userSpaceOnUse" markerWidth="10" markerHeight="10" orient="auto-start-reverse">'
           f'<path d="M0 0L10 5L0 10z" fill="context-stroke" stroke="none"/></marker></defs>' + "".join(body) + "</svg>")
    # context-stroke is not universal: fall back to currentColor
    svg = svg.replace('fill="context-stroke"', 'fill="currentColor"')
    with open(os.path.join(OUT, name + ".svg"), "w") as fh:
        fh.write(svg)
    return svg


figs = {}

# ---------------------------------------------------------------- 1. history: EM -> SR -> GR
b = []
W, H = 1040, 330
xs = [150, 520, 890]
names = ["1865 Maxwell", "1905 Einstein (special)", "1915 Einstein (general)"]
subs = [["waves move at", "c = 1/√(ε₀μ₀)"], ["boosts compose by", "The circle operation (Doppler² lift)"], ["the twist rule made local:", "loop turn = curvature × area"]]
for i, (xc, nm, sb) in enumerate(zip(xs, names, subs)):
    b.append(R(xc - 140, 40, 280, 150, SKY, SKY, 16, 1.5, 0.10))
    b.append(T(xc, 70, nm, 17, weight=700))
    b.append(T(xc, 158, sb[0], 13, fill=SOFT))
    b.append(T(xc, 177, sb[1], 14, weight=700))
# icons
wx = [xs[0] - 70 + k * 2 for k in range(71)]
b.append(P([(xx, 115 - 18 * math.sin((xx - wx[0]) / 70 * 2 * math.pi * 1.5)) for xx in wx], SKY, 2.2))
b.append(L(xs[1] - 30, 135, xs[1] + 30, 85, SKY, 2)); b.append(L(xs[1] + 30, 135, xs[1] - 30, 85, SKY, 2))
b.append(L(xs[1], 138, xs[1], 82, "currentColor", 1, "3 3"))
pts = []
for k in range(41):
    a = k / 40 * 2 * math.pi * 0.86
    pts.append((xs[2] + 34 * math.cos(a), 110 + 22 * math.sin(a)))
b.append(P(pts, SKY, 2.2, arrow=True))
b.append(T(xs[2] + 52, 100, "turned", 12, fill=SKY, anchor="start"))
for i, lab in enumerate(["sets the dial q = −1/c²", "made local, plus a metric"]):
    x1, x2 = xs[i] + 145, xs[i + 1] - 145
    b.append(L(x1, 115, x2, 115, "currentColor", 1.6, arrow=True))
    b.append(T((x1 + x2) / 2, 105, lab.split(",")[0], 12, fill=SOFT))
    if "," in lab:
        b.append(T((x1 + x2) / 2, 134, lab.split(",")[1].strip(), 12, fill=SOFT))
# shared band
b.append(R(10, 225, 1020, 82, AMB, AMB, 14, 1.2, 0.08))
b.append(T(520, 252, "shared by electromagnetism and gravity (both are gauge fields)", 14, weight=700, fill=AMB))
b.append(T(275, 284, "EM: phase around a loop = charge × flux", 14))
b.append(T(780, 284, "gravity: boost around a loop = curvature × area", 14))
figs["history"] = save("history", W, H, "Electromagnetism set the dial, special relativity composed boosts by the circle operation, general relativity made the twist local", b)

# ---------------------------------------------------------------- 2. the dial: Galilean vs Lorentz
b = []
W, H = 1040, 420


def panel_axes(cx, cy, s, title, col):
    o = []
    o.append(T(cx, 34, title, 16, weight=700, fill=col))
    o.append(L(cx - s, cy, cx + s, cy, "currentColor", 1, op=.5)); o.append(L(cx, cy + s, cx, cy - s, "currentColor", 1, op=.5))
    o.append(T(cx + s + 4, cy + 4, "x", 13, anchor="start", fill=SOFT)); o.append(T(cx + 6, cy - s - 4, "t", 13, anchor="start", fill=SOFT))
    o.append(L(cx - s * .9, cy + s * .9, cx + s * .9, cy - s * .9, "currentColor", 1.2, "5 4"))
    o.append(L(cx + s * .9, cy + s * .9, cx - s * .9, cy - s * .9, "currentColor", 1.2, "5 4"))
    return o


s = 150
cxL, cyL = 270, 225
b += panel_axes(cxL, cyL, s, "Galilean boost (q = 0)", PINK)
v = 0.45
# Galilean: x' = x - v t. A light ray x = t (slope 1) becomes x' = (1 - v) t: draw new rays in boosted frame coordinates
for sgn, lab in ((1, "speed c − v"), (-1, "speed c + v")):
    slope = sgn - v  # dx'/dt
    tt = s * .9 / max(abs(slope), 1) if abs(slope) > 1 else s * .9
    x2 = cxL + slope * tt
    y2 = cyL - tt
    b.append(L(cxL, cyL, x2, y2, PINK, 2.4))
    b.append(T(x2 + (8 if sgn > 0 else -8), y2 - 6, lab, 13, anchor="start" if sgn > 0 else "end", fill=PINK))
b.append(T(cxL, 400, "light lines tilt unevenly: the wave equation gains a cross term", 13, fill=SOFT))
cxR, cyR = 770, 225
b += panel_axes(cxR, cyR, s, "Lorentz boost (q = −1/c²)", SKY)
Dd = 1.7
# unit square on light-cone axes, then squeezed
u = 58
sq = [(0, 0), (u, 0), (u, u), (0, u)]  # (forward, backward) light-cone coordinates


def lc(cx, cy, f, bk):
    # forward axis along x=t (up-right), backward along x=-t (up-left)
    k = 1 / math.sqrt(2)
    return (cx + (f - bk) * k, cy - (f + bk) * k)


b.append(P([lc(cxR, cyR, f, bk) for f, bk in sq], "currentColor", 1.4, "3 3", close=True))
b.append(P([lc(cxR, cyR, f * Dd, bk / Dd) for f, bk in sq], SKY, 2.4, fill=SKY, close=True, op=1))
b[-1] = b[-1].replace('fill="#8fd3ff"', 'fill="#8fd3ff" fill-opacity=".18"')
fx, fy = lc(cxR, cyR, u * Dd + 26, 0)
b.append(T(fx + 14, fy + 14, "× D", 15, anchor="start", weight=700, fill=SKY))
bx_, by_ = lc(cxR, cyR, 0, u / Dd + 34)
b.append(T(bx_ - 4, by_, "× 1/D", 15, anchor="end", weight=700, fill=SKY))
b.append(T(cxR, 400, "light lines stay put; areas are kept; D² = (1+v)/(1−v) is the circle-operation lift", 13, fill=SOFT))
figs["dial"] = save("dial", W, H, "A Galilean boost tilts the light lines unevenly; a Lorentz boost squeezes along them and keeps them", b)

# ---------------------------------------------------------------- 3. the field on the light cone
b = []
W, H = 1040, 470
cx, cy, sc = 330, 250, 52
b.append(L(cx - 230, cy, cx + 260, cy, "currentColor", 1, op=.5)); b.append(T(cx + 266, cy + 4, "E_y", 14, anchor="start"))
b.append(L(cx, cy + 210, cx, cy - 210, "currentColor", 1, op=.5)); b.append(T(cx + 6, cy - 214, "B_z", 14, anchor="start"))
k = 1 / math.sqrt(2)
b.append(L(cx - 200 * k * 1.2, cy + 200 * k * 1.2, cx + 230 * k * 1.2, cy - 230 * k * 1.2, PINK, 1.4, "6 4"))
b.append(L(cx - 200 * k * 1.2, cy - 200 * k * 1.2, cx + 230 * k * 1.2, cy + 230 * k * 1.2, PINK, 1.4, "6 4"))
b.append(T(cx + 205, cy - 205, "E_y = B_z (light wave)", 13, anchor="start", fill=PINK))
b.append(T(cx + 205, cy + 215, "E_y = −B_z", 13, anchor="start", fill=PINK))
# hyperbola E^2 - B^2 = const (const = 2^2 in units)
a0 = 2.0
hyp = [(cx + sc * a0 * math.cosh(tt), cy - sc * a0 * math.sinh(tt)) for tt in [i / 30 - 1.05 for i in range(64)]]
b.append(P(hyp, SKY, 2.2))
b.append(T(cx + sc * a0 * math.cosh(1.05) + 8, cy + sc * a0 * math.sinh(1.05) + 4, "E² − B² fixed", 13, anchor="start", fill=SKY))
for j, rap in enumerate([0, math.log(1.7), 2 * math.log(1.7) * 0.5 + math.log(1.7) * 0.5]):
    pass
pts = [0.0, -math.log(1.7) / 2, -math.log(1.7)]  # boosts move the point along the hyperbola
labs = ["field", "after one boost", "after two"]
prev = None
for rap, lab in zip(pts, labs):
    px, py = cx + sc * a0 * math.cosh(rap), cy - sc * a0 * math.sinh(rap)
    b.append(C(px, py, 6, SKY, SKY))
    b.append(T(px + 12, py + 5 if rap < 0 else py - 8, lab, 13, anchor="start"))
    if prev:
        pass
    prev = (px, py)
# wave point on the diagonal, moving inward
for d_, lab, op in ((3.2, "wave", 1), (3.2 / 1.7, "÷ D", 1), (3.2 / 1.7 ** 2, "÷ D²", 1)):
    px, py = cx + sc * d_ * k, cy - sc * d_ * k
    b.append(C(px, py, 5, PINK, PINK))
    b.append(T(px - 10, py - 6, lab, 12, anchor="end", fill=PINK))
# side text
tx = 650
rows = [("E_y + B_z", "× 1/D"), ("E_y − B_z", "× D"), ("their ratio", "× D² = 1 + d(v)"), ("their product", "unchanged"),
        ("two boosts", "D₁ D₂ (so d's the circle operation)")]
b.append(T(tx, 70, "One boost along x does this:", 16, anchor="start", weight=700))
for i, (l1, l2) in enumerate(rows):
    yy = 110 + i * 44
    b.append(T(tx, yy, l1, 15, anchor="start"))
    b.append(T(tx + 150, yy, l2, 15, anchor="start", weight=700, fill=SKY))
b.append(T(tx, 350, "A light wave lives on one light line", 14, anchor="start", fill=PINK, weight=700))
b.append(T(tx, 372, "(null field: E² = B² and E·B = 0). Moving with it,", 13, anchor="start", fill=SOFT))
b.append(T(tx, 392, "its amplitude falls by 1/D and its energy", 13, anchor="start", fill=SOFT))
b.append(T(tx, 412, "density by 1/D², the inverse circle-operation lift.", 13, anchor="start", fill=SOFT))
figs["field"] = save("field", W, H, "Under a boost the light-cone parts of the field scale by D and 1/D; the field slides along a hyperbola; a light wave sits on a light line", b)

# ---------------------------------------------------------------- 4. circle and hyperbola
b = []
W, H = 1040, 430
cx, cy, r = 260, 225, 120
b.append(T(cx, 36, "Gauge phase (q = +1, circle)", 16, weight=700, fill=SKY))
b.append(L(cx - 180, cy, cx + 180, cy, "currentColor", 1, op=.5)); b.append(L(cx, cy + 170, cx, cy - 170, "currentColor", 1, op=.5))
b.append(C(cx, cy, r, SKY, w=2.2))
th = 1.0
px, py = cx + r * math.cos(th), cy - r * math.sin(th)
qx, qy = cx + r * math.cos(th), cy + r * math.sin(th)
b.append(L(cx + r, cy, px, py, AMB, 2.2, arrow=True)); b.append(T((cx + r + px) / 2 + 14, (cy + py) / 2, "a", 15, anchor="start", weight=700, fill=AMB))
b.append(L(cx + r, cy, qx, qy, AMB, 1.6, "5 4", arrow=True)); b.append(T((cx + r + qx) / 2 + 14, (cy + qy) / 2 + 6, "ā", 15, anchor="start", weight=700, fill=AMB))
b.append(C(px, py, 5, SKY, SKY)); b.append(T(px + 8, py - 8, "1 + a = e^(iθ)", 13, anchor="start"))
b.append(C(cx + r, cy, 4, "currentColor", "currentColor")); b.append(T(cx + r + 8, cy + 18, "1", 13, anchor="start"))
b.append(C(cx, cy, 5, PINK, PINK)); b.append(T(cx - 8, cy + 20, "wall: 1 + a = 0", 13, anchor="end", fill=PINK))
b.append(T(cx, 380, "a ∘₁ ā = 0 on the whole circle; it wraps around and never meets the wall", 13, fill=SOFT))
cx2 = 780
b.append(T(cx2, 36, "Boost (q = −1, hyperbola)", 16, weight=700, fill=SKY))
b.append(L(cx2 - 180, cy, cx2 + 180, cy, "currentColor", 1, op=.5)); b.append(L(cx2, cy + 170, cx2, cy - 170, "currentColor", 1, op=.5))
b.append(L(cx2 - 160, cy + 160, cx2 + 160, cy - 160, PINK, 1.5, "6 4")); b.append(L(cx2 - 160, cy - 160, cx2 + 160, cy + 160, PINK, 1.5, "6 4"))
b.append(T(cx2 + 150, cy - 165, "wall (light)", 13, anchor="start", fill=PINK))
hyp = [(cx2 + r * math.cosh(tt), cy - r * math.sinh(tt)) for tt in [i / 25 - 1.0 for i in range(51)]]
b.append(P(hyp, SKY, 2.2))
ph = 0.6
px, py = cx2 + r * math.cosh(ph), cy - r * math.sinh(ph)
qx, qy = cx2 + r * math.cosh(ph), cy + r * math.sinh(ph)
b.append(L(cx2 + r, cy, px, py, AMB, 2.2, arrow=True)); b.append(T((cx2 + r + px) / 2 - 12, (cy + py) / 2 - 4, "a", 15, anchor="end", weight=700, fill=AMB))
b.append(L(cx2 + r, cy, qx, qy, AMB, 1.6, "5 4", arrow=True)); b.append(T((cx2 + r + qx) / 2 - 12, (cy + qy) / 2 + 12, "ā", 15, anchor="end", weight=700, fill=AMB))
b.append(C(px, py, 5, SKY, SKY)); b.append(T(px - 10, py - 12, "1 + a = cosh φ + j sinh φ", 13, anchor="end"))
b.append(T(cx2, 380, "same rule a ∘₁ ā = 0 (with j² = +1); the arms run up to the light walls", 13, fill=SOFT))
figs["circle"] = save("circle", W, H, "The U(1) phase lives on a circle, the boost on a hyperbola; both obey a ∘₁ ā = 0, with the mirror sign of the dial", b)

# ---------------------------------------------------------------- 5. holonomy: plaquettes, solenoid, gravity loop
b = []
W, H = 1040, 400
# lattice: two plaquettes
x0, y0, e = 60, 110, 120
b.append(T(x0 + e, 40, "Lattice: plaquettes combine by the circle operation", 15, weight=700, fill=SKY))
def arrow_seg(p, q, col, w=2.0, off=0):
    (x1, y1), (x2, y2) = p, q
    dx, dy = x2 - x1, y2 - y1
    n = math.hypot(dx, dy); ox, oy = -dy / n * off, dx / n * off
    return L(x1 + dx * .12 + ox, y1 + dy * .12 + oy, x2 - dx * .12 + ox, y2 - dy * .12 + oy, col, w, arrow=True)
A_, B_, Cc, D_ = (x0, y0 + e), (x0 + e, y0 + e), (x0 + e, y0), (x0, y0)
E_, F_ = (x0 + 2 * e, y0 + e), (x0 + 2 * e, y0)
for p, q in ((A_, B_), (B_, Cc), (Cc, D_), (D_, A_)):
    b.append(arrow_seg(p, q, SKY, off=6))
for p, q in ((B_, E_), (E_, F_), (F_, Cc), (Cc, B_)):
    b.append(arrow_seg(p, q, AMB, off=6))
b.append(T(x0 + e / 2, y0 + e / 2 + 5, "a₁", 16, weight=700, fill=SKY))
b.append(T(x0 + 1.5 * e, y0 + e / 2 + 5, "a₂", 16, weight=700, fill=AMB))
b.append(T(x0 + e, y0 + e + 34, "shared link runs both ways: cancels", 12, fill=SOFT))
b.append(T(x0 + e, y0 + e + 56, "big loop = a₁ ∘₁ a₂ (Stokes)", 14, weight=700))
b.append(T(x0 + e, y0 + e + 80, "each plaquette = e^(i ε² F) − 1", 12, fill=SOFT))
# solenoid
cx, cy = 520, 200
b.append(T(cx, 40, "Aharonov–Bohm", 15, weight=700, fill=SKY))
loop = [(cx + 95 * math.cos(a) + 18 * math.cos(3 * a), cy + 75 * math.sin(a) + 10 * math.sin(2 * a)) for a in [i / 60 * 2 * math.pi for i in range(61)]]
b.append(P(loop, SKY, 2.2))
b.append(C(cx, cy, 14, AMB, AMB))
b.append(T(cx, cy + 5, "Φ", 15, weight=700, fill="#0a2540"))
b.append(T(cx, cy + 120, "B = 0 on the path, yet", 13, fill=SOFT))
b.append(T(cx, cy + 140, "phase = charge × flux Φ", 14, weight=700))
b.append(T(cx, cy + 160, "any loop around it, any shape", 12, fill=SOFT))
# gravity loop
gx, gy = 860, 200
b.append(T(gx, 40, "Gravity (GR page, route b)", 15, weight=700, fill=SKY))
b.append(R(gx - 85, gy - 70, 170, 140, SKY, rx=4, sw=2.2))
b.append(L(gx - 85, gy + 70, gx - 40, gy + 70 - 60, "currentColor", 2, arrow=True))
b.append(L(gx - 85, gy + 70, gx - 32, gy + 70 - 36, AMB, 2.2, arrow=True))
b.append(T(gx - 20, gy + 40, "returns boosted", 12, anchor="start", fill=AMB))
b.append(T(gx, gy + 120, "loop in the t–r plane:", 13, fill=SOFT))
b.append(T(gx, gy + 140, "rapidity = curvature × area", 14, weight=700))
b.append(T(gx, gy + 160, "(2M/r³) × area", 12, fill=SOFT))
figs["holonomy"] = save("holonomy", W, H, "Plaquettes compose by the circle operation, the Aharonov-Bohm phase is charge times flux, and a gravity loop returns boosted by curvature times area", b)

# ---------------------------------------------------------------- 6. the corner: a b = a.b + a^b ; nabla A = nabla.A + F
b = []
W, H = 1040, 400
ox, oy, u, av, bv = 70, 330, 150, 110, 95
b.append(T(ox + (u + av) / 2, 40, "The lift square (1 + a)(1 + b)", 15, weight=700))
b.append(R(ox, oy - u, u, u, "currentColor", rx=0, sw=1.4)); b.append(T(ox + u / 2, oy - u / 2 + 5, "1", 16))
b.append(R(ox + u, oy - u, av, u, "currentColor", rx=0, sw=1.4)); b.append(T(ox + u + av / 2, oy - u / 2 + 5, "a", 16))
b.append(R(ox, oy - u - bv, u, bv, "currentColor", rx=0, sw=1.4)); b.append(T(ox + u / 2, oy - u - bv / 2 + 5, "b", 16))
b.append(R(ox + u, oy - u - bv, av, bv, SKY, SKY, 0, 2.2, .18))
cxx, cyy = ox + u, oy - u - bv
b.append(L(cxx, cyy + bv, cxx + av, cyy, SKY, 1.6))
b.append(T(cxx + av * .32, cyy + bv * .72, "a·b", 14, weight=700, fill=SKY))
b.append(T(cxx + av * .70, cyy + bv * .36, "a∧b", 14, weight=700, fill=AMB))
b.append(T(ox + (u + av) / 2, oy + 30, "corner ab = a·b (same both orders) + a∧b (flips)", 12, fill=SOFT))
b.append(L(370, 200, 470, 200, "currentColor", 1.6, arrow=True))
b.append(T(420, 188, "a → ∇, b → A", 13, fill=SOFT))
bx0 = 500
b.append(T(bx0 + 250, 40, "The same split for the potential", 15, weight=700))
b.append(T(bx0 + 250, 110, "∇A  =  ∇·A  +  ∇∧A", 30, weight=700))
b.append(R(bx0 + 40, 150, 200, 100, AMB, AMB, 12, 1.4, .10))
b.append(T(bx0 + 140, 185, "∇·A", 18, weight=700, fill=AMB))
b.append(T(bx0 + 140, 210, "a gauge choice", 13, fill=SOFT))
b.append(T(bx0 + 140, 230, "(Lorenz gauge sets it to 0)", 12, fill=SOFT))
b.append(R(bx0 + 280, 150, 200, 100, SKY, SKY, 12, 1.4, .12))
b.append(T(bx0 + 380, 185, "∇∧A = F", 18, weight=700, fill=SKY))
b.append(T(bx0 + 380, 210, "the field: E and B", 13, fill=SOFT))
b.append(T(bx0 + 380, 230, "in one oriented-area object", 12, fill=SOFT))
b.append(T(bx0 + 250, 300, "the order gap of the multivector circle operation is 2 a∧b (Clifford page);", 13, fill=SOFT))
b.append(T(bx0 + 250, 322, "here the wedge part is exactly the field strength", 13, fill=SOFT))
figs["corner"] = save("corner", W, H, "The corner of the lift square splits into a dot part and a wedge part; for the derivative and the potential the wedge part is the field F", b)

# ---------------------------------------------------------------- 7. one equation
b = []
W, H = 1040, 430
b.append(R(400, 40, 240, 70, SKY, SKY, 16, 2, .14))
b.append(T(520, 86, "∇F = J", 30, weight=700))
b.append(L(470, 112, 250, 170, "currentColor", 1.5, arrow=True)); b.append(L(570, 112, 790, 170, "currentColor", 1.5, arrow=True))
b.append(T(285, 132, "vector part", 13, fill=SOFT)); b.append(T(760, 132, "trivector part", 13, fill=SOFT))
b.append(R(60, 175, 380, 110, "currentColor", rx=12, sw=1.2))
b.append(T(250, 212, "∇·E = ρ   (Gauss)", 16)); b.append(T(250, 248, "∇×B − ∂E/∂t = J   (Ampère–Maxwell)", 16))
b.append(R(600, 175, 380, 110, "currentColor", rx=12, sw=1.2))
b.append(T(790, 212, "∇·B = 0   (no monopoles)", 16)); b.append(T(790, 248, "∇×E + ∂B/∂t = 0   (Faraday)", 16))
b.append(R(60, 320, 440, 84, SKY, SKY, 12, 1.4, .10))
b.append(T(280, 352, "∇·J = ⟨∇∇F⟩₀ = 0", 17, weight=700, fill=SKY))
b.append(T(280, 378, "charge conservation comes free", 13, fill=SOFT))
b.append(R(540, 320, 440, 84, SKY, SKY, 12, 1.4, .10))
b.append(T(760, 352, "∇·G = 0 (Bianchi)", 17, weight=700, fill=SKY))
b.append(T(760, 378, "gravity's energy conservation comes free", 13, fill=SOFT))
figs["one"] = save("one", W, H, "Maxwell's four equations are the vector and trivector parts of one equation; charge conservation follows as gravity's energy conservation does", b)

# ---------------------------------------------------------------- 8. commuting vs not
b = []
W, H = 1040, 330


def loop_fig(cx, cy, title, col, closes, note1, note2):
    o = [T(cx, 36, title, 15, weight=700, fill=col)]
    s = 90
    o.append(L(cx - s, cy + s / 2, cx + s, cy + s / 2, "currentColor", 2, arrow=True))
    o.append(L(cx + s, cy + s / 2, cx + s, cy - s / 2 - 10, "currentColor", 2, arrow=True))
    o.append(L(cx + s, cy - s / 2 - 10, cx - s, cy - s / 2 - 10, "currentColor", 2, arrow=True))
    end_y = cy + s / 2 if closes else cy + s / 2 - 34
    o.append(L(cx - s, cy - s / 2 - 10, cx - s, end_y - 4, "currentColor", 2, arrow=True))
    if closes:
        o.append(C(cx - s, cy + s / 2, 7, SKY, SKY))
        o.append(T(cx - s - 12, cy + s / 2 + 5, "closes", 13, anchor="end", fill=SKY))
    else:
        o.append(L(cx - s, end_y, cx - s, cy + s / 2, PINK, 3))
        o.append(T(cx - s - 10, cy + s / 2 - 12, "gap =", 13, anchor="end", fill=PINK))
        o.append(T(cx - s - 10, cy + s / 2 + 6, "[A_x, A_y]", 13, anchor="end", fill=PINK, weight=700))
    o.append(T(cx, cy + s + 30, note1, 13, fill=SOFT)); o.append(T(cx, cy + s + 50, note2, 13, fill=SOFT))
    return o


b += loop_fig(270, 150, "U(1): phases commute (the circle operation)", SKY, True, "constant potential: no field", "F = ∇∧A, field carries no charge: linear")
b += loop_fig(780, 150, "Boosts or Yang–Mills: they don't", PINK, False, "constant potentials still turn the loop", "field gains A∧A, carries charge: nonlinear")
figs["commute"] = save("commute", W, H, "With commuting phases a loop of constant potentials closes; with non-commuting ones it leaves a gap equal to the commutator", b)

# ---------------------------------------------------------------- 9. walls: Born-Infeld vs velocity
b = []
W, H = 1040, 380
ox, oy, sx, sy = 90, 320, 95, 240
b.append(L(ox, oy, ox + 5.2 * sx, oy, "currentColor", 1, op=.6)); b.append(L(ox, oy, ox, oy - 1.15 * sy, "currentColor", 1, op=.6))
b.append(L(ox, oy - sy, ox + 5.2 * sx, oy - sy, PINK, 1.6, "6 4"))
b.append(T(ox + 5.2 * sx, oy - sy - 8, "wall: c (speed) or b (field)", 13, anchor="end", fill=PINK))
cur = [(ox + p * sx, oy - sy * p / math.sqrt(1 + p * p)) for p in [i / 20 for i in range(101)]]
b.append(P(cur, SKY, 2.6))
b.append(T(ox + 2.6 * sx, oy + 30, "momentum p  (or Born–Infeld's D/b)", 13))
b.append(T(ox - 10, oy - sy / 2, "v or E/b", 13, anchor="end"))
b.append(T(ox + 3.6 * sx, oy - 0.83 * sy, "v = p/√(1+p²)  and  E = D/√(1+D²/b²)", 13, fill=SKY))
tx = 640
b.append(T(tx, 70, "Same wall shape, different rule for adding", 15, anchor="start", weight=700))
b.append(T(tx, 104, "Velocities compose by Einstein's rule (the circle operation).", 13, anchor="start", fill=SOFT))
b.append(T(tx, 126, "Born–Infeld fields compose by adding D (Gauss).", 13, anchor="start", fill=SOFT))
b.append(T(tx, 162, "Example, D₁ = 0.8, D₂ = 1.3 (b = 1):", 13, anchor="start"))
for i, (lab, val, col) in enumerate((("E(D₁ + D₂)", 0.9029, SKY), ("E₁ ⊕ E₂", 0.9479, AMB))):
    yy = 190 + i * 40
    b.append(T(tx, yy + 14, lab, 13, anchor="start"))
    b.append(R(tx + 100, yy, 260 * val, 20, col, col, 4, 1, .5))
    b.append(T(tx + 108 + 260 * val, yy + 15, f"{val:.4f}", 13, anchor="start"))
b.append(T(tx, 300, "Maxwell itself has no wall: its fields add plainly.", 13, anchor="start", fill=SOFT))
b.append(T(tx, 322, "QED bends that only near 1.3 × 10¹⁸ V/m (Schwinger).", 13, anchor="start", fill=SOFT))
figs["walls"] = save("walls", W, H, "Born-Infeld's field cap has the same curve as speed against momentum, but its fields do not add by Einstein's rule", b)

print("wrote", ", ".join(sorted(figs)))
