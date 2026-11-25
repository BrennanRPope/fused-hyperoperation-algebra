"""Writes the SVG pictures for the page (figures/*.svg). Colours: currentColor ink, sky = the circle operation supplies it,
pink = needs something from outside or fails."""
import math, os
OUT = "/mnt/project-files/math-theory/gr-bridges/figures"
SKY, PINK = "#8fd3ff", "#ff5c8a"
HEAD = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{label}" '
        'font-family="Hanken Grotesk, Helvetica, Arial, DejaVu Sans, sans-serif" font-size="13" fill="none" stroke="currentColor">'
        '<defs><marker id="ar{id}" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" markerHeight="11" orient="auto-start-reverse">'
        '<path d="M0 0L10 5L0 10z" fill="currentColor" stroke="none"/></marker></defs>')
def T(x, y, s, size=13, anchor="middle", fill="currentColor", weight=400):
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" font-weight="{weight}" fill="{fill}" stroke="none">{s}</text>'
def save(name, w, h, label, body):
    svg = HEAD.format(w=w, h=h, label=label, id=name) + body.replace("#ar)", f"#ar{name})") + "</svg>"
    open(os.path.join(OUT, name + ".svg"), "w").write(svg)
    return svg

# ---------------------------------------------------------------- 1. the ladder
b = []
W, H = 1040, 330
xs = [20, 280, 540, 800]
titles = ["1. One direction", "2. Two directions", "3. Dial by place", "4. Dial by direction"]
for x0, tt in zip(xs, titles):
    b.append(f'<rect x="{x0}" y="10" width="220" height="300" rx="14" stroke-opacity=".35"/>')
    b.append(T(x0 + 110, 38, tt, 15, weight=700))
# panel 1: line, wall at both ends
x0 = 20
b.append(f'<line x1="{x0+20}" y1="140" x2="{x0+200}" y2="140" stroke-width="2"/>')
for xx in (x0 + 20, x0 + 200):
    b.append(f'<line x1="{xx}" y1="118" x2="{xx}" y2="162" stroke="{PINK}" stroke-width="3"/>')
b.append(T(x0 + 20, 180, "−c", 12, fill=PINK)); b.append(T(x0 + 200, 180, "+c", 12, fill=PINK))
b.append(f'<path d="M{x0+110} 140 Q{x0+130} 100 {x0+150} 140" marker-end="url(#ar)"/>')
b.append(f'<path d="M{x0+150} 140 Q{x0+163} 110 {x0+176} 140" marker-end="url(#ar)"/>')
b.append(f'<circle cx="{x0+110}" cy="140" r="4" fill="currentColor"/>')
b.append(T(x0 + 110, 215, "speeds combine by the circle operation", 13))
b.append(T(x0 + 110, 235, "the wall −1 is light speed", 13))
b.append(T(x0 + 110, 285, "no curvature yet", 13, fill=SKY, weight=700))
# panel 2: curved triangle
x0 = 280
b.append(f'<path d="M{x0+40} 180 Q{x0+110} 160 {x0+180} 180 Q{x0+135} 130 {x0+110} 70 Q{x0+90} 130 {x0+40} 180Z" fill="{SKY}" fill-opacity=".18" stroke="{SKY}" stroke-width="2"/>')
b.append(f'<path d="M{x0+150} 70 A22 22 0 1 1 {x0+128} 92" stroke="{SKY}" stroke-width="2" marker-end="url(#ar)"/>')
b.append(T(x0 + 186, 72, "twist", 12, fill=SKY))
b.append(T(x0 + 110, 215, "a loop of boosts twists", 13))
b.append(T(x0 + 110, 235, "twist = q × area", 13))
b.append(T(x0 + 110, 285, "one fixed curvature", 13, fill=SKY, weight=700))
# panel 3: 2D sheet (ok) and cube (blocked)
x0 = 540
pts = " ".join(f"{x0+20+i*3.6},{92 - 16*math.sin(i/4) * math.exp(-((i-14)/9)**2)}" for i in range(42))
b.append(f'<polyline points="{pts}" stroke="{SKY}" stroke-width="2"/>')
b.append(T(x0 + 20, 128, "2D: any q(x, y) works", 12, anchor="start", fill=SKY))
cx, cy, s = x0 + 140, 205, 40
b.append(f'<path d="M{cx} {cy} h{s} v{-s} h{-s}z M{cx} {cy-s} l16 -12 h{s} l-16 12 M{cx+s} {cy} l16 -12 v{-s}" stroke-width="1.5"/>')
b.append(f'<path d="M{cx} {cy} h{s} v{-s} h{-s}z" fill="{PINK}" fill-opacity=".22" stroke="none"/>')
b.append(f'<path d="M{cx} {cy-s} l16 -12 h{s} l-16 12z" fill="{SKY}" fill-opacity=".22" stroke="none"/>')
b.append(T(x0 + 20, 180, "3D, 4D: each plane", 12, anchor="start", fill=PINK))
b.append(T(x0 + 20, 196, "wants its own q", 12, anchor="start", fill=PINK))
b.append(T(x0 + 110, 245, "Schur: a same-every-way", 13))
b.append(T(x0 + 110, 263, "dial cannot vary", 13))
b.append(T(x0 + 110, 285, "blocked in 3D and 4D", 13, fill=PINK, weight=700))
# panel 4: ellipses (direction dial) on a grid
x0 = 800
for i, (ex, ey, rx, ry, rot) in enumerate([(50, 90, 22, 10, 10), (110, 80, 16, 14, 50), (170, 95, 24, 8, -30),
                                          (60, 160, 12, 20, 0), (120, 150, 26, 12, 70), (180, 165, 14, 14, 0)]):
    b.append(f'<ellipse cx="{x0+ex}" cy="{ey}" rx="{rx}" ry="{ry}" transform="rotate({rot} {x0+ex} {ey})" stroke="{SKY}" stroke-width="2"/>')
b.append(T(x0 + 110, 215, "a dial for each direction", 13))
b.append(T(x0 + 110, 235, "pair: k_ij = the metric", 13))
b.append(T(x0 + 110, 285, "any curvature", 13, fill=SKY, weight=700))
for x0 in xs[:-1]:
    b.append(f'<line x1="{x0+224}" y1="160" x2="{x0+256}" y2="160" stroke-width="2" marker-end="url(#ar)"/>')
fig1 = save("ladder", W, H, "Four steps: the circle operation gives one fixed curvature; a place dial works in 2D but is blocked in 3D and 4D; a direction dial (a metric) gives any curvature", "".join(b))

# ---------------------------------------------------------------- 2. bundle area
b = []
W, H = 620, 330
ox, oy, U = 60, 270, 150
a_, b_ = 0.55, 0.4
b.append(f'<rect x="{ox}" y="{oy-U}" width="{U}" height="{U}" fill="currentColor" fill-opacity=".10" stroke-width="2"/>')
b.append(f'<rect x="{ox+U}" y="{oy-U}" width="{U*a_}" height="{U}" fill="currentColor" fill-opacity=".05"/>')
b.append(f'<rect x="{ox}" y="{oy-U-U*b_}" width="{U}" height="{U*b_}" fill="currentColor" fill-opacity=".05"/>')
b.append(f'<rect x="{ox+U}" y="{oy-U-U*b_}" width="{U*a_}" height="{U*b_}" fill="{SKY}" fill-opacity=".55" stroke="{SKY}" stroke-width="2"/>')
b.append(T(ox + U / 2, oy - U / 2 + 5, "1", 18, weight=700))
b.append(T(ox + U + U * a_ / 2, oy - U / 2 + 5, "a", 16))
b.append(T(ox + U / 2, oy - U - U * b_ / 2 + 5, "b", 16))
b.append(T(ox + U + U * a_ / 2, oy - U - U * b_ / 2 + 5, "ab", 16, weight=700))
b.append(T(ox + (U + U * a_) / 2, oy + 24, "stretch 1 + a across", 13))
b.append(f'<text x="{ox-14}" y="{oy-(U+U*b_)/2}" transform="rotate(-90 {ox-14} {oy-(U+U*b_)/2})" text-anchor="middle" fill="currentColor" stroke="none">stretch 1 + b</text>')
tx = 345
for i, line in enumerate(["A thin beam of light:", "area = (1 + a)(1 + b)", "       = 1 + (a ✳ b)", "",
                          "The ab corner is what makes", "light beams focus", "(Raychaudhuri's θ²/2 + σ² term).",
                          "", "Curvature bends a and b;", "matter's energy does the rest."]):
    col = SKY if "ab corner" in line or "focus" in line else "currentColor"
    b.append(T(tx, 70 + i * 22, line, 14, anchor="start", fill=col, weight=700 if i in (1, 2) else 400))
fig2 = save("bundle", W, H, "The cross-section of a light beam grows as (1+a)(1+b) = 1 + a ∘₁ b; the ab corner is the focusing term", "".join(b))

# ---------------------------------------------------------------- 3. bridges
b = []
rows = [
    ("Jacobson", [("area bookkeeping", 1), ("Unruh temperature", 1), ("entropy = η · area", 0), ("matter's energy flow", 0)]),
    ("Teleparallel", [("flat frame = river", 1), ("torsion = river gradient", 1), ("energy = strain cross terms", 1), ("waves need full frame", 0)]),
    ("Lovelock", [("4 dimensions", 1), ("conserved source", 1), ("metric only (dial k_ij)", 1), ("second order (assumed)", 0)]),
    ("Deser", [("Lorentz rule", 1), ("spin-2 field (dial k_ij)", 1), ("universal coupling", 2), ("scalar dial alone: Nordström, fails", 0)]),
]
W, H = 1040, 60 + 78 * len(rows) + 40
for j, (name, chips) in enumerate(rows):
    y = 60 + 78 * j
    b.append(T(20, y + 30, name, 16, anchor="start", weight=700))
    x = 140
    for lab, kind in chips:
        w = 9 + 7.1 * len(lab)
        if kind == 1:
            b.append(f'<rect x="{x}" y="{y+8}" width="{w}" height="34" rx="17" fill="{SKY}" fill-opacity=".25" stroke="{SKY}" stroke-width="1.5"/>')
            b.append(T(x + w / 2, y + 30, lab, 13))
        elif kind == 0:
            b.append(f'<rect x="{x}" y="{y+8}" width="{w}" height="34" rx="17" stroke="{PINK}" stroke-width="1.5" stroke-dasharray="5 4"/>')
            b.append(T(x + w / 2, y + 30, lab, 13, fill=PINK))
        else:
            b.append(f'<rect x="{x}" y="{y+8}" width="{w}" height="34" rx="17" stroke-opacity=".5" stroke-width="1.5"/>')
            b.append(T(x + w / 2, y + 30, lab, 13))
        x += w + 12
    b.append(f'<line x1="{max(x,880)-4}" y1="{y+25}" x2="905" y2="{y+25}" stroke-width="2" marker-end="url(#ar)"/>')
    b.append(T(970, y + 30, "Einstein + Λ", 14, weight=700))
yl = 60 + 78 * len(rows) + 14
b.append(f'<rect x="140" y="{yl}" width="22" height="16" rx="8" fill="{SKY}" fill-opacity=".25" stroke="{SKY}"/>')
b.append(T(170, yl + 13, "The circle operation, the dial or the river supplies it", 13, anchor="start"))
b.append(f'<rect x="520" y="{yl}" width="22" height="16" rx="8" stroke="{PINK}" stroke-dasharray="4 3"/>')
b.append(T(550, yl + 13, "comes from outside, or fails", 13, anchor="start", fill=PINK))
b.append(f'<rect x="800" y="{yl}" width="22" height="16" rx="8" stroke-opacity=".5"/>')
b.append(T(830, yl + 13, "follows from the others", 13, anchor="start"))
fig3 = save("bridges", W, H, "What each of the four bridges needs, and which parts the circle operation supplies", "".join(b))

# ---------------------------------------------------------------- 4. river and torsion
b = []
W, H = 1040, 270
b.append(f'<circle cx="70" cy="120" r="40" fill="currentColor" fill-opacity=".85" stroke="none"/>')
b.append(T(70, 185, "mass", 13))
for i, rr in enumerate([1.45, 2.0, 2.7, 3.6, 4.8, 6.4, 8.5]):
    x = 70 + rr * 105
    L = 70 / math.sqrt(rr)
    b.append(f'<line x1="{x}" y1="120" x2="{x-L}" y2="120" stroke="{SKY}" stroke-width="{1.5+2/rr}" marker-end="url(#ar)"/>')
b.append(T(560, 70, "the river: space flows inward at escape speed, v² = 2M/r", 14, weight=700))
items = [("speed", "v"), ("torsion", "v′ (how fast the flow speeds up)"), ("Newton's pull", "v v′ = −M/r²"), ("twist (curvature)", "(v v′)′ = 2M/r³")]
for i, (k, v) in enumerate(items):
    b.append(T(300 + 0 * i, 205 + 0, "", 1))
for i, (k, v) in enumerate(items):
    x = 140 + i * 225
    b.append(T(x, 222, k, 13, anchor="start", fill=SKY if i < 3 else "currentColor", weight=700))
    b.append(T(x, 244, v, 13, anchor="start"))
fig4 = save("river", W, H, "The river picture: torsion is the river's speed-up, Newton's pull is speed times torsion, and the curvature twist is its slope", "".join(b))

# ---------------------------------------------------------------- 5. Clifford product
b = []
W, H = 1040, 300
O = (120, 230)
A = (300, 210); B = (190, 90)
b.append(f'<path d="M{O[0]} {O[1]} L{A[0]} {A[1]} L{A[0]+B[0]-O[0]} {A[1]+B[1]-O[1]} L{B[0]} {B[1]}Z" fill="{SKY}" fill-opacity=".22" stroke="none"/>')
b.append(f'<line x1="{O[0]}" y1="{O[1]}" x2="{A[0]}" y2="{A[1]}" stroke-width="2.5" marker-end="url(#ar)"/>')
b.append(f'<line x1="{O[0]}" y1="{O[1]}" x2="{B[0]}" y2="{B[1]}" stroke-width="2.5" marker-end="url(#ar)"/>')
b.append(T(A[0] + 14, A[1] + 6, "a", 16, anchor="start")); b.append(T(B[0] - 14, B[1], "b", 16, anchor="end"))
b.append(T(255, 160, "a ∧ b", 15, fill=SKY, weight=700))
b.append(T(230, 280, "the plane they span = a twist", 13, fill=SKY))
tx = 450
lines = [("ab  =  a·b  +  a ∧ b", 18, "currentColor", 700),
         ("a·b = k_ij aⁱ bʲ : the direction dial (a metric)", 14, "currentColor", 400),
         ("a ∧ b : the twist in the a-b plane", 14, SKY, 400),
         ("", 14, "currentColor", 400),
         ("The circle operation with a Clifford weight k:  a + b + a k b", 16, "currentColor", 700),
         ("1 + k(a ✳ b) = (1 + k a)(1 + k b), associative for every k", 14, "currentColor", 400),
         ("k(x) alone: flat (pure gauge). Add a twist field:", 14, PINK, 400),
         ("gauge-theory gravity = general relativity", 14, "currentColor", 400)]
for i, (s_, sz, col, wt) in enumerate(lines):
    b.append(T(tx, 60 + i * 28, s_, sz, anchor="start", fill=col, weight=wt))
fig5 = save("clifford", W, H, "The Clifford product of two vectors keeps the dot (the direction dial) and the wedge (the twist); the circle operation with a Clifford weight is a + b + a k b", "".join(b))
print("ok")
