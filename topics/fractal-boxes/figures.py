"""Draw the SVG plates for the fractal-boxes page. Run: python3 figures.py"""
from fractions import Fraction as F
import math, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
INK, BLUE, RED, MUTE, RULE, WASH = "#15131a", "#1d2bd8", "#cf3519", "#5d5a66", "#d9d8e0", "#eef0ff"
FONT = "font-family:Inter,system-ui,sans-serif"

def svg(w, h, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{title}">'
            f'<title>{title}</title><rect width="{w}" height="{h}" fill="#fff"/>'
            f'<g style="{FONT}">{body}</g></svg>')
def T(x, y, s, size=13, fill=INK, anchor="start", weight=400, style=""):
    s = s.replace("&", "&amp;").replace("<", "&lt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}" {style}>{s}</text>'
def R(x, y, w, h, fill, stroke="none", sw=1, op=1):
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" fill-opacity="{op}"/>'
def L(x1, y1, x2, y2, c=INK, sw=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{c}" stroke-width="{sw}"{d}/>'

def gbinom(d, k):
    r = F(1)
    for i in range(k): r = r * (d - i) / (i + 1)
    return r

# ---------------------------------------------------------------- plate 1
def plate1():
    b = []
    # left: stick, square, cube corner pieces
    b.append(T(20, 28, "Whole boxes: every corner piece counts +1", 14, INK, weight=600))
    # stick 1 + a
    x0, y0, u = 30, 70, 60
    b.append(R(x0, y0, u, 22, WASH, INK)); b.append(R(x0+u, y0, u*0.8, 22, "#ffd9cf", INK))
    b.append(T(x0+u/2, y0+16, "1", 12, INK, "middle")); b.append(T(x0+u*1.4, y0+16, "a", 12, INK, "middle"))
    b.append(T(x0, y0+42, "stick: 2 pieces", 12, MUTE))
    # square
    x1, y1 = 30, 130
    b.append(R(x1, y1, u, u, WASH, INK)); b.append(R(x1+u, y1, u*.8, u, "#ffd9cf", INK))
    b.append(R(x1, y1+u, u, u*.8, "#ffd9cf", INK)); b.append(R(x1+u, y1+u, u*.8, u*.8, "#ffb3a1", INK))
    for (tx, ty, s) in ((x1+u/2, y1+u/2+4, "1"), (x1+u*1.4, y1+u/2+4, "a"), (x1+u/2, y1+u*1.4+4, "b"), (x1+u*1.4, y1+u*1.4+4, "ab")):
        b.append(T(tx, ty, s, 12, INK, "middle"))
    b.append(T(x1, y1+u*1.8+20, "square: 4 pieces", 12, MUTE))
    # cube (oblique)
    cx, cy = 185, 150
    b.append(T(cx, cy-62, "cube: 8 pieces", 12, MUTE))
    def quad(pts, fill):
        return f'<polygon points="{" ".join(f"{p[0]:.1f},{p[1]:.1f}" for p in pts)}" fill="{fill}" stroke="{INK}" stroke-width="1"/>'
    s, dx, dy = 44, 18, -14
    for i, (fx, fill) in enumerate(((0, WASH), (1, "#ffd9cf"))):
        X = cx + i*s
        b.append(quad([(X, cy), (X+s, cy), (X+s, cy+s), (X, cy+s)], fill))
    for i in range(2):
        X = cx + i*s
        b.append(quad([(X, cy), (X+s, cy), (X+s+dx, cy+dy), (X+dx, cy+dy)], "#ffb3a1" if i else "#ffd9cf"))
        b.append(quad([(X+dx, cy+dy), (X+s+dx, cy+dy), (X+s+2*dx, cy+2*dy), (X+2*dx, cy+2*dy)], "#ff8e75" if i else "#ffb3a1"))
    b.append(quad([(cx+2*s, cy), (cx+2*s+dx, cy+dy), (cx+2*s+dx, cy+s+dy), (cx+2*s, cy+s)], "#ffb3a1"))
    b.append(quad([(cx+2*s+dx, cy+dy), (cx+2*s+2*dx, cy+2*dy), (cx+2*s+2*dx, cy+s+2*dy), (cx+2*s+dx, cy+s+dy)], "#ff8e75"))
    b.append(T(20, 290, "Layer k = pieces made of k sides:", 12, MUTE))
    b.append(T(20, 307, "C(n,k): 1, 3, 3, 1 for the cube (2³ = 8)", 12, MUTE))
    # right: half box layers as signed bars
    X0, Y0, W, H = 360, 50, 380, 250
    b.append(T(X0, 28, "A box of dimension ½: layers C(½,k), some negative", 14, INK, weight=600))
    base = Y0 + 160
    scale = 140
    b.append(L(X0, base, X0+W, base, INK, 1))
    run = F(0)
    pts = []
    for k in range(9):
        v = gbinom(F(1, 2), k); run += v
        bx = X0 + 14 + k*40
        hgt = float(v) * scale
        col = BLUE if v > 0 else RED
        if v > 0: b.append(R(bx, base-hgt, 24, hgt, col, op=.85))
        else: b.append(R(bx, base, 24, -hgt, col, op=.85))
        lab = {0: "1", 1: "½", 2: "−⅛", 3: "1/16", 4: "−5/128", 5: "7/256", 6: "−21/1024", 7: "33/2048", 8: "…"}[k]
        ly = base - hgt - 6 if v > 0 else base - hgt + 14
        b.append(T(bx+12, ly, lab, 10.5 if k < 5 else 9, col, "middle"))
        b.append(T(bx+12, base+52, str(k), 11, MUTE, "middle"))
        pts.append((bx+12, float(run)))
    b.append(T(X0-4, base+56, "layer", 10, MUTE, "end"))
    # running total panel
    py0, ph = base + 80, 0
    b.append(T(X0+W, Y0-8, "", 11))
    # running totals as dots on a small axis above
    ty = base + 78
    b.append(T(X0+W-4, ty-2, "", 11, MUTE))
    b.append(T(X0-4, ty+16, "total", 10, MUTE, "end"))
    for (x, r) in pts:
        b.append(T(x, ty+16, f"{r:.3f}", 9.5, INK, "middle"))
    b.append(T(X0+W-4, base+122, "1 + ½ − ⅛ + 1/16 − 5/128 + … = √2 ≈ 1.414", 13, INK, "end", 600))
    return svg(760, 345, "".join(b), "Corner pieces of whole boxes and of a half box")

# ---------------------------------------------------------------- plate 2: two half boxes join
def plate2():
    b = []
    n = 6
    lay = [gbinom(F(1, 2), k) for k in range(n)]
    lab = ["1", "½", "−⅛", "1/16", "−5/128", "7/256"]
    x0, y0, c = 120, 70, 62
    b.append(T(20, 28, "Joining two half boxes: multiply layer by layer", 14, INK, weight=600))
    b.append(T(20, 48, "Each cell is (layer i) × (layer j). Add each diagonal: the negatives cancel.", 12, MUTE))
    for i in range(n):
        b.append(T(x0 + i*c + c/2, y0-8, lab[i], 11, BLUE if lay[i] > 0 else RED, "middle"))
        b.append(T(x0-10, y0 + i*c + c/2 + 4, lab[i], 11, BLUE if lay[i] > 0 else RED, "end"))
    pal = {0: "#dfe3ff", 1: "#ffe1d9", 2: "#e9f7e7", 3: "#fff4cc", 4: "#efe3ff", 5: "#e3f4fb"}
    for i in range(n):
        for j in range(n):
            v = lay[i]*lay[j]
            k = i + j
            fill = pal[k] if k < 6 else "#f4f4f6"
            b.append(R(x0+j*c, y0+i*c, c, c, fill, "#ffffff", 2))
            s = f"{v}" if v.denominator < 1000 else f"{float(v):.4f}"
            b.append(T(x0+j*c+c/2, y0+i*c+c/2+4, s, 10.5, RED if v < 0 else INK, "middle"))
    sums = [sum(lay[i]*lay[k-i] for i in range(k+1)) for k in range(6)]
    sx = x0 + n*c + 40
    b.append(T(sx, y0-8, "diagonal sums", 12, INK, weight=600))
    names = ["empty piece", "one side", "two sides", "three sides", "four sides", "five sides"]
    for k in range(6):
        y = y0 + 22 + k*42
        b.append(R(sx, y-15, 18, 18, pal[k], INK, .5))
        b.append(T(sx+28, y, f"layer {k}: {sums[k]}", 13, INK, weight=600))
        b.append(T(sx+28, y+15, names[k], 11, MUTE))
    b.append(T(sx, y0+290, "1, 1, 0, 0, … = a stick", 14, BLUE, weight=600))
    b.append(T(sx, y0+308, "(checked exactly through layer 10)", 11, MUTE))
    return svg(760, 445, "".join(b), "Two half boxes joined give a stick")

# ---------------------------------------------------------------- plate 3: half box as fractal
def cantor(ratio, offsets, level):
    iv = [(0.0, 1.0)]
    for _ in range(level):
        iv = [(a + (bb-a)*o, a + (bb-a)*o + (bb-a)*ratio) for (a, bb) in iv for o in offsets]
    return iv
def plate3():
    b = []
    b.append(T(20, 28, "The half box as a fractal: keep the two outer quarters, again and again", 14, INK, weight=600))
    x0, W = 30, 330
    for m in range(6):
        y = 70 + m*34
        for (a, c) in cantor(.25, [0, .75], m):
            b.append(R(x0 + a*W, y, max((c-a)*W, .6), 10, BLUE))
        b.append(T(x0+W+12, y+9, f"{2**m} piece{'s' if m else ''} of 1/{4**m}", 11, MUTE))
    b.append(T(x0, 290, "Each level: twice the pieces, a quarter the size.", 12, INK))
    b.append(T(x0, 307, "Double the size → √2 times as much: dimension ½.", 12, INK))
    # right: four-corner set
    X, Y, S = 500, 60, 230
    b.append(T(X-10, 48, "Half box × half box (levels 1, 2, 3)", 12, INK, weight=600))
    for lev, op in ((1, .15), (2, .35), (3, 1)):
        iv = cantor(.25, [0, .75], lev)
        for (a, c) in iv:
            for (p, q) in iv:
                b.append(R(X + a*S, Y + p*S, (c-a)*S, (q-p)*S, RED, op=op))
    b.append(R(X, Y, S, S, "none", RULE, 1))
    b.append(T(X-10, Y+S+22, "4ᵐ squares of side 1/4ᵐ: dimension 1,", 12, INK))
    b.append(T(X-10, Y+S+39, "like a stick, yet it is dust: no two touch,", 12, INK))
    b.append(T(X-10, Y+S+56, "and its shadow on the side shrinks to nothing.", 12, INK))
    return svg(760, 380, "".join(b), "The half box drawn as a Cantor set and its square")

# ---------------------------------------------------------------- plate 4: many halves
def count_at(iv, eps):
    pts = [a for a, c in iv] + [c for a, c in iv]
    return len({int(p/eps - 1e-12) if p > 0 else 0 for p in pts})
def plate4():
    b = []
    b.append(T(20, 28, "Many different halves: same dimension ½, different shapes, different counts", 14, INK, weight=600))
    shapes = [("two pieces at ¼", .25, [0, .75], 2, BLUE, 7),
              ("three pieces at ⅑", 1/9, [0, 4/9, 8/9], 3, RED, 5),
              ("uneven: 0.09 and 0.49", None, None, 2, "#2a8a3a", 7)]
    def uneven(level):
        iv = [(0.0, 1.0)]
        for _ in range(level):
            iv = [x for (a, c) in iv for x in ((a, a + (c-a)*.09), (c - (c-a)*.49, c))]
        return iv
    x0, W = 30, 320
    ivs = []
    for si, (name, r, offs, n, col, lev) in enumerate(shapes):
        y = 62 + si*78
        b.append(T(x0, y, name, 12, col, weight=600))
        for m in (1, 2, lev if si != 2 else 6):
            iv = cantor(r, offs, m) if r else uneven(m)
            yy = y + 8 + (0 if m == 1 else 14 if m == 2 else 28)
            for (a, c) in iv: b.append(R(x0 + a*W, yy, max((c-a)*W, .6), 8, col))
        ivs.append(cantor(r, offs, lev + 3) if r else uneven(15))
    # right: count / sqrt(1/eps) vs halvings
    X, Y, PW, PH = 420, 60, 310, 200
    b.append(R(X, Y, PW, PH, "none", RULE))
    b.append(T(X, Y-10, "pieces seen at scale 1/2ᵏ, divided by 2^(k/2)", 12, INK, weight=600))
    ks = list(range(4, 21))
    lo, hi = 0.8, 2.6
    def px(k): return X + (k-4)/(20-4)*PW
    def py(v): return Y + PH - (v-lo)/(hi-lo)*PH
    for v in (1.0, 1.5, 2.0, 2.5):
        b.append(L(X, py(v), X+PW, py(v), RULE, .6)); b.append(T(X-6, py(v)+4, f"{v:g}", 10, MUTE, "end"))
    for k in (4, 8, 12, 16, 20):
        b.append(T(px(k), Y+PH+16, str(k), 10, MUTE, "middle"))
    b.append(T(X+PW/2, Y+PH+32, "halvings k", 11, MUTE, "middle"))
    for (name, r, offs, n, col, lev), iv in zip(shapes, ivs):
        pts = [(px(k), py(count_at(iv, 2.0**-k)/2**(k/2))) for k in ks]
        b.append('<polyline fill="none" stroke="%s" stroke-width="2" points="%s"/>' % (col, " ".join(f"{a:.1f},{c:.1f}" for a, c in pts)))
    b.append(T(30, 330, "A whole box has an exact corner count (2ⁿ). A half box only has a growth rate:", 12, INK))
    b.append(T(30, 347, "the count grows like √2 per halving, times a wobble that depends on the shape.", 12, INK))
    b.append(T(30, 364, "So the shape does not pin down one fractional power set; only the number 2^½ survives.", 12, INK))
    return svg(760, 380, "".join(b), "Three different fractals of dimension one half")

# ---------------------------------------------------------------- plate 5: iteration
def build_half(seed):
    h = {2: 3, 3: 4}
    for _ in range(60):
        for k in sorted(list(h)):
            if h[k] not in h and h[k] < 4000:
                h[h[k]] = 2**k
        missing = [x for x in range(2, 60) if x not in h]
        if not missing: break
        x = missing[0]
        lo = max(v for k, v in h.items() if k < x)
        h[x] = lo + 1 + (seed if x > 4 else 0)
    return h
def plate5():
    b = []
    b.append(T(20, 28, "The power set of a shape is the space of all its pieces", 14, INK, weight=600))
    # 3 points -> 7 subsets
    x0, y0 = 30, 70
    b.append(T(x0, y0-12, "At one scale the shape looks like 3 points …", 12, MUTE))
    for i in range(3): b.append(f'<circle cx="{x0+20+i*30}" cy="{y0+14}" r="6" fill="{BLUE}"/>')
    b.append(T(x0, y0+60, "… so its pieces look like 2³ − 1 = 7 points:", 12, MUTE))
    subs = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
    for si, s in enumerate(subs):
        bx = x0 + (si % 4)*78; by = y0 + 80 + (si//4)*40
        b.append(R(bx, by, 70, 28, "#fff", RULE, 1))
        for i in range(3):
            b.append(f'<circle cx="{bx+15+i*20}" cy="{by+14}" r="5" fill="{RED if s[i] else "#eeeef2"}"/>')
    b.append(T(x0, y0+190, "Pieces of a shape seen at one scale = all subsets of what", 12, INK))
    b.append(T(x0, y0+207, "you see there. Geometry's own power set, scale by scale.", 12, INK))
    b.append(T(x0, y0+235, "A half step would need h(h(N)) = 2^N on these counts.", 12, INK, weight=600))
    b.append(T(x0, y0+252, "Products, unions and stretching only change counts", 12, INK))
    b.append(T(x0, y0+269, "linearly in size, so they can never get there twice.", 12, INK))
    # right: log2 counts
    X, Y, PW, PH = 420, 60, 310, 240
    b.append(R(X, Y, PW, PH, "none", RULE))
    b.append(T(X, Y-10, "count N of the shape → count of the new space", 12, INK, weight=600))
    Nmax = 34
    def px(N): return X + (N-2)/(Nmax-2)*PW
    def py(l2): return Y + PH - l2/34*PH
    for l in (0, 8, 16, 24, 32):
        b.append(L(X, py(l), X+PW, py(l), RULE, .6)); b.append(T(X-6, py(l)+4, f"2^{l}", 10, MUTE, "end"))
    for N in (2, 10, 18, 26, 34): b.append(T(px(N), Y+PH+16, str(N), 10, MUTE, "middle"))
    b.append(T(X+PW/2, Y+PH+32, "N", 11, MUTE, "middle"))
    Ns = list(range(2, Nmax+1))
    def line(vals, col, sw=2, dash=""):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        pts = [(px(N), py(min(v, 34))) for N, v in zip(Ns, vals)]
        return '<polyline fill="none" stroke="%s" stroke-width="%s"%s points="%s"/>' % (col, sw, d, " ".join(f"{a:.1f},{c:.1f}" for a, c in pts))
    b.append(line([math.log2(N) for N in Ns], INK))
    b.append(line([N for N in Ns], RED))
    hA, hB = build_half(0), build_half(3)
    b.append(line([math.log2(hA[N]) for N in Ns], BLUE, 2))
    b.append(line([math.log2(hB[N]) for N in Ns], BLUE, 2, "5 4"))
    b.append(T(px(30), py(5.6), "same shape (N)", 11, INK, "middle"))
    b.append(T(px(14), py(17), "all pieces (2^N)", 11, RED, "middle"))
    b.append(T(px(21), py(14.5), "two different", 11, BLUE, "middle"))
    b.append(T(px(21), py(12.7), "half steps", 11, BLUE, "middle"))
    return svg(760, 360, "".join(b), "Pieces of a shape and the half step between")

# ---------------------------------------------------------------- plate 6: infinite sizes
def plate6():
    b = []
    b.append(T(20, 28, "At infinity every fractal has the same number of points", 14, INK, weight=600))
    X = 60
    yN, yC, yP = 290, 120, 40
    b.append(L(X, 60, X, 310, INK, 1.2))
    b.append(T(X-10, yN+4, "ℵ₀", 15, INK, "end", 600)); b.append(T(X-10, yC+4, "2^ℵ₀", 15, INK, "end", 600))
    b.append(T(X+14, yN+4, "the naturals, the rationals: countable sets (dimension 0)", 12, INK))
    for i in range(16): b.append(f'<circle cx="{X+430+i*14}" cy="{yN}" r="2.5" fill="{INK}"/>')
    # gap band
    b.append(R(X+4, yC+22, 680, yN-yC-44, "#fff4cc"))
    b.append(T(X+18, (yC+yN)/2-6, "The gap where a half power set of ℕ would have to live.", 12.5, INK, weight=600))
    b.append(T(X+18, (yC+yN)/2+12, "Open without the axiom of choice (earlier page). No fractal can sit here:", 12, INK))
    b.append(T(X+18, (yC+yN)/2+29, "every closed set is countable or holds a full Cantor set (Cantor–Bendixson).", 12, INK))
    # Cantor sets at continuum
    shapes = [(1/9, [0, 4/9, 8/9], 4, "#2a8a3a", "three at ⅑: dim ½"), (.25, [0, .75], 5, BLUE, "two at ¼: dim ½"),
              (1/3, [0, 2/3], 5, MUTE, "middle thirds: dim 0.63"), (None, None, 0, RED, "the whole stick: dim 1")]
    for i, (r, offs, lev, col, name) in enumerate(shapes):
        x = X + 20 + i*172
        W = 140
        iv = cantor(r, offs, lev) if r else [(0, 1)]
        for (a, c) in iv: b.append(R(x + a*W, yC-6, max((c-a)*W, .6), 9, col))
        b.append(T(x, yC-14, name, 10.5, col))
    b.append(T(X+14, yC+16, "all = P(ℕ) as sets: a point is a left/right path, a path is a subset of ℕ", 12, INK))
    b.append(T(X+14, yP+30, "Dimension measures how spread out the points are, not how many there are.", 12, MUTE))
    return svg(760, 330, "".join(b), "All fractals have the size of the power set of the naturals")

for name, fn in (("p1_half_box_layers", plate1), ("p2_two_halves_join", plate2), ("p3_half_box_fractal", plate3),
                 ("p4_many_halves", plate4), ("p5_pieces_of_a_shape", plate5), ("p6_infinite_sizes", plate6)):
    open(os.path.join(OUT, name + ".svg"), "w").write(fn())
    print("wrote", name)
