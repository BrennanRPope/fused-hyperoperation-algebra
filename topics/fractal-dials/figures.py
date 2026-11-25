"""Draw the SVG plates for the fractal-dials page. Run: python3 figures.py"""
import math, os, json, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractals import HALF, intervals, tube_ratio
OUT = os.path.join(HERE, "figures")
RES = json.load(open(os.path.join(HERE, "results.json")))
INK, BLUE, RED, MUTE, RULE, WASH = "#15131a", "#1d2bd8", "#cf3519", "#5d5a66", "#d9d8e0", "#eef0ff"
AMBER, GREEN, GREY = "#b7791f", "#1f7a4d", "#9a98a3"
FONT = "font-family:Inter,system-ui,sans-serif"

def svg(w, h, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{title}">'
            f'<title>{title}</title><rect width="{w}" height="{h}" fill="#fff"/>'
            f'<g style="{FONT}">{body}</g></svg>')
def T(x, y, s, size=13, fill=INK, anchor="start", weight=400):
    s = s.replace("&", "&amp;").replace("<", "&lt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}">{s}</text>'
def R(x, y, w, h, fill, stroke="none", sw=1, op=1, rx=0):
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" fill-opacity="{op}"/>'
def L(x1, y1, x2, y2, c=INK, sw=1, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{c}" stroke-width="{sw}"{d}/>'
def P(pts, c=INK, sw=1.5, fill="none", dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" fill="{fill}" stroke="{c}" stroke-width="{sw}"{d}/>'
def arrow(x1, y1, x2, y2, c=MUTE):
    a = math.atan2(y2 - y1, x2 - x1); h = 7
    return (L(x1, y1, x2, y2, c, 1.4) +
            f'<polygon points="{x2:.1f},{y2:.1f} {x2-h*math.cos(a-.4):.1f},{y2-h*math.sin(a-.4):.1f} {x2-h*math.cos(a+.4):.1f},{y2-h*math.sin(a+.4):.1f}" fill="{c}"/>')

def bars(maps, x0, y0, w, levels, h=7, gap=5, color=BLUE):
    b = []
    for lv in range(levels + 1):
        for a, c in intervals(maps, lv):
            ww = max(float(c - a) * w, .6)
            b.append(R(x0 + float(a) * w, y0 + lv * (h + gap), ww, h, color, op=1 - .12 * lv if lv < 4 else .55))
    return b

# ---------------------------------------------------------------- plate 1: five half boxes, one dimension
def plate1():
    W, H = 900, 600
    b = [T(20, 30, "Five fractals, all of dimension ½, told apart by hidden dials", 15, INK, weight=600),
         T(20, 50, "Each row: the stick [0,1] cut down level by level (top = level 0). Same dimension, different shapes.", 12, MUTE)]
    cols = [(470, "pieces"), (540, "zoom"), (625, "lattice?"), (705, "content"), (790, "lacunarity")]
    for x, s in cols: b.append(T(x, 82, s, 12, MUTE, "middle", 600))
    levels = {"A": 4, "B": 3, "C": 3, "D": 4, "E": 2}
    y = 100
    for k in "ABCDE":
        v = HALF[k]; r = RES[k]
        b.append(T(20, y + 14, k, 15, RED, weight=700))
        b.append(T(40, y + 14, v["name"], 12, INK))
        b += bars(v["maps"], 40, y + 22, 380, levels[k])
        zoom = ", ".join(sorted({str(z) for z, _ in v["maps"]}))
        wob = (r["minkowski_max"] - r["minkowski_min"]) / r["minkowski_avg"]
        cells = [str(r["pieces"]), zoom, "yes" if r["lattice"] else "no",
                 f'{r["minkowski_avg"]:.2f}' + (" wobbles" if r["lattice"] else " settles"), f'{r["lacunarity"]:.2f}']
        for (x, _), s in zip(cols, cells):
            col = RED if (s == "no" or "settles" in s) else INK
            b.append(T(x, y + 40, s, 12, col, "middle"))
        y += 88
    b.append(L(20, y - 6, W - 20, y - 6, RULE))
    b.append(T(20, y + 12, "Dimension ½ for all five: the square roots of the zooms add to 1 (Lean: half_dimension_certificates).", 12, MUTE))
    b.append(T(20, y + 30, "B and C have the same pieces and zoom; only placement moves, and only lacunarity and content notice.", 12, MUTE))
    return svg(W, H, "".join(b), "Five fractals of dimension one half with their dial values")

# ---------------------------------------------------------------- plate 2: Minkowski wobble
def plate2():
    W, H = 900, 420
    x0, y0, PW, PH = 70, 70, 780, 240
    b = [T(20, 26, "The wobble dial: blur a fractal by ε and measure, then rescale by ε^½", 15, INK, weight=600)]
    lo, hi = 2.55, 3.7
    def X(t): return x0 + PW * t
    def Y(v): return y0 + PH * (hi - v) / (hi - lo)
    b.append(L(x0, y0 + PH, x0 + PW, y0 + PH, INK)); b.append(L(x0, y0, x0, y0 + PH, INK))
    for v in (2.6, 2.8, 3.0, 3.2, 3.4, 3.6):
        b.append(L(x0, Y(v), x0 + PW, Y(v), RULE)); b.append(T(x0 - 8, Y(v) + 4, f"{v:.1f}", 11, MUTE, "end"))
    for d in range(3, 8):
        t = (d - 3) / 4
        b.append(T(X(t), y0 + PH + 18, f"ε = 10^-{d}", 11, MUTE, "middle"))
    series = [("A", BLUE, "A: lattice, wobbles forever"), ("E", AMBER, "E: lattice, wobbles more"),
              ("D", RED, "D: non-lattice, settles to about 3.00")]
    n = 300
    for k, c, lab in series:
        pts = []
        for i in range(n + 1):
            t = i / n; e = 10 ** (-3 - 4 * t)
            pts.append((X(t), Y(tube_ratio(HALF[k]["maps"], .5, e))))
        b.append(P(pts, c, 1.6 if k != "D" else 2.2))
    for i, (k, c, lab) in enumerate(series):
        b.append(R(x0 + 10 + i * 250, y0 + PH + 38, 14, 4, c)); b.append(T(x0 + 30 + i * 250, y0 + PH + 44, lab, 12, INK))
    b.append(T(x0 - 60, y0 - 10, "measure ÷ ε^½", 11, MUTE))
    b.append(T(x0, y0 + PH + 76, "Lattice sets (every zoom a power of one ratio) never settle: they have complex dimensions ½ + i·2πn/ln(1/r),", 12, MUTE))
    b.append(T(x0, y0 + PH + 94, "which act like a phase that turns as you zoom. D's zooms share no power (Lean: uneven_not_lattice), so it settles.", 12, MUTE))
    return svg(W, H, "".join(b), "Minkowski wobble for lattice sets A and E against the settling non-lattice set D")

# ---------------------------------------------------------------- plate 3: weights are t
def plate3():
    W, H = 900, 450
    b = [T(20, 30, "The weight dial is the circle operation's per-element t: same set, same dimension, mass shared unevenly", 15, INK, weight=600)]
    rows = [(F(1, 2), "t = 1   (p = ½, fair)", "all Rényi dimensions D_q equal ½"),
            (F(3, 10), "t = 3/7 (p = 0.3)", "spectrum spreads: α from 0.26 to 0.87"),
            (F(1, 10), "t = 1/9 (p = 0.1)", "spreads further: α from 0.08 to 1.66")]
    lev = 5; x0, w = 230, 520
    y = 60
    for p, lab, note in rows:
        b.append(T(20, y + 50, lab, 13, INK, weight=600)); b.append(T(20, y + 68, note, 11, MUTE))
        masses = []
        for idx in range(2 ** lev):
            m = F(1)
            for j in range(lev): m *= p if (idx >> (lev - 1 - j)) & 1 == 0 else 1 - p
            masses.append(m)
        mx = float(max(masses)); ivs = intervals(HALF["A"]["maps"], lev)
        for (a, c), m in zip(ivs, masses):
            hgt = 80 * float(m) / mx
            b.append(R(x0 + float(a) * w - .6, y + 90 - hgt, max(float(c - a) * w, 1.2), hgt, BLUE))
        b.append(L(x0, y + 90, x0 + w, y + 90, RULE))
        y += 110
    b.append(T(20, y + 14, "Bars: mass on each level-5 piece of set A, scaled to the tallest. The support never changes, so dimension stays ½.", 12, MUTE))
    b.append(T(20, y + 32, "Only the split p : 1−p moves, which is exactly the odds t = p/(1−p) (Lean: weights_odds).", 12, MUTE))
    return svg(W, H, "".join(b), "Binomial measures on the same set with t equal to 1, 3/7 and 1/9")

# ---------------------------------------------------------------- plate 4: topology
def plate4():
    W, H = 900, 420
    b = [T(20, 30, "The topology dial is silent below dimension 1 and speaks once pieces touch", 15, INK, weight=600)]
    # left: B and C dust -> same Cantor set
    b.append(T(20, 62, "Below 1: B and C look different...", 13, INK, weight=600))
    b += bars(HALF["B"]["maps"], 20, 80, 330, 3, color=BLUE)
    b += bars(HALF["C"]["maps"], 20, 140, 330, 3, color=AMBER)
    b.append(T(20, 215, "...but every dust with no isolated points is the", 12, MUTE))
    b.append(T(20, 232, "same Cantor set up to stretching (Brouwer, 1910).", 12, MUTE))
    b.append(T(20, 249, "So shape below 1 is geometry, not topology.", 12, MUTE))
    # right: gasket vs dust at log3/log2
    b.append(T(420, 62, "At dimension log 3 / log 2 ≈ 1.585: same number, different topology", 13, INK, weight=600))
    def tri(x, y, s, lv, out):
        if lv == 0:
            out.append(f'<polygon points="{x:.1f},{y:.1f} {x+s:.1f},{y:.1f} {x+s/2:.1f},{y-s*.866:.1f}" fill="{BLUE}"/>'); return
        tri(x, y, s / 2, lv - 1, out); tri(x + s / 2, y, s / 2, lv - 1, out); tri(x + s / 4, y - s * .433, s / 2, lv - 1, out)
    g = []; tri(430, 330, 210, 4, g); b += g
    b.append(T(535, 360, "gasket: 3 pieces at ½", 12, INK, "middle"))
    b.append(T(535, 377, "pieces share corners at every level", 12, MUTE, "middle"))
    def dust(x, y, s, lv, out):
        if lv == 0: out.append(R(x, y, s, s, AMBER)); return
        for i in range(3):
            for j in range(3): dust(x + i * s * 3 / 8, y + j * s * 3 / 8, s / 4, lv - 1, out)
    g = []; dust(680, 120, 200, 3, g); b += g
    b.append(T(780, 360, "dust: 9 pieces at ¼", 12, INK, "middle"))
    b.append(T(780, 377, "pieces never touch", 12, MUTE, "middle"))
    b.append(T(420, 405, "Lean: same_dimension_gasket_dust, gasket_vertices, gasket_touches.", 11, MUTE))
    return svg(W, H, "".join(b), "Topology: two dusts below dimension one are the same set; gasket and dust differ at dimension 1.585")

# ---------------------------------------------------------------- plate 5: joining two half boxes
def plate5():
    W, H = 900, 400
    b = [T(20, 30, "Join two binary half boxes and you get the binary full box, as a tree", 15, INK, weight=600)]
    # four-corner set
    x0, y0, s = 30, 70, 230
    for lv, op in ((1, .15), (2, .35), (3, 1)):
        iv = intervals(HALF["A"]["maps"], lv)
        for a, c in iv:
            for a2, c2 in iv:
                b.append(R(x0 + float(a) * s, y0 + float(a2) * s, float(c - a) * s, float(c2 - a2) * s, BLUE, op=op))
    b.append(T(x0 + s / 2, y0 + s + 24, "A × A: the four-corner set", 12, INK, "middle"))
    b.append(T(x0 + s / 2, y0 + s + 41, "4 pieces at ¼ each step, dimension ½ + ½ = 1", 11, MUTE, "middle"))
    b.append(arrow(285, 185, 335, 185))
    # binary tree
    tx, ty = 350, 70
    def node(x, y, dx, d):
        if d == 0: return
        for sgn in (-1, 1):
            b.append(L(x, y, x + sgn * dx, y + 42, BLUE, 1.3)); node(x + sgn * dx, y + 42, dx / 2, d - 1)
    node(tx + 130, ty, 64, 5)
    b.append(T(tx + 130, y0 + s + 24, "binary tree: each corner = two bits", 12, INK, "middle"))
    b.append(T(tx + 130, y0 + s + 41, "distances agree within about ×2 either way (check 6b)", 11, MUTE, "middle"))
    b.append(arrow(620, 185, 670, 185))
    # stick with dyadic glue
    sx, sy, sw = 690, 120, 180
    for lv in range(6):
        n = 2 ** lv
        for i in range(n):
            b.append(R(sx + sw * i / n + .8, sy + lv * 14, sw / n - 1.6, 9, BLUE, op=1 - .1 * lv))
    b.append(L(sx + sw / 2, sy - 8, sx + sw / 2, sy + 88, RED, 1.2, "3,3"))
    b.append(T(sx + sw / 2, sy - 14, "0.0111… = 0.1000…", 11, RED, "middle"))
    b.append(T(sx + sw / 2, y0 + s + 24, "the stick [0,1]", 12, INK, "middle"))
    b.append(T(sx + sw / 2, y0 + s + 41, "= the tree with dyadic ends glued", 11, MUTE, "middle"))
    b.append(T(20, 370, "Binary branching + equal zooms + fair weights pins one half box: two pieces at ¼. Joining it with itself is the full box", 12, MUTE))
    b.append(T(20, 388, "as a tree; gluing pairs of ends gives the stick. Lean: four_corner_is_tree, dyadic_glue, pieces_join_closed.", 12, MUTE))
    return svg(W, H, "".join(b), "Half box times half box is the binary tree; gluing dyadic ends gives the stick")

# ---------------------------------------------------------------- plate 6: dial map
def plate6():
    W, H = 900, 560
    b = [T(20, 30, "Which fractal dials the circle operation already has", 15, INK, weight=600)]
    rows = [
        ("dimension D", "fractional circle power (1+a)^D", "has it", BLUE),
        ("weights p : 1−p (multifractal D_q, f(α))", "per-element dial t = p/(1−p)", "has it", BLUE),
        ("branching growth (how fast branching grows)", "rung flow at t = ½ (half-exponential)", "has it, separately", BLUE),
        ("wobble / complex dimensions", "a phase that turns with zoom; like twist or a Clifford rotor", "partial", AMBER),
        ("zoom ratio r (how much each piece shrinks)", "nothing yet: the circle operation has no notion of size", "missing: one candidate", RED),
        ("placement: lacunarity, Minkowski content", "geometry outside the algebra", "outside", GREY),
        ("topology: touching, gluing", "silent below D = 1; only matters at D ≥ 1", "outside", GREY),
    ]
    y = 64
    b.append(T(30, y, "fractal dial", 12, MUTE, weight=600)); b.append(T(380, y, "in the circle operation", 12, MUTE, weight=600))
    b.append(T(740, y, "status", 12, MUTE, weight=600))
    y += 14
    for fd, ma, st, c in rows:
        b.append(R(20, y, 860, 52, c, op=.08, rx=6)); b.append(R(20, y, 6, 52, c, rx=2))
        b.append(T(36, y + 31, fd, 13, INK, weight=600)); b.append(T(380, y + 31, ma, 12, INK))
        b.append(T(740, y + 31, st, 12, c, weight=700))
        y += 62
    b.append(T(20, y + 16, "If the circle operation adds one dial, it is a zoom (size) dial: then dimension stops being an input and comes out of", 12, MUTE))
    b.append(T(20, y + 34, "piece count and zoom, D = log(pieces) / log(1/zoom), via Moran's equation. Everything else is covered or is geometry.", 12, MUTE))
    return svg(W, H, "".join(b), "Map of fractal dials against the dials the circle operation already has")

# ---------------------------------------------------------------- plate 7: why shape dials cannot iterate
def plate7():
    W, H = 900, 430
    x0, y0, PW, PH = 80, 60, 470, 280
    b = [T(20, 30, "Why no shape dial gives the natural half step of the power set", 15, INK, weight=600)]
    nmax, ymax = 6, 64
    def X(n): return x0 + PW * n / nmax
    def Y(v): return y0 + PH * (1 - min(v, ymax) / ymax)
    b.append(L(x0, y0 + PH, x0 + PW, y0 + PH, INK)); b.append(L(x0, y0, x0, y0 + PH, INK))
    for n in range(nmax + 1): b.append(T(X(n), y0 + PH + 18, str(n), 11, MUTE, "middle"))
    for v in (0, 16, 32, 48, 64): b.append(T(x0 - 8, Y(v) + 4, str(v), 11, MUTE, "end"))
    b.append(T(x0 + PW / 2, y0 + PH + 38, "level n", 12, MUTE, "middle"))
    b.append(T(x0 - 50, y0 - 12, "log₂(count)", 12, MUTE))
    for bb, lab in ((2, "branching 2"), (4, "branching 4"), (16, "branching 16")):
        pts = [(X(n), Y(n * math.log2(bb))) for n in range(nmax + 1)]
        b.append(P(pts, BLUE, 1.4)); b.append(T(X(nmax) + 6, Y(nmax * math.log2(bb)) + (10 if bb == 4 else 4), lab, 11, BLUE))
    pts = [(X(n / 10), Y(2 ** (n / 10))) for n in range(0, 61) if 2 ** (n / 10) <= ymax]
    b.append(P(pts, RED, 2)); b.append(T(X(6) - 8, Y(ymax) + 14, "power set 2^(2ⁿ)", 11, RED, "end"))
    # one genuine half step h with h(h(x)) = 2^x: linear on [0, 1/2], then forced by h(x) = 2^(h^-1(x))
    def hinv(y): return y - .5 if y <= 1 else half(math.log2(y))
    def half(x): return .5 + x if x <= .5 else 2 ** hinv(x)
    pts = [(X(n / 20), Y(half(n / 20))) for n in range(0, 121) if half(n / 20) <= ymax]
    b.append(P(pts, AMBER, 2, dash="6,4")); b.append(T(X(6) + 6, Y(half(6)) - 6, "a half step", 11, AMBER))
    tx = 690
    for i, s in enumerate(["Shape dials (zoom, placement,", "weights, gluing) all keep counts", "c·bⁿ: straight lines here, and", "joining keeps them straight.",
                           "", "The power set's counts 2^(2ⁿ)", "curve up and fit no line", "(Lean: powerset_not_selfsimilar).", "",
                           "A half step must bend between", "them: that needs a dial on how", "fast branching grows, which is", "The circle operation's rung flow at t = ½.", "", "The dashed curve is one exact half", "step, h(h(n)) = 2ⁿ; there are many."]):
        b.append(T(tx, 80 + i * 18, s, 12, INK if s and i < 4 else MUTE))
    return svg(W, H, "".join(b), "Self-similar counts are lines; the power set curves up; a half step needs a growth dial")

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for i, f in enumerate((plate1, plate2, plate3, plate4, plate5, plate6, plate7), 1):
        names = {1: "five_half_boxes", 2: "wobble", 3: "weights_t", 4: "topology", 5: "joining_tree", 6: "dial_map", 7: "growth_dial"}
        open(os.path.join(OUT, f"p{i}_{names[i]}.svg"), "w").write(f())
    print("wrote", sorted(os.listdir(OUT)))
