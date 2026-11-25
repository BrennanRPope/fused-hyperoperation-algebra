"""
shapes.py: the shape vocabulary of the fused proof engine, as a renderer.

A "shape spec" is a small JSON-able dict. The same spec is what the Lean side
would emit from a goal (see lean/FusedPictures.sketch.lean); here Python turns
it into a static SVG for pages. One vocabulary, two renderers.

Kinds
  rod      a row of segments laid end to end: [("1",1.0),("a",1.2)]
  rect     a (1+x)(1+y) rectangle split at 1 on each side; the unit square can
           be shown solid, dashed (removed) or hidden; the corner can carry a weight k
  box      a (1+x)(1+y)(1+z) box split at 1 on each side, drawn in oblique view;
           pieces are the 8 corners of the cube (subsets of {x,y,z}); "explode"
           pulls the pieces apart along chosen axes; corner weights shrink pieces
  hyper    the curve xy = 1 with the inverse rectangle (1+a) x 1/(1+a)

Colours follow Byrne's 1847 Euclid: a = blue, b = red, c = yellow, and a corner
piece is the mix of the colours of its sides (ab violet, ac green, bc orange,
abc black). The unit is paper grey. So a piece's colour tells you its product.

Core idea the colours encode: every mult-add expression lifts (x -> 1+x) to a
product of rods, and a product of n rods is a box whose 2^n pieces are the
corners of the n-cube. A fused identity is "the same box, cut another way".
"""

from __future__ import annotations
import math
from itertools import combinations

INK = "#10131a"
PAPER = "#ffffff"
UNIT = "#d9d6cf"

BASE = {"a": "#1d5aa6", "b": "#d0312d", "c": "#e8b417"}
MIX = {
    "": UNIT,
    "a": BASE["a"], "b": BASE["b"], "c": BASE["c"],
    "ab": "#6a3d9a", "ac": "#3f8f4f", "bc": "#e0701f",
    "abc": "#2b2b33",
}
VALUES = {"a": 1.2, "b": 0.8, "c": 0.6}   # drawing lengths only; the proofs are symbolic


def colour(key: str) -> str:
    key = "".join(sorted(key))
    return MIX.get(key, "#777")


def shade(hexcol: str, f: float) -> str:
    """f > 0 lightens toward white, f < 0 darkens toward black."""
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    if f >= 0:
        r, g, b = (int(v + (255 - v) * f) for v in (r, g, b))
    else:
        r, g, b = (int(v * (1 + f)) for v in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"


def text_on(hexcol: str) -> str:
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return INK if (0.299 * r + 0.587 * g + 0.114 * b) > 150 else "#ffffff"


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Canvas:
    def __init__(self, w: int, h: int):
        self.w, self.h = w, h
        self.items: list[str] = []

    def add(self, s: str):
        self.items.append(s)

    def rect(self, x, y, w, h, fill, stroke=INK, sw=1.2, dash=None, opacity=1.0, rx=0):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" '
                 f'fill="{fill}" fill-opacity="{opacity}" stroke="{stroke}" stroke-width="{sw}"{d}/>')

    def poly(self, pts, fill, stroke=INK, sw=1.0, dash=None, opacity=1.0):
        p = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<polygon points="{p}" fill="{fill}" fill-opacity="{opacity}" '
                 f'stroke="{stroke}" stroke-width="{sw}" stroke-linejoin="round"{d}/>')

    def line(self, x1, y1, x2, y2, stroke=INK, sw=1.2, dash=None, arrow=False):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        m = ' marker-end="url(#arr)"' if arrow else ""
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                 f'stroke="{stroke}" stroke-width="{sw}"{d}{m}/>')

    def path(self, d, stroke=INK, sw=1.4, fill="none", dash=None, arrow=False):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        m = ' marker-end="url(#arr)"' if arrow else ""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{da}{m}/>')

    def text(self, x, y, s, size=14, fill=INK, anchor="middle", weight=400, italic=False):
        st = ' font-style="italic"' if italic else ""
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
                 f'text-anchor="{anchor}" font-weight="{weight}"{st} '
                 f'dominant-baseline="middle">{esc(s)}</text>')

    def svg(self) -> str:
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" '
                f'width="{self.w}" height="{self.h}" font-family="Literata, Georgia, serif">'
                f'<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
                f'markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" '
                f'fill="{INK}"/></marker>'
                f'<pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" '
                f'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" '
                f'stroke="{BASE["a"]}" stroke-width="2"/></pattern></defs>'
                f'<rect width="100%" height="100%" fill="{PAPER}"/>')
        return head + "".join(self.items) + "</svg>"


# ---------------------------------------------------------------- rods

def draw_rod(cv: Canvas, x, y, parts, S=60, H=16, labels=True, brace=None):
    """parts: list of (key, length). Returns the end x."""
    cx = x
    for key, L in parts:
        fill = colour(key) if key != "1" else UNIT
        cv.rect(cx, y, L * S, H, fill)
        if labels and L * S > 14:
            lab = key if key else "1"
            cv.text(cx + L * S / 2, y + H + 13, lab, 13, INK, italic=(lab != "1"))
        cx += L * S
    if brace:
        cv.path(f"M{x:.1f},{y - 6:.1f} v-5 H{cx:.1f} v5", sw=1)
        cv.text((x + cx) / 2, y - 20, brace, 13)
    return cx


# ---------------------------------------------------------------- rectangles

def draw_rect(cv: Canvas, ox, oy, xv="a", yv="b", S=70, unit="dashed", k=1.0,
              labels=True, xlen=None, ylen=None, xname=None, yname=None, corner_key=None):
    """Bottom-left corner at (ox, oy) in SVG coordinates (y grows down).
    The rectangle is (1+x) by (1+y). unit: 'solid' | 'dashed' | 'none'.
    k: corner weight; the corner piece is drawn k of its full height."""
    X = VALUES.get(xv, 1.0) if xlen is None else xlen
    Y = VALUES.get(yv, 1.0) if ylen is None else ylen
    u = S
    ck = corner_key if corner_key is not None else xv + yv
    # unit square
    if unit == "solid":
        cv.rect(ox, oy - u, u, u, UNIT)
    elif unit == "dashed":
        cv.rect(ox, oy - u, u, u, PAPER, dash="4 3", sw=1)
        if labels:
            cv.text(ox + u / 2, oy - u / 2, "1", 13, "#777")
    # x strip (right of unit, height 1)
    if X > 0:
        cv.rect(ox + u, oy - u, X * S, u, colour(xv) if len(xv) == 1 else colour(xv))
        if labels and X * S > 16:
            cv.text(ox + u + X * S / 2, oy - u / 2, xname or xv, 15, text_on(colour(xv)), italic=True)
    # y strip (above unit, width 1)
    if Y > 0:
        cv.rect(ox, oy - u - Y * S, u, Y * S, colour(yv))
        if labels and Y * S > 16:
            cv.text(ox + u / 2, oy - u - Y * S / 2, yname or yv, 15, text_on(colour(yv)), italic=True)
    # corner
    if X > 0 and Y > 0:
        cv.rect(ox + u, oy - u - Y * S, X * S, Y * S, PAPER, dash="4 3", sw=1)
        if k > 0:
            h = Y * S * min(k, 1.0)
            cv.rect(ox + u, oy - u - h, X * S, h, colour(ck))
            if labels and h > 16 and X * S > 22:
                lab = (ck if k == 1 else f"{k:g}·{ck}") if not (xname or yname) else (
                    (xname or xv) + "·" + (yname or yv))
                cv.text(ox + u + X * S / 2, oy - u - h / 2, lab, 14, text_on(colour(ck)), italic=True)
    return (ox + u + X * S, oy - u - Y * S)


# ---------------------------------------------------------------- boxes (oblique)

def proj(x, y, z, S, ox, oy, dk=0.5, ang=math.radians(32)):
    return (ox + S * (x + y * dk * math.cos(ang)), oy - S * (z + y * dk * math.sin(ang)))


def draw_cuboid(cv: Canvas, p0, p1, fill, S, ox, oy, dash=None, opacity=1.0, sw=1.0):
    (x0, y0, z0), (x1, y1, z1) = p0, p1
    P = lambda x, y, z: proj(x, y, z, S, ox, oy)
    front = [P(x0, y0, z0), P(x1, y0, z0), P(x1, y0, z1), P(x0, y0, z1)]
    top = [P(x0, y0, z1), P(x1, y0, z1), P(x1, y1, z1), P(x0, y1, z1)]
    right = [P(x1, y0, z0), P(x1, y1, z0), P(x1, y1, z1), P(x1, y0, z1)]
    cv.poly(front, fill, dash=dash, opacity=opacity, sw=sw)
    cv.poly(top, shade(fill, 0.35) if fill != PAPER else PAPER, dash=dash, opacity=opacity, sw=sw)
    cv.poly(right, shade(fill, -0.22) if fill != PAPER else PAPER, dash=dash, opacity=opacity, sw=sw)
    return front


def box_pieces(names=("a", "b", "c"), values=None, explode=(), gap=0.18, weights=None):
    """Yield (subset_key, p0, p1, weight). Axis i carries names[i]; the piece for subset S
    spans [1, 1+v] on the axes in S and [0, 1] on the others. Exploded axes get a gap."""
    values = values or {n: VALUES.get(n, 1.0) for n in names}
    weights = weights or {}
    out = []
    for mask in range(8):
        S = "".join(n for i, n in enumerate(names) if mask >> i & 1)
        p0, p1 = [], []
        for i, n in enumerate(names):
            on = mask >> i & 1
            g = gap if (n in explode and on) else 0.0
            lo, hi = (1 + g, 1 + g + values[n]) if on else (0.0, 1.0)
            p0.append(lo)
            p1.append(hi)
        w = weights.get(S, 1.0)
        out.append((S, p0, p1, w, mask))
    return out


def draw_box(cv: Canvas, ox, oy, names=("a", "b", "c"), S=60, explode=(), gap=0.18,
             unit="dashed", weights=None, labels=True, values=None, keymap=None):
    """Axis order: names[0] -> screen x, names[1] -> depth, names[2] -> up.
    keymap maps a subset string to the colour key (when axes are relabelled)."""
    pieces = box_pieces(names, values, explode, gap, weights)
    # painter: back (depth) first, then low z, then left x
    pieces.sort(key=lambda t: (-t[1][1], t[1][2], t[1][0]))
    for S_, p0, p1, w, mask in pieces:
        key = keymap(S_) if keymap else S_
        if S_ == "":
            if unit == "none":
                continue
            fill = UNIT if unit == "solid" else PAPER
            draw_cuboid(cv, p0, p1, fill, S, ox, oy, dash=("3 3" if unit == "dashed" else None), sw=0.8)
            continue
        q1 = list(p1)
        if w != 1.0:  # weight shrinks the piece along its last axis
            q1[2] = p0[2] + (p1[2] - p0[2]) * w
            draw_cuboid(cv, p0, p1, PAPER, S, ox, oy, dash="2 3", sw=0.7)
        front = draw_cuboid(cv, p0, q1, colour(key), S, ox, oy)
        if labels:
            fx = (front[0][0] + front[1][0]) / 2
            fy = (front[0][1] + front[2][1]) / 2
            wpx = front[1][0] - front[0][0]
            hpx = front[0][1] - front[2][1]
            if wpx > 18 and hpx > 14:
                cv.text(fx, fy, key, 12 if len(key) < 3 else 11, text_on(colour(key)), italic=True)
    return pieces


# ---------------------------------------------------------------- render a spec

def render(spec: dict) -> str:
    cv = Canvas(spec.get("w", 360), spec.get("h", 240))
    k = spec["kind"]
    if k == "rod":
        draw_rod(cv, spec.get("x", 20), spec.get("y", 100), spec["parts"], spec.get("S", 60),
                 brace=spec.get("brace"))
    elif k == "rect":
        draw_rect(cv, spec.get("x", 40), spec.get("y", 200), spec.get("xv", "a"), spec.get("yv", "b"),
                  spec.get("S", 70), spec.get("unit", "dashed"), spec.get("k", 1.0))
    elif k == "box":
        draw_box(cv, spec.get("x", 60), spec.get("y", 200), tuple(spec.get("names", "abc")),
                 spec.get("S", 60), tuple(spec.get("explode", ())), unit=spec.get("unit", "dashed"),
                 weights=spec.get("weights"))
    else:
        raise ValueError(k)
    return cv.svg()
