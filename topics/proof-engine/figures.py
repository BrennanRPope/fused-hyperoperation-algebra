"""
figures.py: draws the vocabulary plates and the worked example for the proof-engine page.
Run:  python3 -I figures.py   (writes figures/*.svg and figures/specs.json)
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapes import (Canvas, draw_rod, draw_rect, draw_box, draw_cuboid, colour, shade,
                    proj, VALUES, UNIT, PAPER, INK, BASE)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
A, B, C = VALUES["a"], VALUES["b"], VALUES["c"]
specs = {}


def save(name, cv, spec=None):
    with open(os.path.join(OUT, name + ".svg"), "w") as f:
        f.write(cv.svg())
    if spec:
        specs[name] = spec


def caption(cv, s, y=None, size=14):
    cv.text(cv.w / 2, y if y else cv.h - 16, s, size)


# ------------------------------------------------------------------ hero
def hero():
    cv = Canvas(760, 330)
    draw_rect(cv, 40, 250, "a", "b", S=78)
    cv.text(150, 290, "a ∘₁ b  =  a + b + ab", 16)
    cv.text(150, 312, "the rectangle (1+a)(1+b), unit square removed", 12, "#555")
    cv.line(285, 160, 345, 160, arrow=True)
    cv.text(315, 142, "one more side", 12, "#555")
    draw_box(cv, 380, 262, ("a", "b", "c"), S=86, explode=("a", "b", "c"), gap=0.12)
    cv.text(560, 300, "(a ∘₁ b) ∘₁ c = box (1+a)(1+b)(1+c) minus the unit cube", 13)
    cv.text(565, 320, "its 7 pieces: a, b, c, ab, ac, bc, abc  (layers 3, 3, 1)", 12, "#555")
    save("hero", cv, {"kind": "box", "names": "abc", "explode": "abc"})


# ------------------------------------------------------------------ vocabulary
def v_lift():
    cv = Canvas(360, 230)
    cv.text(24, 30, "a", 16, anchor="start", italic=True)
    draw_rod(cv, 60, 22, [("a", A)], S=70, labels=False)
    cv.line(180, 70, 180, 98, arrow=True)
    cv.text(190, 84, "lift  x ↦ 1 + x", 13, anchor="start")
    end = draw_rod(cv, 60, 118, [("1", 1.0), ("a", A)], S=70, brace="1 + a")
    caption(cv, "Lift: each number becomes a rod plus the unit", 185, 13)
    caption(cv, "madd_lift:  a ∘₁ b + 1 = (a+1)(b+1)", 208, 12)
    save("v_lift", cv, {"kind": "rod", "parts": [["1", 1.0], ["a", A]]})


def v_rect():
    cv = Canvas(360, 250)
    draw_rect(cv, 60, 195, "a", "b", S=70)
    caption(cv, "The circle operation: two rods make a rectangle; drop the unit", 218, 13)
    caption(cv, "a ∘₁ b = a + b + ab", 238, 13)
    save("v_rect", cv, {"kind": "rect", "xv": "a", "yv": "b", "unit": "dashed"})


def v_swap():
    cv = Canvas(360, 250)
    draw_rect(cv, 22, 190, "a", "b", S=50)
    # transposed: b along x, a along y
    draw_rect(cv, 200, 190, "b", "a", S=50, corner_key="ab")
    cv.path("M150,110 C170,95 185,95 198,110", arrow=True)
    cv.text(175, 86, "flip on the diagonal", 12, "#555")
    caption(cv, "Swap (commutativity): the same pieces, reflected", 214, 13)
    caption(cv, "a ∘₁ b = b ∘₁ a   ·  Nat.mul_comm", 236, 12)
    save("v_swap", cv)


def v_collapse():
    cv = Canvas(360, 250)
    draw_rect(cv, 40, 190, "a", "b", S=55)
    cv.line(185, 140, 215, 140, arrow=True)
    cv.text(200, 124, "b → 0", 12, "#555")
    draw_rect(cv, 225, 190, "a", "b", S=55, ylen=0.0, xlen=A * 0.8)
    caption(cv, "Collapse (identity): a zero rod flattens the box", 214, 13)
    caption(cv, "a ∘₁ 0 = a   ·  Nat.mul_one", 236, 12)
    save("v_collapse", cv)


def v_squash():
    cv = Canvas(360, 270)
    S, ox, oy = 62, 50, 200
    a = A
    w, h = 1 + a, 1 / (1 + a)
    # hyperbola xy = 1
    pts = []
    x = 0.42
    while x <= 3.6:
        pts.append((ox + x * S, oy - (1 / x) * S))
        x += 0.04
    d = "M" + " L".join(f"{px:.1f},{py:.1f}" for px, py in pts if py > 20)
    cv.path(d, stroke="#888", sw=1.2, dash="5 3")
    cv.text(ox + 0.5 * S, 30, "xy = 1", 12, "#777")
    # axes
    cv.line(ox, oy, ox + 3.6 * S, oy, sw=0.8)
    cv.line(ox, oy, ox, 22, sw=0.8)
    # unit square dashed
    cv.rect(ox, oy - S, S, S, PAPER, dash="4 3", sw=1)
    # the inverse rectangle (1+a) x 1/(1+a)
    cv.rect(ox, oy - h * S, w * S, h * S, UNIT, opacity=0.9)
    cv.rect(ox + S, oy - h * S, a * S, h * S, BASE["a"])          # gain: + a/(1+a)
    cv.rect(ox, oy - S, S, (1 - h) * S, "url(#hatch)", stroke=BASE["a"])  # loss: mirror strip
    cv.text(ox + S + a * S / 2, oy - h * S / 2, "+", 18, "#fff", weight=700)
    cv.text(ox + S / 2, oy - S + (1 - h) * S / 2, "−", 20, INK, weight=700)
    cv.text(ox + S * w, oy + 14, "1+a", 12)
    cv.text(ox - 22, oy - h * S, "1/(1+a)", 11, anchor="end") if False else cv.text(ox - 6, oy - h * S, "1/(1+a)", 11, anchor="end")
    cv.add(f'<circle cx="{ox + w * S:.1f}" cy="{oy - h * S:.1f}" r="3.5" fill="{INK}"/>')
    cv.text(ox + 3.0 * S, 120, "a = −1 is the wall:", 12, "#555")
    cv.text(ox + 3.0 * S, 136, "a rod of length 0", 12, "#555")
    cv.text(ox + 3.0 * S, 152, "can't be rescaled", 12, "#555")
    caption(cv, "Squash (inverse): rescale until the area is exactly 1", 238, 13)
    caption(cv, "a ∘₁ ā = 0 with ā = −a/(1+a): gain and mirror cancel", 258, 12)
    save("v_squash", cv)


def v_dial():
    cv = Canvas(360, 250)
    for i, (k, lab) in enumerate([(0.0, "k = 0  plain add"), (0.45, "k = A  soliton"), (1.0, "k = 1  circle operation")]):
        ox = 14 + i * 116
        draw_rect(cv, ox, 170, "a", "b", S=40, k=k, labels=False)
        cv.text(ox + 40 + A * 20, 188, lab, 11)
    caption(cv, "Dial (the weighted k-family): the corner holds k·ab", 214, 13)
    caption(cv, "a ∘ₖ b = a + b + k·ab;  1 + k(a ∘ₖ b) = (1+ka)(1+kb)", 236, 12)
    save("v_dial", cv)


def v_flatten():
    cv = Canvas(360, 250)
    draw_rect(cv, 24, 160, "a", "b", S=40, unit="solid", labels=False)
    cv.text(24 + 20 + 40 * A / 2 + 20, 176, "area (1+a)(1+b)", 11)
    cv.line(150, 110, 190, 110, arrow=True)
    cv.text(170, 94, "log", 13, italic=True)
    la, lb = math.log(1 + A), math.log(1 + B)
    S = 100
    draw_rod(cv, 200, 102, [("a", la), ("b", lb)], S=S, labels=False)
    cv.text(200 + S * la / 2, 132, "L(a)", 12, italic=True)
    cv.text(200 + S * la + S * lb / 2, 132, "L(b)", 12, italic=True)
    cv.path(f"M200,90 v-5 H{200 + S * (la + lb):.1f} v5", sw=1)
    cv.text(200 + S * (la + lb) / 2, 74, "L(a ∘₁ b)", 12, italic=True)
    caption(cv, "Flatten (log): L(x) = log(1+x) turns area into length", 214, 13)
    caption(cv, "L(a ∘₁ b) = L(a) + L(b); on whole numbers: pow_add", 236, 12)
    save("v_flatten", cv)


def v_extrude():
    cv = Canvas(560, 270)
    # rung 1: rods end to end
    draw_rod(cv, 16, 110, [("1", 1.0), ("a", A)], S=34, labels=False)
    draw_rod(cv, 16 + 34 * (1 + A), 110, [("1", 1.0), ("b", B)], S=34, labels=False)
    cv.text(95, 150, "rung 1: length", 13)
    cv.text(95, 168, "(1+a)+(1+b) − 1", 11, "#555")
    # rung 2: rectangle
    draw_rect(cv, 210, 170, "a", "b", S=40, labels=False)
    cv.text(260, 196, "rung 2: area", 13)
    cv.text(260, 214, "(1+a)(1+b) − 1", 11, "#555")
    # rung 3: cube (1+a)^3, i.e. a fused power fexp(a, 2)
    km = lambda S_: {0: "", 1: "a", 2: "ab", 3: "abc"}[len(S_)]
    draw_box(cv, 380, 175, ("x", "y", "z"), S=46, values={"x": A, "y": A, "z": A},
             labels=False, keymap=km, explode=("x", "y", "z"), gap=0.08)
    cv.text(460, 196, "rung 3: volume", 13)
    cv.text(460, 214, "(1+a)³ − 1 = fexp(a, 2)", 11, "#555")
    caption(cv, "Extrude (fused rungs): each rung adds a dimension; rung n+1 iterates rung n", 244, 13)
    caption(cv, "cube pieces by layer 3, 3, 1 (Pascal)  ·  fexp_iter: fexp(a, b+1) = a ∘₁ fexp(a, b)", 262, 12)
    save("v_extrude", cv)


def v_regroup():
    cv = Canvas(560, 270)
    draw_box(cv, 40, 190, ("a", "b", "c"), S=58, explode=("c",), gap=0.3)
    cv.text(140, 222, "slabs: (a ∘₁ b) ∘₁ c", 13)
    cv.text(280, 120, "=", 26)
    draw_box(cv, 330, 190, ("a", "b", "c"), S=58, explode=("a",), gap=0.3)
    cv.text(430, 222, "columns: a ∘₁ (b ∘₁ c)", 13)
    caption(cv, "Regroup (associativity): one box, cut two ways  ·  Nat.mul_assoc", 252, 13)
    save("v_regroup", cv)


def soliton():
    cv = Canvas(560, 300)
    w = {"ab": 0.55, "ac": 0.7, "bc": 0.4}
    w["abc"] = w["ab"] * w["ac"] * w["bc"]
    draw_box(cv, 60, 230, ("a", "b", "c"), S=74, explode=("a", "b", "c"), gap=0.14, weights=w,
             values={"a": 1.0, "b": 1.0, "c": 1.0})
    x = 360
    cv.text(x, 60, "three-wave cube (Hirota τ)", 14, anchor="start", weight=600)
    rows = ["τ = 1 + a + b + c", "  + A₁₂ ab + A₁₃ ac + A₂₃ bc",
            "  + A₁₂A₁₃A₂₃ abc", "", "each pair piece is shrunk by", "its weight; the top corner by",
            "the product of the three", "(dashed: the full cell)"]
    for i, r in enumerate(rows):
        cv.text(x, 88 + i * 20, r, 13 if i < 4 else 12, INK if i < 4 else "#555", anchor="start")
    caption(cv, "Weighted corners: the dial, applied to every corner of the cube at once", 284, 13)
    save("soliton", cv)


# ------------------------------------------------------------------ worked example: circ_assoc
def step_canvas():
    return Canvas(620, 250)


def eqsign(cv, x=310, y=120):
    cv.text(x, y, "=", 26)


def rod_with_label(cv, x, y, length, label, S, unit=False, fill="#5a5a66"):
    if unit:
        cv.rect(x, y, S, 18, UNIT)
        x0 = x + S
    else:
        x0 = x
    cv.rect(x0, y, length * S, 18, fill)
    cv.text(x0 + length * S / 2, y + 9, label, 12, "#fff", italic=True)
    if unit:
        cv.text(x + S / 2, y + 9, "1", 12)
    return x0 + length * S


def worked():
    ab = A + B + A * B
    lhs = ab + C + ab * C
    bc = B + C + B * C
    rhs = A + bc + A * bc
    assert abs(lhs - rhs) < 1e-12
    S = 40

    # step 0: the goal
    cv = step_canvas()
    rod_with_label(cv, 20, 110, lhs / 2, "(a ∘₁ b) ∘₁ c", S)
    eqsign(cv, 310, 119)
    rod_with_label(cv, 340, 110, rhs / 2, "a ∘₁ (b ∘₁ c)", S)
    caption(cv, "Two rods, built in different orders. Why are they the same length?", 200, 13)
    caption(cv, "(half scale)", 222, 11)
    save("w0_goal", cv)

    # step 1: apply lift_inj
    cv = step_canvas()
    rod_with_label(cv, 20, 110, lhs / 2, "(a ∘₁ b) ∘₁ c", S, unit=True)
    eqsign(cv, 336, 119)
    rod_with_label(cv, 360, 110, rhs / 2, "a ∘₁ (b ∘₁ c)", S, unit=True)
    caption(cv, "Attach the unit to both sides (it is safe to take off again).", 200, 13)
    caption(cv, "(half scale)", 222, 11)
    save("w1_lift", cv)

    # step 2: rw [madd_lift]  -- left outer
    cv = step_canvas()
    s2 = 34
    ox, oy = 20, 175
    cv.rect(ox, oy - s2, s2, s2, UNIT)
    cv.rect(ox + s2, oy - s2, ab * s2, s2, "#5a5a66")
    cv.text(ox + s2 + ab * s2 / 2, oy - s2 / 2, "a ∘₁ b", 12, "#fff", italic=True)
    cv.rect(ox, oy - s2 - C * s2, s2, C * s2, colour("c"))
    cv.rect(ox + s2, oy - s2 - C * s2, ab * s2, C * s2, "#9b8f74")
    cv.text(ox + s2 + ab * s2 / 2, oy - s2 - C * s2 / 2, "(a ∘₁ b)·c", 11, "#fff", italic=True)
    cv.text(ox + (1 + ab) * s2 / 2, oy + 16, "1 + (a ∘₁ b)", 12)
    cv.text(ox - 4, oy - s2 - C * s2 / 2, "c", 12, anchor="end", italic=True)
    eqsign(cv, 336, 119)
    rod_with_label(cv, 360, 110, rhs / 2, "a ∘₁ (b ∘₁ c)", S, unit=True)
    caption(cv, "Left side: the outer ∘₁ becomes a rectangle, rod 1+(a ∘₁ b) by rod 1+c.", 222, 13)
    save("w2_left_outer", cv)

    # step 3: rw [madd_lift]  -- left inner -> box (slabs)
    cv = step_canvas()
    draw_box(cv, 40, 185, ("a", "b", "c"), S=56, explode=("c",), gap=0.32)
    cv.text(140, 212, "(1+a)(1+b) · (1+c): a slab of height c", 12)
    eqsign(cv, 336, 119)
    rod_with_label(cv, 360, 110, rhs / 2, "a ∘₁ (b ∘₁ c)", S, unit=True)
    caption(cv, "The rod 1+(a ∘₁ b) opens into its own rectangle: the left side is a box.", 236, 13)
    save("w3_left_box", cv)

    # step 4: rw [madd_lift, madd_lift] -- right side -> box (columns)
    cv = step_canvas()
    draw_box(cv, 40, 185, ("a", "b", "c"), S=56, explode=("c",), gap=0.32)
    cv.text(140, 212, "slabs", 12)
    eqsign(cv, 310, 119)
    draw_box(cv, 350, 185, ("a", "b", "c"), S=56, explode=("a",), gap=0.32)
    cv.text(450, 212, "columns: (1+a) · (1+b)(1+c)", 12)
    caption(cv, "Right side, same two moves in the other order: a column of width a.", 236, 13)
    save("w4_right_box", cv)

    # step 5: rw [Nat.mul_assoc] -- same box; pieces by layer
    cv = Canvas(620, 290)
    draw_box(cv, 40, 200, ("a", "b", "c"), S=62, explode=("a", "b", "c"), gap=0.16)
    x = 330
    cv.text(x, 46, "Same box, so the goal closes.", 14, anchor="start", weight=600)
    layers = [("layer 1", ["a", "b", "c"]), ("layer 2", ["ab", "ac", "bc"]), ("layer 3", ["abc"])]
    y = 80
    for name, keys in layers:
        cv.text(x, y, name, 12, "#555", anchor="start")
        for j, k in enumerate(keys):
            cx = x + 70 + j * 62
            cv.rect(cx, y - 11, 52, 22, colour(k), rx=3)
            cv.text(cx + 26, y, k, 12, "#fff" if k not in ("c",) else INK, italic=True)
        y += 36
    cv.text(x, y + 6, "Take the unit cube off again:", 12, anchor="start")
    cv.text(x, y + 28, "a+b+c + ab+ac+bc + abc  on both sides", 13, anchor="start", italic=True)
    caption(cv, "Nat.mul_assoc says the slabs and the columns are one box: 1 + 3 + 3 + 1 pieces.", 272, 13)
    save("w5_same_box", cv)


if __name__ == "__main__":
    hero()
    v_lift(); v_rect(); v_swap(); v_collapse(); v_squash(); v_dial(); v_flatten()
    v_extrude(); v_regroup(); soliton(); worked()
    with open(os.path.join(OUT, "specs.json"), "w") as f:
        json.dump(specs, f, indent=1)
    print("wrote", sorted(os.listdir(OUT)))
