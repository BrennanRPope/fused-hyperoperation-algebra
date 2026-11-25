"""
film.py: turns the compiled engine's own trace into proof films.

Run Lean on Demo.lean with trace.fused.pictures on, save the output to
build/demo_trace.txt, then:  python3 -I film.py
Every frame here is drawn from a shape that FusedPictures.lean read off a real
goal while the proof ran. Nothing is drawn by hand. Writes figures/film/*.svg
and figures/film/film.json.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from shapes import (Canvas, draw_box, colour as _colour, shade, text_on, VALUES, UNIT, PAPER, INK)

OUT = os.path.join(HERE, "figures", "film")
os.makedirs(OUT, exist_ok=True)
S = 46
DEMOS = [
    ("comm_by_moves", "multAdd a b = multAdd b a", "Commutativity: the rectangle flipped on its diagonal is the same rectangle."),
    ("assoc_by_moves", "multAdd (multAdd a b) c = multAdd a (multAdd b c)", "Associativity: slabs and columns are one box."),
    ("zero_by_moves", "multAdd a 0 = a", "Identity: a side of length 0 flattens the rectangle to a strip."),
    ("wall_by_moves", "multAdd (−1) b = −1", "The wall: a side of length −1 makes the lifted rod 0, and the whole box with it."),
    ("mpow_add_by_moves", "mpow a (m + n) = multAdd (mpow a m) (mpow a n)", "A cube of m + n sides is a cube of m sides times a cube of n."),
    ("wcomm_by_moves", "k · wAdd k a b = k · wAdd k b a", "Dial: scale the weighted corner into a plain box, then flip it."),
    ("fexp_iter_by_moves", "fexp a (b + 1) = multAdd a (fexp a b)", "Rung 2 on whole numbers: one more slab extrudes the fused power."),
    ("inv_by_moves", "multAdd u (multAddInv u) = 0", "Inverse in a field: squash the second side to 1/(1+u); the rectangle has area 1."),
]
CAPTION = {
    "goal": "The goal, as the engine reads it.",
    "lift": "Lift: x ↦ 1+x attaches the unit to both sides. Safe because 1+x = 1+y gives x = y.",
    "multadd": "Mult-add: every ⊕ opens into a product of rods, a box with its unit attached.",
    "swap": "Swap: reflect across the diagonal (CommSemiring.mul_comm).",
    "regroup": "Regroup: re-cut the same box (Semiring.mul_assoc).",
    "collapse": "Collapse: zero lengths and units fall away.",
    "flatten": "Flatten: powers add (pow_add').",
    "dial": "Dial: k times a k-weighted corner is the plain box on k·a and k·b (wAdd_scale).",
    "extrude": "Extrude: rung n+1 is rung n with one more slab (fexp_iter).",
    "squash": "Squash: (1+u) times 1/(1+u) is 1 (inv_one_add, mul_inv_cancel).",
}
MOVE_LEMMA = {
    "lift": "lift_back_ring / lift_inj", "multadd": "multAdd_lift, mpow_lift, …",
    "swap": "CommSemiring.mul_comm", "regroup": "Semiring.mul_assoc",
    "collapse": "mul_one, zero_add, …", "flatten": "pow_add'", "dial": "wAdd_scale",
    "extrude": "fexp_iter", "squash": "inv_one_add",
}


def colour(key):
    c = _colour(key)
    return "#7a7a7a" if len(c) < 7 else c


def length(side):
    n = side["name"].replace(" ", "")
    try:
        return float(n)
    except ValueError:
        pass
    if n.startswith("multAddInv"):
        return "inv"
    return VALUES.get(side["key"], 1.0)


def tidy(name):
    import re
    return name.replace("multAddInv u", "u⁻").replace(" * ", "·")


# ------------------------------------------------------------ sizes and drawing

def size(sh):
    k = sh.get("kind")
    if k == "rod":
        L = length(sh["side"])
        L = L if isinstance(L, float) else 1.0
        tot = (1.0 if sh["lifted"] else 0.0) + max(L, 0.0)
        return (max(tot, 0.6) * S + 10, 52)
    if k == "box":
        b = sh["box"]
        sides = b["sides"]
        Ls = [length(s) for s in sides]
        if len(sides) == 1:
            L = Ls[0] if isinstance(Ls[0], float) else 1.0
            return ((1 + L) * S + 60, (1 + L) * S + 40)
        if len(sides) == 2:
            X, Y = Ls
            X = X if isinstance(X, float) else 1.0
            if Y == "inv":
                return ((1 + X) * S + 60, S + 30)
            lw = max(len(tidy(sides[1]["name"])) * 7 - S, 0) / 2
            return ((1 + max(X, 0)) * S + 12 + lw, (1 + max(Y, 0)) * S + 40)
        if len(sides) == 3:
            v = [L if isinstance(L, float) else 1.0 for L in Ls]
            c, s_ = 0.5 * math.cos(math.radians(32)), 0.5 * math.sin(math.radians(32))
            w = S * (1 + v[0] + 0.25 + (1 + v[1]) * c) + 10
            h = S * (1 + v[2] + 0.25 + (1 + v[1]) * s_) + 24
            return (w, h)
    if k == "piece":
        Ls = [length(s) for s in sh["sides"]]
        Ls = [L if isinstance(L, float) else 1.0 for L in Ls]
        lab = max(len(tidy(x["name"])) for x in sh["sides"]) * 7
        return (max(Ls[0] * S, lab) + 16, Ls[1] * S + 40)
    if k == "flip":
        w, h = size(sh["inner"])
        return (w + 18, h)
    if k == "note":
        return (150, 90)
    return (130, 50)


def draw(cv, sh, x, yb):
    """Draw shape with its bottom-left near (x, yb)."""
    k = sh.get("kind")
    if k == "rod":
        side = sh["side"]
        L = length(side)
        L = L if isinstance(L, float) else 1.0
        y = yb - 30
        cx = x
        if sh["lifted"]:
            if L <= -1:  # 1 + (−1) = 0: the wall
                cv.rect(cx, y, S, 16, PAPER, dash="4 3", sw=1)
                cv.line(cx, y, cx + S, y + 16, stroke=colour("b"), sw=2)
                cv.line(cx, y + 16, cx + S, y, stroke=colour("b"), sw=2)
                cv.text(cx + S / 2, y + 30, "1 + (−1) = 0", 12, INK)
                return
            cv.rect(cx, y, S, 16, UNIT)
            cv.text(cx + S / 2, y + 30, "1", 12, "#666")
            cx += S
        if L > 0:
            cv.rect(cx, y, L * S, 16, UNIT if side["name"].strip() == "1" else colour(side["key"]))
            cv.text(cx + L * S / 2, y + 30, tidy(side["name"]), 12, INK, italic=True)
        elif not sh["lifted"]:
            cv.text(cx + 20, y + 8, tidy(side["name"]), 14, INK, italic=True)
        return
    if k == "box":
        b = sh["box"]
        sides = b["sides"]
        unit = "dashed" if b["dropUnit"] else "solid"
        Ls = [length(s) for s in sides]
        if len(sides) == 1:
            L = Ls[0] if isinstance(Ls[0], float) else 1.0
            side = (1 + L) * S
            col = colour(sides[0]["key"])
            for i in range(3, 0, -1):  # stacked layers suggest "a cube of any size"
                off = i * 5
                cv.rect(x + off, yb - side - off, side, side, shade(col, 0.55 + 0.1 * i), sw=0.8)
            cv.rect(x, yb - side, side, side, shade(col, 0.15))
            cv.rect(x, yb - S, S, S, UNIT if unit == "solid" else PAPER, dash=None if unit == "solid" else "4 3", sw=0.8)
            cv.text(x + side / 2 + S / 2, yb - side / 2 - S / 4, "(1+" + tidy(sides[0]["name"]) + ")", 13,
                    text_on(shade(col, 0.15)), italic=True)
            cv.text(x + side + 20, yb - side - 14, b.get("power") or "n", 13, INK, anchor="start", weight=600)
            return
        if len(sides) == 2:
            X, Y = Ls
            X = X if isinstance(X, float) else 1.0
            kx, ky = sides[0]["key"], sides[1]["key"]
            if X <= -1:
                cv.rect(x, yb - S, S, S, PAPER, dash="4 3", sw=1)
                cv.line(x, yb - S, x + S, yb, stroke=colour("b"), sw=2)
                cv.line(x, yb, x + S, yb - S, stroke=colour("b"), sw=2)
                cv.text(x + S / 2, yb + 12, "side 1 + (−1) = 0", 11, INK)
                return
            if Y == "inv":
                w, h = (1 + X) * S, S / (1 + X)
                cv.rect(x, yb - S, S, S, PAPER, dash="4 3", sw=1)
                cv.rect(x, yb - h, w, h, shade(colour(kx + ky), 0.1), opacity=0.9)
                cv.text(x + w / 2, yb - h / 2, "area 1", 13, text_on(colour(kx + ky)))
                cv.text(x + w / 2, yb + 12, "1+" + tidy(sides[0]["name"]), 12, INK, italic=True)
                cv.text(x + S + 6, yb - S + 8, "height 1/(1+u)", 11, INK, anchor="start", italic=True)
                return
            Y = Y if isinstance(Y, float) else 1.0
            u = S
            if unit == "solid":
                cv.rect(x, yb - u, u, u, UNIT)
            else:
                cv.rect(x, yb - u, u, u, PAPER, dash="4 3", sw=1)
            if X > 0:
                cv.rect(x + u, yb - u, X * S, u, colour(kx))
                nm = tidy(sides[0]["name"])
                if len(nm) * 7 < X * S:
                    cv.text(x + u + X * S / 2, yb - u / 2, nm, 12, text_on(colour(kx)), italic=True)
                else:
                    cv.text(x + u + X * S / 2, yb + 12, nm, 12, INK, italic=True)
            if Y > 0:
                cv.rect(x, yb - u - Y * S, u, Y * S, colour(ky))
                nm = tidy(sides[1]["name"])
                if len(nm) * 7 < u:
                    cv.text(x + u / 2, yb - u - Y * S / 2, nm, 12, text_on(colour(ky)), italic=True)
                else:
                    cv.text(x + u / 2, yb - u - Y * S - 9, nm, 12, INK, italic=True)
            else:
                cv.text(x + u / 2, yb - u - 8, "side 0: flat", 11, INK)
            if X > 0 and Y > 0:
                kk = 0.5 if b.get("weight") else 1.0
                cv.rect(x + u, yb - u - Y * S, X * S, Y * S, PAPER, dash="4 3", sw=1)
                hh = Y * S * kk
                cv.rect(x + u, yb - u - hh, X * S, hh, colour(kx + ky))
            return
        if len(sides) == 3:
            names = tuple(s["key"] for s in sides)
            vals = {s["key"]: (L if isinstance(L, float) else 1.0) for s, L in zip(sides, Ls)}
            cut = b.get("cut", "")
            explode = ()
            if cut.startswith("(("):
                explode = (names[2],)
            elif cut.startswith("(0 ("):
                explode = (names[0],)
            draw_box(cv, x, yb - 4, names, S, explode=explode, gap=0.25, unit=unit, values=vals)
            return
    if k == "piece":
        sides = sh["sides"]
        Ls = [length(s) for s in sides]
        Ls = [L if isinstance(L, float) else 1.0 for L in Ls]
        kk = "".join(sorted({s["key"] for s in sides}))
        lab = max(len(tidy(q["name"])) for q in sides) * 7
        x += max(lab - Ls[0] * S, 0) / 2
        cv.rect(x, yb - Ls[1] * S, Ls[0] * S, Ls[1] * S, shade(colour(kk), 0.25))
        cv.text(x + Ls[0] * S / 2, yb + 12, tidy(sides[0]["name"]), 12, INK, italic=True)
        cv.text(x + Ls[0] * S / 2, yb - Ls[1] * S - 9, tidy(sides[1]["name"]), 12, INK, italic=True)
        return
    if k == "flip":
        cv.text(x + 6, yb - 30, "−", 22, INK)
        draw(cv, sh["inner"], x + 18, yb)
        return
    if k == "note":
        cv.rect(x + 63, yb - 72, 24, 24, INK)
        cv.text(x + 75, yb - 18, "goal closed,", 13, INK)
        cv.text(x + 75, yb - 2, "kernel-checked", 13, INK)
        return
    cv.rect(x, yb - 40, 120, 36, PAPER, dash="3 3", sw=0.8, rx=6)
    cv.text(x + 60, yb - 22, k or "?", 12, "#666")


def frame_svg(sh, label):
    if sh.get("kind") == "rel":
        (wl, hl), (wr, hr) = size(sh["lhs"]), size(sh["rhs"])
        W = max(320, int(wl + wr + 70))
        H = int(max(hl, hr) + 50)
        cv = Canvas(W, H)
        yb = H - 22
        x0 = (W - (wl + wr + 40)) / 2
        draw(cv, sh["lhs"], x0, yb)
        cv.text(x0 + wl + 20, yb - max(hl, hr) / 2, sh.get("rel", "="), 22, INK)
        draw(cv, sh["rhs"], x0 + wl + 40, yb)
    else:
        w, h = size(sh)
        W, H = max(320, int(w + 40)), int(h + 50)
        cv = Canvas(W, H)
        draw(cv, sh, (W - w) / 2, H - 22)
    cv.text(10, 14, label, 13, INK, anchor="start", weight=600)
    return cv.svg()


def main():
    base = os.path.join(HERE, "lean") if os.path.isdir(os.path.join(HERE, "lean", "build")) else HERE
    lines = open(os.path.join(base, "build", "demo_trace.txt")).read().splitlines()
    demo = open(os.path.join(base, "Demo.lean")).read().splitlines()

    def source(name):
        i = next(j for j, l in enumerate(demo) if l.startswith("theorem " + name))
        while not demo[i].rstrip().endswith(":= by"):
            i += 1
        out = []
        for l in demo[i + 1:]:
            if not l.strip():
                break
            out.append(l.strip())
        return out
    moves = [json.loads(l.split("] ", 1)[1]) for l in lines if l.startswith("[fused.pictures]")]
    shapes = [json.loads(l) for l in lines if l.startswith("{")]
    goals = shapes[-len(DEMOS):]
    films, cur = [], []
    for m in moves:
        cur.append(m)
        if m["shape"].get("kind") == "note":
            films.append(cur)
            cur = []
    assert len(films) == len(DEMOS), (len(films), len(DEMOS))
    manifest = []
    for (name, stmt, why), frames, goal in zip(DEMOS, films, goals):
        seq = [{"moves": ["goal"], "shape": goal}]
        for f in frames:
            if seq and json.dumps(seq[-1]["shape"], sort_keys=True) == json.dumps(f["shape"], sort_keys=True):
                seq[-1]["moves"].append(f["move"])
            else:
                seq.append({"moves": [f["move"]], "shape": f["shape"]})
        out = []
        for i, fr in enumerate(seq):
            label = " · ".join(fr["moves"])
            fn = f"{name}_{i}.svg"
            open(os.path.join(OUT, fn), "w").write(frame_svg(fr["shape"], label))
            cap = " ".join(CAPTION.get(m, "") for m in fr["moves"])
            out.append({"file": "figures/film/" + fn, "moves": fr["moves"], "caption": cap})
        tactic = source(name)
        manifest.append({"name": name, "statement": stmt, "why": why, "frames": out, "tactics": tactic})
    json.dump(manifest, open(os.path.join(OUT, "film.json"), "w"), indent=1, ensure_ascii=False)
    print(sum(len(m["frames"]) for m in manifest), "frames for", len(manifest), "proofs")


if __name__ == "__main__":
    main()
