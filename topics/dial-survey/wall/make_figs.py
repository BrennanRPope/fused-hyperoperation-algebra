"""Draw the wall / Lee-Yang pictures as SVG.
Writes standalone files to figs/*.svg (own colours) and figs/*.inline.svg (page classes, no <style>)."""
import math, json, os
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figs")
STYLE = """<style>
.fi{fill:#06222e}.fm{fill:#3d5d69}.fb{fill:#dff1f7}.fk{fill:#ff2e7e}.si{stroke:#06222e}.sl{stroke:#9cc6d4}.sp{stroke:#ff2e7e}.sm{stroke:#3d5d69}
@media (prefers-color-scheme:dark){.fi{fill:#dff1f7}.fm{fill:#9fc3cf}.fb{fill:#06222e}.si{stroke:#dff1f7}.sl{stroke:#24505f}.sm{stroke:#9fc3cf}}
text{font-family:Inter,system-ui,sans-serif}</style><rect width="100%" height="100%" class="fb"/>"""
def wrap(w, h, body, label):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="{label}">'
            + "%STYLE%" + body + "</svg>")
def save(name, svg):
    open(os.path.join(OUT, name + ".svg"), "w").write(svg.replace("%STYLE%", STYLE))
    open(os.path.join(OUT, name + ".inline.svg"), "w").write(svg.replace("%STYLE%", ""))
def T(x, y, s, cls="fi", size=13, anchor="start", weight=None):
    w = f' font-weight="{weight}"' if weight else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" font-size="{size}" text-anchor="{anchor}"{w}>{s}</text>'
def L(x1, y1, x2, y2, cls="si", w=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="{cls}" stroke-width="{w}"{d}/>'
def C(x, y, r, cls="fi", extra=""):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" class="{cls}" {extra}/>'

# ---------- Figure A: the k line ----------
W, H = 640, 360
k0, k1, X0, X1 = -1.35, 3.3, 30, 612
X = lambda k: X0 + (k - k0) / (k1 - k0) * (X1 - X0)
AY = 196
b = []
b.append(T(20, 26, "Where the walls sit, and where the Lee–Yang points sit, for each k", size=15, weight=600))
# legend
b.append(L(24, 58, 52, 58, "sp", 7)); b.append(T(60, 62, "wall at Re s = 0", size=12))
b.append(f'<path d="M{190} 54 l-5 -9 h10 z" class="fk" transform="translate(0,13)"/>'); b.append(T(202, 62, "Lee–Yang points, k = 1 + 2cos(2π/n)", size=12))
b.append('<path d="M450 53 l-5 10 h10 z" class="fi"/>'); b.append(T(462, 62, "at most m copies, k = −1/m", size=12))
# wall band
b.append(L(X(k0), AY, X(k1), AY, "sp", 7))
# axis ticks and numbers
for kk, lab in [(-1, "−1"), (-0.5, "−1/2"), (0, "0"), (1, "1"), (2, "2"), (3, "3")]:
    b.append(T(X(kk), AY + 44, lab, "fm", 12, "middle"))
# copies points below
for m in range(1, 14):
    x = X(-1 / m)
    b.append(f'<path d="M{x:.1f} {AY+14} l-4 9 h8 z" class="fi"/>')
b.append(L(X(-1/3), AY + 54, X(0), AY + 54, "sm", 1))
b.append(T(X(-1/6), AY + 70, "−1/3, −1/4, … crowd to 0", "fm", 11.5, "middle"))
# Lee-Yang points above
ber = [(n, 1 + 2 * math.cos(2 * math.pi / n)) for n in range(2, 41)]
for n, kk in ber:
    b.append(f'<path d="M{X(kk):.1f} {AY-14} l-4 -9 h8 z" class="fk"/>')
for n, kk in ber[:6]:
    b.append(T(X(kk), AY - 30, f"n={n}", "fm", 11, "middle"))
b.append(L(X(3), AY - 60, X(3), AY + 12, "sm", 1, "3 3"))
b.append(T(X(3), AY - 66, "crowd to 3", "fk", 11.5, "middle"))
# no-wall points k = +-1, ring marks where the families meet
for kk in [-1, 0, 1]:
    b.append(f'<circle cx="{X(kk):.1f}" cy="{AY}" r="15" class="sm" fill="none" stroke-width="1.5" stroke-dasharray="3 2"/>')
for kk in [-1, 1]:
    b.append(f'<circle cx="{X(kk):.1f}" cy="{AY}" r="8" class="fb si" stroke-width="3"/>')
b.append(T(X(0), AY - 84, "k = 0: wall known", size=12, anchor="middle"))
b.append(T(X(0), AY - 70, "since 1920", size=12, anchor="middle"))
b.append(L(X(0), AY - 64, X(0), AY - 42, "sm", 1, "2 2"))
b.append(L(X(-1), AY + 50, X(-1), AY + 92, "sm", 1, "2 2"))
b.append(T(X(-1) - 6, AY + 106, "k = −1: no wall (squarefree, ζ(s)/ζ(2s))", size=12))
b.append(L(X(1), AY + 50, X(1), AY + 76, "sm", 1, "2 2"))
b.append(T(X(1) - 6, AY + 90, "k = +1: no wall (all whole numbers, ζ(s))", size=12))
b.append(T(20, H - 26, "Lee–Yang points hit the special k only at −1, 0, 1 (dashed rings). No copies point −1/2, −1/3, … is one.", "fm", 11.5))
b.append(T(20, H - 10, "Past |k| = 1 a few prime factors add branch points further right; the wall stays at Re s = 0.", "fm", 11.5))
figA = wrap(W, H, "".join(b), "A line of k values. A pink wall runs along it except at k = minus 1 and plus 1. Copies points minus 1 over m sit below, crowding to 0; Lee–Yang points 1 + 2cos(2 pi over n) sit above, crowding to 3. They meet only at minus 1, 0 and 1.")
save("walls_by_k", figA)

# ---------- Figure B: singular points rho/N piling onto the imaginary axis ----------
mp.mp.dps = 15
zeros = [float(mp.zetazero(n).imag) for n in range(1, 141)]
W, H = 640, 420
TOP = 25
x0, x1, y0, y1 = 70, 600, 340, 70   # Re s from 0 to 0.55 ; Im s from 0 to TOP
RX = lambda r: x0 + r / 0.55 * (x1 - x0)
IY = lambda t: y0 - t / TOP * (y0 - y1)
b = [T(20, 26, "The singular points ρ/N pile onto the line Re s = 0", size=15, weight=600)]
b.append(L(RX(0), y0, RX(0), y1 - 6, "sp", 5))
b.append(L(RX(0), y0, x1, y0, "sl", 1))
for t in range(0, TOP + 1, 5):
    b.append(T(RX(0) - 10, IY(t) + 4, f"{t}i" if t else "0", "fm", 11, "end"))
for N in [1, 2, 3, 4, 6, 12]:
    b.append(L(RX(0.5 / N), y0, RX(0.5 / N), y0 + 5, "sm", 1))
    b.append(T(RX(0.5 / N), y0 + 19, f"N={N}", "fm", 11, "middle"))
b.append(T(x1, y0 + 38, "column N sits at Re s = 1/(2N) →", "fm", 11, "end"))
pts = []
for N in range(1, 13):
    for g in zeros:
        t = g / N
        if t > TOP: break
        pts.append(C(RX(0.5 / N), IY(t), 2.4, "fi"))
b += pts
b.append(T(RX(0) + 10, y1 - 2, "wall", "fk", 12))
b.append(T(RX(0.5) - 8, IY(14.13) + 4, "ζ's first zero, 1/2 + 14.13i", "fm", 11, "end"))
b.append(T(20, H - 26, "Column N holds the zeros of ζ(Ns): ζ's zeros, squashed by N toward 0. Shown up to N = 12;", "fm", 11.5))
b.append(T(20, H - 10, "as N grows the columns keep closing in, and they fill the line Re s = 0 with no gaps.", "fm", 11.5))
figB = wrap(W, H, "".join(b), "Dots in a strip of the s-plane. Column N sits at real part 1 over 2N and holds the zeta zeros divided by N. The columns crowd onto the imaginary axis, drawn in pink as the wall.")
save("points_pile_up", figB)

# ---------- Figure C: how strong the singularity is, as a function of k ----------
def M(k, N):
    def mu(n):
        r, p, m = 1, 2, n
        while p * p <= m:
            if m % p == 0:
                m //= p
                if m % p == 0: return 0
                r = -r
            p += 1
        return -r if m > 1 else r
    return sum(mu(N // d) * k ** d for d in range(1, N + 1) if N % d == 0) / N
W, H = 640, 390
ka, kb, x0, x1 = -1.6, 1.6, 70, 600
ya, yb, y0, y1 = -1.6, 1.6, 320, 66
KX = lambda k: x0 + (k - ka) / (kb - ka) * (x1 - x0)
VY = lambda v: y0 - (v - ya) / (yb - ya) * (y0 - y1)
b = [T(20, 24, "Strength of the singularity at ρ/N, as k turns", size=15, weight=600)]
for v in [-1, 0, 1]:
    b.append(L(x0, VY(v), x1, VY(v), "sl", 1, None if v == 0 else "3 3"))
    b.append(T(x0 - 8, VY(v) + 4, f"{v}".replace("-", "−"), "fm", 11, "end"))
for kk in [-1, 1]:
    b.append(L(KX(kk), y0, KX(kk), y1, "sm", 1, "4 3"))
for kk in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5]:
    b.append(T(KX(kk), y0 + 18, f"{kk:g}".replace("-", "−"), "fm", 11, "middle"))
b.append(T(x1, y0 + 34, "k →", "fm", 11, "end"))
cls = {2: "sp", 3: "si", 4: "sm", 5: "si", 6: "sm"}
dash = {2: None, 3: None, 4: None, 5: "5 3", 6: "5 3"}
for N in [2, 3, 4, 5, 6]:
    d = []
    for i in range(321):
        kk = ka + (kb - ka) * i / 320
        if abs(kk) < 1e-9: kk = 1e-9
        v = M(kk, N) / kk
        v = max(min(v, 1.75), -1.75)
        d.append(f"{'M' if i == 0 else 'L'}{KX(kk):.1f} {VY(v):.1f}")
    ds = f' stroke-dasharray="{dash[N]}"' if dash[N] else ""
    b.append(f'<path d="{" ".join(d)}" fill="none" class="{cls[N]}" stroke-width="{2.5 if N == 2 else 1.6}"{ds}/>')
    if N in (2, 4, 6):
        lv = M(kb, N) / kb
        b.append(T(KX(kb) + 5, VY(lv) + 4, f"N={N}", "fk" if N == 2 else "fi", 11, "start"))
    else:
        lk = -1.56; lv = M(lk, N) / lk
        b.append(T(KX(lk) + 8, VY(lv) - 2, f"N={N}", "fi", 11, "start"))
for kk, v in [(1, 0), (-1, 0), (-1, -1)]:
    b.append(C(KX(kk), VY(v), 4.5, "fi"))
b.append(T(KX(1), y1 - 12, "k = +1: every curve is 0", "fi", 12, "middle"))
b.append(T(KX(-1), y1 - 12, "k = −1: 0 from N = 3 on", "fi", 12, "middle"))
b.append(T(KX(-1) + 10, VY(-1) + 18, "N = 2 gives −1: a plain pole (the ζ(2s) in ζ(s)/ζ(2s))", "fi", 12))
b.append(T(20, H - 26, "The exponent at ρ/N (one ζ zero there) is M(k,N)/k. A whole number means a pole or a zero;", "fm", 11.5))
b.append(T(20, H - 10, "any other value means a branch point. Only at k = ±1 do all but finitely many vanish.", "fm", 11.5))
figC = wrap(W, H, "".join(b), "Curves of the exponent M of k and N divided by k for N from 2 to 6, plotted against k. At k = 1 every curve is zero; at k = minus 1 every curve from N = 3 on is zero and N = 2 is minus 1.")
save("strength_by_k", figC)
print("wrote", sorted(os.listdir(OUT)))
