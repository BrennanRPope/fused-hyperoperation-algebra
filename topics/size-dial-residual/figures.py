"""Draws the three SVG figures for the page (light/dark aware via CSS variables when inlined;
standalone files use fixed colours)."""
import mpmath as mp
BLUE, ORANGE, PURPLE, GRAY, INK, SUB = "#1d4fd8", "#f08a24", "#8a3ffc", "#b9bcc4", "#14161c", "#5b5f6b"
FONT = 'font-family="Inter,system-ui,sans-serif"'
def share(beta, t=0.5):
    if beta == 0: return float(mp.log(1+t)/mp.log(2))
    return float(mp.re(mp.polylog(beta+1, -t)/mp.polylog(beta+1, -1)))

# ---------- Figure 1: the dimension dial ----------
W, H, L, R, T, B = 760, 380, 70, 30, 40, 70
x = lambda b: L + (W-L-R)*(b/(b+1))           # beta in [0,inf) squeezed to [0,1)
y = lambda f: T + (H-T-B)*(0.59 - f)/(0.59-0.495)
pts = [i/400*60 for i in range(401)]
path = " ".join(f"{'M' if i==0 else 'L'}{x(b):.1f},{y(share(b)):.1f}" for i, b in enumerate(pts))
marks = [(0, "finite set", "sharp count"), (1, "naturals by size", ""), (2, "lattice in a plane", ""),
         (3, "3-D radiation", "")]
s = [f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="f1t"><title id="f1t">Turning every dial to one half gives a share that depends only on the dimension of the size measure: 58.5% for a sharp count, 54.5% for the naturals by size, falling to exactly 50% for the primes and every other set measured by logarithmic size.</title>']
for f in [0.50, 0.52, 0.54, 0.56, 0.58]:
    s.append(f'<line x1="{L}" x2="{W-R}" y1="{y(f):.1f}" y2="{y(f):.1f}" stroke="{GRAY}" stroke-width="0.6"/>'
             f'<text x="{L-8}" y="{y(f)+4:.1f}" {FONT} font-size="13" fill="{SUB}" text-anchor="end">{f*100:.0f}%</text>')
s.append(f'<line x1="{L}" x2="{W-R}" y1="{y(0.5):.1f}" y2="{y(0.5):.1f}" stroke="{INK}" stroke-width="1" stroke-dasharray="4 4"/>')
s.append(f'<path d="{path}" fill="none" stroke="{BLUE}" stroke-width="3.5"/>')
for b, lab, sub in marks:
    cx, cy = x(b), y(share(b))
    dx = 14 if b == 0 else 10
    s.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="6" fill="{BLUE}"/>'
             f'<text x="{cx+dx:.1f}" y="{cy-10-(16 if b==0 else 0):.1f}" {FONT} font-size="14" font-weight="600" fill="{INK}">{lab}</text>'
             f'<text x="{cx+dx:.1f}" y="{cy+8-(16 if b==0 else 0):.1f}" {FONT} font-size="13" fill="{SUB}">{share(b)*100:.1f}%{"  " + sub if sub else ""}</text>')
s.append(f'<circle cx="{W-R:.1f}" cy="{y(0.5):.1f}" r="6" fill="{ORANGE}"/>'
         f'<text x="{W-R-8}" y="{y(0.5)+22:.1f}" {FONT} font-size="14" font-weight="600" fill="{INK}" text-anchor="end">primes, integers, rationals, number fields, random reals: exactly 50%</text>')
for b, lab in [(0, "0"), (1, "1"), (2, "2"), (3, "3"), (5, "5"), (10, "10")]:
    s.append(f'<text x="{x(b):.1f}" y="{H-B+22}" {FONT} font-size="13" fill="{SUB}" text-anchor="middle">{lab}</text>')
s.append(f'<text x="{W-R}" y="{H-B+22}" {FONT} font-size="13" fill="{SUB}" text-anchor="middle">∞</text>')
s.append(f'<text x="{(L+W-R)/2}" y="{H-14}" {FONT} font-size="14" fill="{INK}" text-anchor="middle">dimension β of the size measure (how fast the count of elements grows with size)</text>')
s.append(f'<text x="{W-R}" y="{T-16}" {FONT} font-size="14" fill="{INK}" text-anchor="end">share of the total when every dial is ½</text></svg>')
open("fig1-dimension-dial.svg", "w").write("\n".join(s))

# ---------- Figure 2: explained vs left over, per set ----------
rows = [  # label, measure, [(share, kind)], note
 ("Primes, reals, algebraic numbers", "as sets", [(1, "e")], "nothing left over"),
 ("Naturals by divisibility (ζ)", "growth", [(1, "e")], "all of it"),
 ("", "count", [(0.6079, "e"), (0.3921, "r")], "6/π² explained"),
 ("Integers", "count", [(0.6079, "e"), (0.3921, "r")], "same as the naturals"),
 ("Gaussian integers Z[i]", "count", [(0.6637, "e"), (0.3363, "r")], "1/ζ_K(2)"),
 ("Rationals by height", "count", [(0.4717, "e"), (0.5283, "r")], "repeats + p·(1/p) = 1"),
 ("Z[√−5] numbers", "count", [(0.5315*0.5, "e"), (0.4685*0.5, "r"), (0.5, "c")], "class group takes half"),
 ("Naturals by size (partitions)", "log-total", [(0.5, "e"), (0.5, "r")], "β = 1"),
 ("Lattice in a plane", "log-total", [(0.75, "e"), (0.25, "r")], "β = 2"),
 ("3-D radiation", "log-total", [(0.875, "e"), (0.125, "r")], "β = 3: the 7/8 law"),
 ("Algebraic numbers, multiplied", "", [(1, "x")], "no atoms: can't start"),
]
col = {"e": BLUE, "r": ORANGE, "c": PURPLE, "x": GRAY}
W2, rowh, top, lab_w, bar_w = 900, 34, 84, 400, 300
H2 = top + rowh*len(rows) + 50
s = [f'<svg viewBox="0 0 {W2} {H2}" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="f2t"><title id="f2t">For each set, the share that mult-add with dials explains (blue) and what is left over: repetition (orange), the class group (purple), or no foothold at all (grey).</title>']
legend = [("explained by mult-add", BLUE), ("repetition: the set again at double size", ORANGE), ("class group (sign dial)", PURPLE), ("out of reach", GRAY)]
for j, (t_, c) in enumerate(legend):
    lx = 10 + (j % 2)*420; ly = 14 + (j // 2)*24
    s.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{c}"/><text x="{lx+20}" y="{ly+12}" {FONT} font-size="13" fill="{INK}">{t_}</text>')
for i, (lab, meas, segs, note) in enumerate(rows):
    yy = top + i*rowh
    s.append(f'<text x="{lab_w-84}" y="{yy+17}" {FONT} font-size="14" font-weight="600" fill="{INK}" text-anchor="end">{lab}</text>'
             f'<text x="{lab_w-8}" y="{yy+17}" {FONT} font-size="13" fill="{SUB}" text-anchor="end">{meas}</text>')
    xx = lab_w
    for f, k in segs:
        s.append(f'<rect x="{xx:.1f}" y="{yy+4}" width="{bar_w*f:.1f}" height="20" fill="{col[k]}"/>')
        xx += bar_w*f
    e = sum(f for f, k in segs if k == "e")
    if segs[0][1] == "e" and e < 1:
        s.append(f'<text x="{lab_w+6}" y="{yy+19}" {FONT} font-size="13" font-weight="600" fill="#fff">{e*100:.1f}%</text>')
    s.append(f'<text x="{lab_w+bar_w+10}" y="{yy+19}" {FONT} font-size="13" fill="{SUB}">{note}</text>')
s.append(f'<text x="10" y="{H2-14}" {FONT} font-size="13" fill="{SUB}">growth: exponent of how fast the count grows. count: the leading constant. log-total: log of the weighted total at a fixed size scale.</text></svg>')
open("fig2-explained-residual.svg", "w").write("\n".join(s))

# ---------- Figure 3: the doubling ladder ----------
W3, H3 = 760, 286
s = [f'<svg viewBox="0 0 {W3} {H3}" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="f3t"><title id="f3t">The leftover is the same set at double size, which is itself mostly explained, leaving a smaller leftover at four times the size, and so on. In dimension 1 the pieces are one half, one quarter, one eighth; in dimension 2 three quarters, three sixteenths and so on.</title>']
for r_, (beta, lab) in enumerate([(1, "dimension 1 (naturals by size)"), (2, "dimension 2 (lattice in a plane)")]):
    y0 = 40 + r_*130
    s.append(f'<text x="10" y="{y0-10}" {FONT} font-size="14" font-weight="600" fill="{INK}">{lab}</text>')
    xx = 10; total = 740; rem = 1.0; j = 0
    while rem > 0.004 and j < 8:
        piece = rem*(1-2**-beta)
        wpx = total*piece
        shade = ["#1d4fd8", "#3f6fe6", "#6a90ee", "#97b3f4", "#c0d1f8", "#dbe5fb", "#e9effc", "#f2f6fd"][j]
        s.append(f'<rect x="{xx:.1f}" y="{y0}" width="{wpx:.1f}" height="44" fill="{shade}" stroke="#fff" stroke-width="1"/>')
        if wpx > 60:
            s.append(f'<text x="{xx+8:.1f}" y="{y0+20}" {FONT} font-size="13" font-weight="600" fill="{"#fff" if j<3 else INK}">size ×{2**j}</text>'
                     f'<text x="{xx+8:.1f}" y="{y0+37}" {FONT} font-size="13" fill="{"#fff" if j<3 else INK}">{piece*100:.1f}%</text>')
        xx += wpx; rem -= piece; j += 1
    s.append(f'<text x="10" y="{y0+66}" {FONT} font-size="13" fill="{SUB}">each rung is mult-add on dials t, t², t⁴, …; the first rung alone explains 1 − 1/{2**beta} = {(1-2**-beta)*100:.1f}%</text>')
s.append('</svg>')
open("fig3-doubling-ladder.svg", "w").write("\n".join(s))
print("wrote 3 figures")
