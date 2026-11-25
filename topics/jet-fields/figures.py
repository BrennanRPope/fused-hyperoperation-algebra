"""Draws the two SVG figures for the jet-fields page (writes page/cube.svg, page/secant.svg)."""
import math, os
OUT = os.path.join(os.path.dirname(__file__), 'page')
CS = 0.62*math.cos(math.pi/4)
def P(x, y, z, S=150, ox=110, oy=440):
    return (ox + (x + CS*y)*S, oy - (z + CS*y)*S)
ha, hb, hc, g = 0.46, 0.40, 0.34, 0.30
cls = {0:'l0', 1:'l1', 2:'l2', 3:'l3'}
boxes = []
for bx in (0,1):
    for by in (0,1):
        for bz in (0,1):
            x0 = 0 if not bx else 1 + g; x1 = 1 if not bx else 1 + g + ha
            y0 = 0 if not by else 1 + g; y1 = 1 if not by else 1 + g + hb
            z0 = 0 if not bz else 1 + g; z1 = 1 if not bz else 1 + g + hc
            boxes.append((bx+by+bz, (x0,x1,y0,y1,z0,z1)))
# farthest first: viewer at (+x, -y, +z)
boxes.sort(key=lambda b: -(b[1][2] - b[1][0] - b[1][4]))
def poly(pts, c):
    return '<polygon class="%s" points="%s"/>' % (c, ' '.join('%.1f,%.1f' % p for p in pts))
parts = []
for layer, (x0,x1,y0,y1,z0,z1) in boxes:
    c = cls[layer]
    parts.append(poly([P(x0,y0,z1),P(x1,y0,z1),P(x1,y1,z1),P(x0,y1,z1)], c+' top'))
    parts.append(poly([P(x1,y0,z0),P(x1,y1,z0),P(x1,y1,z1),P(x1,y0,z1)], c+' rt'))
    parts.append(poly([P(x0,y0,z0),P(x1,y0,z0),P(x1,y0,z1),P(x0,y0,z1)], c+' lf'))
def lab(x, y, z, text, anchor='start', dx=0, dy=0, c='lab'):
    px, py = P(x, y, z)
    return '<text class="%s" x="%.1f" y="%.1f" text-anchor="%s">%s</text>' % (c, px+dx, py+dy, anchor, text)
G = 1 + g
RX = P(G + ha, G + hb, 0)[0] + 40
def lead(pt, text, sub, cl, dy=0):
    px, py = P(*pt)
    ty = py + dy
    return ('<path class="ldr" d="M%.1f,%.1f L%.1f,%.1f L%.1f,%.1f"/>' % (px + 4, py, RX - 30, ty, RX - 8, ty)
            + '<circle class="dot" cx="%.1f" cy="%.1f" r="3.5"/>' % (px + 4, py)
            + '<text class="lab %s" x="%.1f" y="%.1f">%s</text>' % (cl, RX, ty + 2, text)
            + '<text class="sub" x="%.1f" y="%.1f">%s</text>' % (RX, ty + 20, sub))
labels = [
  lab(0.5, 0, 0.5, '1', 'middle', 0, 2, 'lab big'),
  lab(0.5, 0, 0.5, 'old volume', 'middle', 0, 24, 'sub'),
  lead((G + ha, G + hb/2, G + hc/2), 'R · h³', 'far corner: stretching', 'l3t', -40),
  lead((G + ha, 0.5, G + hc/2), 'Q · h²', 'three rods: spin vs strain', 'l2t', 10),
  lead((G + ha, 0.5, 0.5), 'P · h', 'three slabs: divergence', 'l1t', 60),
]
W = RX + 210
svg = ('<svg viewBox="95 20 %d 440"' % (W - 95) + '  xmlns="http://www.w3.org/2000/svg" role="img" '
       'aria-label="Exploded box: unit cube, three slabs, three rods and one far corner">'
       + ''.join(parts) + ''.join(labels) + '</svg>')
open(os.path.join(OUT, 'cube.svg'), 'w').write(svg)

# secant -> tangent
W, H = 560, 300
f = lambda x: 0.35*x*x + 0.25
X = lambda x: 40 + x*120; Y = lambda y: H - 30 - y*80
x0 = 1.0
path = 'M' + ' L'.join('%.1f,%.1f' % (X(t/50), Y(f(t/50))) for t in range(0, 141))
els = ['<path class="curve" d="%s"/>' % path]
for k, op in ((1.2, .35), (0.7, .6), (0.3, .85)):
    s = (f(x0+k) - f(x0))/k
    xa, xb = x0 - 0.7, x0 + k + 0.5
    els.append('<line class="sec" style="opacity:%s" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>' % (op, X(xa), Y(f(x0)+s*(xa-x0)), X(xb), Y(f(x0)+s*(xb-x0))))
    els.append('<circle class="pt2" style="opacity:%s" cx="%.1f" cy="%.1f" r="5"/>' % (op, X(x0+k), Y(f(x0+k))))
    els.append('<text class="sub" style="opacity:%s" x="%.1f" y="%.1f">k = %s</text>' % (op, X(x0+k)+9, Y(f(x0+k))+16, k))
s = 2*0.35*x0
xa, xb = x0 - 0.8, x0 + 1.3
els.append('<line class="tan" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>' % (X(xa), Y(f(x0)+s*(xa-x0)), X(xb), Y(f(x0)+s*(xb-x0))))
els.append('<circle class="pt" cx="%.1f" cy="%.1f" r="6"/>' % (X(x0), Y(f(x0))))
els.append('<text class="lab" x="%.1f" y="%.1f">x</text>' % (X(x0)-6, Y(f(x0))+26))
els.append('<text class="lab tant" x="%.1f" y="%.1f">k = 0: the tangent</text>' % (X(2.25), Y(0.62)))
els.append('<text class="sub" x="%.1f" y="%.1f">ε² = 0, corner dropped</text>' % (X(2.25), Y(0.62)+20))
svg2 = ('<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" '
        'aria-label="Secants through x and x+k turning into the tangent as k shrinks">' % (W, H) + ''.join(els) + '</svg>')
open(os.path.join(OUT, 'secant.svg'), 'w').write(svg2)
print('ok')
