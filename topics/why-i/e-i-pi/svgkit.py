"""Tiny SVG plotting kit (numpy only), for pages where matplotlib is not installed."""
import numpy as np
from xml.sax.saxutils import escape
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY, DARK, AMBER = "#2a6fdb", "#e07b1a", "#2e9e5b", "#d0443a", "#7b4fc9", "#8a8a8a", "#222222", "#d9a400"

class Fig:
    def __init__(self, w, h):
        self.w, self.h, self.el, self.ops = w, h, [], []
        self.el.append(f'<rect x="0" y="0" width="{w}" height="{h}" fill="white"/>')
    def add(self, s): self.el.append(s)
    def text(self, x, y, s, size=14, anchor="start", color=DARK, weight="normal", italic=False, rotate=None):
        st = ' font-style="italic"' if italic else ''
        rt = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate is not None else ''
        self.ops.append(("text", x, y, s, size, anchor, color, weight, rotate))
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="{size}" '
                 f'text-anchor="{anchor}" fill="{color}" font-weight="{weight}"{st}{rt}>{escape(s)}</text>')
    def lines(self, x, y, s, size=13, dy=None, **kw):
        dy = dy or size*1.3
        for i, line in enumerate(s.split("\n")): self.text(x, y + i*dy, line, size=size, **kw)
    def line(self, x1, y1, x2, y2, color=DARK, width=1.5, dash=None, arrow=False, opacity=1):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        m = f' marker-end="url(#ah{color[1:]})"' if arrow else ''
        if arrow: self._marker(color)
        self.ops.append(("poly", [(x1, y1), (x2, y2)], color, width, dash, "none", opacity, False, arrow))
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}"{d}{m} opacity="{opacity}"/>')
    def _marker(self, color):
        key = f'id="ah{color[1:]}"'
        if any(key in e for e in self.el): return
        self.el.insert(1, f'<defs><marker {key} viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
                          f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker></defs>')
    def poly(self, pts, color=DARK, width=2, dash=None, fill="none", opacity=1, close=False, arrow=False):
        if len(pts) < 2: return
        d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + (" Z" if close else "")
        ds = f' stroke-dasharray="{dash}"' if dash else ''
        m = f' marker-end="url(#ah{color[1:]})"' if arrow else ''
        if arrow: self._marker(color)
        self.ops.append(("poly", list(pts), color, width, dash, fill, opacity, close, arrow))
        self.add(f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="{fill}" opacity="{opacity}"{ds}{m} stroke-linejoin="round"/>')
    def circle(self, x, y, r, color=DARK, fill="none", width=1.5, opacity=1, dash=None):
        ds = f' stroke-dasharray="{dash}"' if dash else ''
        self.ops.append(("ellipse", x, y, r, r, color, fill, width, opacity))
        self.add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" stroke="{color}" fill="{fill}" stroke-width="{width}" opacity="{opacity}"{ds}/>')
    def ellipse(self, x, y, rx, ry, color=DARK, fill="none", width=1.5, opacity=1):
        self.ops.append(("ellipse", x, y, rx, ry, color, fill, width, opacity))
        self.add(f'<ellipse cx="{x:.1f}" cy="{y:.1f}" rx="{rx:.1f}" ry="{ry:.1f}" stroke="{color}" fill="{fill}" stroke-width="{width}" opacity="{opacity}"/>')
    def rect(self, x, y, w, h, color=DARK, fill="none", width=1.2, rx=0, opacity=1):
        self.ops.append(("rect", x, y, w, h, color, fill, width, rx, opacity))
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" stroke="{color}" fill="{fill}" stroke-width="{width}" opacity="{opacity}"/>')
    def title(self, s, size=18): self.text(self.w/2, 28, s, size=size, anchor="middle", weight="bold")
    def save(self, path):
        with open(path, "w") as f:
            f.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" height="{self.h}" viewBox="0 0 {self.w} {self.h}">\n')
            f.write("\n".join(self.el)); f.write("\n</svg>\n")
        self.png(path[:-4] + ".png")
    def png(self, path, S=2):
        from PIL import Image, ImageDraw, ImageFont
        img = Image.new("RGBA", (self.w*S, self.h*S), "white")
        FD = "/usr/share/fonts/truetype/dejavu/"
        fonts = {}
        def font(size, bold):
            key = (size, bold)
            if key not in fonts: fonts[key] = ImageFont.truetype(FD + ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), int(size*S*0.92))
            return fonts[key]
        def rgba(c, op):
            if c in (None, "none"): return None
            c = c.lstrip("#")
            if len(c) == 3: c = "".join(ch*2 for ch in c)
            r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
            if c == "ffffff" or c == "white": pass
            return (r, g, b, int(255*op))
        named = {"white": "#ffffff"}
        def layer(): return Image.new("RGBA", img.size, (0, 0, 0, 0))
        for op in self.ops:
            L = layer(); d = ImageDraw.Draw(L)
            if op[0] == "poly":
                _, pts, color, width, dash, fill, opac, close, arrow = op
                color = named.get(color, color); fill = named.get(fill, fill)
                P = [(x*S, y*S) for x, y in pts]
                if fill not in (None, "none"): d.polygon(P, fill=rgba(fill, opac))
                if color not in (None, "none"):
                    cc = rgba(color, opac); ww = max(1, int(round(width*S)))
                    seq = P + ([P[0]] if close else [])
                    if dash:
                        on, off = [float(v)*S for v in dash.split(",")]
                        for (xa, ya), (xb, yb) in zip(seq[:-1], seq[1:]):
                            Ls = np.hypot(xb - xa, yb - ya); t = 0
                            while t < Ls:
                                t2 = min(t + on, Ls)
                                d.line([(xa + (xb-xa)*t/Ls, ya + (yb-ya)*t/Ls), (xa + (xb-xa)*t2/Ls, ya + (yb-ya)*t2/Ls)], fill=cc, width=ww)
                                t = t2 + off
                    else:
                        d.line(seq, fill=cc, width=ww, joint="curve")
                    if arrow and len(seq) >= 2:
                        (xa, ya), (xb, yb) = seq[-2], seq[-1]; ang = np.arctan2(yb - ya, xb - xa); hl = 9*S
                        d.polygon([(xb, yb), (xb - hl*np.cos(ang - 0.4), yb - hl*np.sin(ang - 0.4)), (xb - hl*np.cos(ang + 0.4), yb - hl*np.sin(ang + 0.4))], fill=cc)
            elif op[0] == "ellipse":
                _, x, y, rx, ry, color, fill, width, opac = op
                color = named.get(color, color); fill = named.get(fill, fill)
                box = [(x - rx)*S, (y - ry)*S, (x + rx)*S, (y + ry)*S]
                d.ellipse(box, fill=rgba(fill, opac) if fill not in (None, "none") else None,
                          outline=rgba(color, opac) if color not in (None, "none") else None, width=max(1, int(round(width*S))))
            elif op[0] == "rect":
                _, x, y, w, h, color, fill, width, rx, opac = op
                color = named.get(color, color); fill = named.get(fill, fill)
                box = [x*S, y*S, (x + w)*S, (y + h)*S]
                d.rounded_rectangle(box, radius=rx*S, fill=rgba(fill, opac) if fill not in (None, "none") else None,
                                    outline=rgba(color, opac) if color not in (None, "none") else None, width=max(1, int(round(width*S))))
            elif op[0] == "text":
                _, x, y, s_, size, anchor, color, weight, rot = op
                fnt = font(size, weight == "bold"); tw = fnt.getlength(s_)
                asc = fnt.getmetrics()[0]
                if rot is None:
                    xx = x*S - {"start": 0, "middle": tw/2, "end": tw}[anchor]
                    d.text((xx, y*S - asc), s_, font=fnt, fill=rgba(named.get(color, color), 1))
                else:
                    T = Image.new("RGBA", (int(tw) + 4, int(size*S*1.5)), (0, 0, 0, 0))
                    ImageDraw.Draw(T).text((2, 0), s_, font=fnt, fill=rgba(named.get(color, color), 1))
                    T = T.rotate(-rot, expand=True)
                    L.paste(T, (int(x*S - T.width/2), int(y*S - T.height/2)), T)
            img = Image.alpha_composite(img, L)
        img.convert("RGB").save(path)

class Axes:
    def __init__(self, fig, x, y, w, h, xlim, ylim, equal=False):
        self.f, self.x, self.y, self.w, self.h = fig, x, y, w, h
        self.xlim, self.ylim = xlim, ylim
    def X(self, v): return self.x + (np.asarray(v) - self.xlim[0])/(self.xlim[1] - self.xlim[0])*self.w
    def Y(self, v): return self.y + self.h - (np.asarray(v) - self.ylim[0])/(self.ylim[1] - self.ylim[0])*self.h
    def frame(self, xticks=(), yticks=(), xlabel="", ylabel="", grid=False, fmt="{:g}", box=True):
        f = self.f
        if box: f.rect(self.x, self.y, self.w, self.h, color="#999", width=1)
        for t in xticks:
            lab = t[1] if isinstance(t, tuple) else fmt.format(t); t = t[0] if isinstance(t, tuple) else t
            X = self.X(t); f.line(X, self.y + self.h, X, self.y + self.h + 5, color="#999", width=1)
            if grid: f.line(X, self.y, X, self.y + self.h, color="#eee", width=1)
            f.text(X, self.y + self.h + 18, lab, size=12, anchor="middle", color="#444")
        for t in yticks:
            lab = t[1] if isinstance(t, tuple) else fmt.format(t); t = t[0] if isinstance(t, tuple) else t
            Y = self.Y(t); f.line(self.x - 5, Y, self.x, Y, color="#999", width=1)
            if grid: f.line(self.x, Y, self.x + self.w, Y, color="#eee", width=1)
            f.text(self.x - 8, Y + 4, lab, size=12, anchor="end", color="#444")
        if xlabel: f.text(self.x + self.w/2, self.y + self.h + 38, xlabel, size=13, anchor="middle")
        if ylabel: f.text(self.x - 46, self.y + self.h/2, ylabel, size=13, anchor="middle", rotate=-90)
    def plot(self, xs, ys, color=BLUE, width=2, dash=None, clip=True, fill="none", opacity=1, arrow=False, close=False):
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        if clip:
            ys = np.clip(ys, self.ylim[0] - 0.02*(self.ylim[1]-self.ylim[0]), self.ylim[1] + 0.02*(self.ylim[1]-self.ylim[0]))
        self.f.poly(list(zip(self.X(xs), self.Y(ys))), color=color, width=width, dash=dash, fill=fill, opacity=opacity, arrow=arrow, close=close)
    def dot(self, x, y, r=4, color=BLUE, fill=None): self.f.circle(self.X(x), self.Y(y), r, color=color, fill=fill or color)
    def text(self, x, y, s, **kw): self.f.text(self.X(x), self.Y(y), s, **kw)
    def hline(self, y, **kw): self.f.line(self.x, self.Y(y), self.x + self.w, self.Y(y), **kw)
    def vline(self, x, **kw): self.f.line(self.X(x), self.y, self.X(x), self.y + self.h, **kw)
    def arrow(self, x1, y1, x2, y2, color=DARK, width=1.5): self.f.line(self.X(x1), self.Y(y1), self.X(x2), self.Y(y2), color=color, width=width, arrow=True)
