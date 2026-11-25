"""Shared page builder for the chemistry pages (one look for the whole series).
Each subject's make_page.py calls build(...) and writes index.html beside it."""
import html, json

CSS = r"""
/* VOICE work: a technical write-up read by one curious person, often on a phone
   IDEA a page from a chemistry handbook: each subject opens with its element-style tile and one finding
   PICTURE one large drawn chart or diagram per page, carrying the subject's mult-add finding
   COLOUR white ground, near-black ink, sodium-flame yellow slab behind the verdict, copper-blue marks, flame red accent
   TYPE Archivo 800 titles + Source Sans 3 text
   LAYOUT title and tile across the top, picture full width beneath, figures in a row; phone: stacked */
:root{color-scheme:light dark;--bg:light-dark(#ffffff,#0f1216);--ink:light-dark(#14161a,#e9ecef);--mute:light-dark(#555c66,#a3abb5);--rule:light-dark(#d9dde2,#2a3038);--blue:light-dark(#1f5fbf,#6ea3ff);--red:light-dark(#d23c2a,#ff7a66);--yel:light-dark(#ffd23f,#8a6d00);--yelink:light-dark(#14161a,#fff4c7);--panel:light-dark(#f4f6f8,#171b21);--green:light-dark(#1d7f4e,#5fd39a)}
:root[data-theme=light]{color-scheme:light}:root[data-theme=dark]{color-scheme:dark}
*{box-sizing:border-box;min-width:0}
html{scroll-padding-top:72px}
body{margin:0;background:var(--bg);color:var(--ink);font:17px/1.55 "Source Sans 3",system-ui,sans-serif;overflow-x:clip}
.band{padding-block:clamp(24px,4vw,56px);padding-left:max(5vw,16px);padding-right:max(5vw,16px)}
.wrap{max-width:1240px;margin:0 auto}
h1{font:800 clamp(30px,4.4vw,48px)/1.06 Archivo,system-ui,sans-serif;letter-spacing:-.02em;margin:0 0 14px;text-wrap:balance}
h2{font:800 clamp(22px,2.6vw,28px)/1.15 Archivo,system-ui,sans-serif;letter-spacing:-.01em;margin:0 0 10px;text-wrap:balance}
h3{font:700 18px/1.3 Archivo,system-ui,sans-serif;margin:18px 0 6px}
p{margin:0 0 12px;max-width:72ch}
.top{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:24px;align-items:start}
.tile{width:118px;height:132px;border:3px solid var(--ink);padding:8px 10px;display:flex;flex-direction:column;justify-content:space-between;background:var(--bg)}
.tile .z{font:700 15px Archivo,sans-serif}.tile .sym{font:800 46px/1 Archivo,sans-serif;letter-spacing:-.03em}.tile .nm{font-size:13px;line-height:1.2;color:var(--mute)}
.lead{font-size:19px;max-width:70ch}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(10em,1fr));gap:14px;margin:18px 0 6px}
.figs div{border-top:3px solid var(--ink);padding-top:6px}
.figs b{display:block;font:800 clamp(24px,3vw,34px)/1.1 Archivo,sans-serif;white-space:nowrap}
.figs span{font-size:14px;color:var(--mute)}
.fig{position:relative;min-height:360px;margin:12px 0 4px}
.fig.tall{min-height:440px}.fig.short{min-height:280px}
.fig>svg{position:absolute;inset:0;width:100%;height:100%;overflow:visible}
.cap{font-size:14px;color:var(--mute);margin:4px 0 0;max-width:none}
svg text{font:14px "Source Sans 3",sans-serif;fill:var(--ink);paint-order:stroke;stroke:var(--bg);stroke-width:4px;stroke-linejoin:round}
svg .ax{stroke:var(--rule);stroke-width:1}svg .axk{stroke:var(--ink);stroke-width:1.2}
.verdict{background:var(--yel);color:var(--yelink);padding:clamp(18px,3vw,30px) max(5vw,16px)}
.verdict h2{margin-bottom:6px}.verdict p{max-width:80ch;font-size:18px}
.sec{border-top:1px solid var(--rule)}
.two{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr);gap:clamp(18px,3vw,40px);align-items:start}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:15px}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--rule);vertical-align:top}
th{font-weight:700;border-bottom:2px solid var(--ink)}
.st{white-space:nowrap;font-weight:700}
.st.known{color:var(--mute)}.st.tidy{color:var(--blue)}.st.new{color:var(--red)}
code{font:14px ui-monospace,Menlo,monospace;overflow-wrap:anywhere;background:var(--panel);padding:1px 4px;border-radius:3px}
.lean{columns:2 22em;column-gap:32px;font-size:15px;padding-left:18px}.lean li{break-inside:avoid;margin-bottom:6px}
.eq{font:17px/1.4 "Source Sans 3",sans-serif;background:var(--panel);padding:10px 14px;margin:8px 0 14px;overflow-x:auto;white-space:nowrap}
a{color:var(--blue)}
.nav{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:15px;margin-top:10px}
footer{font-size:14px;color:var(--mute)}
@media(max-width:700px){.two{grid-template-columns:minmax(0,1fr)}.top{grid-template-columns:minmax(0,1fr)}.tile{width:96px;height:108px}.tile .sym{font-size:36px}}
@media(max-width:600px){table{font-size:13px}th,td{padding:6px 4px}.stack thead{display:none}.stack tr{display:grid;padding-block:10px}.stack td{display:block;border:0;padding:2px 0}.stack td::before{content:attr(data-l) ": ";font-weight:700}.lean{columns:1}}
"""

JS_HELPERS = r"""
const NS='http://www.w3.org/2000/svg';
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function nice(lo,hi,n){const span=hi-lo,step0=span/n,mag=Math.pow(10,Math.floor(Math.log10(step0))),f=step0/mag;const st=(f<1.5?1:f<3?2:f<7?5:10)*mag;const out=[];for(let v=Math.ceil(lo/st-1e-9)*st;v<=hi+1e-9;v+=st)out.push(+v.toPrecision(12));return out}
function fmt(v){const a=Math.abs(v);if(a!==0&&(a<1e-3||a>=1e5))return v.toExponential(0).replace('e+','e');return String(+v.toPrecision(4))}
function plot(el,o){
  const draw=()=>{const W=el.clientWidth,H=el.clientHeight,small=W<520;
    const m={l:o.ml||(small?46:60),r:o.mr||(small?14:150),t:24,b:o.xlab?52:30};
    const tx=v=>o.xlog?Math.log10(v):v,ty=v=>o.ylog?Math.log10(v):v;
    const X=v=>m.l+(tx(v)-tx(o.x[0]))/(tx(o.x[1])-tx(o.x[0]))*(W-m.l-m.r);
    const Y=v=>H-m.b-(ty(v)-ty(o.y[0]))/(ty(o.y[1])-ty(o.y[0]))*(H-m.t-m.b);
    let s='';
    const xt=o.xticks||(o.xlog?[]:nice(o.x[0],o.x[1],small?4:7)),yt=o.yticks||(o.ylog?[]:nice(o.y[0],o.y[1],5));
    yt.forEach(v=>{s+=`<line class="ax" x1="${m.l}" x2="${W-m.r}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${m.l-8}" y="${Y(v)+5}" text-anchor="end">${esc(o.yfmt?o.yfmt(v):fmt(v))}</text>`});
    xt.forEach((v,i)=>{if(small&&i%2&&xt.length>5)return;s+=`<text x="${X(v)}" y="${H-m.b+20}" text-anchor="middle">${esc(o.xfmt?o.xfmt(v):fmt(v))}</text>`});
    s+=`<line class="axk" x1="${m.l}" x2="${W-m.r}" y1="${H-m.b}" y2="${H-m.b}"/>`;
    if(o.xlab)s+=`<text x="${(m.l+W-m.r)/2}" y="${H-8}" text-anchor="middle">${esc(o.xlab)}</text>`;
    if(o.ylab)s+=`<text x="${m.l}" y="14">${esc(o.ylab)}</text>`;
    (o.bands||[]).forEach(b=>{s+=`<rect x="${X(b.x0)}" y="${m.t}" width="${X(b.x1)-X(b.x0)}" height="${H-m.t-m.b}" fill="${b.fill}" opacity="${b.op||.14}"/>`+(b.label?`<text x="${(X(b.x0)+X(b.x1))/2}" y="${m.t+16}" text-anchor="middle">${esc(b.label)}</text>`:'')});
    (o.hlines||[]).forEach(h=>{s+=`<line x1="${m.l}" x2="${W-m.r}" y1="${Y(h.y)}" y2="${Y(h.y)}" stroke="${h.c||'var(--mute)'}" stroke-dasharray="5 4"/>`+(h.label?`<text x="${W-m.r+6}" y="${Y(h.y)+5}">${esc(small?(h.short||''):h.label)}</text>`:'')});
    (o.series||[]).forEach(se=>{const c=se.c||'var(--blue)';
      if(se.line!==false&&se.pts.length>1)s+=`<polyline fill="none" stroke="${c}" stroke-width="${se.w||3}" ${se.dash?`stroke-dasharray="${se.dash}"`:''} points="${se.pts.map(p=>X(p[0])+','+Y(p[1])).join(' ')}"/>`;
      if(se.dots)se.pts.forEach(p=>{s+=`<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="${se.r||5}" fill="${c}"/>`});
      if(se.name){const p=se.pts[se.pts.length-1];const lx=small?X(p[0])-4:X(p[0])+8;s+=`<text x="${lx}" y="${Y(p[1])+(se.dy||0)+5}" text-anchor="${small?'end':'start'}" fill="${c}">${esc(small?(se.short||se.name):se.name)}</text>`}
      (se.labels||[]).forEach(L=>{s+=`<text x="${X(L.x)+(L.dx||8)}" y="${Y(L.y)+(L.dy||5)}" text-anchor="${L.a||'start'}">${esc(L.t)}</text>`});
    });
    if(o.extra)s+=o.extra({X,Y,W,H,m,small});
    el.innerHTML=s};
  const ro=new ResizeObserver(draw);ro.observe(el);draw();
}
"""

def table(rows, heads=("Claim", "Status", "Checked by")):
    th = "".join(f"<th>{html.escape(h)}</th>" for h in heads)
    body = ""
    for r in rows:
        cells = []
        for i, c in enumerate(r):
            l = html.escape(heads[i])
            if i == 1:
                cls = {"known": "known", "tidy reading": "tidy", "new here": "new"}.get(c, "known")
                cells.append(f'<td data-l="{l}"><span class="st {cls}">{html.escape(c)}</span></td>')
            else:
                cells.append(f'<td data-l="{l}">{c}</td>')
        body += "<tr>" + "".join(cells) + "</tr>"
    return f'<div class="scroll"><table class="stack"><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'

def build(path, *, title, tile, h1, lead, figs, picture, caption, verdict_h, verdict, sections, claims,
          lean_file, lean_items, files, sources, series_nav):
    figs_html = "".join(f"<div><b>{html.escape(v)}</b><span>{html.escape(t)}</span></div>" for v, t in figs)
    sec_html = ""
    for i, s in enumerate(sections):
        sid = s.get("id", f"s{i}")
        if s.get("fig"):
            inner = f'<div class="two"><div>{s["body"]}</div><div><div class="fig {s.get("figcls","short")}"><svg id="{sid}-fig"></svg></div><p class="cap">{s.get("cap","")}</p></div></div>'
        else:
            inner = s["body"]
        sec_html += f'<section class="band sec" id="{sid}"><div class="wrap"><h2>{s["h"]}</h2>{inner}</div></section>'
    nav = "".join(f'<a href="#{s.get("id", f"s{i}")}">{html.escape(s.get("nav", s["h"]))}</a>' for i, s in enumerate(sections))
    lean_html = "".join(f"<li><code>{html.escape(n)}</code>: {d}</li>" for n, d in lean_items)
    js = JS_HELPERS + "\n{\n" + picture + "\n" + "".join(s.get("js", "") for s in sections) + "\n}"
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Archivo:wght@700;800&family=Source+Sans+3:wght@400;600;700&display=swap" rel="stylesheet">
<style>{CSS}</style></head><body>
<header class="band"><div class="wrap">
<div class="top"><div><h1>{h1}</h1><p class="lead">{lead}</p></div>
<div class="tile" aria-hidden="true"><span class="z">{tile[0]}</span><span class="sym">{tile[1]}</span><span class="nm">{html.escape(tile[2])}</span></div></div>
<div class="figs">{figs_html}</div>
<div class="fig tall"><svg id="pic" role="img" aria-label="{html.escape(caption)}"></svg></div>
<p class="cap">{caption}</p>
<nav class="nav">{nav}<a href="#claims">Every claim</a><a href="#lean">Lean</a><a href="#files">Files</a></nav>
</div></header>
<section class="verdict"><div class="wrap"><h2>{verdict_h}</h2><p>{verdict}</p></div></section>
{sec_html}
<section class="band sec" id="claims"><div class="wrap"><h2>Every claim, what it is, and how it was checked</h2>
<p>Status: <b>known</b> is standard science restated; <b>tidy reading</b> is a known fact read in mult-add terms; <b>new here</b> is new as far as I could find.</p>{table(claims)}</div></section>
<section class="band sec" id="lean"><div class="wrap"><h2>Proved in Lean</h2>
<p><code>{html.escape(lean_file)}</code>: core Lean 4.37.0-rc1, no mathlib, no <code>sorry</code>, compiled against the shared library; every theorem uses only Lean's standard axioms.</p>
<ul class="lean">{lean_html}</ul></div></section>
<section class="band sec" id="files"><div class="wrap"><h2>Files and sources</h2><p>{files}</p><p>{sources}</p>
<footer><p>{series_nav}</p></footer></div></section>
<script>{js}</script></body></html>"""
    open(path, "w").write(doc)
    return len(doc)
