"""Build one biology page from a spec module. usage: python3 build_page.py <subject-dir>
The spec (page.py in the subject dir) defines SPEC = dict(title, h1, lede, stats[(num,label)],
hero_fig, hero_cap, known (html), lead (html), sections[(id, nav, h2, html)], table[(rows)], lean_file, sources[(text,url)])"""
import sys, os, html, importlib.util, csv
d=os.path.abspath(sys.argv[1])
spec=importlib.util.spec_from_file_location("page",os.path.join(d,"page.py")); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
S=m.SPEC
def fig(src,cap,alt):
    return f'<figure class="plate"><img src="{src}" alt="{html.escape(alt)}" loading="lazy"><figcaption>{cap}</figcaption></figure>'
# Lean theorem list from theorems.tsv
names=[]
tsv=os.path.join(d,"lean","theorems.tsv")
if os.path.exists(tsv):
    for row in csv.reader(open(tsv),delimiter="\t"):
        if row: names.append(row[1].split(".")[-1])
nav="".join(f'<a href="#{i}">{html.escape(n)}</a>' for i,n,_,_ in S["sections"])+'<a href="#table">Verdicts</a><a href="#lean">Lean</a><a href="#sources">Sources</a>'
secs="".join(f'<section class="band" id="{i}"><div class="wrap"><h2>{h2}</h2>{body}</div></section>' for i,_,h2,body in S["sections"])
rows="".join("<tr>"+"".join(f"<td>{c}</td>" for c in r)+"</tr>" for r in S["table"])
stats="".join(f"<li><b>{a}</b>{b}</li>" for a,b in S["stats"])
src="".join(f'<li><a href="{u}">{html.escape(t)}</a></li>' if u else f"<li>{html.escape(t)}</li>" for t,u in S["sources"])
leanlist=", ".join(f"<code>{html.escape(n)}</code>" for n in names)
page=f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(S['title'])}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Inter:wght@400;600&display=swap" rel="stylesheet">
<style>
/* VOICE a careful reading document in the biology series, honest about what is known and what might be new
   IDEA a naturalist's field plate pinned over a proof notebook: one living system, one dial, the corner cell in view
   PICTURE the lead plate in the first screen, then one plate per finding, a verdict table and the Lean list
   COLOUR white ground, near-black ink, leaf green for the living thing and the lead, ultramarine for the circle operation, vermilion for warnings; night is deep ink with plates kept white
   TYPE Newsreader for titles, Inter for text
   LAYOUT finding and lead plate first; known/new and lead-verdict boxes; one band per finding; verdict table; Lean; sources */
:root{{color-scheme:light dark;--bg:light-dark(#ffffff,#0e1414);--ink:light-dark(#15171a,#e8eeea);--mute:light-dark(#5b615d,#a6b2ab);--green:light-dark(#1d6b3c,#7fd6a0);--blue:light-dark(#1d2bd8,#8f9bff);--red:light-dark(#c2361c,#ff8a6b);--wash:light-dark(#eef6f0,#16241c);--rule:light-dark(#d8ddd9,#2a3a31)}}
:root[data-theme=light]{{color-scheme:light}}:root[data-theme=dark]{{color-scheme:dark}}
*{{box-sizing:border-box;min-width:0}}html{{scroll-padding-top:64px}}
body{{margin:0;background:var(--bg);color:var(--ink);font:16.5px/1.6 Inter,system-ui,sans-serif;overflow-x:clip}}
h1,h2,h3{{font-family:Newsreader,Georgia,serif;font-weight:600;line-height:1.08;text-wrap:balance;margin:0 0 .4em}}
h1{{font-size:clamp(32px,5vw,50px)}}h2{{font-size:clamp(25px,3.2vw,34px)}}h3{{font-size:20px;margin-top:1.2em}}
p,li{{max-width:70ch}}a{{color:var(--blue)}}code{{overflow-wrap:anywhere;font-size:.9em}}
.band{{padding:clamp(26px,4vw,56px) max(5vw,16px)}}.wrap{{max-width:1240px;margin:0 auto}}
.hero{{display:grid;grid-template-columns:minmax(0,5fr) minmax(0,6fr);gap:clamp(20px,3vw,44px);align-items:center}}
.kicker{{color:var(--green);font-weight:600;margin:0 0 .6em}}
.lede{{font-size:18.5px}}
.plate{{background:#fff;border:1px solid var(--rule);border-radius:6px;padding:10px;margin:18px 0}}
.plate img{{display:block;width:100%;height:auto}}
.plate figcaption{{color:#5b615d;font-size:13.5px;padding:8px 4px 2px}}
.stats{{display:flex;flex-wrap:wrap;gap:12px 28px;margin:18px 0 0;padding:0;list-style:none}}
.stats li{{font-size:14px;color:var(--mute)}}.stats b{{display:block;font:600 28px/1.1 Newsreader,serif;color:var(--ink)}}
.box{{background:var(--wash);border-radius:6px;padding:18px 22px;margin:22px 0}}.box p{{margin:.4em 0}}
.box.lead{{border-left:4px solid var(--green)}}
nav.scroll{{display:flex;gap:18px;white-space:nowrap;padding:12px max(5vw,16px);border-block:1px solid var(--rule);font-size:14.5px;background:var(--bg)}}
nav.scroll a{{flex:none;text-decoration:none}}
@media(min-width:900px){{nav.scroll{{position:sticky;top:0;z-index:2}}}}
.eq{{font-family:Newsreader,Georgia,serif;font-size:21px;background:var(--wash);border-radius:6px;padding:12px 18px;overflow-x:auto;white-space:nowrap}}
table{{border-collapse:collapse;width:100%;font-size:15px}}th,td{{text-align:left;vertical-align:top;padding:9px 10px;border-bottom:1px solid var(--rule)}}th{{font-weight:600}}
.scroll{{overflow-x:auto}}
@media(max-width:760px){{.hero{{grid-template-columns:1fr}}}}
@media(max-width:600px){{table{{font-size:13px}}.stack thead{{display:none}}.stack tr{{display:grid;padding-block:10px}}.stack td{{display:flow-root;border:0;padding:3px 0;padding-inline-start:8em;text-align:start}}.stack td::before{{content:attr(data-l);float:inline-start;width:7.5em;margin-inline-start:-8em;white-space:normal;color:var(--mute)}}.stack td:empty{{display:none}}}}
</style><script>new MutationObserver(()=>{{for(const c of document.querySelectorAll('.stack td:not([data-l])')){{let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}}}}).observe(document,{{childList:true,subtree:true}})</script></head><body>
<header class="band"><div class="wrap hero"><div>
<p class="kicker">{S.get('kicker','Biology through the circle operation')}</p>
<h1>{S['h1']}</h1><p class="lede">{S['lede']}</p><ul class="stats">{stats}</ul></div>
{fig(S['hero_fig'],S['hero_cap'],S.get('hero_alt',S['hero_cap']))}
</div></header>
<nav class="scroll">{nav}</nav>
<section class="band"><div class="wrap"><div class="box"><p><b>Known or new.</b> {S['known']}</p></div><div class="box lead"><p><b>Is there a lead worth following?</b> {S['lead']}</p></div></div></section>
{secs}
<section class="band" id="table"><div class="wrap"><h2>{S.get('table_h2','Where the circle operation shows up, and how much it gives')}</h2>
<table class="stack"><thead><tr><th>Place</th><th>What the circle operation does there</th><th>Known since</th><th>Verdict</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section class="band" id="lean"><div class="wrap"><h2>Lean: {len(names)} theorems, compiled, no sorry</h2>
<p>File <code>{S['lean_file']}</code>, core Lean 4.37 with no mathlib, built on the shared library's base and using ∘₁ for the circle operation. Every theorem uses only Lean's standard axioms. It goes into the shared proof library as <code>FusedAlgebra/Topics/Biology/{os.path.basename(S['lean_file'])}</code>.</p>
<p class="scroll" style="font-size:14px">{leanlist}</p>
<p>Numeric checks: {S['checks']}, from <code>checks.py</code> beside this page.</p></div></section>
<section class="band" id="sources"><div class="wrap"><h2>Sources</h2><ul>{src}</ul>
<p style="color:var(--mute);font-size:14px">{S.get('footnote','')}</p></div></section>
</body></html>"""
open(os.path.join(d,"index.html"),"w").write(page)
print("wrote",os.path.join(d,"index.html"),len(page),"bytes;",len(names),"theorems")
