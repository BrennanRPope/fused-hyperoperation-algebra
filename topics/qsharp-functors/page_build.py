import base64, re, sys, markdown
src = sys.argv[1]; out = sys.argv[2]
md = open(src + "/README.md").read()
def img(m):
    alt, path = m.group(1), m.group(2)
    b = base64.b64encode(open(src + "/" + path, "rb").read()).decode()
    return f'![{alt}](data:image/png;base64,{b})'
import os
if False:
    om = open(src + "/orbits/README.md").read()
    om = re.sub(r'!\[([^\]]*)\]\(figures/', r'![\1](orbits/figures/', om)
    om = om.replace("\n## ", "\n### ")
    om = re.sub(r'^# DRAFT · ', '## Orbits section · ', om)
    md = md + "\n\n---\n\n" + om
md = re.sub(r'!\[([^\]]*)\]\(([^)]+\.png)\)', img, md)
body = markdown.markdown(md, extensions=["tables"])
body = body.replace("<table>", '<div class="scroll"><table>').replace("</table>", "</table></div>")
html = f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Q# forms on our shapes</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Red+Hat+Display:wght@700&family=Red+Hat+Text:wght@400;500;700&display=swap" rel="stylesheet">
<style>
/* VOICE: a research report, calm and exact, conclusion first
   IDEA: a lab poster built around one sphere and its corners
   PICTURE: the labeled figures carry the argument
   COLOUR and TYPE: oxford blue band, white panels, magenta for the verdicts; Red Hat
   LAYOUT: one readable column, figures full width */
:root{{color-scheme:light dark;--bg:light-dark(#f4f6fb,#070f26);--panel:light-dark(#ffffff,#0f1b3d);--ink:light-dark(#0b1f4b,#e8ecf7);--mute:light-dark(#4a5778,#a9b3cf);--acc:#e0115f;--line:light-dark(#d9deea,#24325c)}}
:root[data-theme=light]{{color-scheme:light}}:root[data-theme=dark]{{color-scheme:dark}}
*{{box-sizing:border-box;min-width:0}}
body{{margin:0;background:var(--bg);color:var(--ink);font:16.5px/1.6 "Red Hat Text",system-ui,sans-serif;overflow-x:clip}}
header{{background:#0b1f4b;color:#fff;padding:clamp(28px,4vw,56px) max(5vw,16px)}}
header p{{max-width:68ch;color:#cfd8f0;margin:.6em 0 0}}
.tag{{display:inline-block;background:var(--acc);color:#fff;font-weight:700;border-radius:4px;padding:2px 8px;font-size:13px;margin-bottom:10px}}
main{{max-width:1040px;margin:0 auto;padding:24px max(4vw,14px) 64px}}
article{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:clamp(16px,3vw,40px)}}
h1,h2,h3{{font-family:"Red Hat Display",sans-serif;line-height:1.1;text-wrap:balance}}
header h1{{font-size:clamp(30px,4.4vw,48px);margin:0}}
article h1{{display:none}} h2{{font-size:clamp(22px,2.6vw,30px);margin-top:1.8em;padding-top:.6em;border-top:1px solid var(--line)}} h3{{font-size:20px;margin-top:1.5em}}
p,li{{max-width:72ch}} em{{color:var(--mute)}}
img{{width:100%;height:auto;display:block;background:#fff;border-radius:8px;border:1px solid var(--line);margin:18px 0 6px}}
table{{border-collapse:collapse;font-size:14.5px;width:100%}} th,td{{border-bottom:1px solid var(--line);padding:7px 8px;text-align:left;vertical-align:top}}
.scroll{{overflow-x:auto}} code{{overflow-wrap:anywhere;font-size:.92em}} hr{{border:0;border-top:1px solid var(--line);margin:2em 0}}
strong{{color:inherit}}
</style></head><body>
<header><span class="tag">DRAFT · Nov 26, 2026</span><h1>Q#’s body, Adjoint and Controlled on our shapes</h1>
<p>Hypercube Composition Theory. Gates on the 3-sphere in our coordinates: the body is our core law, Adjoint is its inverse (a mirror), Controlled multiplies by a corner of the control cube. Orange in the pictures is known math; blue is our reading; red does not hold.</p></header>
<main><article>{body}</article></main></body></html>'''
open(out, "w").write(html)
print(len(html))
