import re
svg = {k: open(f).read() for k, f in [("f1","fig1-dimension-dial.svg"),("f2","fig2-explained-residual.svg"),("f3","fig3-doubling-ladder.svg")]}
html = open("page_template.html").read()
for k, v in svg.items(): html = html.replace("{{"+k+"}}", v)
open("index.html","w").write(html)
print(len(html))
