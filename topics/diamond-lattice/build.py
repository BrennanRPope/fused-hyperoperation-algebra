"""Builds index.html from template.html: figure data from figures/figdata.json, counts from checks and Lean."""
import re, json
out = open('checks_output.txt').read()
n = re.search(r'(\d+)/(\d+) checks passed', out)
thms = len(re.findall(r'^theorem ', open('lean/DiamondLattice.lean').read(), re.M))
d = json.load(open('figures/figdata.json'))
d['julia'] = [{'name': j['name'], 'Q': j['Q'], 'note': j['note'], 'png': j['png']} for j in d['julia']]
d['nu'] = [[round(a, 4), round(b, 4)] for a, b in d['nu']]
html = open('template.html').read().replace('__CHECKS__', n.group(1)).replace('__TOTAL__', n.group(2)).replace('__THMS__', str(thms)).replace('__DATA__', json.dumps(d, separators=(',', ':')))
open('index.html', 'w').write(html); print('built', n.group(0), thms, 'theorems', len(html.splitlines()), 'lines', len(html)//1024, 'KB')
