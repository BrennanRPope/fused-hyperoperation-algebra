"""Builds index.html from template.html + figures.json, with counts from checks_output.txt and the Lean file."""
import json, re
d = json.load(open('figures.json'))
data = {'coef': d['coef'], 'grid': [r[:81] for r in d['grid']]}
out = open('checks_output.txt').read()
n = re.search(r'(\d+) of (\d+) checks pass', out)
thms = len(re.findall(r'^theorem ', open('lean/Combinatorics.lean').read(), re.M))
html = open('template.html').read().replace('__DATA__', json.dumps(data, separators=(',', ':'))).replace('__CHECKS__', n.group(1), 1).replace('__CHECKS__', n.group(2), 1).replace('__THMS__', str(thms))
open('index.html', 'w').write(html); print('built', n.group(0), thms, 'theorems', len(html.splitlines()), 'lines')
