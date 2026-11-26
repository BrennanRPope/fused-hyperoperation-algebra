"""Merge the three raw bridge tables, map fields to map nodes, tag the framework hub, flag likely duplicates."""
import csv, re, json, difflib, collections, os
D = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in ['physics', 'math', 'shapes-applied']:
    for x in csv.DictReader(open(f'{D}/raw/{f}.tsv'), delimiter='\t'):
        x['src'] = f; rows.append(x)

def node(x):
    f = (x['field'] + ' ' + x['folders']).lower()
    fl = x['field'].lower()
    rules = [
        ('Solitons & integrable', ['integrable', 'soliton', 'kdv']),
        ('Chemistry', ['chem', 'kinetics', 'enzym', 'pharmac', 'biochem']),
        ('Biology', ['biolog', 'ecolog', 'evolution', 'neuro', 'genetic', 'physiolog', 'population', 'behav']),
        ('Economics, actuarial & engineering', ['actuar', 'econ', 'finance', 'engineering', 'datacenter', 'distributed', 'consensus', 'rulecraft', 'online learning', 'likeness']),
        ('Fluids', ['fluid', 'navier']),
        ('Quantum gravity & holography', ['quantum gravity', 'holograph', 'string', 'ads']),
        ('Gravity & GR', ['general relativity', 'gravity', 'cosmolog']),
        ('Relativity & EM', ['special relativity', 'electromag', 'relativity', 'particle']),
        ('Quantum mechanics', ['quantum mechanics', 'quantum foundations', 'quantum computing', 'quantum field', 'quantum chem', 'atomic']),
        ('Statistical physics & thermo', ['statistical', 'thermo', 'condensed', 'mathematical physics']),
        ('Information & entropy', ['information', 'evidence theory', 'fuzzy']),
        ('Probability & statistics', ['probab', 'statistic']),
        ('Classical mechanics', ['classical mechanics', 'physics']),
        ('Knots & topology', ['knot', 'topolog']),
        ('Number theory', ['number theory', 'arithmetic']),
        ('Formal groups & algebra', ['formal group', 'algebra ', 'algebra/', 'algebra', 'clifford']),
        ('Geometry', ['geometr', 'convex', 'tropical']),
        ('Analysis & dynamics', ['analysis', 'dynamic', 'numerical', 'complex']),
        ('Logic, sets & computation', ['set theory', 'logic', 'computation', 'category']),
        ('Combinatorics & games', ['combinator', 'game', 'optimization', 'meta']),
    ]
    for name, keys in rules:
        if any(k in fl for k in keys): return name
    for name, keys in rules:
        if any(k in f for k in keys): return name
    return 'Combinatorics & games'

def hub(x):
    t = (x['our_side'] + ' ' + x['note'] + ' ' + x['their_side']).lower()
    rules = [
        ('Swap rule (turn = d x area)', ['swap', 'precess', 'thomas', 'berry', 'foucault', 'phase', 'holonomy', 'area']),
        ('Cube corners & weights', ['corner', 'cube', 'lucas', 'minor', 'hirota', 'tau', 'boolean', 'inclusion', 'subset', 'tesseract', 'mobius', 'möbius', 'permutohedr', 'associahedr']),
        ('Rungs & height flow', ['rung', 'height', 'e^x', 'iterat', 'abel', 'branching', 'tetration', 'super', 'hyperop']),
        ('Shapes: simplex, ball, sphere', ['simplex', 'sphere', 'ball', 'hilbert', 'polar', 'fock', 'elliptope', 'roundness', 'norm']),
        ('Chains & path sums', ['chain', 'path', 'product integral', 'checkerboard', 'telegraph', 'propagator', 'semiring', 'tropical', 'maslov']),
        ('Plane of laws (k, d)', ['plane', '(k, d)', 'wall', 'parabola', 'landscape', 'formal group', 'velocity', 'tan', 'möbius group', 'riccati', 'logistic', 'd =', 'q =', 'curvature']),
    ]
    for name, keys in rules:
        if any(k in t for k in keys): return name
    return 'Dial k (one law a+b+kab)'

def verdict(x):
    v = x['verdict'].upper()
    return 'POSSIBLY-NEW' if v.startswith('POSS') else v.split()[0] if v else 'FRAMING'

for x in rows:
    x['node'] = node(x); x['hub'] = hub(x); x['v'] = verdict(x)
    x['key'] = re.sub(r'[^a-z0-9 ]', ' ', (x['our_side'] + ' ' + x['their_side']).lower())

# duplicates across the three tables (same bridge filed twice)
dup = set()
for i, a in enumerate(rows):
    for b in rows[i+1:]:
        if a['src'] == b['src'] or b['id'] in dup: continue
        r = difflib.SequenceMatcher(None, a['key'], b['key']).ratio()
        if r > 0.62: dup.add(b['id'])
# hand-checked cross-table repeats (second copy of the same bridge)
dup |= set('S128 S11 S126 S129 A08 S127 S130 M228 M225 M229 M222 S122 M132 M146 S121 S119 S101 S104 S103 S105 S102 S106 S17 S107 S110 S112 S113 S114 S115 S117 A89 S135 S133 S97 S98 S88 S86 S139 S99 S109 P109'.split())
dup -= {'S112b'}
for x in rows: x['dup'] = x['id'] in dup
json.dump(rows, open(f'{D}/bridges.json', 'w'), ensure_ascii=False, indent=0)
u = [x for x in rows if not x['dup']]
print('rows', len(rows), 'distinct', len(u))
print(collections.Counter(x['v'] for x in u))
print(collections.Counter(x['node'] for x in u).most_common())
print(collections.Counter(x['hub'] for x in u).most_common())


# export the merged table
VL = {'KNOWN': 'known', 'FRAMING': 'our reading of known pieces', 'POSSIBLY-NEW': 'possibly new', 'NEGATIVE': 'tested and failed'}
with open(f'{D}/bridges.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['id', 'our_side', 'their_side', 'map_field', 'framework_part', 'verdict', 'source_cited', 'folders', 'note', 'repeat_of_another_row'])
    for x in rows:
        w.writerow([x['id'], x['our_side'], x['their_side'], x['node'], x['hub'], VL.get(x['v'], x['v']), x['source_cited'], x['folders'], x['note'], 'yes' if x['dup'] else ''])
