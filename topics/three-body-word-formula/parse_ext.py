"""DRAFT. Parse the two outside lists that carry words (both from the sjtu-liao/three-body GitHub repository):
Li, Jing & Liao 2018 unequal masses (1, 1, m3) and Li & Liao 2019 free fall (masses in the label). Writes
data/unequal-mass.json and data/free-fall.json. Files are read as data only."""
import re, json, os, sys
H = os.path.dirname(os.path.abspath(__file__)); E = os.path.join(H, 'data', 'external')
def rows(path):
    t = open(path).read()
    for tr in t.split('<tr>')[1:]: yield re.findall(r'<t[hd][^>]*>(.*?)</t[hd]>', tr, re.S)
def clean(c): return re.sub(r'\s+', '', re.sub(r'<[^>]+>', '', re.sub(r'<sup>.*?</sup>', '', c, flags=re.S)))
def num(c):
    try: return float(clean(c))
    except ValueError: return None
# unequal mass
P, W = {}, {}
for c in rows(os.path.join(E, 'liao-unequal-mass', 'three-body-unequal-mass-pictures.md')):
    if len(c) >= 6 and num(c[1]) is not None: P[clean(c[0])] = dict(v1=num(c[1]), v2=num(c[2]), T=num(c[3]), Tstar=num(c[4]), Lf=int(num(c[5])))
for c in rows(os.path.join(E, 'liao-unequal-mass', 'three-body-unequal-mass-free-group-word.md')):
    if len(c) >= 2:
        w = clean(c[-1])
        if re.fullmatch(r'[abAB]+', w): W[clean(c[0])] = w
U = []
for k, v in P.items():
    if k in W:
        m = re.search(r'\(([\d.]+)\)$', k); v = dict(v, cls=k, word=W[k], m3=float(m.group(1))); U.append(v)
print('unequal', len(P), len(W), len(U), 'Lf mismatch', sum(len(o['word']) != o['Lf'] for o in U), file=sys.stderr)
json.dump(U, open(os.path.join(H, 'data', 'unequal-mass.json'), 'w'))
# free fall
P, W = {}, {}
for c in rows(os.path.join(E, 'liao-free-fall', 'free-fall-3b-pictures.md')):
    if len(c) >= 4 and num(c[1]) is not None: P[clean(c[0])] = dict(x=num(c[1]), y=num(c[2]), T=num(c[3]))
for c in rows(os.path.join(E, 'liao-free-fall', 'free-fall-3b-free-group-word.md')):
    if len(c) >= 2:
        w = clean(c[-1])
        if re.fullmatch(r'[abAB]+', w): W[clean(c[0])] = w
Fo = []
for k, v in P.items():
    if k in W:
        m = re.search(r'\(([\d.]+),([\d.]+),([\d.]+)\)', k); ms = [float(x) for x in m.groups()]
        x, y = v['x'], v['y']; r12 = 1.0; r13 = ((x + .5)**2 + y**2)**.5; r23 = ((x - .5)**2 + y**2)**.5
        En = -(ms[0] * ms[1] / r12 + ms[0] * ms[2] / r13 + ms[1] * ms[2] / r23)
        Fo.append(dict(v, cls=k, word=W[k], masses=ms, E=En, Tstar=v['T'] * abs(En)**1.5))
print('free fall', len(P), len(W), len(Fo), file=sys.stderr)
json.dump(Fo, open(os.path.join(H, 'data', 'free-fall.json'), 'w'))
