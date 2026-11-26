# DRAFT. Parse Li & Liao 2017 (695 families): class, v1, v2, T, T*, Lf, word.
import re, json, sys, os
H = os.path.dirname(os.path.abspath(__file__))
def rows(path):
    t = open(path).read()
    for tr in t.split('<tr>')[1:]:
        cells = re.findall(r'<t[hd][^>]*>(.*?)</t[hd]>', tr, re.S)
        yield cells
def name(c):
    c = re.sub(r'<sup>.*?</sup>', '', c, flags=re.S)
    c = re.sub(r'<[^>]+>', '', c)
    return c.replace(' ', '').strip()
P = {}
for cells in rows(os.path.join(H,'data','li-liao-2017-pictures.md')):
    if len(cells) >= 6:
        try:
            P[name(cells[0])] = dict(v1=float(cells[1]), v2=float(cells[2]), T=float(cells[3]), Tstar=float(cells[4]), Lf=int(re.sub(r'<[^>]+>','',cells[5])))
        except ValueError: pass
W = {}
for cells in rows(os.path.join(H,'data','li-liao-2017-words.md')):
    if len(cells) == 2:
        w = re.sub(r'<[^>]+>','',cells[1]).strip()
        if re.fullmatch(r'[abAB]+', w): W[name(cells[0])] = w
out = []
for k,v in P.items():
    if k in W:
        v = dict(v); v['cls']=k; v['word']=W[k]; out.append(v)
print(len(P), len(W), len(out), file=sys.stderr)
bad = [o for o in out if len(o['word'])!=o['Lf']]
print('Lf mismatch', len(bad), file=sys.stderr)
json.dump(out, open(os.path.join(H,'data','li-liao-695.json'),'w'), indent=0)
