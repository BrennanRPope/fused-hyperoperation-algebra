import json,numpy as np
from wind_ode import windings
from words import walk_area
d=json.load(open('../three-body-word-formula/data/li-liao-695.json'))
d=sorted(d,key=lambda o:len(o['word']))[:90]
out=[]
for o in d:
    w,c=windings(o['v1'],o['v2'],o['T'],rtol=1e-11); x,y,_=walk_area(o['word'])
    out.append(dict(cls=o['cls'],word=o['word'],w=w,close=c,netA=x,netB=y))
json.dump(out,open('wind_cat.json','w'))
import collections
tab=collections.defaultdict(collections.Counter)
for e in out: tab[(e['netA'],e['netB'])][tuple(int(round(v)) for v in e['w'].values())]+=1
for k,v in sorted(tab.items()): print(k,dict(v))
print('max |w - round(w)|',max(abs(v-round(v)) for e in out for v in e['w'].values()),'sum nonzero',sum(1 for e in out if abs(sum(e['w'].values()))>0.5))
