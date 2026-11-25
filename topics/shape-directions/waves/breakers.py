import sympy as sp, itertools, json
X=sp.symbols('A01 A02 A03 A12 A13 A23'); k=sp.symbols('k')
R=sp.sympify(open('relation4.txt').read())
pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
out={}
def arr(edges):  # weight k on edges, 1 elsewhere
    return {X[n]:(k if p in edges else 1) for n,p in enumerate(pairs)}
graphs={'K4 (every pair k)':set(pairs),
 'K4 minus one edge':set(pairs)-{(2,3)},
 'cycle C4':{(0,1),(1,2),(2,3),(0,3)},
 'path P4':{(0,1),(1,2),(2,3)},
 'star':{(0,1),(0,2),(0,3)},
 'triangle + lone wave':{(0,1),(0,2),(1,2)},
 'two separate pairs':{(0,1),(2,3)},
 'one pair':{(0,1)}}
for name,E in graphs.items():
    f=sp.factor(R.subs(arr(E)))
    out[name]=str(f); print(name,':',f)
json.dump(out,open('breakers_graphs.json','w'),indent=1)
