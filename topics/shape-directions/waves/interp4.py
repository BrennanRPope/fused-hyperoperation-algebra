import numpy as np, itertools, random, sympy as sp, json
random.seed(5); np.random.seed(5)
pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
def sample():
    while True:
        a=np.random.uniform(-3,3,4); b=np.random.uniform(-3,3,4)
        pts=np.concatenate([a,b])
        if min(abs(x-y) for x,y in itertools.combinations(pts,2))>0.3: break
    return [(a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j])) for i,j in pairs]
for maxdeg in (1,2):
    mons=list(itertools.product(range(maxdeg+1),repeat=6))
    rows=[]
    for _ in range(len(mons)+150):
        A=sample(); rows.append([np.prod([A[k]**e[k] for k in range(6)]) for e in mons])
    M=np.array(rows); M/=np.linalg.norm(M,axis=1,keepdims=True)
    U,S,Vt=np.linalg.svd(M)
    print('maxdeg',maxdeg,'smallest sv',S[-4:])
    if S[-1]<1e-9*S[0]:
        v=Vt[-1]; v/=v[np.argmax(abs(v))]
        nz=[abs(x) for x in v if abs(x)>1e-6]; v/=min(nz)
        print('nullity',sum(S<1e-9*S[0]))
        A01,A02,A03,A12,A13,A23=sp.symbols('A01 A02 A03 A12 A13 A23'); X=[A01,A02,A03,A12,A13,A23]
        expr=sum(sp.Integer(int(round(c)))*sp.prod([X[k]**e[k] for k in range(6)]) for c,e in zip(v,mons) if abs(c)>1e-6)
        print('max rounding err',max(abs(c-round(c)) for c in v))
        print(sp.factor(expr)); print('terms',len(sp.Add.make_args(sp.expand(expr))))
        open('relation4.txt','w').write(str(sp.expand(expr)))
        break
