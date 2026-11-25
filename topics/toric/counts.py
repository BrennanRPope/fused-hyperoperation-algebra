"""Brute-force point counts over F_q for toric and non-toric shapes."""
import itertools
def gf(q):
    # returns (elements, add, mul) for GF(q), q = p^k, k<=3
    p=[x for x in range(2,q+1) if q%x==0][0]; k=0; t=q
    while t>1: t//=p; k+=1
    if k==1:
        E=list(range(q)); return E,(lambda a,b:(a+b)%p),(lambda a,b:(a*b)%p)
    irr={(2,2):[1,1,1],(2,3):[1,1,0,1],(3,2):[2,2,1]}[(p,k)]  # monic irreducible, low->high
    E=list(itertools.product(range(p),repeat=k))
    def add(a,b): return tuple((x+y)%p for x,y in zip(a,b))
    def mul(a,b):
        r=[0]*(2*k-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
        for d in range(2*k-2,k-1,-1):
            c=r[d]
            if c:
                for i in range(k+1): r[d-k+i]=(r[d-k+i]-c*irr[i])%p
        return tuple(r[:k])
    return E,add,mul
def proj_points(E,n,zero,one):
    for v in itertools.product(E,repeat=n+1):
        nz=[x for x in v if x!=zero]
        if not nz: continue
        # normalise: first nonzero coord = one
        i=next(i for i,x in enumerate(v) if x!=zero)
        if v[i]==one: yield v
def counts(q):
    E,add,mul=gf(q); zero=E[0]; one=E[1] if q==len(range(q)) else None
    if isinstance(E[0],tuple): zero=tuple([0]*len(E[0])); one=tuple([1]+[0]*(len(E[0])-1))
    else: zero,one=0,1
    pw=lambda a,n: (lambda r:r)(__import__('functools').reduce(mul,[a]*n,one))
    res={}
    res['torus G_m^2']=sum(1 for x in E for y in E if x!=zero and y!=zero)
    res['A^2']=q*q
    P2=list(proj_points(E,2,zero,one)); P3=list(proj_points(E,3,zero,one))
    res['P^2']=len(P2)
    res['xy=zw (P1xP1)']=sum(1 for x,y,z,w in P3 if mul(x,y)==mul(z,w))
    res['xy=z^2 (cone)']=sum(1 for x,y,z,w in P3 if mul(x,y)==mul(z,z))
    res['xyz=w^3']=sum(1 for x,y,z,w in P3 if mul(mul(x,y),z)==pw(w,3))
    res['Fermat cubic']=sum(1 for x,y,z,w in P3 if add(add(pw(x,3),pw(y,3)),add(pw(z,3),pw(w,3)))==zero)
    # y^2 z = x^3 - x z^2  ->  y^2 z + x z^2 = x^3
    res['elliptic y^2=x^3-x']=sum(1 for x,y,z in P2 if add(mul(mul(y,y),z),mul(x,mul(z,z)))==pw(x,3))
    return res
pred={'torus G_m^2':lambda q:(q-1)**2,'A^2':lambda q:q*q,'P^2':lambda q:q*q+q+1,
 'xy=zw (P1xP1)':lambda q:(q+1)**2,'xy=z^2 (cone)':lambda q:q*q+q+1,'xyz=w^3':lambda q:q*q+q+1}
if __name__=='__main__':
    import json
    out={}
    for q in [2,3,4,5,7,8,9,11,13]:
        c=counts(q); out[q]=c
        ok=all(c[k]==f(q) for k,f in pred.items())
        print(q,c,'toric formulas match' if ok else 'MISMATCH')
    json.dump(out,open('counts.json','w'))
