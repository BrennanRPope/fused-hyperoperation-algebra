"""Mirror of the Lean algorithm in Knots.lean, step for step, to produce the exact lists the
theorems state. Bracket <D> = A^(-3n) * P(A); f = (-A^3)^(-w) <D>."""
import json
def padd(p,q):
    L=max(len(p),len(q)); p=p+[0]*(L-len(p)); q=q+[0]*(L-len(q)); return [a+b for a,b in zip(p,q)]
def pmul(p,q):
    if not p: return []
    return padd([p[0]*c for c in q],[0]+pmul(p[1:],q))
def ppow(p,n): return [1] if n==0 else pmul(p,ppow(p,n-1))
def trim(p):
    while p and p[-1]==0: p=p[:-1]
    return p
def loops(m,pairs):
    comp=list(range(m+1))
    for a,b in pairs:
        ra,rb=comp[a],comp[b]
        if ra!=rb: comp=[ra if c==rb else c for c in comp]
    return sum(1 for i,c in enumerate(comp) if i>=1 and c==i)
def states(n):
    if n==0: return [[]]
    return [[False]+s for s in states(n-1)]+[[True]+s for s in states(n-1)]
def smooth(x,b):  # b False = A-smoothing
    a,bb,c,d=x
    return [(a,d),(bb,c)] if b else [(a,bb),(c,d)]
def term(n,m,pd,st):
    pairs=[pr for x,b in zip(pd,st) for pr in smooth(x,b)]
    L=loops(m,pairs); na=sum(1 for b in st if not b)
    e=2*n+2*na+2-2*L; assert e>=0
    sign=(-1)**(L-1)
    return [0]*e+[sign*c for c in ppow([1,0,0,0,1],L-1)], L
def bracketP(n,m,pd):
    tot=[]
    for st in states(n): tot=padd(tot,term(n,m,pd,st)[0])
    return trim(tot)
def normLP(lo,cs):
    cs=trim(cs)
    while cs and cs[0]==0: cs=cs[1:]; lo+=1
    return (lo if cs else 0, cs)
def fpoly(n,m,pd,w):
    P=bracketP(n,m,pd); sign=(-1)**(w%2)
    return normLP(-3*w-3*n,[sign*c for c in P])
def mirror(lp):
    lo,cs=lp; return (-(lo+len(cs)-1), cs[::-1])
K={
 'kinkUnknot':(1,2,[(1,1,2,2)],1),
 'trefoilL':(3,6,[(1,4,2,5),(3,6,4,1),(5,2,6,3)],-3),
 'trefoilR':(3,6,[(4,2,5,1),(6,4,1,3),(2,6,3,5)],3),
 'figEight':(4,8,[(4,2,5,1),(8,6,1,5),(6,3,7,4),(2,7,3,8)],0),
 'hopf':(2,4,[(4,1,3,2),(2,3,1,4)],2),
 'cinquefoilL':(5,10,[(1,6,2,7),(3,8,4,9),(5,10,6,1),(7,2,8,3),(9,4,10,5)],-5),
}
def main():
    res={}
    for k,(n,m,pd,w) in K.items():
        P=bracketP(n,m,pd); f=fpoly(n,m,pd,w)
        lo,cs=f
        t={}
        for i,c in enumerate(cs):
            if c: t[(lo+i)/-4]=c
        res[k]={'P':P,'f':f,'jones_t':{str(e):c for e,c in sorted(t.items())}}
        print(k,'P=',P,'\n   f=',f,' V(t)=',t)
        if all((lo+i)%4==0 for i,c in enumerate(cs) if c):
            print('   V(1)=',sum(cs),' V(-1)=',sum(c*(-1)**((lo+i)//4) for i,c in enumerate(cs)))
    print('mirror check',mirror(res['trefoilL']['f'])==res['trefoilR']['f'], mirror(res['figEight']['f'])==res['figEight']['f'])
    # trefoil state cube
    n,m,pd,w=K['trefoilL']
    cube=[(sum(st),term(n,m,pd,st)[1]) for st in states(n)]
    print('trefoil states (number of B-smoothings, loops):',cube)
    res['trefoilL_cube']=cube
    json.dump(res,open('knots_results.json','w'),indent=1)

if __name__=="__main__":
    main()
