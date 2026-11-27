"""Round 4 checks: the true dimension off the tangent line. Run: python3 checks4.py"""
import random, math, json, itertools
from fractions import Fraction as F
import numpy as np
import planelaws as P
from offline import bounds, zc, rho
random.seed(4)
res=[]
def ok(name,cond,extra=""):
    res.append(bool(cond)); print(("PASS " if cond else "FAIL ")+name+("  "+extra if extra else ""))
s1,s2,s3=F(1,5),F(3,20),F(1,4)
def gens(q):
    return [(0,P.map_with(F(1),F(q),s)[0],s) for s in (s1,s2)]+[(1,P.map_with(F(-1),F(-q),s3)[0],s3)]
walls=lambda q: [(F(1),F(q)),(F(-1),F(-q))]

# U1 block derivative formula, exact
good=True
for q in (2,3,12):
    G=gens(q)
    for _ in range(40):
        law=random.choice([0,1]); mats=[g for g in G if g[0]==law]
        word=[random.choice(mats) for _ in range(random.randint(1,4))]
        M=((F(1),F(0)),(F(0),F(1))); lam=F(1)
        for _,A,s in word: M=P.mm(A,M); lam*=s
        y=F(random.randint(-99,99),100); p,qq=walls(q)[law]
        z=(y-p)/(y-qq); a,b=M[0]; c,d=M[1]
        exact=(a*d-b*c)/(c*y+d)**2
        good&= exact==lam*((1-z)/(1-lam*z))**2
ok("U1 a block's derivative is lam*((1-z)/(1-lam z))^2 exactly (120 random blocks, fractions)",good)

# U2 orientation preserving, strictly invariant interval V=[-1.05,1.05] (uniform hyperbolicity for Solomyak-Takahashi)
good=True
for q in (1.5,2,3,12,100):
    G=gens(F(q).limit_denominator())
    for _,A,s in G:
        a,b=A[0]; c,d=A[1]; good&= a*d-b*c>0
        for x in (F(-105,100),F(105,100)):
            y=P.app(A,x); good&= F(-105,100)<y<F(105,100)
ok("U2 every map has positive determinant and sends V=[-1.05,1.05] strictly inside itself (q=1.5..100)",good)

# U3 no accidental collisions for the first-return pieces S_2 (u in law1, degree<=2; v = law-2 power <=2): words of length 3
q=3; G=gens(q); a,b,c=[g[1] for g in G]
def prod(ms):
    M=((F(1),F(0)),(F(0),F(1)))
    for A in ms: M=P.mm(M,A)
    return P.norm(M)
U=[[a],[b],[a,a],[a,b],[b,b]]; V=[[c],[c,c]]
S=[u+v for u in U for v in V]
words=[sum(w,[]) for w in itertools.product(S,repeat=3)]
distinct=len({prod(w) for w in words})
ok("U3 the 10 first-return pieces S_2 give 1000 different maps for their 1000 words of length 3 (free, walls ±3)",distinct==1000,f"{distinct}")

# U4 soundness of the bounds: along random monoid words, every block's input lies in the region the bound assumed,
#    and the word's true derivative lies between the products of the lower and upper block factors
good=True; tested=0
for q in (2,3,12):
    cl=[(1.0,float(q),[0.2,0.15]),(-1.0,-float(q),[0.25])]
    lo,up,ds,info=bounds(cl,nb=3)
    out=info['outputs']; inp=info['inputs']
    for _ in range(300):
        nblk=random.randint(2,8); law=random.choice([0,1]); x=random.uniform(*inp[law])   # start where the bound's states start
        y=x; prev=None; dlo=1.0; dup=1.0; dtrue=1.0
        for k in range(nblk):
            p,qq,ms=cl[law]; deg=random.choice([1,1,2,3,5,9]); lam=math.prod(random.choice(ms) for _ in range(deg))
            z=zc(y,p,qq); dtrue*=lam*rho(z,lam)
            if prev is not None:
                j,bj=prev; ylo,yhi=out[j][bj]
                good&= ylo-1e-12<=y<=yhi+1e-12 and inp[law][0]-1e-12<=y<=inp[law][1]+1e-12
                zs=sorted([zc(ylo,p,qq),zc(yhi,p,qq)])
                dlo*=lam*rho(zs[1],lam); dup*=lam*rho(zs[0],lam)
            else:
                zs=sorted([zc(inp[law][0],p,qq),zc(inp[law][1],p,qq)]); dlo*=lam*rho(zs[1],lam); dup*=lam*rho(zs[0],lam)
            y=(p-qq*(lam*z))/(1-lam*z)
            prev=(law,min(deg,3)-1); law=1-law; tested+=1
        good&= dlo<=dtrue*(1+1e-9) and dtrue<=dup*(1+1e-9)
ok("U4 the bounds are sound on 900 random words started inside the attractor's input region (every input in the assumed region, true derivative between the bounds)",good,f"{tested} blocks")

# U5 brackets, and strictly below D*
PL=json.load(open("planelaws.json"))
B={}
for name,q in [("close",2),("near",3),("far",12)]:
    lo,up,ds,_=bounds([(1.0,float(q),[0.2,0.15]),(-1.0,-float(q),[0.25])],nb=6)
    B[q]=(lo,up,ds)
    ok(f"U5 walls ±{q}: {lo:.4f} <= delta <= {up:.4f}, width < 0.002, and strictly below D* = {ds:.4f}", up-lo<0.002 and up<ds-0.04)

# U6 Richardson check: the level-n estimates D_n (n = 7 and 14, from dn_levels.json) extrapolate into the bracket
DN=json.load(open("dn_levels.json"))
for q in (2,3,12):
    d=dict((int(n),v) for n,c,v in DN[str(q)]); ext=2*d[14]-d[7]
    lo,up,_=B[q]
    ok(f"U6 walls ±{q}: D_n falls slowly (D_7 {d[7]:.4f}, D_14 {d[14]:.4f}); the 1/n extrapolation {ext:.4f} is within 0.002 of the bracket",
       lo-0.002<=ext<=up+0.002)

# U7 limits: the bound reproduces D* when the repelling walls are very far, and the two-law threshold identity
lo,up,ds,_=bounds([(1.0,1e7,[0.2,0.15]),(-1.0,-1e7,[0.25])],nb=3)
ok("U7 repelling walls at ±10^7: bracket closes on D* (line case)",abs(lo-ds)<1e-4 and abs(up-ds)<1e-4,f"{lo:.5f} {up:.5f} {ds:.5f}")
x1,x2=3.7,1/3.7
ok("U8 two laws: spectral radius of [[0,x1],[x2,0]] is 1 exactly when x1/(1+x1)+x2/(1+x2)=1",
   abs(max(abs(np.linalg.eigvals(np.array([[0,x1],[x2,0]]))))-1)<1e-12 and abs(x1/(1+x1)+x2/(1+x2)-1)<1e-12)

# U9 sweep: brackets rise toward D* as the walls move out
SW=[]
for q in (1.5,2,3,5,8,12,20,40,100):
    lo,up,ds,_=bounds([(1.0,q,[0.2,0.15]),(-1.0,-q,[0.25])],nb=4); SW.append((q,lo,up,ds))
ok("U9 sweep: the bracket rises with the distance of the repelling walls and stays below D*",
   all(SW[i][1]<SW[i+1][1] and SW[i][2]<SW[i+1][2] for i in range(len(SW)-1)) and all(u<d for q,l,u,d in SW),
   " ".join(f"{q}:[{l:.3f},{u:.3f}]" for q,l,u,d in SW))
json.dump(dict(brackets={str(k):v for k,v in B.items()},sweep=SW),open("offline.json","w"),indent=1)
print(f"{sum(res)}/{len(res)} pass")
