"""Round 5 (DRAFT): no accidental collisions, proved for a family by two-prime ping-pong (README section 8.2).
Criterion (ours): plane laws with rational walls and multipliers, and primes p_1..p_r such that
  (i)  every multiplier has negative p-adic valuation for every chosen prime (so every block contracts p-adically
       toward its REPELLING wall);
  (ii) within each law the valuation vectors (v_p1, ..., v_pr) of its multipliers are linearly independent;
  (iii) for every chosen prime there are radii making the ping-pong sets U_i (p-adic balls around the repelling walls)
       disjoint, with every law-i block sending all U_j and the base point x0 into U_i (checked exactly below).
Then different words of the trace monoid give different maps."""
from fractions import Fraction as F
import itertools, json
def v(x,p):
    x=F(x); n,d=x.numerator,x.denominator; k=0
    while n%p==0: n//=p; k+=1
    while d%p==0: d//=p; k-=1
    return k
def ab(x,p): return F(p)**(-v(x,p))
def pingpong_sets(walls,x0,p):
    """exact check of (iii) for one prime; returns (ok, radii rho_i in the w_i = (y-q_i)/(y-p_i) coordinate, ball radii)"""
    L=len(walls); w=lambda i,y: (y-walls[i][1])/(y-walls[i][0])
    rho=[max([ab(w(i,walls[j][1]),p) for j in range(L) if j!=i]+[ab(w(i,x0),p)])/p for i in range(L)]
    if any(r>=1 for r in rho): return False,rho,None
    r=[rho[i]*ab(walls[i][1]-walls[i][0],p) for i in range(L)]          # U_i = ball(q_i, r_i) because rho_i < 1
    ok=True
    for i in range(L):
        for j in range(L):
            if i!=j:
                ok&= r[j]<ab(walls[j][1]-walls[i][1],p) and r[j]<ab(walls[j][1]-walls[i][0],p)   # |w_i| constant on U_j
                ok&= max(r[i],r[j])<ab(walls[i][1]-walls[j][1],p)                                 # U_i, U_j disjoint
        ok&= ab(x0-walls[i][1],p)>r[i]                                                            # x0 outside U_i
    return ok,rho,r
def criterion(classes,primes,x0):
    walls=[(c[0],c[1]) for c in classes]
    c1=all(v(m,p)<0 for c in classes for m in c[2] for p in primes)
    c2=True
    for c in classes:
        V=[[v(m,p) for p in primes] for m in c[2]]
        if len(V)>len(primes): c2=False; continue
        # rank check by fraction Gaussian elimination
        Mx=[[F(x) for x in row] for row in V]; rk=0; cols=len(primes)
        for col in range(cols):
            piv=next((r for r in range(rk,len(Mx)) if Mx[r][col]!=0),None)
            if piv is None: continue
            Mx[rk],Mx[piv]=Mx[piv],Mx[rk]
            for r in range(len(Mx)):
                if r!=rk and Mx[r][col]!=0:
                    f=Mx[r][col]/Mx[rk][col]; Mx[r]=[a-f*b for a,b in zip(Mx[r],Mx[rk])]
            rk+=1
        c2&= rk==len(V)
    c3=[pingpong_sets(walls,x0,p) for p in primes]
    return c1,c2,all(x[0] for x in c3),c3
EXAMPLE=[(F(1),F(21,5),[F(5,18),F(5,24)]),(F(-1),F(-21,5),[F(7,24)])]
PRIMES=(2,3); X0=F(0)
if __name__=="__main__":
    c1,c2,c3,info=criterion(EXAMPLE,PRIMES,X0)
    print("criterion (i), (ii), (iii):",c1,c2,c3)
    for p,(ok,rho,r) in zip(PRIMES,info): print(f"  p={p}: rho={rho}, ball radii={r}")
