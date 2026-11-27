"""Round 4 (DRAFT): the true dimension off the tangent line.
Proved two-sided bounds for delta, the critical exponent of the Poincare series over the trace monoid,
for plane-law Mobius maps. Theorem P3 (README section 7) says dim K = min(1, delta).

Coordinates: for law i with attracting wall p and repelling wall q, z = (y - p)/(y - q). A block B of law i
(any composite of law-i maps, multiplier lam_B = product of multipliers) acts as z -> lam_B z, and its derivative at y is
   B'(y) = lam_B * ((1 - z)/(1 - lam_B z))^2 .                                   (exact, checked in checks4.py)
Words of the trace monoid = sequences of nonempty blocks with consecutive blocks from different laws (unique).
The derivative of a word at a point is the product of its blocks' derivatives at the points the word passes through.
Bounding each block's factor above / below by the extreme z over the set of points it can receive gives
transfer matrices T_up(t) >= true weights >= T_low(t), so  D_low <= delta <= D_up  where rho(T(D)) = 1.
States remember the previous block's law and a multiplier bucket, which narrows where the next block's input lies."""
import itertools, math, json
import numpy as np
from scipy.optimize import brentq

def zc(y,p,q): return (y-p)/(y-q)
def yc(z,p,q): return (p-q*z)/(1-z)
def rho(z,lam): return ((1-z)/(1-lam*z))**2

def blocks(mults, N):
    """nonempty multisets of the generators of one law up to total degree N: list of (degree, lam)"""
    out=[]; m=len(mults)
    for deg in range(1,N+1):
        for combo in itertools.combinations_with_replacement(range(m),deg):
            out.append((deg, math.prod(mults[c] for c in combo)))
    return out

def bounds(classes, nb=4, N=40):
    """classes: list of (p, q, [multipliers]); all repelling walls outside the hull of attracting walls.
    nb: multiplier buckets by block degree 1..nb-1 and >= nb. Returns (D_low, D_up, D_star)."""
    L=len(classes); ps=[c[0] for c in classes]; I=(min(ps),max(ps))
    for p,q,_ in classes: assert not (I[0] < q < I[1]) and q not in I, "repelling wall inside the hull"
    blk=[blocks(c[2],N) for c in classes]
    def bucket(deg): return min(deg,nb)-1
    # lam range per (law, bucket) among ALL blocks (degree unbounded): bucket nb-1 has lam in (0, max at degree nb]
    lamrng=[[None]*nb for _ in range(L)]
    for i in range(L):
        for b in range(nb):
            ls=[lam for deg,lam in blk[i] if bucket(deg)==b]
            lamrng[i][b]=(0.0 if b==nb-1 else min(ls), max(ls))
    # where can the input of a law-i block lie?  y-interval, outer approximation, iterate to a fixed point
    def yint_from_z(i, zlo, zhi):
        p,q,_=classes[i]; ys=[yc(zlo,p,q), yc(zhi,p,q)]; return (min(ys),max(ys))
    def zint_from_y(i, ylo, yhi):
        p,q,_=classes[i]; zs=[zc(ylo,p,q), zc(yhi,p,q)]; return (min(zs),max(zs))
    inp=[I for _ in range(L)]       # inputs to law-i blocks (y-interval)
    for _ in range(60):
        out=[[None]*nb for _ in range(L)]
        for j in range(L):
            zlo,zhi=zint_from_y(j,*inp[j])
            for b in range(nb):
                lo,hi=lamrng[j][b]
                cand=[lo*zlo,lo*zhi,hi*zlo,hi*zhi]
                out[j][b]=yint_from_z(j,min(cand),max(cand))
        new=[]
        for i in range(L):
            ivs=[out[j][b] for j in range(L) if j!=i for b in range(nb)]
            new.append((min(v[0] for v in ivs), max(v[1] for v in ivs)))
        inp=new
    # transfer matrices: state (j,b) = previous block's law and bucket; next block B of law i != j
    def T(t, which):
        n=L*nb; M=np.zeros((n,n))
        for j in range(L):
            for b in range(nb):
                ylo,yhi=out[j][b]
                for i in range(L):
                    if i==j: continue
                    z1,z2=zint_from_y(i,ylo,yhi)
                    z=z1 if which=='up' else z2          # rho decreasing in z
                    for deg,lam in blk[i]:
                        M[i*nb+bucket(deg), j*nb+b]+= (lam*rho(z,lam))**t
                    if which=='up':                       # tail beyond degree N, bounded by the pure-multiplier tail
                        cap=max(1.0,(1-z)**2)
                        full=np.prod([1/(1-s**t) for s in classes[i][2]])-1
                        part=sum(lam**t for deg,lam in blk[i])
                        M[i*nb+nb-1, j*nb+b]+= cap*max(full-part,0.0)
        return M
    def sr(t,which): return max(abs(np.linalg.eigvals(T(t,which))))
    Dlow=brentq(lambda t: sr(t,'low')-1, 1e-3, 3)
    Dup=brentq(lambda t: sr(t,'up')-1, 1e-3, 3)
    # D*: same with rho = 1 (pure multipliers), two-class formula check
    def srstar(t):
        x=[np.prod([1/(1-s**t) for s in c[2]])-1 for c in classes]
        return sum(xi/(1+xi) for xi in x)-1
    Dstar=brentq(srstar,1e-3,3)
    return Dlow, Dup, Dstar, dict(inputs=inp, outputs=out, lamrng=lamrng, T=T, sr=sr)

if __name__=="__main__":
    s1,s2,s3=1/5,3/20,1/4
    res={}
    PL=json.load(open("planelaws.json")); SW=json.load(open("planelaws_sweep.json"))
    cases=[("no shared wall, close",2),("no shared wall, near",3),("no shared wall, far",12)]
    for name,q in cases:
        for nb in (1,3,6):
            lo,up,ds,_=bounds([(1.0,float(q),[s1,s2]),(-1.0,-float(q),[s3])],nb=nb)
            print(f"q=±{q} nb={nb}: {lo:.4f} <= delta <= {up:.4f}   D*={ds:.4f}   (measured box {PL[name]['box']:.4f}, D_n {PL[name]['Dn_wall1']:.4f})")
