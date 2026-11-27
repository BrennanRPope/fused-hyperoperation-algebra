"""Round 5 (DRAFT): airtight brackets for delta (README section 8).
Same bound as offline.py (section 7.3), redone so that every number is certified:
  - walls, multipliers, block multipliers, input regions and z-ranges are exact fractions (regions rounded outward);
  - powers and sums are done in ball arithmetic (python-flint arb, 128 bits), which returns guaranteed enclosures;
  - the spectral-radius claims are certified with a Collatz-Wielandt vector v > 0:
      T_up(t_up) v < v componentwise   =>  spectral radius < 1  =>  the Poincare series converges at t_up  =>  delta <= t_up
      T_low(t_low) v > v componentwise =>  spectral radius > 1  =>  it diverges at t_low                =>  delta >= t_low
No floating-point number enters a certified inequality; floats are only used to *guess* t and v."""
import itertools, math
from fractions import Fraction as F
import numpy as np
import flint
from flint import arb, ctx
ctx.prec=128

def zc(y,p,q): return (y-p)/(y-q)
def yc(z,p,q): return (p-q*z)/(1-z)
def down(x,D=10**12): return F(math.floor(x*D),D)
def up(x,D=10**12): return F(math.ceil(x*D),D)
def A(x): return arb(x.numerator)/arb(x.denominator)

def blocks(mults,N):
    out=[]
    for deg in range(1,N+1):
        for combo in itertools.combinations_with_replacement(range(len(mults)),deg):
            lam=F(1)
            for c in combo: lam*=mults[c]
            out.append((deg,lam))
    return out

class Setup:
    def __init__(self, classes, nb=4, N=40, iters=40):
        """classes: list of (p, q, [multipliers]) as Fractions."""
        self.cl=classes; self.nb=nb; self.L=L=len(classes)
        ps=[c[0] for c in classes]; I=(min(ps),max(ps))
        for p,q,ms in classes:
            assert not (I[0]<=q<=I[1]) and all(0<m<1 for m in ms)
        self.blk=[blocks(c[2],N) for c in classes]
        b=lambda deg: min(deg,nb)-1; self.bucket=b
        self.lamrng=[[(F(0) if k==nb-1 else min(l for d,l in self.blk[i] if b(d)==k), max(l for d,l in self.blk[i] if b(d)==k))
                      for k in range(nb)] for i in range(L)]
        def zint(i,lo,hi):
            p,q,_=classes[i]; zs=[zc(lo,p,q),zc(hi,p,q)]; return min(zs),max(zs)
        def yint(i,lo,hi):
            p,q,_=classes[i]; ys=[yc(lo,p,q),yc(hi,p,q)]; return min(ys),max(ys)
        self.zint=zint
        inp=[I]*L
        for _ in range(iters):
            out=[[None]*nb for _ in range(L)]
            for j in range(L):
                zlo,zhi=zint(j,*inp[j]); assert zhi<1
                for k in range(nb):
                    lo,hi=self.lamrng[j][k]; c=[lo*zlo,lo*zhi,hi*zlo,hi*zhi]
                    a,bb=yint(j,min(c),max(c)); out[j][k]=(down(a),up(bb))   # outward rounding keeps it an outer bound
            inp=[(min(out[j][k][0] for j in range(L) if j!=i for k in range(nb)),
                  max(out[j][k][1] for j in range(L) if j!=i for k in range(nb))) for i in range(L)]
            inp=[(max(v[0],I[0]),min(v[1],I[1])) for v in inp]
        self.inp=inp; self.out=out
        # z-ranges of each law's input given the previous state, exact
        self.zr={}
        for j in range(L):
            for k in range(nb):
                for i in range(L):
                    if i!=j:
                        z1,z2=zint(i,*out[j][k]); assert z2<1
                        self.zr[(j,k,i)]=(z1,z2)
    def T(self,t,which):
        """which='up': upper weights (smallest z) plus a tail bound; 'low': lower weights (largest z), truncated."""
        L,nb=self.L,self.nb; n=L*nb; M=[[arb(0)]*n for _ in range(n)]
        tt=A(t) if isinstance(t,F) else arb(t)
        for (j,k,i),(z1,z2) in self.zr.items():
            z=z1 if which=='up' else z2
            acc=arb(0); accb=[arb(0)]*nb
            for deg,lam in self.blk[i]:
                base=lam*((1-z)/(1-lam*z))**2; assert base>0 and (which=='low' or base<1)
                accb[self.bucket(deg)]+=A(base)**tt
            if which=='up':
                cap=max(F(1),(1-z)**2)
                full=arb(1)
                for s in self.cl[i][2]: full*=1/(1-A(s)**tt)
                full-=1
                part=arb(0)
                for deg,lam in self.blk[i]: part+=A(lam)**tt
                tail=full-part
                accb[nb-1]+=A(cap)*arb(tail.upper()) if tail.upper()>0 else arb(0)
            for b2 in range(nb): M[i*nb+b2][j*nb+k]=accb[b2]
        return M
    def floatT(self,t,which):
        M=self.T(F(t).limit_denominator(10**6),which); return np.array([[float(x.mid()) for x in r] for r in M])
    def sr(self,t,which): return max(abs(np.linalg.eigvals(self.floatT(t,which))))
    def certify(self,t,which):
        """returns (ok, min margin). v from the float Perron vector, made rational and positive."""
        Mf=self.floatT(t,which); w,V=np.linalg.eig(Mf); v=np.abs(np.real(V[:,np.argmax(abs(w))]))
        v=[F(max(x,1e-9)).limit_denominator(10**9) for x in v]
        M=self.T(t,which); n=len(v); ok=True; margins=[]
        for r in range(n):
            s=arb(0)
            for c in range(n): s+=M[r][c]*A(v[c])
            if which=='up':  ok&= bool(s<A(v[r])); margins.append(float((A(v[r])-s).lower()/A(v[r])))
            else:            ok&= bool(s>A(v[r])); margins.append(float((s-A(v[r])).lower()/A(v[r])))
        return ok, min(margins)
    def bracket(self, digits=4):
        from scipy.optimize import brentq
        lo=brentq(lambda t: self.sr(t,'low')-1,1e-3,3); hi=brentq(lambda t: self.sr(t,'up')-1,1e-3,3)
        step=F(1,10**digits)
        tl=F(math.floor(lo*10**digits),10**digits)-step; tu=F(math.ceil(hi*10**digits),10**digits)+step
        okl,ml=self.certify(tl,'low'); oku,mu=self.certify(tu,'up')
        return dict(t_low=tl,t_up=tu,low_ok=okl,up_ok=oku,low_margin=ml,up_margin=mu,float_low=lo,float_up=hi)

def mu_star(classes,t):
    """the shared-centre function 1 + sum[prod(1 - s^t) - 1] in ball arithmetic (D* is its root)"""
    tt=A(t); tot=arb(1)
    for p,q,ms in classes:
        pr=arb(1)
        for s in ms: pr*=1-A(s)**tt
        tot+=pr-1
    return tot

if __name__=="__main__":
    res={}
    for q in (2,3,12):
        cl=[(F(1),F(q),[F(1,5),F(3,20)]),(F(-1),F(-q),[F(1,4)])]
        S=Setup(cl,nb=6); r=S.bracket()
        # strictly below D*: mu_star < 0 means t < D*? mu is increasing through its root; check sign at t_up
        m=mu_star(cl,r['t_up'])
        print(f"walls ±{q}: certified {r['t_low']} <= delta <= {r['t_up']}  (low ok {r['low_ok']}, up ok {r['up_ok']}, margins {r['low_margin']:.2e}, {r['up_margin']:.2e});"
              f" shared-centre function at t_up: {m}")
