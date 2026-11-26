"""DRAFT: order readings of three-body words (no dynamics). Convention (checked from bigrams): lower-case a,b turn one way,
capitals the other; a run = alternating letters of one case. Letters as matrices:
a -> P=[[1,2],[0,1]], b -> Q=[[1,0],[-2,1]] at t=-1 (pure braid generators sigma1^2, sigma2^2 under reduced Burau)."""
import numpy as np, json, math, sys
sys.path.insert(0,'.')
def runs(w):
    n=len(w); cyc=w
    brk=[i for i in range(n) if cyc[i].islower()!=cyc[(i+1)%n].islower()]
    if not brk: return [n]
    out=[]
    for j in range(len(brk)):
        out.append((brk[j]-brk[j-1])%n or n)
    return out
def burau(t):
    s1=np.array([[-t,1],[0,1]],complex); s2=np.array([[1,0],[t,-t]],complex)
    P=s1@s1; Q=s2@s2
    return {'a':P,'b':Q,'A':np.linalg.inv(P),'B':np.linalg.inv(Q)}
def mat(w,t):
    G=burau(t); M=np.eye(2,dtype=complex)
    for ch in w: M=M@G[ch]
    return M
def rho(w,t):
    return max(abs(np.linalg.eigvals(mat(w,t))))
def trace_m1(w):
    # exact integer trace at t=-1
    P=np.array([[1,2],[0,1]],dtype=object); Q=np.array([[1,0],[-2,1]],dtype=object)
    G={'a':P,'b':Q,'A':np.array([[1,-2],[0,1]],dtype=object),'B':np.array([[1,0],[2,1]],dtype=object)}
    M=np.array([[1,0],[0,1]],dtype=object)
    for ch in w: M=M.dot(G[ch])
    return int(M[0,0]+M[1,1])
def walk_area(w):
    step={'a':(1,0),'A':(-1,0),'b':(0,1),'B':(0,-1)}
    x=y=0; area=0
    for ch in w:
        dx,dy=step[ch]; area+=x*dy-y*dx; x+=dx; y+=dy
    return x,y,area/2
def readings(w):
    tr=trace_m1(w); lam=(abs(tr)+math.sqrt(max(tr*tr-4,0)))/2
    x,y,ar=walk_area(w)
    return dict(L=len(w),runs=runs(w),tr=tr,ell=2*math.log(lam) if lam>0 else 0.0,netA=x,netB=y,area=ar)
if __name__=='__main__':
    d=json.load(open('../three-body-word-formula/data/li-liao-695.json'))
    out=[]
    for o in d:
        r=readings(o['word']); r.update(cls=o['cls'],word=o['word'],Tstar=o['Tstar'],v1=o['v1'],v2=o['v2'],T=o['T'])
        out.append(r)
    json.dump(out,open('words.json','w'))
    import collections
    print('net letters distribution',collections.Counter((r['netA'],r['netB']) for r in out).most_common(5))
    print('area distribution',collections.Counter(r['area'] for r in out).most_common(8))
    print('figure-eight',[r for r in out if r['cls']=='I.A1'][0])
    print('traces sample',[(r['L'],r['tr']) for r in out[:12]])
