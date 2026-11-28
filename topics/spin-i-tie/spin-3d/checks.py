# DRAFT. Round 16: do two of the law's flips that don't commute glue into full 3D spin?
# Each axis uses the round-14 step at k = -2: U_A(t) = (1 - x) I + x A, x = (1 - e^{it})/2, for a flip A with A^2 = I.
import json, numpy as np
from scipy.linalg import expm
rng=np.random.default_rng(16); res={}
I=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],complex); Z=np.diag([1,-1]).astype(complex)
xs=lambda t:(1-np.exp(1j*t))/2
U=lambda A,t:(1-xs(t))*np.eye(len(A))+xs(t)*A
# S1: each flip undoes itself (the XOR law), and the two flips anticommute: XZ = -ZX
res['S1']=np.allclose(X@X,I) and np.allclose(Z@Z,I) and np.allclose(X@Z,-Z@X)
# S2: the flip made from both, ZX, squares to -1, so it is a quarter turn, not a flip. Multiplying by i makes it a flip again (the usual Y).
Y=-1j*Z@X
res['S2']=np.allclose((Z@X)@(Z@X),-I) and np.allclose(Y@Y,I) and np.allclose(Y,np.array([[0,-1j],[1j,0]]))
# S3: law steps along Z, then X, then Z reach every 2x2 turn (up to an overall phase): 300 random turns rebuilt
def haar():
    q,r=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))); return q*(np.diag(r)/abs(np.diag(r)))
def euler(V):
    V=V/np.sqrt(np.linalg.det(V)); a,b=V[0,0],V[0,1]    # V = [[a, b], [-b*, a*]] = e^{-i al Z/2} e^{-i th X/2} e^{-i ga Z/2}
    s1=-2*np.angle(a); s2=-2*np.angle(1j*b); th=2*np.arctan2(abs(b),abs(a))
    return (s1+s2)/2, th, (s1-s2)/2, V
err=[]
for _ in range(300):
    V=haar(); al,th,ga,Vs=euler(V)
    W=expm(-1j*al*Z/2)@expm(-1j*th*X/2)@expm(-1j*ga*Z/2)
    L=U(Z,al)@U(X,th)@U(Z,ga)            # law steps differ from W only by an overall phase
    ph=np.vdot(L.flatten(),V.flatten()); ph/=abs(ph)
    err.append(max(min(np.abs(W-Vs).max(),np.abs(W+Vs).max()), np.abs(L*ph-V).max()))
res['S3']=max(err)<1e-9
# S4: if the two flips commute instead (two separate bits, X on one and X on the other), every product commutes: no 3D turns
A=np.kron(X,I); B=np.kron(I,X)
els=[U(A,t)@U(B,s) for t,s in rng.uniform(-7,7,(20,2))]
res['S4']=np.allclose(A@B,B@A) and max(np.abs(e@f-f@e).max() for e in els for f in els)<1e-12
# S5: turning by t about X turns the spin's direction by t (Z goes to cos t Z - sin t Y), while the state turned by t/2: the double cover
d=[]
for t in rng.uniform(-7,7,100):
    R=expm(-1j*t*X/2); d.append(np.abs(R@Z@R.conj().T-(np.cos(t)*Z-np.sin(t)*Y)).max())
R=expm(-1j*2*np.pi*X/2)
res['S5']=max(d)<1e-12 and np.allclose(R,-I) and np.allclose(R@R,I)
# S6: the two law generators X, Z, by their commutator, produce the third axis with an i: [X, Z] = -2i Y. Together they span all 3 directions.
comm=X@Z-Z@X
gens=[1j*X,1j*Z,comm]; M=np.array([g.flatten() for g in gens])
res['S6']=np.allclose(comm,-2j*Y) and np.linalg.matrix_rank(np.hstack([M.real,M.imag]))==3
# S7: for qubit flips n.sigma, anticommuting is the same as having mutually unbiased readings
#     (every outcome of one leaves the other at 50/50): checked on 500 random pairs plus 200 perpendicular ones
def flip(n): n=n/np.linalg.norm(n); return n[0]*X+n[1]*Y+n[2]*Z
ok=True
for i in range(700):
    n=rng.normal(size=3); m=rng.normal(size=3)
    if i>=500: m=np.cross(n,m)
    A,B=flip(n),flip(m); va=np.linalg.eigh(A)[1]; vb=np.linalg.eigh(B)[1]
    mub=np.allclose(np.abs(va.conj().T@vb)**2,0.5,atol=1e-9); anti=np.allclose(A@B,-B@A,atol=1e-9)
    ok&= (mub==anti); ok&= (anti if i>=500 else True)
res['S7']=ok
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({k:bool(v) for k,v in res.items()},open('results.json','w'),indent=1)
