# DRAFT. Round 17: does the uncertainty rule for 50/50 readings bring out the flips' sign (XZ = -ZX)?
# Known: Maassen & Uffink, PRL 60, 1103 (1988): H(A) + H(B) >= -2 log2 c, c = largest overlap between the two sets of readings.
# Flips A = n.sigma, B = m.sigma with angle phi between n and m.
import json, numpy as np
rng=np.random.default_rng(17); res={}
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
flip=lambda n: n[0]*X+n[1]*Y+n[2]*Z
def H(p,q=1.0):
    p=np.clip(p,1e-300,1); return -np.sum(p*np.log2(p)) if q==1.0 else (1-np.sum(p**q))/(q-1)
def probs(A,psi): v=np.linalg.eigh(A)[1]; return np.abs(v.conj().T@psi)**2
def min_sum(A,B,q=1.0,N=4000):
    th=np.arccos(rng.uniform(-1,1,N)); ph=rng.uniform(0,2*np.pi,N)
    best=min(H(probs(A,s),q)+H(probs(B,s),q) for s in np.stack([np.cos(th/2),np.exp(1j*ph)*np.sin(th/2)],1))
    return best
phis=np.radians(np.arange(0,181,10)); n=np.array([0,0,1.])
def m_at(phi): return np.array([np.sin(phi),0,np.cos(phi)])
MU=lambda phi: -np.log2(max(np.cos(phi/2)**2,np.sin(phi/2)**2))
curve={q:[min_sum(flip(n),flip(m_at(p)),q) for p in phis] for q in (1.0,2.0)}
# U1: the Maassen-Uffink bound holds at every angle, and is reached exactly at 90 degrees (1 bit).
#     Between 0 and 90 it is loose for qubits, as known (Ghirardi, Marinatto & Romano 2003 give the tight curve).
res['U1']=all(c>=MU(p)-1e-9 for c,p in zip(curve[1.0],phis)) and abs(curve[1.0][9]-1)<2e-3
# U2: the bound reaches its top value, 1 full bit, only at 90 degrees
mus=[MU(p) for p in phis]; res['U2']=abs(max(mus)-1)<1e-12 and [round(np.degrees(p)) for p,u in zip(phis,mus) if u>1-1e-9]==[90]
# U3: exactly: AB + BA = 2 cos(phi) I, so the top uncertainty is the same thing as the anticommuting sign
res['U3']=max(np.abs(flip(n)@flip(m_at(p))+flip(m_at(p))@flip(n)-2*np.cos(p)*np.eye(2)).max() for p in phis)<1e-12
# U4: the law's dial doesn't move it: with Tsallis q = 2 (k = -1) the smallest total uncertainty also peaks at 90 degrees
res['U4']=round(np.degrees(phis[int(np.argmax(curve[2.0]))]))==90 and round(np.degrees(phis[int(np.argmax(curve[1.0]))]))==90
# U5: the sign alone does not need i: X and Z are real, anticommute, and reach the full 1 bit
res['U5']=np.allclose(X.imag,0) and np.allclose(Z.imag,0) and np.allclose(X@Z,-Z@X) and abs(min_sum(X,Z)-1)<2e-3
# U6: real flips (2x2 real symmetric, trace 0) span only 2 directions, so at most 2 can be pairwise 50/50;
#     complex flips span 3. A third 50/50 question on a bit needs i.
real_basis=[np.array([[1,0],[0,-1]]),np.array([[0,1],[1,0]])]
def dim_span(mats): return np.linalg.matrix_rank(np.array([np.r_[m.real.flatten(),m.imag.flatten()] for m in mats]))
rr=[]
for _ in range(200):     # random real flips: any three are linearly dependent, so no three pairwise anticommute
    v=[a*real_basis[0]+b*real_basis[1] for a,b in rng.normal(size=(3,2))]; rr.append(dim_span(v))
res['U6']=max(rr)==2 and dim_span([X,Y,Z])==3 and np.allclose(X@Y,-Y@X) and np.allclose(Y@Z,-Z@Y)
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({'passed':{k:bool(v) for k,v in res.items()},'angles_deg':[int(round(np.degrees(p))) for p in phis],
           'min_shannon':[float(c) for c in curve[1.0]],'min_tsallis2':[float(c) for c in curve[2.0]],'MU_bound':[float(u) for u in mus]},open('results.json','w'),indent=1)
