# DRAFT. Round 14: can the law's flip (k = -2, XOR, lift 1 - 2x) be built from smooth turns?
# Setup: the lift L = 1 + k x turns the law into multiplication: L(a o b) = L(a) L(b). At k = -2, x = 1 (flip) has L = -1.
import json, numpy as np
from scipy.linalg import expm
rng=np.random.default_rng(14); res={}
law=lambda a,b,k: a+b+k*a*b
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex)
xc=lambda t: (1-np.exp(1j*t))/2        # complex weight with lift e^{i t}
# T1: with real weights, any path from x = 0 to x = 1 crosses x = 1/2, where the lift is 0 (the step erases; it can't be undone)
lift_real=1-2*np.linspace(0,1,100001)
res['T1']=np.min(np.abs(lift_real))<1e-4 and np.all(np.diff(np.sign(lift_real))<=0)
# T2: complex weights x(t) compose by the law exactly like turn angles add
res['T2']=max(abs(law(xc(a),xc(b),-2)-xc(a+b)) for a,b in rng.uniform(-7,7,(500,2)))<1e-12
# T3: the law's step (1 - x) I + x X is unitary and equals e^{i t/2} exp(-i t X / 2), a qubit turn by t about one axis
d=[]
for t in rng.uniform(-7,7,300):
    U=(1-xc(t))*I2+xc(t)*X; d.append(max(np.abs(U.conj().T@U-I2).max(),np.abs(U-np.exp(1j*t/2)*expm(-1j*t*X/2)).max()))
res['T3']=max(d)<1e-12
# T4: flip odds |x|^2 = sin^2(t/2), the textbook Rabi formula; t = pi is the full flip x = 1
ts=np.linspace(0,2*np.pi,721); res['T4']=np.allclose(np.abs(xc(ts))**2,np.sin(ts/2)**2) and abs(xc(np.pi)-1)<1e-15
# T5: unitarity for every t pins the dial. With x = (e^{i t} - 1)/k (lift = e^{i t}), the eigenvalue 1 - 2x has size 1 for all t only at k = -2
def dev(k): x=(np.exp(1j*ts)-1)/k; return np.max(np.abs(np.abs(1-2*x)-1))
grid=[complex(a,b) for a in np.linspace(-4,4,161) for b in np.linspace(-4,4,161) if abs(complex(a,b))>1e-9]
devs={k:dev(k) for k in grid}; best=min(devs,key=devs.get)
res['T5']=abs(best-(-2))<1e-9 and devs[best]<1e-12 and sorted(devs.values())[1]>1e-2
# T6: at one full turn (t = 2 pi) the law's step is the identity, but the spin part exp(-i t X/2) is -I:
#     the 720-degree sign sits in the factor e^{i t/2} that the law's form absorbs
t=2*np.pi; U=(1-xc(t))*I2+xc(t)*X
res['T6']=np.allclose(U,I2) and np.allclose(expm(-1j*t*X/2),-I2)
# T7: the law is commutative, so every group of turns it makes is commutative; turns about two axes do not commute
Z=np.diag([1,-1]).astype(complex); A=expm(-1j*.7*X/2); B=expm(-1j*.7*Z/2)
pairs=rng.normal(size=(100,4)); pairs=[(complex(p[0],p[1]),complex(p[2],p[3])) for p in pairs]
res['T7']=max(abs(law(a,b,-2)-law(b,a,-2)) for a,b in pairs)<1e-14 and np.abs(A@B-B@A).max()>0.1
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
print(f"T5: best dial on a complex grid = {best}, size-1 error {devs[best]:.1e}; next best error {sorted(devs.values())[1]:.3f}")
json.dump({'passed':{k:bool(v) for k,v in res.items()},'T5_best_k':str(best)},open('results.json','w'),indent=1)
