"""Self-feeding (rung-2 squared) growth: the restricted Euler model of the velocity gradient A
(Vieillefosse 1982):  dA/dt = -A^2 + (tr(A^2)/d) I,  trace-free A, d = 2 or 3.
In 2D the right side is identically 0 (A^2 = -det(A) I for trace-free 2x2): no self-feeding.
In 3D it blows up in finite time like 1/(T-t), the Riccati pole of Y' = Y^2."""
import numpy as np
from scipy.integrate import solve_ivp
rng=np.random.default_rng(3)
def rhs(t,y,d):
    A=y.reshape(d,d); A2=A@A; return (-A2+np.trace(A2)/d*np.eye(d)).ravel()
out=[]
A=rng.standard_normal((2,2)); A-=np.trace(A)/2*np.eye(2)
r=rhs(0,A.ravel(),2); out.append(f"2D: |dA/dt| at a random trace-free A = {np.abs(r).max():.1e}  (identically zero)")
A=rng.standard_normal((3,3)); A-=np.trace(A)/3*np.eye(3)
ev=lambda t,y,d: np.abs(y).max()-1e6; ev.terminal=True
s=solve_ivp(rhs,(0,50),A.ravel(),args=(3,),rtol=1e-11,atol=1e-12,events=ev,dense_output=True)
T=s.t[-1]; out.append(f"3D: |A| reaches 1e6 at t = {T:.6f}")
tt=T-np.logspace(-1,-5,9); nrm=[np.linalg.norm(s.sol(t)) for t in tt]
slope=np.polyfit(np.log(T-tt),np.log(nrm),1)[0]
out.append(f"3D: log|A| vs log(T-t) slope = {slope:.4f}  (Riccati pole Y'=Y^2 gives -1)")
print("\n".join(out)); open('restricted_euler_output.txt','w').write("\n".join(out)+"\n")
np.savez('restricted_euler.npz',t=s.t,n=np.linalg.norm(s.y,axis=0),T=T)
