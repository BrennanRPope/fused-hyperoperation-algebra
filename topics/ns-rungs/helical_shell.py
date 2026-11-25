"""Helical shell model with the cross-handed (heterochiral) triads on a dial eps.
Shells k_n = phi^n (golden ratio), two helicities u_n^+ and u_n^-.  For each triad of neighbouring shells
(n, n+1, n+2) and each sign triple (s0,s1,s2) -- the 8 corners of the triad-sign cube --
  dX/dt += i g k_n c0 conj(Y) Z,  dY/dt += i g k_n c1 conj(X) Z,  dZ/dt += i g k_n c2 X Y
  c0 = (s1*l - s2*l^2) G,  c1 = (s2*l^2 - s0) G,  c2 = c0 + c1,  G = s0 + s1*l + s2*l^2   (l = phi)
These follow Waleffe's helical triad form (s_p p - s_q q)(s_k k + s_p p + s_q q); they conserve
energy E = sum |u|^2 and helicity H = sum s k |u^s|^2 for every corner separately.
Weight g = 1 on the two homochiral corners (+++), (---) [cube layers 0 and 3] and
g = eps on the six heterochiral corners [cube layers 1 and 2]. Inviscid, free decay.
Arrival time t_n: first time shell n's energy reaches theta * k_n^(-2/3) (the Kolmogorov level).
Converging t_n = energy reaches infinitely small scales in finite time (blow-up signature)."""
import numpy as np, itertools, sys, json
from scipy.integrate import solve_ivp
lam=(1+5**0.5)/2   # golden ratio: k_n + k_{n+1} = k_{n+2}, shells form (flat) triangles
def build(N, eps):
    corners=[]
    for s in itertools.product((1,-1),repeat=3):
        s0,s1,s2=s
        G=s0+s1*lam+s2*lam**2            # Waleffe geometric factor (s_k k + s_p p + s_q q)/k_n
        c0=(s1*lam-s2*lam**2)*G; c1=(s2*lam**2-s0)*G; c2=c0+c1
        g=(1.0 if s0==s1==s2 else eps) if eps!=float('inf') else (0.0 if s0==s1==s2 else 1.0)
        if g: corners.append((s,c0*g,c1*g,c2*g))
    return corners
def check_conservation():
    res=[]
    for s in itertools.product((1,-1),repeat=3):
        s0,s1,s2=s; c0=s1*lam-s2*lam**2; c1=s2*lam**2-s0; c2=c0+c1
        res.append((s, c0+c1-c2, s0*1*c0+s1*lam*c1-s2*lam**2*c2))
    return res
def make_rhs(N,eps):
    k=lam**np.arange(N); corners=build(N,eps); idx={1:0,-1:1}
    def rhs(t,y):
        u=(y[:2*N]+1j*y[2*N:]).reshape(2,N); du=np.zeros_like(u)
        for (s0,s1,s2),c0,c1,c2 in corners:
            X=u[idx[s0],:-2]; Y=u[idx[s1],1:-1]; Z=u[idx[s2],2:]; kk=k[:-2]
            du[idx[s0],:-2]+=1j*kk*c0*np.conj(Y)*Z
            du[idx[s1],1:-1]+=1j*kk*c1*np.conj(X)*Z
            du[idx[s2],2:]+=1j*kk*c2*X*Y
        du=du.ravel(); return np.concatenate([du.real,du.imag])
    return rhs,k
def init(N,seed=0):
    r=np.random.default_rng(seed); u=np.zeros((2,N),complex)
    amp=np.zeros(N); amp[:3]=[1.0,0.7,0.5]
    u=amp*np.exp(2j*np.pi*r.random((2,N)))+1e-12*(r.standard_normal((2,N))+1j*r.standard_normal((2,N)))*(np.arange(N)>=3)
    return np.concatenate([u.ravel().real,u.ravel().imag])
def run(eps,N=30,tmax=12.0,theta=1e-2,seed=0):
    rhs,k=make_rhs(N,eps); y0=init(N,seed)
    E=lambda y: np.sum(y**2)
    def Hel(y):
        u=(y[:2*N]+1j*y[2*N:]).reshape(2,N); return np.sum(k*(abs(u[0])**2-abs(u[1])**2))
    def Hplus(y):
        u=(y[:2*N]+1j*y[2*N:]).reshape(2,N); return np.sum(k*abs(u[0])**2)
    target=N-3
    def ev(t,y):
        u=(y[:2*N]+1j*y[2*N:]).reshape(2,N); return np.sum(abs(u[:,target])**2)-theta*k[target]**(-2/3)
    ev.terminal=True; ev.direction=1
    sol=solve_ivp(rhs,(0,tmax),y0,method='DOP853',rtol=1e-9,atol=1e-14,events=ev,dense_output=False,max_step=0.02)
    U=(sol.y[:2*N]+1j*sol.y[2*N:]).reshape(2,N,-1); En=np.sum(abs(U)**2,axis=0)  # (N, T)
    tn=[]
    for n in range(N):
        hit=np.where(En[n]>=theta*k[n]**(-2/3))[0]; tn.append(float(sol.t[hit[0]]) if len(hit) else None)
    Z=np.sum((k[:,None]**2)*En,axis=0)
    return dict(eps=eps,t=sol.t,tn=tn,Z=Z,E0=E(y0),E1=E(sol.y[:,-1]),H0=Hel(y0),H1=Hel(sol.y[:,-1]),
                Hp0=Hplus(y0),Hp1=Hplus(sol.y[:,-1]),tend=float(sol.t[-1]),status=sol.status,En=En,k=k)
if __name__=="__main__":
    for s,dE,dH in check_conservation(): print("corner",s," energy coeff sum",dE," helicity coeff sum",dH)
    eps_list=[float(a) for a in sys.argv[1:]] or [0.0,0.1,0.3,1.0]
    results={}
    for eps in eps_list:
        r=run(eps)
        print(f"eps={eps:5.2f}  t_end={r['tend']:.4f} status={r['status']}  E {r['E0']:.6f}->{r['E1']:.6f}  H {r['H0']:.5f}->{r['H1']:.5f}  H+ {r['Hp0']:.4f}->{r['Hp1']:.4f}  Zmax={r['Z'].max():.3e}")
        print("   arrival t_n:", [None if v is None else round(v,4) for v in r['tn']])
        results[eps]={k_:v for k_,v in r.items() if k_ in('tn','tend','E0','E1','H0','H1','Hp0','Hp1')}
        np.savez(f"shell_eps{eps:.2f}.npz" if eps!=float('inf') else 'shell_epsinf.npz',t=r['t'],Z=r['Z'],En=r['En'],k=r['k'])
    json.dump({str(a):b for a,b in results.items()},open('shell_results_'+'_'.join(f'{a:.2f}' for a in eps_list)+'.json','w'),indent=1)
