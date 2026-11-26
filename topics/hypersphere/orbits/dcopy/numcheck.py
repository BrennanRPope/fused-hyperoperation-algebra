import mpmath as mp, json
mp.mp.dps=40
def advance(a,b,c,M=1,L=200,e=0.3):
    # pick E from Newtonian-ish turning points, then find exact roots of W
    p=L**2/(a*M); r1=p/(1+e); r2=p/(1-e)
    def W(u,E): return (1-2*c*M*u)*((1+2*a*M*u)*E**2+4*b**2*M**2*u**2*E**2/(1-2*c*M*u)-1-L**2*u**2)
    # choose E so that W(1/r1)=0
    E=mp.re(mp.findroot(lambda E: W(mp.mpf(1)/r1,E), mp.mpf(1)-a*M/(2*p)))
    u1=mp.mpf(1)/r1
    u2=mp.re(mp.findroot(lambda u: W(u,E), mp.mpf(1)/r2))
    lo,hi=min(u1,u2),max(u1,u2)
    mid=(lo+hi)/2; half=(hi-lo)/2
    f=lambda th: L*half*mp.sin(th)/mp.sqrt(abs(W(mid-half*mp.cos(th),E)))
    half_lap=mp.quad(f,[0,mp.pi])
    return float(2*half_lap-2*mp.pi)
out={}
for tag in [(1,1,1),(1,1,0),(1,0,1),(1,0,0),(1,0.5,1),(1,0.5,0),(2,1,1)]:
    d=advance(*tag); pred=mp.pi*(4*tag[1]**2/tag[0]+2*tag[2])/(200**2)*tag[0]**0
    # with a!=1, x~aM^2/L^2 so prediction pi*(4b^2+2ac)M^2/L^2
    pred=float(mp.pi*(4*tag[1]**2+2*tag[0]*tag[2])/200**2)
    out[str(tag)]=(d,pred)
    print(tag,"advance/lap=%.6e  predicted pi(4b^2+2ac)M^2/L^2=%.6e  ratio=%.4f"%(d,pred,d/pred if pred else float('nan')))
json.dump(out,open('numcheck.json','w'),indent=1)
