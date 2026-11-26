"""DRAFT: largest Floquet multiplier by segmented finite differences (monodromy = product of short-segment Jacobians)."""
import sys,json,numpy as np
sys.path.insert(0,'../three-body'); import tbpy
d=json.load(open('../three-body-word-formula/data/li-liao-695.json'))
d=sorted(d,key=lambda o:len(o['word']))[:240]
out=[]
for o in d:
    X0=tbpy.suvakov(o['v1'],o['v2']); T=o['T']; nseg=max(20,6*len(o['word'])); dt=T/nseg; h=1e-7
    M=np.eye(12); X=X0.copy(); ok=True; logdet=0
    try:
        for s in range(nseg):
            J=np.zeros((12,12))
            for j in range(12):
                e=np.zeros(12); e[j]=h; J[:,j]=(tbpy.flow(X+e,dt)-tbpy.flow(X-e,dt))/(2*h)
            M=J@M; nrm=np.abs(M).max(); M/=nrm; logdet+=np.log(nrm)*12 - 0  # keep scale
            X=tbpy.flow(X,dt)
            if s==0: scale=[]
            scale.append(np.log(nrm))
    except RuntimeError: continue
    lam=np.sort(np.abs(np.linalg.eigvals(M)))[::-1]; lmax=np.log(lam[0])+sum(scale)
    # symplectic check: product of all |eigs| should be 1 -> sum log = 0
    sumlog=float(np.sum(np.log(lam))+12*sum(scale))
    out.append(dict(cls=o['cls'],word=o['word'],Tstar=o['Tstar'],T=T,logmaxlam=float(lmax),sumlog=sumlog))
json.dump(out,open('stability.json','w')); print(len(out),'orbits; |sum log eig| max',max(abs(o['sumlog']) for o in out))
