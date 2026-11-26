"""DRAFT: Floquet multipliers (finite-difference monodromy) of the 9 orbits not in the 2017 list."""
import sys, json, numpy as np
sys.path.insert(0, '/mnt/project-files/math-theory/three-body')
import tbpy
r = json.load(open('/mnt/project-files/math-theory/three-body/results.json'))['E']
rows = [e for e in r if not e['match']]
out=[]
for e in rows:
    X0 = tbpy.suvakov(e['p1'], e['p2']); T = e['T']
    M = np.zeros((12,12)); h=1e-6
    for j in range(12):
        d=np.zeros(12); d[j]=h
        M[:,j]=(tbpy.flow(X0+d,T)-tbpy.flow(X0-d,T))/(2*h)
    lam=np.linalg.eigvals(M); a=np.sort(np.abs(lam))[::-1]
    out.append(dict(p1=e['p1'],p2=e['p2'],T=T,Tstar=e['Tstar'],maxlam=float(a[0]),mods=[float(x) for x in a]))
    print(f"{e['p1']:.6f} {e['p2']:.6f} T={T:.4f} T*={e['Tstar']:.2f} |lam| top: {a[:4]}")
json.dump(out,open('/mnt/project-files/math-theory/three-body/stability_unlisted.json','w'),indent=1)
