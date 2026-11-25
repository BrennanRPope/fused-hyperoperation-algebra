import sympy as sp, random
X=sp.symbols('A01 A02 A03 A12 A13 A23')
R=sp.Poly(sp.sympify(open('relation4.txt').read()),*X)
names=['w01','w02','w03','w12','w13','w23']
terms=[]
for mon,c in R.terms():
    fac=[]
    for n,e in zip(names,mon): fac+= [n]*e
    body=' * '.join(fac) if fac else '1'
    c=int(c)
    terms.append((c,body))
lines=[]
for k,(c,body) in enumerate(terms):
    sgn='-' if c<0 else '+'
    mag=abs(c)
    t=body if mag==1 else (f'{mag} * {body}' if body!='1' else f'{mag}')
    lines.append((sgn,t))
out='def rule4 (w01 w02 w03 w12 w13 w23 : Rat) : Rat :=\n'
first=True; cur='  '
for sgn,t in lines:
    piece=(('-' if sgn=='-' else '')+t) if first else f' {sgn} {t}'
    first=False
    if len(cur)+len(piece)>98: out+=cur.rstrip()+'\n'; cur='    '+piece.strip()
    else: cur+=piece
out+=cur+'\n'
open('rule4.lean.txt','w').write(out)
print(out[:600]); print(len(out.splitlines()),'lines')
# concrete points
def cr(a,b,i,j): return sp.Rational((a[i]-a[j])*(b[i]-b[j]),(a[i]-b[j])*(b[i]-a[j]))
a=[0,1,3,7]; b=[2,5,-1,4]
vals={(i,j):cr(a,b,i,j) for i in range(4) for j in range(i+1,4)}
print(vals)
sub=dict(zip(X,[vals[(0,1)],vals[(0,2)],vals[(0,3)],vals[(1,2)],vals[(1,3)],vals[(2,3)]]))
print('R at points',R.as_expr().subs(sub))
print('R at all 2', R.as_expr().subs({x:2 for x in X}))
# sixth weight's other root at these five
A23=X[5]; q=sp.Poly(R.as_expr().subs({k:v for k,v in sub.items() if k!=A23}),A23)
print('quadratic roots',sp.roots(q))
