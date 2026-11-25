# x^3+y^3+z^3+w^3=0 over F_p: affine cone count, then projective = (cone-1)/(p-1)
from sympy import primerange
for p in primerange(5,80):
    c={}
    for x in range(p): c[pow(x,3,p)]=c.get(pow(x,3,p),0)+1
    pair={}
    for a,u in c.items():
        for b,v in c.items(): pair[(a+b)%p]=pair.get((a+b)%p,0)+u*v
    cone=sum(u*pair.get((-k)%p,0) for k,u in pair.items())
    N=(cone-1)//(p-1); t=(N-p*p-1)//p
    print(p,p%3,N,'t=',t)
