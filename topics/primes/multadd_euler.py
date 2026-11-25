import mpmath
from sympy import primerange
def madd(a,b): return a+b+a*b
for s in [2,3]:
    acc=0; acc2=0
    for p in primerange(2,200000):
        acc=madd(acc,1/(p**s-1)); acc2=madd(acc2,-p**-s)
    print(s, acc, float(mpmath.zeta(s))-1, acc2, 1/float(mpmath.zeta(s))-1)
# twin constant as mult-add of local corrections
acc=0
for p in primerange(3,10**6): acc=madd(acc,(p*(p-2))/((p-1)**2)-1)
print('C2', 1+acc)
