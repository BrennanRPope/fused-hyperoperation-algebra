import time, random, mpmath as mp
from common import *
import flint
random.seed(1)
# speed/correctness test on a 4-dim de Weger lattice with C=10^260
mp.mp.dps=330
K=strip(T2(100,80,'N'))
al=[mp.log(K),mp.log(2),mp.log(5),mp.log(19)]
C=10**260
a=[int(mp.nint(C*x)) for x in al]
B=[[1,0,0,a[0]],[0,1,0,a[1]],[0,0,1,a[2]],[0,0,0,a[3]]]
t=time.time(); R=lll_int(B); t1=time.time()-t
t=time.time(); ok,nrm=check_lll_reduced(R); t2=time.time()-t
print('own LLL %.4fs, check %.4fs, ok=%s'%(t1,t2,ok))
M=flint.fmpz_mat(B); t=time.time(); Rf=M.lll(delta=0.99); t3=time.time()-t
print('flint %.4fs'%t3)
import math
print('b1 own', math.log10(math.sqrt(sum(x*x for x in R[0]))), ' flint', math.log10(math.sqrt(sum(int(Rf[0,j])**2 for j in range(4)))))
# det check
print(abs(flint.fmpz_mat(R).det())==abs(a[3]))
