from dialgas import *
from fractions import Fraction as F
import mpmath as mp
# 1. exact wall check: Z(-1/k) for cycles and tori, exact rationals
print("Wall check Z(z=-1/k, k) exact:")
for k in (F(3),F(2),F(1,2),F(-1),F(-2),F(5,3)):
    row=[]
    for n in range(3,10):
        row.append(str(evalZ(cycle_counts(n),-1/k,k)))
    print(" cycle k=",k, row)
Ct={}
for L,M in ((3,3),(3,4),(3,5),(4,4),(5,5),(4,5)):
    Ct[(L,M)]=torus_counts(L,M)
for k in (F(3),F(1,2),F(-1),F(-2)):
    print(" torus k=",k,[ (LM, str(evalZ(C,-1/k,k))) for LM,C in Ct.items()])
