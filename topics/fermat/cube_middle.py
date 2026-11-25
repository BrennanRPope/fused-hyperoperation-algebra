# Mult-add p-th power vs ordinary p-th power: [p](s) = (1+s)^p - 1, middle of Pascal cube = [p](s) - s^p.
# Height 1: middle always divisible by p. Wieferich/Mirimanoff: a case-I Fermat counterexample forces it mod p^2.
from sympy import primerange
def agree2(p, s):
    return ((1+s)**p - 1 - s**p) % (p*p) == 0
for p in [1093, 3511]:
    print(p, [s for s in range(1, 6) if agree2(p, s)])
hits = [(p, s) for p in primerange(3, 20000) for s in range(1, 4) if agree2(p, s)]
print("primes<20000 where middle of cube divisible by p^2, s in 1..3:", hits)
# all p: middle divisible by p
print(all(((1+s)**p-1-s**p) % p == 0 for p in primerange(3,200) for s in range(1,30)))
