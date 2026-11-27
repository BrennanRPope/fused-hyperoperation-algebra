import mpmath as mp
from epstein import Z
mp.mp.dps = 50
f = lambda y: Z(mp.mpf(0.5)+mp.mpf(10)**-25, y, K=30).real
for y in (8, 10, 15, 20, 30, 50): print(y, mp.nstr(f(y), 8))
