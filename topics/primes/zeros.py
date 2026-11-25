import mpmath, numpy as np, time
mpmath.mp.dps = 15
t=time.time(); g=[]
for n in range(1,1001):
    g.append(float(mpmath.zetazero(n).imag))
np.save('zeta_zeros_1000.npy', np.array(g)); print(len(g), g[:3], g[-1], time.time()-t)
