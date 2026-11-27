import mpmath as mp
mp.mp.dps = 30
def Z(s, y, K=60):
    """Epstein zeta sum' (m^2 + y n^2)^-s via Chowla-Selberg."""
    s = mp.mpc(s); r = mp.sqrt(y)
    t1 = 2*mp.zeta(2*s)
    t2 = 2*mp.sqrt(mp.pi)*mp.gamma(s-0.5)/mp.gamma(s)*mp.zeta(2*s-1)*y**(0.5-s)
    acc = 0
    for k in range(1, K+1):
        sig = sum(mp.mpf(d)**(1-2*s) for d in range(1, k+1) if k % d == 0)
        acc += k**(s-0.5)*sig*mp.besselk(s-0.5, 2*mp.pi*k*r)
    t3 = 8*mp.pi**s/mp.gamma(s)*y**(0.25-0.5*s)*acc
    return t1+t2+t3
def brute(s, y, N=400):
    tot = mp.mpf(0)
    for m in range(-N, N+1):
        for n in range(-N, N+1):
            if m or n: tot += (m*m+y*n*n)**(-s)
    return tot
if __name__ == "__main__":
    for y in (1, 2, 5.3):
        print(y, mp.nstr(Z(3, y), 12), mp.nstr(brute(3, y, 150), 12))
    print("y=1 vs 4 zeta beta:", mp.nstr(Z(0.5+14j, 1), 12), mp.nstr(4*mp.zeta(0.5+14j)*mp.dirichlet(0.5+14j, [0,1,0,-1]), 12))
