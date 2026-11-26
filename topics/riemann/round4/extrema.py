"""Round 4: the local fused discriminant D_e = -8 Xi/Xi'' at every extremum of Xi(gamma)
on the critical line, computed straight from Xi (via Z and the gamma factor), no zeros used
except to bracket the extremum.  Also kappa_e = (log Xi'')'' at the extremum, which sets the
exact rewind rate  dD_e/dt = 8 + D_e*kappa_e  (gamma units: same formula, see README).
Xi(1/2+i g) = 1/2 (g^2+1/4) pi^(-1/4) |Gamma(1/4+ig/2)| Z(g) up to sign."""
import mpmath as mp, json, sys
from multiprocessing import Pool
mp.mp.dps = 20
G = [mp.mpf(s) for s in json.load(open('../elliptic-object/zeros.json'))['gamma']]
L = json.load(open('../round3/zeros_lehmer.json')); GL = [mp.mpf(s) for s in L['gamma']]; NL = L['n']

def lF(g, k):
    """k-th derivative of log|F|, F = (g^2+1/4)|Gamma(1/4+ig/2)|"""
    w = mp.mpc(0.25, g/2)
    r = 2*mp.re((-1)**(k-1)*mp.factorial(k-1)/mp.mpc(g,0.5)**k)
    r += mp.re((1j/2)**k*mp.psi(k-1, w))
    return r

def h(g):  # Xi'/F up to the factor: Z' + a Z
    return mp.siegelz(g,derivative=1) + lF(g,1)*mp.siegelz(g)

def extremum(pair):
    g1, g2 = pair
    ge = mp.findroot(h, (g1+mp.mpf('1e-9')*(g2-g1), g2-mp.mpf('1e-9')*(g2-g1)), solver='illinois', tol=mp.mpf(10)**-16)
    Zk = [mp.siegelz(ge,derivative=k) for k in range(5)]
    l = [None]+[lF(ge,k) for k in range(1,5)]
    f = [1, l[1], l[2]+l[1]**2, l[3]+3*l[1]*l[2]+l[1]**3,
         l[4]+4*l[1]*l[3]+3*l[2]**2+6*l[1]**2*l[2]+l[1]**4]
    X = [sum(mp.binomial(k,j)*f[j]*Zk[k-j] for j in range(k+1)) for k in range(5)]  # Xi^(k)/F
    De = -8*X[0]/X[2]
    kap = X[4]/X[2]-(X[3]/X[2])**2
    return dict(g1=float(g1), g2=float(g2), ge=float(ge), gap=float(g2-g1), De=float(De),
                kappa=float(kap), rate=float(8+De*kap), xi1_rel=float(X[1]/X[2]))

if __name__ == '__main__':
    nmain = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    pairs = [(G[i], G[i+1]) for i in range(nmain-1)]
    lp = [(GL[i], GL[i+1]) for i in range(len(GL)-1)]
    with Pool(4) as P:
        main = P.map(extremum, pairs, chunksize=8)
        leh = P.map(extremum, lp, chunksize=8)
    for i, r in enumerate(main): r['n'] = i+1
    for i, r in enumerate(leh): r['n'] = NL[i] if isinstance(NL, list) else NL+i
    json.dump(dict(main=main, lehmer=leh), open('extrema.json','w'))
    print('done', len(main), len(leh), min(r['De'] for r in main+leh))
