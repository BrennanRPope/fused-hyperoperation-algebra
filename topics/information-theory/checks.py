"""Checks for the information-theory page (Fused Hyperoperation Algebra).

Part A: every ring identity stated in lean/InfoFused.lean, checked symbolically (sympy).
Part B: the list lemmas and the laws on random distributions, channels and codes (numpy).
Run: python3 checks.py > checks_output.txt
"""
import itertools, math
import numpy as np
import sympy as sp

rng = np.random.default_rng(20261125)
ok_count = 0
def ok(name, cond):
    global ok_count
    assert cond, name
    ok_count += 1
    print(f"  ok  {name}")

a, b, c, k, p, q, r, P, R, N, M, o1, o2, e1, e2, cc = sp.symbols('a b c k p q r P R N M o1 o2 e1 e2 cc')
w = lambda k, a, b: a + b + k*a*b
madd = lambda a, b: a + b + a*b
z = lambda e: sp.expand(e) == 0

print("A. Symbolic identities (the Lean file's ring facts)")
ok("K1 lift", z(1 + k*w(k,a,b) - (1+k*a)*(1+k*b)))
ok("scale to plain mult-add", z(k*w(k,a,b) - madd(k*a,k*b)))
ok("K3 flip", z(w(k,-a,-b) + w(-k,a,b)))
ok("assoc", z(w(k,w(k,a,b),c) - w(k,a,w(k,b,c))))
ok("union rule", z(w(-1,p,q) - (1-(1-p)*(1-q))))
ok("union = flipped mult-add", z(w(-1,p,q) + madd(-p,-q)))
ok("inclusion-exclusion 3", z(w(-1,w(-1,p,q),r) - (p+q+r-p*q-p*r-q*r+p*q*r)))
ok("gini-simpson", z(1-P*R - w(-1,1-P,1-R)))
ok("support count", z(N*M-1 - madd(N-1,M-1)))
ok("tsallis via lift", z((1+(P-1))*(1+(R-1)) - P*R))
ok("kraft slack", z(1-P*R - w(-1,1-P,1-R)))
ok("campbell", z(P*R-1 - madd(P-1,R-1)))
bscM = lambda p: sp.Matrix([[1-p,p],[p,1-p]])
prod = bscM(p)*bscM(q)
ok("bsc cascade = matrix product off-diagonal", z(prod[0,1] - w(-2,p,q)))
ok("bsc cascade diagonal", z(prod[0,0] - (1-w(-2,p,q))))
ok("polarization", z(1-2*w(-2,p,q) - (1-2*p)*(1-2*q)))
ok("piling-up", z(sp.Rational(1,2)-w(-2,p,q) - 2*(sp.Rational(1,2)-p)*(sp.Rational(1,2)-q)))
ok("q-ary cascade", z(1-((1-e1)*(1-e2)+cc*e1*e2) - w(-(1+cc),e1,e2)))
for m in [2,3,4,5,8]:
    T = lambda e: sp.Matrix(m, m, lambda i,j: 1-e if i==j else e/(m-1))
    TT = T(e1)*T(e2)
    ok(f"q-ary formula matches {m}x{m} matrix product, k = -{m}/{m-1}",
       z((1-TT[0,0]) - w(-sp.Rational(m,m-1),e1,e2)))
ok("parity of 3 = odd layers", z(w(-2,w(-2,w(-2,0,p),p),p) - (3*p*(1-p)**2 + p**3)))
ok("check node", z(-2*w(-2,p,q) - madd(-2*p,-2*q)))
ok("variable node", z(o1*o2-1 - madd(o1-1,o2-1)))
ok("odds to velocity numerator", z((o1-1)*(o2+1)+(o2-1)*(o1+1) - 2*(o1*o2-1)))
ok("odds to velocity denominator", z((o1+1)*(o2+1)+(o1-1)*(o2-1) - 2*(o1*o2+1)))
ok("velocity corners +", z((1+a*b)+(a+b) - (1+a)*(1+b)))
ok("velocity corners -", z((1+a*b)-(a+b) - (1-a)*(1-b)))
x1,x2,y1,y2,u1,u2 = sp.symbols('x1 x2 y1 y2 u1 u2')
vadd = lambda X,Y: (X[0]*Y[1]+Y[0]*X[1], X[1]*Y[1]+X[0]*Y[0])
L1 = vadd(vadd((x1,x2),(y1,y2)),(u1,u2)); R1 = vadd((x1,x2),vadd((y1,y2),(u1,u2)))
ok("velocity addition associative (pairs)", z(L1[0]-R1[0]) and z(L1[1]-R1[1]))
n = sp.symbols('n', integer=True)
for nn in range(0,8):
    acc = 0
    for _ in range(nn): acc = w(-2, acc, p)
    ok(f"n flips, n={nn}", z(1-2*acc - (1-2*p)**nn))

print("\nB. Laws on numbers")
def rand_dist(n):
    x = rng.random(n); return x/x.sum()
def H(x): return -np.sum(x*np.log2(x))
def tsallis(x,qq): return (np.sum(x**qq)-1)/(1-qq)
def renyi(x,al): return np.log2(np.sum(x**al))/(1-al)
for trial in range(5):
    A, B = rand_dist(rng.integers(2,7)), rand_dist(rng.integers(2,7))
    J = np.outer(A,B).ravel()
    ok(f"[{trial}] joint sums to 1", abs(J.sum()-1)<1e-12)
    ok(f"[{trial}] Shannon additive", abs(H(J)-H(A)-H(B))<1e-12)
    for qq in [0.0, 0.5, 2.0, 3.0]:
        kk = 1-qq
        ok(f"[{trial}] Tsallis q={qq}: S_AB = S_A (+)_k S_B, k={kk}",
           abs(tsallis(J,qq) - (tsallis(A,qq)+tsallis(B,qq)+kk*tsallis(A,qq)*tsallis(B,qq)))<1e-10)
        ok(f"[{trial}] Renyi a={qq} additive", abs(renyi(J,qq)-renyi(A,qq)-renyi(B,qq))<1e-10)
        ok(f"[{trial}] Renyi = log2(1+k S_T)/k, q={qq}",
           abs(renyi(A,qq) - np.log2(1+kk*tsallis(A,qq))/kk)<1e-10)
    G = lambda x: 1-np.sum(x**2)
    ok(f"[{trial}] Gini-Simpson by union rule", abs(G(J) - (G(A)+G(B)-G(A)*G(B)))<1e-12)
    # 'pair differs' by simulation-free enumeration
    differs = sum(J[i]*J[j] for i in range(len(J)) for j in range(len(J)) if i!=j)
    ok(f"[{trial}] Gini = chance two draws differ", abs(differs-G(J))<1e-12)

for Nn in [2,3,5]:
    for mm in [1,2,3]:
        qq = 1-mm
        U = np.full(Nn, 1/Nn)
        ok(f"uniform N={Nn}, k={mm}: k*S_T = (N-1) fused power m = N^m - 1",
           abs(mm*tsallis(U,qq) - (Nn**mm-1))<1e-9)

print("\n  channels")
def bsc(p): return np.array([[1-p,p],[p,1-p]])
def h2(p): return 0 if p in (0,1) else -p*np.log2(p)-(1-p)*np.log2(1-p)
for p1,p2 in [(0.1,0.2),(0.05,0.3),(0.25,0.25),(0.01,0.4)]:
    T = bsc(p1)@bsc(p2); pc = T[0,1]
    ok(f"BSC {p1},{p2}: cascade crossover {pc:.4f} = p (+)_-2 q", abs(pc-(p1+p2-2*p1*p2))<1e-12)
    C1,C2,C = 1-h2(p1),1-h2(p2),1-h2(pc)
    print(f"      capacities {C1:.4f}, {C2:.4f} -> cascade {C:.4f}")
    ok("  DPI: cascade capacity <= min", C <= min(C1,C2)+1e-12)
    ok("  SDPI: (1-2p)^2 multiplies", abs((1-2*pc)**2-(1-2*p1)**2*(1-2*p2)**2)<1e-12)
for m in [3,4,6]:
    e1v,e2v = 0.1,0.15
    T = lambda e: np.where(np.eye(m)==1, 1-e, e/(m-1))
    TT = T(e1v)@T(e2v); cv = 1/(m-1)
    ok(f"{m}-ary symmetric: error {1-TT[0,0]:.5f} = weight -(1+1/{m-1})",
       abs((1-TT[0,0]) - (e1v+e2v-(1+cv)*e1v*e2v))<1e-12)
def bec(e):  # inputs 0,1 ; outputs 0,?,1
    return np.array([[1-e,e,0],[0,e,1-e]])
e1v,e2v = 0.2,0.35
# cascade of BEC: erasure stays erasure; second channel sees ? as input that stays ?
M2 = np.array([[1-e2v,e2v,0],[0,1,0],[0,e2v,1-e2v]])
TT = bec(e1v)@M2
ok(f"BEC cascade erasure {TT[0,1]:.4f} = union rule", abs(TT[0,1]-(e1v+e2v-e1v*e2v))<1e-12)
ok("BEC capacity multiplies", abs((1-TT[0,1])-(1-e1v)*(1-e2v))<1e-12)

print("\n  compression")
def huffman_lengths(pr):
    import heapq
    heap = [(x,[i]) for i,x in enumerate(pr)]; L = [0]*len(pr); heapq.heapify(heap)
    while len(heap)>1:
        x1,i1 = heapq.heappop(heap); x2,i2 = heapq.heappop(heap)
        for i in i1+i2: L[i]+=1
        heapq.heappush(heap,(x1+x2,i1+i2))
    return np.array(L)
for trial in range(3):
    A, B = rand_dist(5), rand_dist(4)
    la, lb = huffman_lengths(A), huffman_lengths(B)
    KA, KB = np.sum(2.0**-la), np.sum(2.0**-lb)
    lj = (la[:,None]+lb[None,:]).ravel()
    KJ = np.sum(2.0**-lj)
    ok(f"[{trial}] Kraft sums multiply ({KA:.3f}*{KB:.3f})", abs(KJ-KA*KB)<1e-12)
    ok(f"[{trial}] Kraft slack composes by union rule", abs((1-KJ)-((1-KA)+(1-KB)-(1-KA)*(1-KB)))<1e-12)
    J = np.outer(A,B).ravel()
    ok(f"[{trial}] expected lengths add", abs(np.sum(J*lj)-np.sum(A*la)-np.sum(B*lb))<1e-12)
    t = 0.7
    MA, MB, MJ = np.sum(A*2**(t*la)), np.sum(B*2**(t*lb)), np.sum(J*2**(t*lj))
    ok(f"[{trial}] Campbell M multiplies, M-1 mult-adds", abs((MJ-1)-((MA-1)+(MB-1)+(MA-1)*(MB-1)))<1e-12)
    al = 1/(1+t)
    ok(f"[{trial}] Campbell cost >= Renyi(1/(1+t)) bound (Campbell 1965)",
       np.log2(MA)/t >= renyi(A,al)-1e-12)

print("\n  codes and decoding")
for rr in range(2,7):
    nn = 2**rr-1; mm = nn-rr
    ok(f"Hamming r={rr}: 2^{mm}*(1+{nn}) = 2^{nn}", 2**mm*(1+nn)==2**nn)
for pv in [0.1, 0.3]:
    for nn in [3,5]:
        podd = sum(math.comb(nn,i)*pv**i*(1-pv)**(nn-i) for i in range(1,nn+1,2))
        ok(f"parity of {nn} bits, p={pv}: odd cube layers = (1-(1-2p)^n)/2", abs(podd-(1-(1-2*pv)**nn)/2)<1e-12)
for pv,qv in [(0.1,0.3),(0.2,0.45)]:
    # check node by enumeration: chance parity of two independent bits is 1
    par = pv*(1-qv)+qv*(1-pv)
    ok(f"check node p={pv},q={qv}: rho multiplies", abs((1-2*par)-(1-2*pv)*(1-2*qv))<1e-12)
    # variable node: two independent looks at one bit (uniform prior), each says 'chance of 1'
    post1 = pv*qv/(pv*qv+(1-pv)*(1-qv))
    rho = lambda x: 1-2*x
    vel = (rho(pv)+rho(qv))/(1+rho(pv)*rho(qv))
    ok(f"variable node p={pv},q={qv}: rho adds like velocities", abs(rho(post1)-vel)<1e-12)
    L = lambda x: math.log((1-x)/x)
    ok(f"variable node: LLRs add", abs(L(post1)-(L(pv)+L(qv)))<1e-12)
    ok(f"rho = tanh(LLR/2)", abs(rho(pv)-math.tanh(L(pv)/2))<1e-12)
    phi = lambda x: -math.log(math.tanh(x/2))
    ok("Gallager phi: check node in LLRs, phi(phi(x)) = x", abs(phi(phi(1.3))-1.3)<1e-9)

print(f"\nAll {ok_count} checks passed.")
