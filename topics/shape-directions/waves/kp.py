import sympy as sp, random, itertools, json
random.seed(3)
x,y,t=sp.symbols('x y t')
def tau(g,p,q,Aov=None):
    eta=[(p[i]-q[i])*x+(p[i]**2-q[i]**2)*y+(p[i]**3-q[i]**3)*t for i in range(g)]
    def A(i,j):
        if Aov and (i,j) in Aov: return Aov[(i,j)]
        return (p[i]-p[j])*(q[i]-q[j])/((p[i]-q[j])*(q[i]-p[j]))
    T=0
    for S in itertools.product([0,1],repeat=g):
        idx=[i for i in range(g) if S[i]]
        w=sp.Integer(1)
        for i,j in itertools.combinations(idx,2): w*=A(i,j)
        T+=w*sp.exp(sum(eta[i] for i in idx))
    return T
def kp(T):
    d=lambda f,*v: sp.diff(f,*v)
    H=2*(T*d(T,x,4)-4*d(T,x)*d(T,x,3)+3*d(T,x,2)**2) - 4*2*(T*d(T,x,t)-d(T,x)*d(T,t)) + 3*2*(T*d(T,y,2)-d(T,y)**2)
    return H
res={}
def rnd(): return sp.Rational(random.randint(-30,30),random.randint(1,7))
for g in (1,2,3,4):
    while True:
        p=[rnd() for _ in range(g)]; q=[rnd() for _ in range(g)]
        if len(set(p+q))==2*g: break
    H=sp.expand(kp(tau(g,p,q)))
    ok= H==0
    # perturb one pair weight
    bad=None
    if g>=2:
        A01=(p[0]-p[1])*(q[0]-q[1])/((p[0]-q[1])*(q[0]-p[1]))
        Hb=sp.expand(kp(tau(g,p,q,{(0,1):A01*sp.Rational(3,2)})))
        bad= Hb!=0
    res[g]=dict(solves=ok,perturbed_fails=bad); print(g,res[g])
json.dump({str(k):v for k,v in res.items()},open('kp.json','w'))
