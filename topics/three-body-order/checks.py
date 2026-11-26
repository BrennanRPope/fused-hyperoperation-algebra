"""DRAFT checks for three-body-order. Run: python3 checks.py > checks_output.txt  (needs numpy, sympy; reads the JSON results made by
words.py, physics.py/polish_swap.py, wind_ode.py, wind_cat.py, stability.py, pairs_ctrl.py, figures.py, dial_nested.json)."""
import json, math, collections, numpy as np, sympy as sp
from words import runs, trace_m1, walk_area, readings
ok=0; bad=0
def check(name,cond,detail=''):
    global ok,bad
    if cond: ok+=1; print('PASS',name,detail)
    else: bad+=1; print('FAIL',name,detail)
D=json.load(open('../three-body-word-formula/data/li-liao-695.json')); R=json.load(open('words.json'))
# A. the letters
bg=collections.Counter(w[i]+w[(i+1)%len(w)] for w in (o['word'] for o in D) for i in range(len(w)))
check('A1 same-case pairs (ab, ba, AB, BA) are the common ones; case changes are the reversals',min(bg['ab'],bg['AB'])>max(bg['aB'],bg['Ab']),str(dict(bg)))
rc=collections.Counter(k for o in D for k in runs(o['word']))
check('A2 runs as same-case blocks reproduce the earlier run table 24/8052/10196/7830/56',[rc[k] for k in range(1,6)]==[24,8052,10196,7830,56],str(dict(rc)))
# B. order dial on letters (Burau), exact algebra
t=sp.symbols('t'); s1=sp.Matrix([[-t,1],[0,1]]); s2=sp.Matrix([[1,0],[t,-t]])
check('B1 t = 1: a pure twist sigma1^2 becomes the identity (the dial forgets pure braids)',(s1*s1).subs(t,1)==sp.eye(2))
P=(s1*s1).subs(t,-1); Q=(s2*s2).subs(t,-1)
check('B2 t = -1: letters become the Mobius matrices [[1,2],[0,1]] and [[1,0],[-2,1]]',P==sp.Matrix([[1,2],[0,1]]) and Q==sp.Matrix([[1,0],[-2,1]]))
check('B3 ab is parabolic (trace -2): a run circles the third collision point',(P*Q).trace()==-2)
phi=(1+sp.sqrt(5))/2
check('B4 figure-eight BabA: trace 18 = phi^6 + phi^-6',trace_m1('BabA')==18 and sp.simplify(phi**6+phi**-6-18)==0)
check('B5 figure-eight Mobius length = 12 ln(phi) (braid entropy 6 ln phi per period)',abs(readings('BabA')['ell']-12*math.log((1+5**.5)/2))<1e-12)
check('B6 full twist is a scalar under the dial: (sigma1 sigma2)^3 = t^3 I',sp.simplify((s1*s2)**3-t**3*sp.eye(2))==sp.zeros(2))
# definite (unitary) range: for t = exp(i x pi), x < 2/3, every word's matrix has spectral radius 1
from words import rho
xs=[0.2,0.5,0.6]; mx=max(rho(r['word'],np.exp(1j*np.pi*x)) for r in R[:200] for x in xs)
check('B7 below x = 2/3 the dial shows no growth at all (spectral radius 1)',abs(mx-1)<1e-6,'max %.9f'%mx)
# C. hyperbolic length adds run by run
K=np.zeros((len(R),5)); ell=np.array([r['ell'] for r in R]); T=np.array([r['Tstar'] for r in R]); L=np.array([r['L'] for r in R])
for i,r in enumerate(R):
    for k in r['runs']: K[i,min(k,5)-1]+=1
c=np.linalg.lstsq(K,ell,rcond=None)[0]; r2=1-np.var(ell-K@c)/np.var(ell)
check('C1 Mobius length is a sum of per-run prices (R2 > 0.999)',r2>0.999,'R2 %.5f prices %s'%(r2,np.round(c,3)))
pl=c/np.arange(1,6); ct=np.linalg.lstsq(K,T,rcond=None)[0]; tl=ct/np.arange(1,6)
check('C2 per letter: Mobius price falls from run 2 to 5, time price rises',pl[1]>pl[2]>pl[3]>pl[4] and tl[1]<tl[2]<tl[3],'mobius %s time %s'%(np.round(pl,3),np.round(tl,3)))
check('C3 Mobius length per letter vs period per letter: strongly anti-correlated',np.corrcoef(ell/L,T/L)[0,1]<-0.9,'r = %.3f'%np.corrcoef(ell/L,T/L)[0,1])
# D. the word's wedge (Heisenberg area)
inv=lambda w: ''.join(ch.swapcase() for ch in reversed(w))
closed=[o['word'] for o in D if walk_area(o['word'])[:2]==(0,0)]
check('D1 time reversal (reverse and invert the word) flips the walk area, on all closed words',all(walk_area(inv(w))[2]==-walk_area(w)[2] for w in closed),'%d words'%len(closed))
def loo(X,y):
    H=X@np.linalg.pinv(X); e=(y-H@y)/(1-np.diag(H))/y; return np.median(abs(e))*100
cm=np.array([r['netA']==0 and r['netB']==0 for r in R]); ar=np.abs([r['area'] for r in R])
a,b=loo(K[cm],T[cm]),loo(np.column_stack([K[cm],ar[cm]]),T[cm])
check('D2 |area| does not sharpen the run law (held-out change under 10%)',b>0.9*a,'%.4f%% -> %.4f%%'%(a,b))
# E. swap rule
S=json.load(open('swap_polished.json')); worst=max(e['worst'] for e in S.values()); npairs=sum(len(e['rows']) for e in S.values())
check('E1 swap rule on the shape sphere: line-up turn = -1/2 enclosed area (mod pi), after polish',worst<2e-4,'%d arcs, worst %.1e rad'%(npairs,worst))
# F. windings
W=json.load(open('wind_cat.json'))
check('F1 windings are whole numbers',max(abs(v-round(v)) for e in W for v in e['w'].values())<1e-4)
check('F2 windings = (0, n, -n) with n = net count of a, on the 90 shortest orbits',all(tuple(round(e['w'][k]) for k in ['12','13','23'])==(0,e['netA'],-e['netA']) and e['netA']==e['netB'] for e in W))
TW=json.load(open('twin_windings.json'))
check('F3 all 7 same-word twins have identical windings',all(all(round(v['cat'][k])==round(v['twin'][k]) for k in v['cat']) for v in TW.values()))
# G. period and order
C=json.load(open('pairs_ctrl.json'))
check('G1 which runs sit next to which sharpens the period (held out)',C['real'][0]<0.5*C['counts'][0],'%.4f%% -> %.4f%%'%(C['counts'][0],C['real'][0]))
check('G2 real run order beats all 20 shuffles',C['real'][0]<min(C['shuffled']),'real %.4f%%, shuffles %.4f..%.4f%%'%(C['real'][0],min(C['shuffled']),max(C['shuffled'])))
dn=json.load(open('dial_nested.json'))
check('G3 dial setting chosen on training data alone stays near x = 0.7',all(0.66<=x<=0.72 for x in dn['chosen']),str(dn['chosen']))
# H. stability
St=json.load(open('stability.json')); ly=np.array([s['logmaxlam'] for s in St]); el=np.array([readings(s['word'])['ell']/2 for s in St])
Ks=np.zeros((len(St),5))
for i,s in enumerate(St):
    for k in runs(s['word']): Ks[i,min(k,5)-1]+=1
def r2h(X,y):
    X=np.column_stack([np.ones(len(y))]+X); H=X@np.linalg.pinv(X); e=(y-H@y)/(1-np.diag(H)); return 1-np.mean(e**2)/np.var(y)
a,b=r2h([Ks],ly),r2h([Ks,el],ly)
check('H1 braid entropy adds under 0.05 held-out R2 to runs for instability',b-a<0.05,'%.3f -> %.3f'%(a,b))
check('H2 the figure-eight stays nearly stable (largest multiplier < 1.001)',[s for s in St if s['cls']=='I.A1'][0]['logmaxlam']<1e-3)
print(f'\n{ok} passed, {bad} failed')
