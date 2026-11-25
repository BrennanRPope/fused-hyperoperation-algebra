"""Checks for 'the landscape from the ground up' (simple-laws). Exact (sympy) and numeric.
Run: python3 checks.py   (needs sympy, numpy, scipy)"""
import sympy as sp, numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
ok=0; bad=0
def check(name, cond):
    global ok,bad
    if bool(cond): ok+=1; print('PASS', name)
    else: bad+=1; print('FAIL', name)

x,y,a,b,c,d,k,q,u,t,P1,P2,D,lam,th,al,be=sp.symbols('x y a b c d k q u t P1 P2 D lam th al be')
law=lambda K,Q,X,Y:(X+Y+K*X*Y)/(1-Q*X*Y)

# 1. number system: (1+xu)(1+yu) with u^2 = k u - d, read as slope = landscape law
e=sp.expand((1+x*u)*(1+y*u)); e=sp.rem(sp.Poly(e,u),sp.Poly(u**2-k*u+q,u)).as_expr()
re_=e.coeff(u,0); im_=e.coeff(u,1)
check('number system slope = landscape law', sp.simplify(im_/re_-law(k,q,x,y))==0)
# norm multiplicative
N=lambda A,B:A**2+k*A*B+q*B**2
pr=sp.rem(sp.Poly(sp.expand((a+b*u)*(c+d*u)),u),sp.Poly(u**2-k*u+q,u)).as_expr()
check('norm multiplies', sp.expand(N(pr.coeff(u,0),pr.coeff(u,1))-N(a,b)*N(c,d))==0)
# 2. regions of each placed point
def region(K,Q):
    D_=K*K-4*Q
    return 'two walls' if D_>0 else ('border' if D_==0 else 'turning')
pts={'addition':(0,0),'The circle operation':(1,0),'union of failures / cooling':(-1,0),'tangent / Pythagoras':(0,1),
     'velocity addition (Einstein)':(0,-1),'lens pair ℓ=0.1 m':(-0.1,0),'rung 2':(2,sp.Rational(1,2)),
     'law of cosines C=60deg':(-1,1),'law of cosines C=90deg':(0,1),'flat triangle':(2,1)}
exp={'addition':'border','The circle operation':'two walls','union of failures / cooling':'two walls','tangent / Pythagoras':'turning',
     'velocity addition (Einstein)':'two walls','lens pair ℓ=0.1 m':'two walls','rung 2':'two walls',
     'law of cosines C=60deg':'turning','law of cosines C=90deg':'turning','flat triangle':'border'}
for n,(K,Q) in pts.items(): check('region '+n+' = '+exp[n], region(K,Q)==exp[n])
# 3. rotation seen in a frame tilted to angle theta is the law at (2 cos theta, 1)
def tilt_law(thv):
    e1=np.array([1,0.]); e2=np.array([np.cos(thv),np.sin(thv)]); Bi=np.linalg.inv(np.column_stack([e1,e2]))
    def T(phi):
        R=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]]); X,Y=Bi@R@e1; return Y/X
    worst=0
    for A in [0.1,0.25,-0.15]:
        for Bv in [0.05,-0.2,0.3]:
            pa=brentq(lambda p:T(p)-A,-0.6,0.6); pb=brentq(lambda p:T(p)-Bv,-0.6,0.6)
            worst=max(worst,abs(T(pa+pb)-law(2*np.cos(thv),1,A,Bv)))
    return worst
for deg in [90,60,120,75]:
    w=tilt_law(np.radians(deg)); check(f'rotation in a {deg} deg frame = law at (2cos,1), err {w:.1e}', w<1e-10)
# 4. geometry identities
check('gnomon (1+a)^2-1 = 2a+a^2', sp.expand((1+a)**2-1-(2*a+a**2))==0)
check('cube layers', sp.Poly(sp.expand((1+a)**3),a).all_coeffs()[::-1]==[1,3,3,1])
check('triangles compose (Brahmagupta)', sp.expand((1+x**2)*(1+y**2)-((1-x*y)**2+(x+y)**2))==0)
check('Euclid triples', sp.expand((q**2-x**2)**2+(2*x*q)**2-(q**2+x**2)**2)==0)
z=sp.expand((1+sp.I*x)**3)
check('triple angle: re 1-3x^2, im 3x-x^3', sp.re(z.subs(x,sp.Symbol('r',real=True))).equals(1-3*sp.Symbol('r',real=True)**2))
for tv in [0.2,0.7,-0.4]:
    check(f'tan 3theta from the cube, t={tv}', abs(np.tan(3*np.arctan(tv))-(3*tv-tv**3)/(1-3*tv**2))<1e-12)
for Cdeg in [30,60,90,120,150]:
    A,B=3.0,5.0; C=np.radians(Cdeg); P=np.array([A,0]); Q=B*np.array([np.cos(C),np.sin(C)])
    check(f'law of cosines as norm at k=-2cosC, C={Cdeg}', abs(np.sum((P-Q)**2)-(A*A+B*B-2*np.cos(C)*A*B))<1e-12)
# 5. Kepler III with real orbits (semi-major axis AU, period yr)
orb={'Mercury':(0.3871,0.2408),'Venus':(0.7233,0.6152),'Mars':(1.5237,1.8808),'Jupiter':(5.2029,11.862),'Saturn':(9.5371,29.457)}
for n,(A,T) in orb.items():
    al_=A-1; t_=T-1
    check(f'Kepler III {n}: t(+)t vs al(+)al(+)al within 0.5%', abs(((1+t_)**2)/((1+al_)**3)-1)<0.005)
# 6. lenses: ray trace two thin lenses
for p1,dd,p2 in [(5.,0.1,3.),(10.,0.05,-4.),(2.,0.5,7.)]:
    M=np.array([[1,0],[-p2,1]])@np.array([[1,dd],[0,1]])@np.array([[1,0],[-p1,1]])
    check(f'two lenses apart P={p1},{p2} ℓ={dd}', abs(-M[1,0]-(p1+p2-dd*p1*p2))<1e-12)
M=np.array([[1,0],[-9.,1]])@np.array([[1,0.2],[0,1]])@np.array([[1,0],[-5.,1]])
check('field lens: P1 ℓ = 1 leaves total power = P1', abs(-M[1,0]-5.)<1e-12)
# 7. Boyle: pressure gain a, volume factor 1/(1+a); hyperbola asymptotes at the walls
for av in [0.5,1.0,3.0,-0.5]:
    bv=-av/(1+av); check(f'Boyle partner of {av} is its circle inverse', abs(av+bv+av*bv)<1e-14)
check('Boyle hyperbola: b -> -1 as a -> infinity (the wall)', abs(-1e9/(1+1e9)+1)<1e-8)
# 8. ladder: F(2n)/F(2n-1) -> golden ratio, invariant n^2-nd-d^2=-1
nn,dd_=0,1; phi=(1+5**0.5)/2
for i in range(12): nn,dd_=2*nn+dd_,nn+dd_
check('ladder invariant after 12 sections', nn*nn-nn*dd_-dd_*dd_==-1)
check('ladder -> golden ratio', abs(nn/dd_-phi)<1e-8)
# 9. transmission line Z' = r - g Z^2 composes like velocity addition in units Z0, length sqrt(rg)
r_,g_=4.0,1.0; Z0=(r_/g_)**0.5
def run(z0,L):
    s=solve_ivp(lambda X,Z:[r_-g_*Z[0]**2],[0,L],[z0],rtol=1e-12,atol=1e-14); return s.y[0,-1]
z1=run(0,0.2); z2=run(0,0.35); z12=run(0,0.55)
check('line sections compose by velocity addition (units of Z0)', abs(z12/Z0-law(0,-1,z1/Z0,z2/Z0))<1e-9)
check('matched load (Z = Z0, the wall) is unchanged by any length', abs(run(Z0,1.3)-Z0)<1e-9)
# 10. cooling: gap fraction closes by union of failures
check('cooling: 1-(1-a)^2 = a + a - a*a', sp.expand(1-(1-a)**2-(a+a-a*a))==0)
# 11. Galilean limit: Einstein (0,-1/c^2) slides to (0,0) as c grows
cc=sp.Symbol('cc',positive=True)
check('Galilean is the c -> infinity limit', sp.limit((x+y)/(1+x*y/cc**2),cc,sp.oo)==x+y)
check('units: Einstein in m/s is the point (0,-1/c^2)', sp.simplify((x+y)/(1+x*y/cc**2)-law(0,-1/cc**2,x,y))==0)
# 12. relative error of a product is the circle operation exactly
check('relative errors multiply-add', sp.expand((1+a)*(1+b)-1-(a+b+a*b))==0)
# 13. nonempty subsets
for m_,n_ in [(2,3),(5,1),(4,4)]:
    A=2**m_-1;B=2**n_-1; check(f'nonempty subsets {m_},{n_}', A+B+A*B==2**(m_+n_)-1)
print(f'\n{ok} passed, {bad} failed')
