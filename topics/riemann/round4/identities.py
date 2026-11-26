"""Round 4: exact identities, checked symbolically with sympy.
Convention (de Bruijn-Newman): dH/dt = -H_xx.  Rewind = t decreasing.
Jet variables h_k stand for the k-th x-derivative of H; the flow sends h_k -> -h_{k+2}."""
import sympy as sp, json
x,t,a,c = sp.symbols('x t a c', real=True)
h = sp.symbols('h0:9')
Dx = lambda e: sum(sp.diff(e,h[k])*h[k+1] for k in range(8))
Dt = lambda e: sum(sp.diff(e,h[k])*(-h[k+2]) for k in range(7))
shift = lambda e: e.subs({h[k]:h[k+1] for k in reversed(range(8))}, simultaneous=True)
Lag = lambda: h[1]**2 - h[0]*h[2]
out = {}
# 1. ladder PDE: L_t = -L_xx + 2 L(H')
L = Lag()
out['ladder_pde'] = sp.expand(Dt(L) - (-Dx(Dx(L)) + 2*shift(L))) == 0
# 2. local fused discriminant at an extremum (h1=0): D_e=-8H/H'', dD_e/dt = 8 + D_e*(log H'')''
De = -8*h[0]/h[2]
xdot = h[3]/h[2]                         # from H'(x_e(t),t)=0
kappa = h[4]/h[2] - (h[3]/h[2])**2       # (log H'')''
out['De_rate'] = sp.simplify((Dt(De) + Dx(De)*xdot - (8 + De*kappa)).subs(h[1],0)) == 0
# 3. quadratic: D_e equals the fused discriminant s^2-4p (=4c), min L = Delta/2
q = (x-a)**2 - c
out['quad_De_equals_Delta'] = sp.simplify((-8*q/sp.diff(q,x,2)).subs(x,a) - 4*c) == 0
Lq = sp.expand(sp.diff(q,x)**2 - q*sp.diff(q,x,2)); xm = sp.solve(sp.diff(Lq,x),x)[0]
out['quad_minL_half_Delta'] = sp.simplify(Lq.subs(x,xm) - 2*c) == 0
# 4. sin is a steady state: dD_e/dt = 0 at its extrema
s = sp.sin(x); d = lambda k: sp.diff(s,x,k).subs(x,sp.pi/2)
out['sin_steady'] = sp.simplify(8 + (-8*d(0)/d(2))*(d(4)/d(2)-(d(3)/d(2))**2)) == 0
# 5. heat flow of a quadratic: q_t = -q_xx gives D(t) = D0 + 8t
qt = (x-a)**2 - c - 2*t
out['quad_flow'] = bool(sp.simplify(sp.diff(qt,t)+sp.diff(qt,x,2))==0 and
                        sp.simplify((-8*qt/sp.diff(qt,x,2)).subs(x,a)-(4*c+8*t))==0)
# 6. rewind form: L_s = L_xx - 2 L(H') with s=-t ; lone quadratic: L_xx=4, 2L(H')=8
out['quad_rewind_rates'] = (sp.expand(sp.diff(Lq,x,2))==4 and
     sp.expand(sp.diff(q,x,2)**2 - sp.diff(q,x)*sp.diff(q,x,3))==4)
print(out); json.dump({k:bool(v) for k,v in out.items()}, open('identities.json','w'), indent=1)
