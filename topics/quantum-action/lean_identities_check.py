"""Sympy check of every identity in lean/QuantumAction.lean (the Lean file itself is uncompiled here)."""
import sympy as sp
q1,p1,q2,p2,s1,s2,q3,p3,s3,a,b,c,d,e,f,u,v,w,t,x,y,k,r,D = sp.symbols("q1 p1 q2 p2 s1 s2 q3 p3 s3 a b c d e f u v w t x y k r D")
def mul(A,B): return (A[0]+B[0], A[1]+B[1], A[2]+B[2]+(A[0]*B[1]-A[1]*B[0])/2)
z = lambda E: sp.simplify(sp.expand(E)) == 0
res = {}
A,B,C = (q1,p1,s1),(q2,p2,s2),(q3,p3,s3)
res["mul_assoc"] = all(z(i-j) for i,j in zip(mul(mul(A,B),C), mul(A,mul(B,C))))
res["mul_inv"] = all(z(i) for i in mul(A,(-q1,-p1,-s1)))
res["swap_costs_area"] = z(mul(A,B)[2]-mul(B,A)[2]-(q1*p2-p1*q2))
res["parallel_no_cost"] = z(q1*(r*p1)-p1*(r*q1))
tri = mul(mul((q1,p1,0),(q2,p2,0)),(-(q1+q2),-(p1+p2),0))
res["triangle_loop"] = z(tri[0]) and z(tri[1]) and z(tri[2]-(q1*p2-p1*q2)/2)
quad = mul(mul(mul((a,b,0),(c-a,d-b,0)),(e-c,f-d,0)),(-e,-f,0))
res["quad_loop"] = z(quad[0]) and z(quad[1]) and z(quad[2]-((a*d-b*c)+(c*f-d*e))/2)
res["boost_keeps_form"] = z((t+D*u*x)**2 + D*(x-u*t)**2 - (1+D*u*u)*(t*t+D*x*x))
M = lambda u: sp.Matrix([[1, D*u],[-u, 1]])
prod = M(u)*M(v); wv = (u+v)/(1-D*u*v)
res["boost_compose"] = all(z(i) for i in (prod - (1-D*u*v)*M(wv)))
res["k_is_a_slide"] = z((x+y+k*x*y)*((1+k/2*x)*(1+k/2*y)-(D-k*k/4)*x*y) - (x*(1+k/2*y)+y*(1+k/2*x))*((1-D*x*y)+k/2*(x+y+k*x*y)))
V,E,F = sp.symbols("V E F")
res["regge2d_total"] = z((6*V-2*E-6*(V-E+F)).subs(E, sp.Rational(3,2)*F))
for kk,vv in res.items(): print(("PASS" if vv else "FAIL"), kk)
print(sum(res.values()), "of", len(res))
