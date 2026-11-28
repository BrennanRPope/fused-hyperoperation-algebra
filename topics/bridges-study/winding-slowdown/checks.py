"""DRAFT. Phase winds at E/hbar along each clock's line; does gravity's slowdown come out of the law as slower winding?
Units G = c = hbar = 1. Checks W1-W7."""
import sympy as sp, json
M,r,a,v,al = sp.symbols('M r a v alpha', positive=True)
res = {}
def chk(n, ok, msg): res[n]=bool(ok); print(n, "PASS" if ok else "FAIL", " ", msg)
def Z(e):
    """zero by exact simplify, else by high-precision evaluation at 9 points (0 < a < M)"""
    e = sp.sympify(e)
    if sp.simplify(e) == 0: return True
    syms = sorted(e.free_symbols, key=str)
    import itertools, random
    random.seed(1)
    for _ in range(9):
        sub = {}
        for s_ in syms:
            sub[s_] = sp.Rational(random.randint(11, 97), 100)
        if M in sub and a in sub: sub[a] = sub[M]*sp.Rational(random.randint(5, 95), 100)
        if M in sub and r in sub: sub[r] = sub[M]*sp.Rational(random.randint(250, 900), 100)
        if abs(sp.N(e.subs(sub), 40)) > 1e-25: return False
    return True
law = lambda x,y,k=1: x+y+k*x*y                       # the law as it is (c = 0)
w2 = lambda vv: 1 + law(sp.I*vv, sp.I*vv)              # 1 + x∘x with x = i v  ( = (1+iv)^2 )

# W1 Schwarzschild: river speed v^2 = 2M/r; winding factor N^2 = 1 - v^2 = Re(1 + x∘x)
vS = sp.sqrt(2*M/r)
chk("W1", Z(sp.re(sp.expand(w2(vS))) - (1-2*M/r)) and Z((1-2*M/r).subs(r,2*M)),
    "far-seen winding rate = (E/ħ)·N, N² = Re(1 + x∘x) with x = i·(river speed); winding stops where the river reaches light speed (r = 2M)")
# W2 Kerr: shifted source 1/(r - i a). Re part is the potential, Im part the shell turn rate
zinv = 1/(r - sp.I*a)
vK2 = 2*M*sp.re(sp.expand_complex(zinv))
Delta = r**2-2*M*r+a**2
rp = M+sp.sqrt(M**2-a**2)
chk("W2", Z(1 - vK2 - Delta/(r**2+a**2)) and Z(vK2.subs(r,rp)-1) and Z(sp.im(sp.expand_complex(zinv)) - a/(r**2+a**2)),
    "Kerr: 1 − 2M·Re[1/(r−ia)] = Δ/(r²+a²) (shell clock factor), = 0 at r₊; Im[1/(r−ia)] = a/(r²+a²) = shell turn rate")
# W3 surface gravity = slope of the clock factor = pull of the shifted source -> 1/(r - ia)^2 brings the square
kap = sp.simplify(sp.diff(Delta/(r**2+a**2), r).subs(r,rp)/2)
pull = sp.simplify(M*sp.re(sp.expand_complex(1/(rp - sp.I*a)**2)))
vr = a/rp
chk("W3", Z(kap - pull) and Z(kap - (1-vr**2)/(4*M)),
    "κ = ½·slope of the clock factor = M·Re[1/(r₊−ia)²] = (1 − v²)/4M, v = a/r₊: the square comes from the pull being the slope of 1/(r−ia)")
# W4 the slowdown factor is the law composing the rim with itself
chk("W4", Z(sp.re(sp.expand(w2(v))) - (1-v**2)) and Z(sp.simplify(pull*4*M) - sp.re(sp.expand(w2(vr)))),
    "1 − v² = Re(1 + x∘x), x = i·(rim speed): the rim's slowdown, squared, is the law composing the rim with itself")
# W5 helix = Im/Re of the same composed number; Euclidean continuation
Om = sp.im(sp.expand_complex(1/(rp - sp.I*a)))
helix = sp.simplify(Om/kap)
ratio = sp.im(sp.expand(w2(vr)))/sp.re(sp.expand(w2(vr)))
t = sp.symbols('t', positive=True)
helixE = sp.simplify((2*sp.I*t/(1+t**2))/sp.I)
rpE = M+sp.sqrt(M**2+al**2)
chk("W5", Z(helix - ratio) and Z(helixE.subs(t, al/rpE) - al/sp.sqrt(M**2+al**2)),
    "helix Ω_H/κ = Im/Re of (1 + x∘x): slowdown and helix are the real and imaginary parts of one composition; a = iα gives tanh u")
# W6 winding per imaginary-time loop counts entropy
beta_S = 8*sp.pi*M; S_S = 4*sp.pi*M**2
J = M*a; S_K = sp.pi*(rp**2+a**2); beta_K = 2*sp.pi/kap
chk("W6", Z(M*beta_S - 2*S_S) and Z(sp.simplify(beta_K*(M - 2*Om*J) - 2*S_K)),
    "mass's own phase per imaginary-time loop: M·β = 2S (Schwarzschild); with spin β(M − 2Ω_H J) = 2S (Smarr). Bigger mass: faster winding, longer loop, phase ∝ M²")
# W7 where G enters: river speed from the source's own winding rate needs the Planck length
G,hb,c0,m,R = sp.symbols('G hbar c m R', positive=True)
omega_m = m*c0**2/hb; lP2 = hb*G/c0**3
chk("W7", Z(2*G*m/(R*c0**2) - 2*omega_m*lP2/(c0*R)),
    "v² = 2GM/(rc²) = 2·(source winding rate)·l_P²/(c·r): the source's winding gives the river speed only with l_P² (so G) put in")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump(res, open("results.json","w"), indent=1)
