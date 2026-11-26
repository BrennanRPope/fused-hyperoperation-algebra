"""Polishing step (cheap, ~seconds): rerun the load-bearing identities at 30 digits, and give the
Monte Carlo checks an error bar. Run: python3 polish.py > polish_output.txt"""
import mpmath as mp, numpy as np, math, time
mp.mp.dps = 30
t0 = time.time()
out = []
def rep(name, val): out.append(name); print(f"{name}: {val}")
# F4: averages over the loop phase, by quadrature of the exact cosine form
u, v = mp.mpf('0.4225'), mp.mpf('0.16')
A = u + v + 2*u*v; B = 2*mp.sqrt(u*v*(1+u)*(1+v))
rho = lambda ph: A - B*mp.cos(ph)
m1 = mp.quad(rho, [0, 2*mp.pi])/(2*mp.pi)
m2 = mp.quad(lambda ph: mp.log(1 + rho(ph)), [0, 2*mp.pi])/(2*mp.pi)
m3 = mp.quad(lambda ph: 1/(1 + rho(ph)), [0, 2*mp.pi])/(2*mp.pi)
rep("F4.3 polish |<rho> - u∘₂v|", mp.nstr(abs(m1 - A), 3))
rep("F4.4 polish |<ln(1+rho)> - ln(1+u) - ln(1+v)|", mp.nstr(abs(m2 - mp.log(1+u) - mp.log(1+v)), 3))
rep("F4.5 polish |1/<T> - 1 - (u+v)|", mp.nstr(abs(1/m3 - 1 - (u + v)), 3))
for nu in (2, mp.mpf(1)/2, mp.mpc(-0.5, 1.7)):
    lhs = mp.quad(lambda ph: mp.legenp(nu, 0, 1 + 2*rho(ph), type=3), [0, 2*mp.pi])/(2*mp.pi)
    rhs = mp.legenp(nu, 0, 1 + 2*u, type=3)*mp.legenp(nu, 0, 1 + 2*v, type=3)
    rep(f"F4.6 polish Legendre product, nu={mp.nstr(nu, 3)}, relative error", mp.nstr(abs(lhs - rhs)/abs(rhs), 3))
# F4.7: exact mean for N barriers is the ∘₂ power; compare Monte Carlo with its standard error
def delta_M(g, kw):
    be = g/(2*kw); return np.array([[1 - 1j*be, -1j*be], [1j*be, 1 + 1j*be]])
def prop(kw, Lx): return np.array([[np.exp(1j*kw*Lx), 0], [0, np.exp(-1j*kw*Lx)]])
def odds(M): T = 1/abs(M[1, 1])**2; return (1 - T)/T
rng = np.random.default_rng(11); Mg = delta_M(1.0, 1.0); u1 = odds(Mg); N = 6; n = 60000
vals = np.empty(n); logs = np.empty(n)
for i in range(n):
    M = Mg.copy()
    for _ in range(1, N): M = Mg @ prop(1.0, rng.uniform(0, np.pi)) @ M
    vals[i] = odds(M)
logs = np.log1p(vals)
pred = ((1 + 2*u1)**N - 1)/2
rep("F4.7 polish mean odds, 6 barriers (Monte Carlo ± 1 s.e. vs exact ∘₂ power)",
    f"{vals.mean():.4f} ± {vals.std(ddof=1)/math.sqrt(n):.4f} vs {pred:.4f}  ({(vals.mean()-pred)/(vals.std(ddof=1)/math.sqrt(n)):+.2f} s.e.)")
rep("F4.8 polish log-mean, 6 barriers (± 1 s.e. vs N ln(1+u))",
    f"{logs.mean():.5f} ± {logs.std(ddof=1)/math.sqrt(n):.5f} vs {N*math.log1p(u1):.5f}  ({(logs.mean()-N*math.log1p(u1))/(logs.std(ddof=1)/math.sqrt(n)):+.2f} s.e.)")
# F3.5: Jacobi addition at 30 digits
m = mp.mpf('0.8'); a, b = mp.mpf('0.31'), mp.mpf('-0.77')
sn = lambda x: mp.ellipfun('sn', x, m=m); cn = lambda x: mp.ellipfun('cn', x, m=m); dn = lambda x: mp.ellipfun('dn', x, m=m)
x1, y1, x2, y2 = sn(a), cn(a)/dn(a), sn(b), cn(b)/dn(b)
den = m*x1*x2*y1*y2
rep("F3.5 polish |Edwards sum - sn(a+b)| at m = 0.8, 30 digits", mp.nstr(abs((x1*y2 + y1*x2)/(1 + den) - sn(a + b)), 3))
# F1.5 at 30 digits
kk, dd = mp.mpf('0.7'), mp.mpf('0.3'); kap = dd - kk**2/4
T = lambda t: mp.tan(mp.sqrt(kap)*t)/mp.sqrt(kap); G = lambda t: T(t)/(1 - kk*T(t)/2)
Lw = lambda x, y: (x + y + kk*x*y)/(1 - dd*x*y)
rep("F1.5 polish |G(u1+u2) - L(G(u1),G(u2))|", mp.nstr(abs(G(mp.mpf('0.21') + mp.mpf('0.13')) - Lw(G(mp.mpf('0.21')), G(mp.mpf('0.13')))), 3))
print(f"polish done in {time.time()-t0:.1f} s")
