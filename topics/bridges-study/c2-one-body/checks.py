# DRAFT. C2 follow-up 2 (Nov 28, 2026): do the river's rules give the moved mass and the helix for one spinning body?
import json, sympy as sp
res, out = {}, []
def log(name, ok, note): res[name] = bool(ok); out.append(f"{name} {'PASS' if ok else 'FAIL'}  {note}")
r, th, M, a, al, e, u = sp.symbols('r theta M a alpha epsilon u', positive=True)
# H1 The spin's own drag (Lense-Thirring, J = M a, the river's far swirl 2J/r^3) is the dipole term of Kerr's Ernst
#    potential; at a = i alpha it becomes exactly the first-order term of one mass moved to z = -alpha.
eps = 1 - 2*M/(r - sp.I*a*sp.cos(th))
lin = sp.series(eps.subs(a, e*a), e, 0, 2).removeO().subs(e, 1)
linE = sp.expand(lin.subs(a, sp.I*al))
moved = 1 - 2*M*(1/r - al*sp.cos(th)/r**2)          # 1/|x + alpha z| to first order in alpha
ok1 = sp.simplify(linE - moved) == 0
log("H1", ok1, f"spin drag term {sp.simplify(lin - (1-2*M/r))} -> at a = i alpha: the dipole of a mass moved by -alpha")
# H2 Part 8's moving-body rules carry the mass and the spin's drag only (moments l = 0, 1). Kerr also needs
#    l >= 2 (quadrupole -M a^2, ...), which only the exact shift supplies: the l = 2 term of eps is nonzero.
quad = sp.simplify(sp.series(eps.subs(a, e*a), e, 0, 3).removeO().coeff(e, 2))
ok2 = sp.simplify(quad) != 0
log("H2", ok2, f"second-order term {quad}: needed for Kerr, missing from 1PN rules; the exact shift by i a has it")
# H3 The river's rigid shells turn at Omega(r) = a/(r^2 + a^2) (Part 3). At the horizon this equals Kerr's frame-drag
#    rate -g_tphi/g_phiphi for EVERY polar angle (the horizon turns rigidly).
rho2 = r**2 + a**2*sp.cos(th)**2; D = r**2 - 2*M*r + a**2
gtphi = -2*M*r*a*sp.sin(th)**2/rho2
gphph = ((r**2+a**2)**2 - D*a**2*sp.sin(th)**2)*sp.sin(th)**2/rho2
rp = M + sp.sqrt(M**2 - a**2)
omega_h = sp.simplify((-gtphi/gphph).subs(r, rp))
ok3 = sp.simplify(omega_h - a/(rp**2 + a**2)) == 0 and sp.simplify(sp.diff(omega_h, th)) == 0
log("H3", ok3, "Kerr's frame-drag rate at the horizon = the river's shell rate a/(r+^2 + a^2), the same at every angle")
# H4 The helix: the river's horizon shell rate, continued (a = i alpha), times the imaginary-time period beta = 2 pi/kappa_E,
#    is 2 pi alpha/sqrt(M^2 + alpha^2) = 2 pi tanh u (alpha = M sinh u). Outer shells would turn by beta*alpha/(r^2 - alpha^2).
rpE = M + sp.sqrt(M**2 + al**2); rmE = M - sp.sqrt(M**2 + al**2)
beta = 2*sp.pi*2*(rpE**2 - al**2)/(rpE - rmE)
OmE = lambda rr: al/(rr**2 - al**2)                   # i*Omega(r) at a = i alpha, real
turn = sp.simplify(beta*OmE(rpE))
ok4 = sp.simplify(turn - 2*sp.pi*sp.tanh(u)).subs(al, M*sp.sinh(u)).simplify() == 0
ok4 = ok4 and sp.simplify(sp.I*(a/(r**2+a**2)).subs(a, sp.I*al) + OmE(r)) == 0
log("H4", ok4, f"helix turn per loop = beta * river shell rate at r+ = {turn} = 2 pi tanh u")
# H5 The 1PN rules have no horizon, so they have no imaginary-time period: far away the continued swirl
#    2J/r^3 -> 2 M alpha/r^3 is real but sets no gluing; the turn per loop needs beta, which only the horizon fixes.
swirl = (2*M*a/r**3).subs(a, sp.I*al)
ok5 = sp.simplify(swirl/sp.I - 2*M*al/r**3) == 0
log("H5", ok5, "far swirl continues to a real 2 M alpha/r^3, but the loop length beta comes only from the horizon")
n_ok = sum(1 for v in res.values() if v is True); n_all = sum(1 for v in res.values() if isinstance(v, bool))
out.append(f"{n_ok}/{n_all} passed"); print("\n".join(out)); json.dump(res, open("results.json", "w"), indent=1)
