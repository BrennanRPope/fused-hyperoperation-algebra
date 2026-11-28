# DRAFT. C2 follow-up (Nov 28, 2026): does the river's two-body rule (Part 7: for bodies at rest the logs of the
# clock rates add, i.e. static Weyl superposition) on two half-masses at +-alpha rebuild Euclidean Kerr?
import json, sympy as sp
res, out = {}, []
def log(name, ok, note): res[name] = bool(ok); out.append(f"{name} {'PASS' if ok else 'FAIL'}  {note}")
r, th, M, a, al, z = sp.symbols('r theta M a alpha z', positive=True)
# G1 Euclidean Kerr (t = -i tau, a = i alpha) keeps a real tau-phi cross term, so it is not a body at rest.
rho2L = r**2 + a**2*sp.cos(th)**2
gtphi = -2*M*r*a*sp.sin(th)**2/rho2L            # Lorentzian g_{t phi}
gtauphi = sp.simplify((gtphi*(-sp.I)).subs(a, sp.I*al))   # dt = -i dtau
ok1 = sp.simplify(sp.im(gtauphi.subs({r: 3, th: 1, M: 1, al: sp.Rational(1, 2)}))) == 0 and gtauphi != 0
log("G1", ok1, f"Euclidean g_tau_phi = {gtauphi}: real and nonzero, so not static")
# G2 Kerr's Ernst potential eps = 1 - 2M/(r - i a cos theta) gives the clock rate f = -g_tt (check on the axis).
eps_axis = 1 - 2*M/(r - sp.I*a)
f_axis = 1 - 2*M*r/(r**2 + a**2)
ok2 = sp.simplify(sp.re(sp.expand_complex(eps_axis)) - f_axis) == 0
log("G2", ok2, "Re(Ernst) on the axis = Kerr's clock rate squared")
# G3 On the axis (Weyl z = r - M) xi = (1 - eps)/(1 + eps) = M/(z - i a): one mass at z = i a, moments M (i a)^l
#    (Hansen 1974). Euclidean a = i alpha: xi = M/(z + alpha), ONE real mass at z = -alpha, moments M(-alpha)^l.
xi = sp.simplify(((1-eps_axis)/(1+eps_axis)).subs(r, z+M))
ok3 = sp.simplify(xi - M/(z - sp.I*a)) == 0
xiE = sp.simplify(xi.subs(a, sp.I*al)); w = sp.symbols('w', positive=True)
ser = sp.series(xiE.subs(z, 1/w), w, 0, 7).removeO()
momE = [sp.simplify(ser.coeff(w, l+1)) for l in range(6)]
ok3 = ok3 and all(sp.simplify(momE[l] - M*(-al)**l) == 0 for l in range(6))
log("G3", ok3, f"Euclidean axis potential {xiE}: moments {momE}")
# G4 The two-body rule on two half-masses at +-alpha: logs of clock rates add, so the axis potentials add:
#    moments (M/2)(alpha^l + (-alpha)^l): every odd moment is 0, but Euclidean Kerr's dipole is -M alpha.
ser2 = sp.series((M/2/(z-al) + M/2/(z+al)).subs(z, 1/w), w, 0, 7).removeO()
mom2 = [sp.simplify(ser2.coeff(w, l+1)) for l in range(6)]
ok4 = sp.simplify(mom2[1]) == 0 and sp.simplify(momE[1] + M*al) == 0
log("G4", ok4, f"two half-masses: moments {mom2}; dipole 0 vs Kerr's -M*alpha, so the rule fails")
res["moments_kerr"] = [str(m) for m in momE]; res["moments_two"] = [str(m) for m in mom2]
# G5 What carries the spin: smoothness identifies (tau, phi) ~ (tau + beta, phi + beta*Omega_E), a helix in imaginary
#    time. Turn per period = beta*Omega_E = 2 pi alpha / sqrt(M^2 + alpha^2); with alpha = M sinh u it is 2 pi tanh u.
rp = M + sp.sqrt(M**2+al**2); rm = M - sp.sqrt(M**2+al**2)
kapE = (rp-rm)/(2*(rp**2-al**2)); beta = 2*sp.pi/kapE; OmE = al/(rp**2-al**2)
turn = sp.simplify(beta*OmE)
u = sp.symbols('u', positive=True)
ok5 = sp.simplify(turn - 2*sp.pi*al/sp.sqrt(M**2+al**2)) == 0 and sp.simplify(turn.subs(al, M*sp.sinh(u)) - 2*sp.pi*sp.tanh(u)) == 0
log("G5", ok5, f"turn per imaginary-time period = {turn} = 2 pi tanh u for alpha = M sinh u")
n_ok = sum(1 for v in res.values() if v is True); n_all = sum(1 for v in res.values() if isinstance(v, bool))
out.append(f"{n_ok}/{n_all} passed"); print("\n".join(out)); json.dump(res, open("results.json", "w"), indent=1)
