"""DRAFT, round 18 (Nov 27, 2026): does Kepler's law in the shift's distance survive orbits tilted out of the equator?

Round 17 found that in the equator, the shifted source's pull is exactly -M/r in the shift's own radius r (the
distance to Kerr's ring), so Kepler's law Omega^2 = M/r^3 holds exactly there and the whole tower of shapes leaves
one -3. This file asks whether that is an equatorial accident or holds for every orbit.

1. Newton around the shifted source, in its own (oblate) coordinates: the pull is -M r/Sigma with
   Sigma = r^2 + a^2 cos^2(theta), and the motion separates (Euler's two fixed centres M/2 at +-ia). Spherical
   orbits (fixed r) exist at every tilt; their up-and-down and around-the-axis rates are found by quadrature.
2. Answer: Kepler's period law in r holds only in the equator. Tilted orbits feel the tower again. The polar orbit
   is exact in closed form: its period is Kepler's times (2/pi)(2K(k) - E(k)), k = a/r (complete elliptic
   integrals), an endless series 1 - 3k^2/2 + 33k^4/32 - ...
3. What survives exactly: (a) Kepler's energy law E = -M/(2r) for every orbit with no spin about the axis;
   (b) in the equator, the up-and-down rate nu^2/Omega^2 = 1 + 3a^2/r^2, so radial plus vertical gives exactly 2
   (Laplace's rule: the shifted source is empty space everywhere off its ring).
4. Kerr: the equatorial up-and-down rate is (Omega_theta/Omega_phi)^2 = 1 - 4 chi u^(3/2) + 3 chi^2 u^2 (known), and
   radial plus vertical is exactly 2 - 6u + 4 chi u^(3/2): Kerr's shape term cancels, exactly as Laplace's rule says
   it must for a Newtonian body in empty space. Both of Kerr's shape terms equal the shifted source's, term for term.
Run: python3 tilt.py  (writes tilt.json, about 20 seconds)
"""
import json
import os

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30
R, z, a, M, r, th = sp.symbols("R z a M r theta", positive=True)


# ---------------------------------------------------------------- 1. Newton, the shifted source in cylindrical form
def phi_cyl():
    """Real part of the pull of a point mass M at z = i a, written as the mean of the sources at +-ia."""
    return -M/2*((R**2 + (z - sp.I*a)**2)**sp.Rational(-1, 2) + (R**2 + (z + sp.I*a)**2)**sp.Rational(-1, 2))


def equator_rates():
    """Omega^2, kappa^2/Omega^2 (= T) and nu^2/Omega^2 in the equator, in flat R and in r = sqrt(R^2 - a^2)."""
    P = phi_cyl()
    PR = sp.simplify(sp.diff(P, R).subs(z, 0))
    PRR = sp.simplify(sp.diff(P, R, 2).subs(z, 0))
    Pzz = sp.simplify(sp.diff(P, z, 2).subs(z, 0))
    W2 = sp.simplify(PR/R)
    T = sp.simplify(3 + R*PRR/PR)
    V = sp.simplify(Pzz/W2)
    to_r = {R: sp.sqrt(r**2 + a**2)}
    return dict(omega2=W2, T=T, V=V, omega2_r=sp.simplify(W2.subs(to_r)), T_r=sp.simplify(T.subs(to_r)),
                V_r=sp.simplify(V.subs(to_r)), laplace=sp.simplify(T + V))


def oblate_potential_ok():
    """In oblate coordinates (rho = sqrt(r^2+a^2) sin, z = r cos) the complex distance is r - i a cos, so the pull
    is -M r/Sigma. Checked by squaring."""
    c = sp.cos(th); s = sp.sin(th)
    d2 = (r**2 + a**2)*s**2 + (r*c - sp.I*a)**2
    return sp.simplify(sp.expand(d2 - (r - sp.I*a*c)**2)) == 0


# ---------------------------------------------------------------- 2. spherical orbits at every tilt (quadrature)
def EK(Mv, av, rv, L):
    """Energy and separation constant of the spherical orbit at radius r with axial spin L (Newton, per unit mass).
    Radial: R(r) = 2E r^2 + 2M r + a^2 L^2/(r^2+a^2) - K, with R = R' = 0."""
    E = (-2*Mv + 2*av**2*L**2*rv/(rv**2 + av**2)**2)/(4*rv)
    K = 2*E*rv**2 + 2*Mv*rv + av**2*L**2/(rv**2 + av**2)
    return E, K


def sphere(Mv, av, rv, L):
    """Up-and-down and around-the-axis rates of a spherical orbit, by quadrature in Mino time
    (dt = Sigma dlambda, p_theta^2 = K - L^2/sin^2 + 2E a^2 cos^2, dphi/dlambda = L (1/sin^2 - a^2/(r^2+a^2)))."""
    E, K = EK(Mv, av, rv, L)
    A = -2*E*av**2; B = 2*E*av**2 - K; C = K - L**2          # in w = cos^2: (dw)^2 ~ A w^2 + B w + C
    disc = mp.sqrt(B**2 - 4*A*C)
    w1, w2 = sorted([(-B - disc)/(2*A), (-B + disc)/(2*A)])
    f = lambda p: 1/(mp.sqrt(A)*mp.sqrt(w2 - w1*mp.sin(p)**2))
    Tt = 4*mp.quad(lambda p: (rv**2 + av**2*w1*mp.sin(p)**2)*f(p), [0, mp.pi/2])
    Ph = 4*mp.quad(lambda p: L*(1/(1 - w1*mp.sin(p)**2) - av**2/(rv**2 + av**2))*f(p), [0, mp.pi/2])
    return dict(E=E, K=K, tilt_deg=mp.degrees(mp.asin(mp.sqrt(w1))), Om_th=2*mp.pi/Tt, Om_ph=Ph/Tt)


def L_equator(Mv, av, rv):
    """Axial spin of the circular equatorial orbit: the turning latitude closes, K = L^2."""
    return mp.findroot(lambda L: EK(Mv, av, rv, L)[1] - L**2, mp.sqrt(Mv*rv))


def polar_exact(k):
    """Polar orbit (L = 0): E = -M/2r, K = M r, and Omega_theta^2 r^3/M = (pi/2)^2/(2K(k) - E(k))^2, k = a/r."""
    return (mp.pi/2)**2/(2*mp.ellipk(k**2) - mp.ellipe(k**2))**2


def polar_series(n=5):
    m = sp.symbols("m", positive=True)
    ser = sp.series((sp.pi/2)**2/(2*sp.elliptic_k(m) - sp.elliptic_e(m))**2, m, 0, n).removeO()
    return [sp.nsimplify(ser.coeff(m, j)) for j in range(n)]


def polar_symbolic():
    """L = 0 in the formulas above gives E = -M/(2r) and K = M r: Kepler's energy law in r, exactly."""
    Lz = sp.symbols("L")
    E = (-2*M + 2*a**2*Lz**2*r/(r**2 + a**2)**2)/(4*r)
    K = 2*E*r**2 + 2*M*r + a**2*Lz**2/(r**2 + a**2)
    return sp.simplify(E.subs(Lz, 0) + M/(2*r)), sp.simplify(K.subs(Lz, 0) - M*r)


def spherical_ok(Mv, av, rv, L):
    """The constants really give a spherical orbit: the radial function and its slope vanish at r."""
    E, K = EK(Mv, av, rv, L)
    Rf = lambda x: 2*E*x**2 + 2*Mv*x + av**2*L**2/(x**2 + av**2) - K
    return abs(Rf(rv)), abs(mp.diff(Rf, rv))


# ---------------------------------------------------------------- 3. Kerr, equatorial up-and-down rate
def kerr_vertical():
    rr, MM, aa, E, L = sp.symbols("r M a E L", positive=True)
    s, chi = sp.symbols("s chi", positive=True)
    D = rr**2 - 2*MM*rr + aa**2
    P = E*(rr**2 + aa**2) - aa*L
    den = sp.sqrt(1 - 3*s**2 + 2*chi*s**3)
    sub = {rr: MM/s**2, aa: chi*MM, E: (1 - 2*s**2 + chi*s**3)/den,           # Bardeen, Press & Teukolsky 1972
           L: (MM/s)*(1 - 2*chi*s**3 + chi**2*s**4)/den}
    phid = -(aa*E - L) + aa*P/D                                               # dphi/dlambda in the equator
    V = sp.expand(sp.simplify(((aa**2*(1 - E**2) + L**2)/phid**2).subs(sub)))  # Carter's theta equation, small z
    T = 1 - 6*s**2 + 8*chi*s**3 - 3*chi**2*s**4
    return dict(V=V, T=T, sum=sp.expand(V + T), s=s, chi=chi)


def main():
    out = {}
    eq = equator_rates()
    out["newton_equator"] = {k: str(v) for k, v in eq.items()}
    out["oblate_potential_ok"] = bool(oblate_potential_ok())
    print("equator", out["newton_equator"], flush=True)
    Mv = mp.mpf(1)
    curves = {}
    for av, rv in (("0.8", "3"), ("0.8", "5")):
        av, rv = mp.mpf(av), mp.mpf(rv)
        Leq = L_equator(Mv, av, rv)
        rows = []
        for fr in ("0.9999", "0.99", "0.95", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1",
                   "0.05", "0.01", "0.001"):
            L = mp.mpf(fr)*Leq
            s = sphere(Mv, av, rv, L)
            k = Mv/rv**3
            res = spherical_ok(Mv, av, rv, L)
            rows.append(dict(L_frac=fr, tilt_deg=float(s["tilt_deg"]), Om_th2=float(s["Om_th"]**2/k),
                             Om_ph2=float(s["Om_ph"]**2/k), energy=float(-2*rv*s["E"]/Mv),
                             residual=float(max(res))))
        curves[f"a={mp.nstr(av, 2)}, r={mp.nstr(rv, 2)}"] = dict(
            L_equator=float(Leq), rows=rows, equator_vertical=float(1 + 3*av**2/rv**2),
            polar_exact=float(polar_exact(av/rv)))
        print(f"a={av} r={rv}", [(round(x["tilt_deg"], 1), round(x["Om_th2"], 5), round(x["Om_ph2"], 5))
                                 for x in rows], "polar", polar_exact(av/rv), flush=True)
    out["spheres"] = curves
    out["polar_series"] = [str(c) for c in polar_series()]
    pe, pk = polar_symbolic()
    out["polar_energy_minus_kepler"] = str(pe); out["polar_K_minus_Mr"] = str(pk)
    kv = kerr_vertical()
    out["kerr"] = {"vertical": str(kv["V"]), "T": str(kv["T"]), "sum": str(kv["sum"])}
    print("kerr", out["kerr"], "polar series", out["polar_series"], flush=True)
    json.dump(out, open(os.path.join(HERE, "tilt.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
