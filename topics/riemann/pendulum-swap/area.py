"""The '1/2 is an area' finding, tested against the critical line. DRAFT.
g(x) = 1/(4 sinh^2(x/2)) - 1/x^2  (quantum minus classical oscillator curve).
Checks: area = -1/2 = zeta(0); Mellin of g = Gamma(s) zeta(s-1) on 0<Re s<2 (so the line Re s = 3/2 is inside);
partial fractions g = sum_{k!=0} 1/(x+2 pi i k)^2 (Matsubara / Poisson swap); zero-point energy reading;
kernel K(u) = e^{3u/2} g(e^u) is one-signed and its Fourier transform vanishes at zeta's zeros."""
import mpmath as mp, json
mp.mp.dps = 30
def g(x):
    x = mp.mpf(x) if not isinstance(x, mp.mpc) else x
    if abs(x) < mp.mpf('1e-6'): return -mp.mpf(1)/12 + x**2/240
    return 1/(4*mp.sinh(x/2)**2) - 1/x**2
def main():
    out = {}
    area = mp.quad(g, [0, 1, 10, mp.inf]); out["area"] = str(area); out["zeta0"] = str(mp.zeta(0))
    # zero-point energy: E_quantum(x) - E_classical(x) = coth(x/2)/2 - 1/x -> 1/2 as x -> inf
    out["E_q_minus_E_c_at_x=40"] = str(mp.coth(20)/2 - mp.mpf(1)/40)
    mel = {}
    for s in [mp.mpf('0.5'), mp.mpf(1), mp.mpf('1.5'), mp.mpc(1.5, 5), mp.mpc(1.5, 14.134725141734693)]:
        lhs = mp.quad(lambda x: x**(s-1)*g(x), [0, 1, 10, 40, mp.inf])
        rhs = mp.gamma(s)*mp.zeta(s-1)
        mel[str(s)] = [str(lhs), str(rhs), float(abs(lhs-rhs))]
    out["mellin"] = mel
    pf = {}
    for x in [0.3, 1.0, 4.0]:
        S = mp.nsum(lambda k: 2*(x**2 - 4*mp.pi**2*k**2)/(x**2 + 4*mp.pi**2*k**2)**2, [1, mp.inf])
        pf[str(x)] = [float(g(x)), float(S)]
    out["partial_fractions"] = pf
    # one sign: g < 0 on a grid
    out["g_negative_on_grid"] = all(g(mp.mpf(10)**(e/10.)) < 0 for e in range(-40, 31))
    # Fourier transform of K at t: = Gamma(3/2+it) zeta(1/2+it); at first zeros
    z = [float(mp.im(mp.zetazero(k))) for k in (1, 2, 3)]
    out["FT_K_at_zeros"] = [float(abs(mp.gamma(mp.mpc(1.5, t))*mp.zeta(mp.mpc(0.5, t)))) for t in z]
    out["first_zeros"] = z
    # xi(0) = xi(1) = 1/2 and zeta(0) = -1/2 as mirror of pole residue 1
    xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    out["xi_near_0_and_1"] = [str(xi(mp.mpf('1e-12'))), str(xi(1+mp.mpf('1e-12')))]
    json.dump(out, open("area.json", "w"), indent=1); print(json.dumps(out, indent=1))
if __name__ == "__main__": main()
