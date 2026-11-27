"""DRAFT, round 28 (Nov 27, 2026): the black hole's shadow, split into the river's time part, push and twist.

Rounds 25-27 found that, per unit energy, Kerr's orbits are Newton's spinning ring (the time part) plus the push plus
the twist. Light is the cleanest case: for a photon "per unit energy" is the only way to count, so the time part is
exactly Newton's pull on a corpuscle moving at c. Brennan asked for the shadow, the outline the Event Horizon Telescope
photographed, split the same way.

Findings:
1. Light bending, no spin (exact). The time part alone bends a ray by exactly 2 arctan(M/b), Newton's hyperbola for
   light (Soldner 1801; Einstein 1911): 2M/b far out, half of Einstein's 4M/b. The push gives the other half and all
   of the (15 pi/4)(M/b)^2.
2. No time part shadow (exact). With the time part alone, R/r^2 = r^2 + 2Mr - b^2 rises with r: there is no photon
   sphere and no horizon, every ray but the head-on one comes back out, and the hole casts no shadow at all. The
   horizon is the push's: A = r^2 + a^2 - 2Mr wPP vanishes only with P.P on. With the push the photon sphere is at 3M
   and the shadow's radius is sqrt(27) M (5.196 M).
3. With spin. Time part + twist still has no horizon and no shadow. Time part + push (no twist) has the horizon and a
   shadow centred on the hole, slightly squashed (a = 0.9, edge-on: 10.21 M wide, 9.65 M tall). The twist then shifts
   it sideways by about 2a sin i (frame dragging: prograde rays get closer), flattens the prograde side into Kerr's
   D shape, and restores the height to exactly 2 sqrt(27) M edge-on. All three together match Bardeen's (1973)
   closed-form Kerr outline.
4. Sgr A*: from the stellar orbits' mass and distance (no EHT input), GM/(c^2 D) is about 4.9-5.0 microarcseconds, so
   the push's shadow is about 2 sqrt(27) x 5.0 = 51-52 microarcseconds across, within about +-4% for any spin
   (inclination below 50 degrees). The EHT (2022) measured a ring of 51.8 +- 2.3 microarcseconds and inferred a shadow
   of 48.7 +- 7.0. With the time part alone the prediction is no shadow.
Run: python3 shadow_split.py  (writes shadow_split.json, about a minute)
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
r, M, a, xi, eta, b = sp.symbols("r M a xi eta b", real=True)
NAMES = ("EE", "PP", "SS", "EP", "ES", "PS")
CASES = {"time part (Newton's ring)": ("EE",), "ring + push": ("EE", "PP", "EP"),
         "ring + twist": ("EE", "SS", "ES"), "Kerr (all)": NAMES}
D = r**2 + a**2


def Rnull(on):
    """Null rays per unit energy in the Kerr-Schild shell with the given pieces on; xi = L/E, eta = Carter's Q/E^2."""
    d = {n: (1 if n in on else 0) for n in NAMES}
    f = 2*M*r
    S = -a*xi/D
    C = eta + xi**2
    A = D - f*d["PP"]
    B = -f*(d["EP"] + d["PS"]*S)
    Cr = -a**2*xi**2/D - r**2 + C - f*(d["EE"] + d["SS"]*S**2 + 2*d["ES"]*S)
    return A, B, Cr, sp.together(B**2 - A*Cr)


def identities():
    zero = lambda x: sp.simplify(x) == 0
    RK = (D - a*xi)**2 - (r**2 - 2*M*r + a**2)*(eta + (xi - a)**2)
    out = {"all_pieces_is_kerr_null": zero(Rnull(NAMES)[3] - RK)}
    for n, on in CASES.items():
        A = Rnull(on)[0]
        out[f"horizon_factor[{n}]"] = str(sp.factor(A))
    out["horizon_only_with_push"] = all(sp.Poly(Rnull(on)[0], r).degree() == 2 and
                                        (sp.diff(Rnull(on)[0], M) != 0) == ("PP" in on) for on in CASES.values())
    # no spin
    R_time = sp.factor(Rnull(("EE",))[3].subs({a: 0, eta: b**2 - xi**2}))
    R_push = sp.factor(Rnull(("EE", "PP", "EP"))[3].subs({a: 0, eta: b**2 - xi**2}))
    out["a0_time_part"] = str(R_time)
    out["a0_ring_plus_push"] = str(R_push)
    out["a0_time_part_is_newton_light"] = zero(R_time - r**2*(r**2 + 2*M*r - b**2))
    out["a0_push_is_schwarzschild"] = zero(R_push - r**2*(r**4 - b**2*r*(r - 2*M))/r**2 * 1) or \
        zero(sp.cancel(R_push) - (r**4 - b**2*r*(r - 2*M)))
    # photon sphere: time part has none (R/r^2 strictly increasing for r > 0), push: r = 3M, b^2 = 27M^2
    g = r**2 + 2*M*r - b**2
    out["a0_time_part_no_photon_sphere"] = sp.solve([g, sp.diff(g, r)], [r, b], dict=True) == [] or all(
        s_[r].subs(M, 1) <= 0 for s_ in sp.solve([g, sp.diff(g, r)], [r, b], dict=True))
    h = r**3 - b**2*(r - 2*M)
    sol = [s_ for s_ in sp.solve([h, sp.diff(h, r)], [r, b], dict=True) if s_[b].subs(M, 1) > 0]
    out["a0_push_photon_sphere"] = str(sol)
    out["a0_push_photon_sphere_3M_sqrt27M"] = len(sol) == 1 and sp.simplify(sol[0][r] - 3*M) == 0 and \
        sp.simplify(sol[0][b] - 3*sp.sqrt(3)*M) == 0
    return out


def bending():
    """Deflection of light, no spin: time part exactly 2 arctan(M/b); with the push, Schwarzschild."""
    rows = []
    for bv in (8.0, 20.0, 100.0, 1000.0):
        u0 = 1/(-1 + np.sqrt(1 + bv**2))
        g = lambda t: bv*2*t*u0/np.sqrt(1 + 2*u0*(1 - t**2) - bv**2*u0**2*(1 - t**2)**2)
        d_time = 2*quad(g, 0, 1, limit=400, epsabs=1e-14, epsrel=1e-13)[0] - np.pi
        # Schwarzschild: (du/dphi)^2 = 1/b^2 - u^2 + 2 u^3, turning point u0 = largest root below 1/3
        rts = np.roots([2, -1, 0, 1/bv**2]); u0s = min(x.real for x in rts if abs(x.imag) < 1e-12 and 0 < x.real)
        gs = lambda t: 2*t*u0s/np.sqrt(1/bv**2 - (u0s*(1 - t**2))**2 + 2*(u0s*(1 - t**2))**3)
        d_full = 2*quad(gs, 0, 1, limit=400, epsabs=1e-14, epsrel=1e-13)[0] - np.pi
        rows.append({"b": bv, "time_part": d_time, "newton_2_arctan": 2*np.arctan(1/bv), "full": d_full,
                     "full_series_4/b+15pi/4b^2": 4/bv + 15*np.pi/4/bv**2, "push_share": d_full - d_time,
                     "time_over_full": d_time/d_full})
    return rows


def photon_orbits(on, av):
    A, B, Cr, Rt = Rnull(on)
    N = sp.numer(Rt).subs({M: 1, a: av})
    R0, R1 = N.subs(eta, 0), sp.diff(N, eta)
    etas = -R0/R1
    eq = sp.numer(sp.together(sp.diff(N, r).subs(eta, etas)))
    P = sp.Poly(sp.expand(eq), xi)
    cf = [sp.lambdify(r, c) for c in P.all_coeffs()]
    fe = sp.lambdify((r, xi), etas)
    R2 = sp.lambdify((r, xi, eta), sp.diff(N, r, 2))
    rh = 1 + np.sqrt(1 - av**2) if "PP" in on else None

    def scan(rs):
        out = []
        for rv in rs:
            for x in np.roots([c(rv) for c in cf]):
                if abs(x.imag) < 1e-9:
                    e = fe(rv, x.real)
                    out.append((rv, x.real, e, R2(rv, x.real, e)))
        return np.array(out).reshape(-1, 4)
    return scan, rh


def outline(on, av, inc, n=24000):
    scan, rh = photon_orbits(on, av)
    lo = (rh + 1e-7) if rh else 0.05

    def keep(orb):
        pts = []
        for rv, x, e, r2 in orb:
            if r2 <= 0:
                continue
            b2 = e + av**2*np.cos(inc)**2 - x**2/np.tan(inc)**2
            if b2 >= 0:
                pts.append((-x/np.sin(inc), np.sqrt(b2), rv))
        return np.array(pts).reshape(-1, 3)
    o = keep(scan(np.linspace(lo, 6, 30000)))
    if len(o):
        r1, r2 = o[:, 2].min(), o[:, 2].max()
        dr = 6/30000
        o = keep(scan(np.linspace(max(lo, r1 - 2*dr), r2 + 2*dr, n)))
        o = o[np.argsort(o[:, 0])]
    return o, rh


def bardeen(av, inc, n=40000):
    rp = 2*(1 + np.cos(2/3*np.arccos(-av)))
    rm = 2*(1 + np.cos(2/3*np.arccos(av)))
    rr = np.linspace(rp, rm, n)
    x = (rr**2*(3 - rr) - av**2*(rr + 1))/(av*(rr - 1))
    e = rr**3*(4*av**2 - rr*(rr - 3)**2)/(av**2*(rr - 1)**2)
    b2 = e + av**2*np.cos(inc)**2 - x**2/np.tan(inc)**2
    ok = b2 >= 0
    return np.c_[-x[ok]/np.sin(inc), np.sqrt(b2[ok])]


def shadows():
    res = {}
    for av in (0.5, 0.9, 0.998):
        for incd in (90, 50, 17):
            inc = np.radians(incd)
            for name, on in CASES.items():
                o, rh = outline(on, av, inc)
                key = f"a={av}, i={incd}, {name}"
                if rh is None:
                    res[key] = {"horizon": False, "shadow": False,
                                "light_rings_r": [float(o[:, 2].min()), float(o[:, 2].max())] if len(o) else None}
                    continue
                d = {"horizon": True, "shadow": True, "alpha_min": float(o[:, 0].min()),
                     "alpha_max": float(o[:, 0].max()), "centre": float((o[:, 0].min() + o[:, 0].max())/2),
                     "width": float(o[:, 0].max() - o[:, 0].min()), "height": float(2*o[:, 1].max())}
                if on == NAMES:
                    bd = bardeen(av, inc)
                    d.update({"bardeen_alpha_min": float(bd[:, 0].min()), "bardeen_alpha_max": float(bd[:, 0].max()),
                              "bardeen_height": float(2*bd[:, 1].max())})
                if (av == 0.9 and incd in (90, 17)) or (av == 0.998 and incd == 90):
                    step = max(1, len(o)//400)
                    d["outline"] = o[::step, :2].tolist()
                res[key] = d
                print(key, {k_: v_ for k_, v_ in d.items() if k_ != "outline"}, flush=True)
    return res


def eht(sh):
    G, c, Msun, pc = 6.67430e-11, 2.99792458e8, 1.98847e30, 3.0857e16
    uas = np.degrees(1)*3600e6
    priors = {"GRAVITY 2019 (VLTI)": (4.152e6, 8178.0), "Do et al. 2019 (Keck)": (3.975e6, 7946.0)}
    out = {}
    # spread of Kerr's shadow size over spin and inclination below 50 degrees: mean of width and height
    sizes = [0.5*(v_["width"] + v_["height"]) for k_, v_ in sh.items()
             if k_.endswith("Kerr (all)") and (", i=50," in k_ or ", i=17," in k_)]
    for nm, (m, dist) in priors.items():
        tg = G*m*Msun/c**2/(dist*pc)*uas
        out[nm] = {"theta_g_uas": tg, "push_shadow_a0_uas": 2*np.sqrt(27)*tg,
                   "kerr_range_uas": [min(sizes)*tg, max(sizes)*tg], "time_part_alone_uas": 0.0}
    out["EHT_2022_ring_diameter_uas"] = [51.8, 2.3]
    out["EHT_2022_shadow_diameter_uas"] = [48.7, 7.0]
    return out


def main():
    res = {"identities": identities()}
    print(res, flush=True)
    res["bending"] = bending()
    print(res["bending"], flush=True)
    res["shadows"] = shadows()
    res["eht"] = eht(res["shadows"])
    print(res["eht"], flush=True)
    json.dump(res, open(os.path.join(HERE, "shadow_split.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
