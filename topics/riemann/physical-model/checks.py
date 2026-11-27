"""Checks for the physical-model page. Run after phys.py."""
import json, numpy as np
O = json.load(open("phys.json")); res = []
def ck(n, ok, d=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + n + ("  " + d if d else ""))
# field rules around one zero at the origin
coul = lambda z: 1 / np.conj(z); vort = lambda z: 1j / np.conj(z); ours = lambda z: 1 / z
x = np.array([0.5, 1.3, -2.0]); y = 1j * np.array([0.5, 1.3, -2.0]); d = 0.9 * np.exp(1j * np.pi / 4)
ck("Q1 on the wire, ζ's rule = charges pushing apart", np.allclose(ours(x), coul(x)))
ck("Q2 across the wire, ζ's rule = charges pulling together (Coulomb mirrored)", np.allclose(ours(y), -coul(y)))
ck("Q3 at 45°, ζ's rule is perpendicular to the push (a swirl)", abs(np.real(ours(d) * np.conj(coul(d)))) < 1e-12)
ck("Q4 vortices are the Coulomb field turned 90°, not mirrored: on the wire they move across it",
   np.allclose(np.real(vort(x)), 0) and not np.allclose(vort(y), ours(y)))
ck("Q5 any real potential's force is anti-holomorphic, -conj(g'), so no Euclidean gradient system matches off the line",
   O["gradient_is_antiholomorphic"] < 1e-8, f"{O['gradient_is_antiholomorphic']:.1e}")
ck("Q6 ζ's rule is the gradient of 2 log r in the metric dx² − dy² (ascent along, descent across)",
   O["lorentzian_gradient"] < 1e-8, f"{O['lorentzian_gradient']:.1e}")
ck("Q7 Burgers (Cole–Hopf): poles of the gas velocity = zeros of the heat solution, exactly",
   O["burgers_poles"]["max_err"] < 1e-10, f"{O['burgers_poles']['max_err']:.1e}")
sl = O["burgers_front_slope"]
ck("Q8 a pole pair at distance y off the line makes a front of slope 4/y² (within 5%)",
   all(abs(sl[k] * float(k)**2 / 4 - 1) < 0.05 for k in sl), str({k: round(v * float(k)**2 / 4, 3) for k, v in sl.items()}))
cm = O["calogero_moser"]
ck("Q9 second-order form: z'' = −8 Σ 1/(z−w)³ (Calogero–Moser), exact", cm["max_rel_err_accel"] < 1e-12)
ck("Q10 the heat flow is the zero-energy surface of that system, on and off the line",
   cm["energy"] < 1e-9 and cm["energy_offline"] < 1e-9 and cm["identity_sumF2"] < 1e-9)
ck("Q11 2D dipoles pointing along the wire have exactly that pair energy (−cos 2φ / r²) on the line", O["dipole_match"] < 1e-12)
c = O["collapse"]
ck("Q12 rewound beads: 1/r push collapses as (time left)^½, 1/r² push as (time left)^⅓",
   abs(c["one_over_r"]["exponent"] - 0.5) < 0.01 and abs(c["one_over_r2"]["exponent"] - 1 / 3) < 0.01,
   f"{c['one_over_r']['exponent']:.3f}, {c['one_over_r2']['exponent']:.3f}")
h = O["hill"]
ck("Q13 our zeros sit at the Hill radius: median gap 1.0 ± 0.05 Hill radii, about half inside",
   abs(h["median_dstar_over_gap"] - 1) < 0.05 and 0.45 < h["share_inside"] < 0.55,
   f"median ratio {h['median_dstar_over_gap']:.3f}, inside {h['share_inside']:.3f}")
print(f"\n{sum(res)} of {len(res)} checks pass")
