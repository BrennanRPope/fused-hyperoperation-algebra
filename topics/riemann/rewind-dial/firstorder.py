"""First-order motion of the zeros as the rewind angle leaves 0. F = sum w e^{z X}, X = cos(t) a M + sin(t) H.
At t = 0: dz/dt = -F_t/F_z with F_t = z * sum w H e^{z a M}. On the line z = iy: F_t is imaginary and F_z is
imaginary, so dz/dt is REAL: the zeros leave the line at first order, at rate proportional to height y. DRAFT."""
import numpy as np, json
from rotate import setup, F, start_zeros, newton
out = {}
for n in (8, 10, 12, 14):
    S = setup(n); z0 = start_zeros(S, 36.0)
    G = (S["w"]*S["H"]*np.exp(np.multiply.outer(z0, S["a"]*S["M"]))).sum(-1)
    _, Fz = F(S, z0, 0.0); zp = -(z0*G)/Fz
    # numerical check with a tiny angle
    t = 1e-5; z = newton(S, z0, t); num = (z - z0)/t
    rate = zp.real/z0.imag
    out[n] = dict(y=z0.imag.tolist(), dRe=zp.real.tolist(), dIm=zp.imag.tolist())
    print("n=%d: |Im dz/dt|max %.1e (should be ~0); formula vs numeric max diff %.1e; Re dz/dt / height: median %.3f, range [%.3f, %.3f]"
          % (n, np.abs(zp.imag).max(), np.abs(num - zp).max(), np.median(np.abs(rate)), rate.min(), rate.max()))
json.dump(out, open("firstorder.json", "w"))
