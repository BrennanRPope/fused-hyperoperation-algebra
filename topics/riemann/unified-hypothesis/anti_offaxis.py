"""Find the non-real zeros of the antiderivative parent F(z) = P(z) - L, P(z) = int_0^z H_0.
F(-z) = -(P(z) + L), so zeros of F in the left half are zeros of P + L in the right half."""
import numpy as np, json, kern
R = kern.riemann(600)
L = np.pi / 2 * kern.Phi(np.array([0.0]))[0]
def P(z):
    z = np.asarray(z, complex); sh = z.shape; zf = z.reshape(-1); out = np.empty(zf.shape, complex)
    for i in range(0, zf.size, 4000):
        out[i:i+4000] = np.sin(np.outer(zf[i:i+4000], R.u)) @ (R.w / R.u)
    return out.reshape(sh)
X, Y = np.meshgrid(np.linspace(0.5, 200, 1600), np.linspace(0.05, 140, 1100))
for sgn, name in [(+1, "P+L (F zeros mirrored from left half)")]:
    A = np.abs(P(X + 1j * Y) + sgn * L)
    m = (A < np.roll(A, 1, 0)) & (A < np.roll(A, -1, 0)) & (A < np.roll(A, 1, 1)) & (A < np.roll(A, -1, 1))
    m[0, :] = m[-1, :] = m[:, 0] = m[:, -1] = False
    z = (X + 1j * Y)[m]
    for _ in range(60):
        z = z - (P(z) + sgn * L) / R.H(z)
    ok = (np.abs(P(z) + sgn * L) < 1e-9) & (z.real > 0) & (z.real < 200) & (z.imag > -1e-6) & (z.imag < 140)
    z = np.unique(np.round(z[ok], 5))
    print(name, len(z)); print(z)
    json.dump({'zeros': [[float(a.real), float(a.imag)] for a in z]}, open(f'anti_offaxis_{"plus" if sgn>0 else "minus"}.json', 'w'))
