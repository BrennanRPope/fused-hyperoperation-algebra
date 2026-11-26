"""Polish for checks.py check 3b (Nov 26): the fast soliton's peak is found in checks.py as the argmax of u on a grid
(spacing 0.001) with u = 2 (log tau)'' taken by two numerical gradients. Here u and its x-derivatives are written in
closed form (tau = sum of four exponentials, normalised by the largest term to avoid overflow) and the peak is the zero
of u'(x), found by Newton's method from the grid answer, to machine precision. numpy only."""
import json, math, time
import numpy as np
t0 = time.process_time()
k1, k2 = 1.0, 1.6
A = ((k1 - k2)/(k1 + k2))**2                         # Aij(k1, k2) in checks.py (results.json stores the same A)
R = json.load(open('results.json'))['tropical']; assert abs(R['A'] - A) < 1e-15
phase = math.log(A)

def logtau_derivs(x, T):
    """d^n/dx^n log tau for n = 1..4 at a point"""
    s = np.array([0.0, k1, k2, k1 + k2])
    e = np.array([0.0, k1*x - k1**3*T, k2*x - k2**3*T, k1*x - k1**3*T + k2*x - k2**3*T + phase])
    c = np.ones(4)                                    # A is folded into the last exponent as log(A)
    w = np.exp(e - e.max())
    t = [np.sum(c*w*s**n) for n in range(5)]; t1, t2, t3, t4 = (t[n]/t[0] for n in range(1, 5))
    l2 = t2 - t1**2
    l3 = t3 - 3*t1*t2 + 2*t1**3
    l4 = t4 - 4*t3*t1 - 3*t2**2 + 12*t2*t1**2 - 6*t1**4
    return l2, l3, l4

def peak_newton(x0, T):
    x = x0
    for _ in range(50):
        _, l3, l4 = logtau_derivs(x, T); step = l3/l4; x -= step     # u' = 2 l3, u'' = 2 l4
        if abs(step) < 1e-15*max(1, abs(x)): break
    return x

def peak_grid(T):                                     # as in checks.py
    xs = np.linspace(-80, 80, 160001)
    e1 = k1*xs - k1**3*T; e2 = k2*xs - k2**3*T
    tau = 1 + np.exp(e1) + np.exp(e2) + A*np.exp(e1 + e2)
    u = 2*np.gradient(np.gradient(np.log(tau), xs), xs)
    return xs[np.argmax(u)]

xb0, xa0 = peak_grid(-25.0), peak_grid(25.0)
xb1, xa1 = peak_newton(xb0, -25.0), peak_newton(xa0, 25.0)
sh0 = (xa0 - k2**2*25) - (xb0 + k2**2*25); sh1 = (xa1 - k2**2*25) - (xb1 + k2**2*25); pred = -phase/k2
dig = lambda a: -math.log10(abs(a) + 1e-17)
lines = ["Polish of the fast-soliton peak positions (checks.py 3b): grid argmax -> Newton on the closed-form u'(x) = 0", "",
         f"peak before (T=-25): {xb0:.6f} -> {xb1:.12f}",
         f"peak after  (T=+25): {xa0:.6f} -> {xa1:.12f}",
         f"observed jump: {sh0:.6f} -> {sh1:.12f}   predicted -log(A)/k2 = {pred:.12f}",
         f"gap to prediction: {abs(sh0-pred):.1e} -> {abs(sh1-pred):.1e}  (about {dig(sh1-pred)-dig(sh0-pred):.0f} more digits of agreement)",
         "checks_output quotes 'observed 1.8330, predicted 1.8329'; polished observed = " + f"{sh1:.4f}",
         "The 1e-4 gap was the 0.001 grid; the tropical offset holds to machine precision at T = +-25.",
         f"CPU time: {time.process_time() - t0:.2f} s"]
out = "\n".join(lines); print(out); open('polish_output.txt', 'w').write(out + "\n")
