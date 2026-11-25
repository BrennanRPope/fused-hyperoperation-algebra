"""Light and sound in a moving medium, through the circle operation.
1) Fizeau: light in moving water. Relativistic speed addition equals
   the circle operation of (Doppler factor - 1). 2) Draining-bathtub acoustic metric:
   horizon and ergosphere radii."""
import math
c = 299_792_458.0; n = 1.333; v = 7.06   # Fizeau 1851 water speed (approx)
w_rel = (c/n + v)/(1 + v/(n*c))
drag = (w_rel - c/n)/v
print(f"Fizeau: light speed in still water {c/n:,.1f} m/s; in water moving at {v} m/s: {w_rel:,.4f} m/s")
print(f"  drag fraction = {drag:.5f}; Fresnel 1 - 1/n^2 = {1-1/n**2:.5f}")
D = lambda beta: math.sqrt((1+beta)/(1-beta))
mult_add = lambda a, b: a + b + a*b
a, b = D(1/n) - 1, D(v/c) - 1
beta_from = lambda Dm1: ((Dm1+1)**2 - 1)/((Dm1+1)**2 + 1)
print(f"  via the circle operation of Doppler factors: {beta_from(mult_add(a, b))*c:,.4f} m/s (matches: {abs(beta_from(mult_add(a,b))*c - w_rel) < 1e-6})")
# Draining bathtub: v = (-A r_hat + B theta_hat)/r, sound speed cs
A, B, cs = 1.0, 2.0, 1.0
print(f"Bathtub: horizon r = A/c = {A/cs:.3f}; ergosphere r = sqrt(A^2+B^2)/c = {math.hypot(A,B)/cs:.3f}")
print("  inside the ergosphere the swirl (spin) alone outruns sound: no wave can stay still there")
