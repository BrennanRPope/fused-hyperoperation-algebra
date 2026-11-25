# Binary pulsar spin precession through the twist bookkeeping (velocity twist + space twist + companion's drag twist)
import numpy as np, sympy as sp, math
ok=0; tot=0
def check(name, cond):
    global ok, tot
    tot+=1; ok+=bool(cond); print(('PASS ' if cond else 'FAIL ')+name)
m1,m2 = sp.symbols('m1 m2', positive=True); M = m1+m2
# relative velocity v = v1 - v2; centre of mass: v1 = (m2/M) v, v2 = -(m1/M) v
w1 = m2/M; w2 = -m1/M
vel_twist = sp.Rational(1,2)*w1          # Thomas-type: own velocity through companion's field
space_twist = 1*w1                        # curved space, own velocity
drag_twist = -2*w2                        # companion's moving mass drags frames (gravitomagnetic)
coef = sp.simplify(m2*(vel_twist+space_twist+drag_twist))
check("three twists sum to the standard spin-orbit coefficient m₂(4m₁ + 3m₂)/(2M)", sp.simplify(coef - m2*(4*m1+3*m2)/(2*M))==0)
check("test-gyroscope limit m₁ → 0: coefficient → (3/2)m₂, the Gravity Probe B geodetic law", sp.limit(coef, m1, 0)==sp.Rational(3,2)*m2)
check("in that limit the drag twist vanishes (a gyroscope can't drag its heavy partner)", sp.limit(m2*drag_twist, m1, 0)==0)

# orbit average for an eccentric Kepler orbit: <|r×v|/r³> = n·... ; check the closed form numerically
def kepler_avg(e, N=200000):
    # units: G M = 1, a = 1, n = 1
    Mm = np.linspace(0, 2*np.pi, N, endpoint=False)
    E = Mm.copy()
    for _ in range(50): E = E - (E - e*np.sin(E) - Mm)/(1 - e*np.cos(E))
    r = 1 - e*np.cos(E); h = math.sqrt(1-e*e)       # |r×v| = sqrt(GM a(1-e²))
    return np.mean(h/r**3)
for e in (0.0, 0.0877775, 0.6171):
    got = kepler_avg(e); want = 1/(1-e*e)
    check(f"orbit average of |r×v|/r³ = n²·... → 1/(1−e²) in units GM=a=1 (e={e}): {got:.8f} vs {want:.8f}", abs(got-want) < 1e-8)

# Double pulsar J0737-3039, precession of B (p = B, companion c = A)
Tsun = 4.925490947e-6                     # GM_sun/c³ in seconds
mA, mB = 1.3381, 1.2489; Pb = 0.10225156248*86400; e = 0.0877775
n = 2*math.pi/Pb; Mt = mA+mB
def rate(mp, mc):  # rad/s
    return n**(5/3)*Tsun**(2/3)*mc*(4*mp+3*mc)/(2*Mt**(4/3))/(1-e*e)
degyr = 180/math.pi*365.25*86400
OmB = rate(mB, mA)*degyr; OmA = rate(mA, mB)*degyr
print(f"   predicted: B {OmB:.4f} deg/yr (period {360/OmB:.1f} yr), A {OmA:.4f} deg/yr (period {360/OmA:.1f} yr)")
check("double pulsar B: twist bookkeeping gives ≈ 5.07 °/yr (precession period ≈ 71 yr)", abs(OmB-5.07) < 0.01)
check("measured (Breton et al. 2008): 4.77 +0.66 −0.65 °/yr: within 1σ", abs(OmB-4.77) < 0.66)
parts = {'velocity': 1.5*0+0.5*mA/Mt*mA, 'space': mA/Mt*mA, 'drag': 2*mB/Mt*mA}
S = sum(parts.values())
for k,v in parts.items(): print(f"   {k} twist: {v/S*100:.1f}% = {OmB*v/S:.2f} deg/yr")
no_drag = OmB*(parts['velocity']+parts['space'])/S; no_vel = OmB*(parts['space']+parts['drag'])/S
check(f"without the companion's drag twist: {no_drag:.2f} °/yr, off by {(4.77-no_drag)/0.65:.1f}σ (needed)", (4.77-no_drag)/0.65 > 3)
check(f"without the velocity twist: {no_vel:.2f} °/yr, still within 1σ (pulsar can't isolate it)", abs(no_vel-4.77) < 0.66)
# strength of field: M/r
a = (Tsun*Mt*299792458.0*(Pb/(2*math.pi))**2*299792458.0**2)**(1/3)   # GM = Tsun*c³*Mt
x = Tsun*Mt*299792458.0/a
print(f"   separation a = {a/1e3:.0f} km, M/r = {x:.2e} (Gravity Probe B: {3.986004418e14/299792458.0**2/7027.4e3:.2e})")
check("field strength M/r ≈ 4e-6: 10,000× Earth orbit's, but the ✳ corner (second order) is still ~1e-6 of the twist", 3e-6 < x < 6e-6)
check("so the exact fusion's corner is far below the 14% measurement error: not testable here", x < 0.14/1000)
# Hulse–Taylor B1913+16 prediction (no precise measured rate; profile changes seen)
mp, mc, Pb2, e2 = 1.438, 1.390, 0.322997448918*86400, 0.6171340
n2 = 2*math.pi/Pb2; Mt2 = mp+mc
Om_HT = n2**(5/3)*Tsun**(2/3)*mc*(4*mp+3*mc)/(2*Mt2**(4/3))/(1-e2*e2)*degyr
print(f"   Hulse–Taylor predicted {Om_HT:.3f} deg/yr")
check("Hulse–Taylor pulsar: same bookkeeping predicts ≈ 1.21 °/yr (profile changes observed, consistent)", abs(Om_HT-1.21) < 0.01)
print(f"\n{ok}/{tot} checks pass")
