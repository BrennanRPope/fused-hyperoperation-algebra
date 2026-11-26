"""Round 7 checks: every number quoted in the README, re-read from the saved runs or recomputed."""
import json, math
import numpy as np, mpmath as mp
from sympy import divisors, mobius
PG = json.load(open('prime_gas.json')); FI = json.load(open('fisher.json')); HP = json.load(open('hilbert_polya.json'))
RM = json.load(open('rmt.json')); TR = json.load(open('trace.json')); SS = json.load(open('solved_shapes.json'))
CC = json.load(open('count_check.json'))
res = []
def check(name, ok, detail=''):
    res.append(ok); print(('PASS ' if ok else 'FAIL ') + name + ('  ' + detail if detail else ''))

# prime gas
for q, v in PG['identity'].items():
    gaps = [abs(r[2]) for r in v['partial']]
    check(f'identity at beta=2, q={q}: partial-sum gap halves per doubling', all(0.4 < gaps[i+1]/gaps[i] < 0.6 for i in range(3)), str([f'{g:.1e}' for g in gaps]))
g = [abs(r[2]) for r in PG['q_minus1']['partial']]
check('q=-1 closed form (1+P(b))/(1-P(2b)): gap halves', all(0.4 < g[i+1]/g[i] < 0.6 for i in range(3)))
check('q=-1 coefficient rule, n < 3000', PG['q_minus1']['coef_rule_mismatches'] == 0)
check('Hagedorn wall P(b) = 1 at 1.39943', abs(PG['hagedorn'] - 1.3994333) < 1e-6, f"{PG['hagedorn']:.7f}")
for q, v in PG['exponent'].items():
    check(f'exponent 1/(1-q) at q={q}', abs(v['slope'] - v['pred']) < 2e-4, f"{v['slope']:.6f} vs {v['pred']:.6f}")
mp.mp.dps = 20
for q in [0.3, 0.5, 0.7]:
    qq = mp.mpf(q)
    c0 = sum(sum(mobius(N//m)/(1-qq**m) for m in divisors(N))/N*mp.log(mp.zeta(N)) for N in range(2, 80))
    amp = PG['exponent'][str(q)]['qpoch']*float(mp.e**c0)
    check(f'amplitude (q;q)_inf exp(c0) at q={q}', abs(amp/PG['exponent'][str(q)]['amp'] - 1) < 1e-3, f"{amp:.5f} vs {PG['exponent'][str(q)]['amp']:.5f}")
for q in ['0.2', '0.5', '0.8']:
    check(f'no Fisher zeros for q={q} in 0.505<Re<1.6, 1<Im<50', len(FI[q]['zeros']) == 0)
for q in ['-0.5', '-0.8']:
    z = FI[q]['zeros']
    check(f'q={q}: one Fisher zero beside each of the first 10 Riemann zeros, all Re > 1/2',
          len(z) == 10 and all(r[0] > 0.55 for r in z) and all(abs(r[1]-g) < 0.5 for r, g in zip(z, FI['riemann'])))
# Hilbert-Polya
check('BBM intertwining for zeta', HP['bbm_zeta']['intertwine_err'] < 1e-25)
check('DH functional equation', HP['dh_fe'] < 1e-25)
for d in HP['dh']:
    check(f"DH zero {d['zero'][0]:.6f}+{d['zero'][1]:.4f}i: psi(0)=0, intertwining, E non-real",
          d['psi0'] < 1e-25 and d['inter'] < 1e-25 and d['E'][1] > 0.1)
check('x^2+5y^2 genus formula = direct lattice sum at s=2', abs(HP['x2+5y2']['genus_s2'] - HP['x2+5y2']['direct_s2']) < 1e-12)
check('x^2+5y^2 zero matches round 6', abs(HP['x2+5y2']['zero'][0]-0.93297) < 1e-5 and abs(HP['x2+5y2']['zero'][1]-15.66825) < 1e-5)
# RMT
check('zeta: KS to GUE below 0.035 and far below two-GUE', RM['zeta']['ks'][0] < 0.035 and RM['zeta']['ks'][1] > 0.15)
check('x^2+y^2: KS to two-GUE below 0.035, to GUE above 0.12', RM['x2+y2']['ks'][1] < 0.035 and RM['x2+y2']['ks'][0] > 0.12)
check('x^2+5y^2 on-line zeros fit neither (rescaled KS > 0.08 both)', min(RM['x2+5y2']['ks_rescaled'][:2]) > 0.08)
check('DH on-line zeros fit neither (rescaled KS > 0.08 both)', min(RM['DH']['ks_rescaled'][:2]) > 0.08)
check('Euler-product cases lose no zeros (|deficit| < 2)', all(abs(RM[k]['deficit']) < 2 for k in ['zeta', 'L-4', 'x2+y2']))
check('x^2+5y^2: 498 of 1855 off the line in [500,1400]', round(RM['x2+5y2']['deficit']) == 498)
check('DH: 176 of 1987 off the line in [500,2300]', round(RM['DH']['deficit']) == 176)
offs = {k: v[0]-v[1] for k, v in CC.items()}
check('argument-principle windows: off-line counts even, x^2+y^2 none',
      all(abs(o - round(o)) < 0.01 and round(o) % 2 == 0 for o in offs.values()) and round(offs['x2+y2 500-520']) == 0, str({k: round(v) for k, v in offs.items()}))
# fast Hurwitz
from fastl import zeta as fz
pts = [0.5+1000j, 0.7+15j, 2+0j]
err = max(abs(complex(fz(np.array([p]))[0]) - complex(mp.zeta(p))) for p in pts)
check('vectorised Hurwitz zeta vs mpmath', err < 1e-8, f'{err:.1e}')
# trace
check('primes < 3000 rebuild the first 10 zeros to 0.05', all(abs(a-b) < 0.05 for a, b in zip(TR['peaks'][:10], TR['rz'][:10])))
lp = [math.log(n) for n in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23)]
check('zeros rebuild log of prime powers to 0.01', len(TR['dips']) == 13 and all(abs(a-b) < 0.01 for a, b in zip(TR['dips'], lp)))
check('cyclotomic identity at s = 3', all(abs(r[1]-r[2]) < 1e-10 for r in TR['cyclotomic']))
check('Witt numbers t=2: 2,1,2,3,6,9,18,30', [round(x) for x in TR['mult_t']['2']] == [2, 1, 2, 3, 6, 9, 18, 30])
check('whole multiplicities exactly at whole t (grid -2..4)', TR['integral_t'] == [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
# solved shapes
r = np.array(SS['graphs']['petersen']['roots']); m = np.abs(r[:, 0]+1j*r[:, 1])
check('Petersen: Ihara roots on |u| = 1/sqrt2', np.allclose(m, 1/math.sqrt(2)) and SS['graphs']['petersen']['ramanujan'])
r = np.array(SS['graphs']['necklace10']['roots']); m = np.abs(r[:, 0]+1j*r[:, 1])
check('necklace: not Ramanujan, roots leave the circle', (not SS['graphs']['necklace10']['ramanujan']) and m.max() > 0.9)
check('Hasse for y^2=x^3+x+1, p < 2000', SS['elliptic']['all_hasse'])
z = np.array(SS['lee_yang']['petersen ferro J=0.5']); check('Lee-Yang ferro zeros on |z|=1', np.allclose(np.abs(z[:, 0]+1j*z[:, 1]), 1, atol=1e-6))
z = np.array(SS['lee_yang']['petersen antiferro J=-0.5']); check('Lee-Yang antiferro zeros off the circle', np.abs(np.abs(z[:, 0]+1j*z[:, 1])-1).max() > 0.5)
print(f"\n{sum(res)} of {len(res)} checks pass")
