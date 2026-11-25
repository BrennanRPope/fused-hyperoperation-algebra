import json
r = json.load(open('/mnt/project-files/math-theory/pi/pi_results.json'))
c = json.load(open('/mnt/project-files/math-theory/pi/constants_results.json'))
half = json.load(open('/mnt/project-files/math-theory/pi/half_curve.json'))
v = c['values']
def short(x, n=22): return x[:n] + '…'
constants = [
 dict(name='π, and every turning period 2π/√(4d − p²)', where='the trip from 0 to −2 round the wall; period of each turning law', value='3.14159265358979…', search='found X = 2π for the period at rung −1 (detector check)', status='p proved', why='Lindemann, 1882. Known.'),
 dict(name='2πi', where="The circle operation's own period: log(1 + a) once round the wall", value='6.28318530717958… i', search='(same as π)', status='p proved', why='Lindemann. Known.'),
 dict(name='e^(−2π)', where='gluing number q of the most symmetric elliptic curve, y² = x³ − x', value=short(v['gluing_q_tau_i']), search='found log q = −2π', status='p proved', why='Gelfond–Schneider, 1934: e^π = (−1)^(−i). Known.'),
 dict(name='−e^(−π√3)', where='gluing number q of y² = x³ + 1, the other symmetric curve', value=short(v['gluing_q_tau_rho']), search='found log|q| = −π√3', status='p proved', why='Gelfond–Schneider: e^(π√3) = (−1)^(−i√3). Known.'),
 dict(name='Ω = W(1)', where='fixed point of the fused power with base 1/e (y = (1/e)^y)', value=short(v['omega_fused_fixed']), search='found Ω + log Ω = 0', status='p proved', why='Lindemann–Weierstrass: if Ω were algebraic, e^(−Ω) = Ω would be too. Known.'),
 dict(name='e^(1/e)', where='rung-3 kiss (fused-constants page)', value=short(v['kiss_rung3']), search='no relation', status='o open', why='Believed transcendental; no proof I know of. Known number.'),
 dict(name='1.6353244967…', where='rung-4 kiss (Kneser heights, fused-constants page)', value='1.63532449671527…', search='not rerun here (published to 24 digits)', status='o open', why='Known number (Tetration Forum). Status unknown.'),
 dict(name='E½(1)', where='half-rung: the half-iterate of eˣ − 1, at 1', value=short(v['E_half_at_1']), search='no relation; not algebraic up to degree 12', status='o open', why='Known value: Finch 2025, Morris and Szekeres 1961.'),
 dict(name='B½(1,1)', where='halfway between 1 + 1 = 2 and 1 ∘₁ 1 = 3', value=short(v['B_half_1_1']), search='no relation; not algebraic up to degree 12', status='o open', why='No published value found. The new one.'),
 dict(name='Feigenbaum δ', where='the road from circle doubling (c = 0) to turning doubling (c = −2)', value='4.6692016091…', search='not rerun (known constant)', status='o open', why='Feigenbaum, 1975. Not even proved irrational.'),
]
D = dict(halving=r['halving'], steps=r['step_count'], periods=r['turning_periods'], digits=r['digits_1000'], constants=constants, half=half)
t = open('/mnt/project-files/math-theory/pi/page/template.html').read()
t = t.replace("Six halvings would reach Archimedes\\u2019 96-gon (48 steps) and beyond.", "Four halvings reach 48 steps, exactly Archimedes\\u2019 96-gon.")
assert '__DATA__' in t
js = json.dumps(D, ensure_ascii=False).replace('</', '<\\/')
t = t.replace('__DATA__', js)
open('/mnt/project-files/math-theory/pi/page/index.html', 'w').write(t)
print(len(t))
