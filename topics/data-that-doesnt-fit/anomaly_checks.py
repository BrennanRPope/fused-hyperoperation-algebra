import math
print("=== Muon puzzle: Heitler-Matthews cascade as k-fold circle operation ===")
xi=10.0; Nm=50.0
for E in [1e7,1e8,1e9,1e10]:  # GeV
    k=math.log(E/xi)/math.log(Nm)
    print(f"E={E:.0e} GeV  generations k={k:.2f}")
k10=math.log(1e10/xi)/math.log(Nm)
f=1/3; a=1-f
print("hadronic fraction (1-f)^k at 10 EeV:", round(a**k10,4), " linear (no cross terms) 1-k f:", round(1-k10*f,3))
for ex in [1.3,1.6]:
    r=ex**(1/k10)
    print(f"excess {ex}: alpha'/alpha={r:.4f}, f' (from 1/3)={1-a*r:.4f}")
# only generations above 8 TeV cms
mp=0.938; E=1e10; j=0
while math.sqrt(2*mp*E)>8000: j+=1; E/=Nm
print("generations above 8 TeV cms in 10 EeV shower:", j)
for ex in [1.3,1.6]:
    r=ex**(1/j); print(f"  only those changed, excess {ex}: f'={1-a*r:.4f}")
# Pascal layers of the cascade cube, k=5
k=5
print("layer weights C(k,j)(-f)^j:", [round(math.comb(k,jj)*(-f)**jj,4) for jj in range(k+1)], "sum", round(sum(math.comb(k,jj)*(-f)**jj for jj in range(k+1)),4))
print("energy lost at generation j: f(1-f)^(j-1):", [round(f*a**(jj-1),4) for jj in range(1,k+1)])
print()
print("=== Redshift composition vs Hubble tension ===")
c=299792.458
gap=73.2/67.4-1; print("H0 gap SH0ES/Planck:", round(gap*100,2),"%")
for zc in [0.023,0.05,0.15]:
    for v in [300,369]:
        zp=v/c; zobs=(1+zc)*(1+zp)-1
        zadd=zobs-zp
        print(f"z_cos={zc} v={v}: additive z err rel={ (zadd-zc)/zc*100:.3f}%")
print()
print("=== Neutron lifetime ===")
bottle=878.4; beam=887.7
b=1-bottle/beam; print("dark branch needed:", round(b*100,3),"%")
e=b/2; print("two equal hidden losses e each, cross term e^2:", round(e,5), round(e*e,7), "s-equivalent:", round(e*e*bottle,4))
print()
print("=== Cosmic dipole ===")
beta=369.8/c; x=1.0; al=0.75
D=(2+x*(1+al))*beta; print("kinematic D:",round(D,5)," beta^2:",f"{beta**2:.2e}", " needed extra factor 2.7-3.7")
print()
print("=== Lithium: two destruction channels ===")
s=1/3.5; d=1-math.sqrt(s); print("total destruction",round(1-s,3)," each",round(d,3)," cross term",round(d*d,3))
print()
print("=== Hubble: two partial fixes ===")
for aa,bb in [(0.04,0.045),(0.06,0.025)]:
    print(aa,bb,"The circle operation",round(aa+bb+aa*bb,4),"add",aa+bb,"cross",round(aa*bb,4))
