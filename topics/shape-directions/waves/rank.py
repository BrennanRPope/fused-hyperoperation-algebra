import sympy as sp, random, json
random.seed(7)
def A(a,b,i,j): return (a[i]-a[j])*(b[i]-b[j])/((a[i]-b[j])*(b[i]-a[j]))
out={}
for g in range(2,9):
    a=sp.symbols(f'a0:{g}'); b=sp.symbols(f'b0:{g}')
    pairs=[(i,j) for i in range(g) for j in range(i+1,g)]
    F=sp.Matrix([A(a,b,i,j) for i,j in pairs])
    vars_=list(a)+list(b)
    J=F.jacobian(vars_)
    ranks=[]
    for trial in range(3):
        while True:
            vals={v:sp.Rational(random.randint(-500,500),random.randint(1,50)) for v in vars_}
            M=J.subs(vals)
            if len(set(vals.values()))==len(vals) and not M.has(sp.zoo,sp.nan,sp.oo): break
        ranks.append(M.rank())
    # KdV: b=-a
    Fk=F.subs({b[i]:-a[i] for i in range(g)}); Jk=Fk.jacobian(list(a))
    rk=[Jk.subs({v:sp.Rational(random.randint(1,900),random.randint(1,50)) for v in a}).rank() for _ in range(3)]
    exp_kp=min(len(pairs),2*g-3); exp_kdv=min(len(pairs),g-1)
    out[g]=dict(pairs=len(pairs),kp_rank=ranks,expect_kp=exp_kp,kdv_rank=rk,expect_kdv=exp_kdv)
    print(g,out[g],'OK' if all(r==exp_kp for r in ranks) and all(r==exp_kdv for r in rk) else 'FAIL')
json.dump(out,open('rank.json','w'),indent=1)
