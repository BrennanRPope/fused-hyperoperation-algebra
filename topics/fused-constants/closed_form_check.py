# Checks two claims on the page: the closed form on the second unit of height, and the n-copy tower of 1+1/n
def Hf(k,b,x):
    if x>200: return float('inf')
    if k==3: return b**x
    if x<=0: return 1+x
    return Hf(k-1,b,Hf(k,b,x-1))
worst=0
for b in (1.2,1.5,1.8,1.99):
    s=lambda x,b=b: b**(x-1)
    for k in range(3,15):
        for i in range(0,101):
            x=1+i/100; y=x
            for _ in range(k-3): y=s(y)
            worst=max(worst,abs(Hf(k,b,x)-b**y))
print("closed form H_k(b,x) = b^(s^(k-3)(x)) on [1,2], worst error:",worst)
for n in (10,100,1000,10000):
    b=1+1/n; t=b
    for _ in range(n-1): t=b**t
    print(f"tower of {n} copies of 1+1/{n}: {t:.6f}")
