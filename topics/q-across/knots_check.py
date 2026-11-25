"""Kauffman bracket from PD codes, checked against the standard Jones polynomials."""
import sympy as sp
A,t=sp.symbols('A t')
def loops(pairs):
    # pairs: list of 2-tuples of labels; count cycles
    parent={}
    def f(x):
        parent.setdefault(x,x)
        while parent[x]!=x: x=parent[x]
        return x
    labs=set()
    for a,b in pairs:
        labs|={a,b}; ra,rb=f(a),f(b)
        if ra!=rb: parent[ra]=rb
    return len({f(x) for x in labs})
def bracket(pd,conv):
    n=len(pd); tot=0
    for s in range(2**n):
        pairs=[];na=0
        for i,(a,b,c,d) in enumerate(pd):
            if (s>>i)&1==0:
                na+=1; pairs+= [(a,b),(c,d)] if conv==0 else [(a,d),(b,c)]
            else:
                pairs+= [(a,d),(b,c)] if conv==0 else [(a,b),(c,d)]
        L=loops(pairs)
        tot+=A**(na-(n-na))*(-A**2-A**-2)**(L-1)
    return sp.expand(tot)
knots={
 'unknot_kink':([(1,1,2,2)],None),
 'trefoil_KA':([(1,4,2,5),(3,6,4,1),(5,2,6,3)],None),
 'fig8':([(4,2,5,1),(8,6,1,5),(6,3,7,4),(2,7,3,8)],None),
 'hopf':([(4,1,3,2),(2,3,1,4)],None),
 'cinq':([(1,6,2,7),(3,8,4,9),(5,10,6,1),(7,2,8,3),(9,4,10,5)],None),
}
for k,(pd,_) in knots.items():
    for conv in (0,1):
        b=bracket(pd,conv)
        print(k,conv,b)
